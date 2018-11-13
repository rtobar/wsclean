#include "lofarbeamterm.h"

#include "../banddata.h"
#include "../lane.h"
#include "../matrix2x2.h"
#include "../system.h"

#include "../units/imagecoordinates.h"

#include "../wsclean/logger.h"

#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>
#include <casacore/measures/Measures/MEpoch.h>

#include <algorithm>

LofarBeamTerm::LofarBeamTerm(casacore::MeasurementSet& ms, size_t width, size_t height, double dl, double dm, double phaseCentreDL, double phaseCentreDM) :
	_width(width),
	_height(height),
	_dl(dl), _dm(dm),
	_phaseCentreDL(phaseCentreDL),
	_phaseCentreDM(phaseCentreDM),
	_useDifferentialBeam(false),
	_useChannelFrequency(true)
{
	casacore::MSAntenna aTable(ms.antenna());

	size_t nStations = aTable.nrow();
	size_t nCPUs = System::ProcessorCount();
	_nThreads = std::min(nCPUs, nStations);
	_threads.resize(_nThreads);

	casacore::MPosition::ROScalarColumn antPosColumn(aTable, aTable.columnName(casacore::MSAntennaEnums::POSITION));
	_arrayPos = antPosColumn(0);
	_stations.resize(aTable.nrow());
	readStations(ms, _stations.begin());
	
	BandData band(ms.spectralWindow());
	_subbandFrequency = band.CentreFrequency();
	
	casacore::MSField fieldTable(ms.field());
	
	casacore::MEpoch::ROScalarColumn timeColumn(ms, ms.columnName(casacore::MSMainEnums::TIME));
	casacore::MDirection::ROScalarColumn phaseDirColumn(fieldTable, fieldTable.columnName(casacore::MSFieldEnums::PHASE_DIR));
	casacore::MDirection phaseDir = phaseDirColumn(0);
	casacore::MEpoch curtime = timeColumn(0);
	casacore::MeasFrame frame(_arrayPos, curtime);
	casacore::MDirection::Ref j2000Ref(casacore::MDirection::J2000, frame);
	casacore::MDirection j2000 = casacore::MDirection::Convert(phaseDir, j2000Ref)();
	casacore::Vector<casacore::Double> j2000Val = j2000.getValue().get();
	_phaseCentreRA = j2000Val[0];
	_phaseCentreDec = j2000Val[1];

	casacore::ScalarMeasColumn<casacore::MDirection> delayDirColumn(fieldTable, casacore::MSField::columnName(casacore::MSFieldEnums::DELAY_DIR));
	if(fieldTable.nrow() != 1)
		throw std::runtime_error("Set has multiple fields");
	_delayDir = delayDirColumn(0);
	
	casacore::ScalarMeasColumn<casacore::MDirection> referenceDirColumn(fieldTable, casacore::MSField::columnName(casacore::MSFieldEnums::REFERENCE_DIR));
	_referenceDir = referenceDirColumn(0);

	if(fieldTable.tableDesc().isColumn("LOFAR_TILE_BEAM_DIR")) {
		casacore::ArrayMeasColumn<casacore::MDirection> tileBeamDirColumn(fieldTable, "LOFAR_TILE_BEAM_DIR");
		_tileBeamDir = *(tileBeamDirColumn(0).data());
	} else {
		throw std::runtime_error("LOFAR_TILE_BEAM_DIR column not found");
	}
}

void setITRF(const casacore::MDirection& itrfDir, LOFAR::StationResponse::vector3r_t& itrf)
{
	const casacore::Vector<double>& itrfVal = itrfDir.getValue().getValue();
	itrf[0] = itrfVal[0];
	itrf[1] = itrfVal[1];
	itrf[2] = itrfVal[2];
}

bool LofarBeamTerm::calculateBeam(std::complex<float>* buffer, double time, double frequency)
{
	ao::lane<size_t> lane(_nThreads);
	_lane = &lane;

	casacore::MEpoch timeEpoch(casacore::Quantity(time, "s"));
	casacore::MeasFrame frame(_arrayPos, timeEpoch);
	casacore::MDirection::Ref j2000Ref(casacore::MDirection::J2000, frame);
	casacore::MDirection::Ref itrfRef(casacore::MDirection::ITRF, frame);
	casacore::MDirection::Convert j2000ToITRFRef(j2000Ref, itrfRef);

	setITRF(j2000ToITRFRef(_delayDir), _station0);
	setITRF(j2000ToITRFRef(_tileBeamDir), _tile0);

	const casacore::Unit radUnit("rad");

	casacore::MDirection lDir(casacore::MVDirection(
		casacore::Quantity(_phaseCentreRA + M_PI/2, radUnit),
		casacore::Quantity(0, radUnit)),
		j2000Ref);
	setITRF(j2000ToITRFRef(lDir), _l_vector_itrf);

	casacore::MDirection mDir(casacore::MVDirection(
		casacore::Quantity(_phaseCentreRA, radUnit),
		casacore::Quantity(_phaseCentreDec + M_PI/2, radUnit)),
		j2000Ref);
	setITRF(j2000ToITRFRef(mDir), _m_vector_itrf);

	casacore::MDirection nDir(casacore::MVDirection(
		casacore::Quantity(_phaseCentreRA, radUnit),
		casacore::Quantity(_phaseCentreDec, radUnit)),
		j2000Ref);
	setITRF(j2000ToITRFRef(nDir), _n_vector_itrf);

	vector3r_t diffBeamCentre;
	setITRF(j2000ToITRFRef(_referenceDir), diffBeamCentre);
	_inverseCentralGain.resize(_stations.size());
	for(size_t a=0; a!=_stations.size(); ++a)
	{
		double sbFrequency = _useChannelFrequency ? frequency : _subbandFrequency;
		matrix22c_t gainMatrix = _stations[a]->response(time, frequency, diffBeamCentre, sbFrequency, _station0, _tile0);
		if(_useDifferentialBeam)
		{
			_inverseCentralGain[a][0] = gainMatrix[0][0];
			_inverseCentralGain[a][1] = gainMatrix[0][1];
			_inverseCentralGain[a][2] = gainMatrix[1][0];
			_inverseCentralGain[a][3] = gainMatrix[1][1];
			if(!_inverseCentralGain[a].Invert())
			{
				_inverseCentralGain[a] = MC2x2F::Zero();
			}
		}
	}

	for(size_t i=0; i!=_nThreads; ++i)
	{
		_threads[i] = std::thread(&LofarBeamTerm::calcThread, this, buffer, time, frequency);
	}

	for(size_t y=0; y!=_height; ++y)
	{
		for(size_t antennaIndex=0; antennaIndex!=_stations.size(); ++antennaIndex)
		{
			size_t job_id = y*_stations.size() + antennaIndex;
			lane.write(job_id);
		}
	}

	lane.write_end();
	for(size_t i=0; i!=_nThreads; ++i)
		_threads[i].join();
	
	saveATermsIfNecessary(buffer, _stations.size(), _width, _height);

	return true;
}

void LofarBeamTerm::calcThread(std::complex<float>* buffer, double time, double frequency)
{
	const size_t valuesPerAntenna = _width * _height * 4;
	const casacore::Unit radUnit("rad");

	size_t job_id;
	while(_lane->read(job_id))
	{
		size_t antennaIndex = job_id % _stations.size();
		size_t y = job_id / _stations.size();
		for(size_t x=0; x!=_width; ++x)
		{
			double l, m, n, ra, dec;
			ImageCoordinates::XYToLM(x, y, _dl, _dm, _width, _height, l, m);
			l += _phaseCentreDL;
			m += _phaseCentreDM;
			n = sqrt(1.0 - l*l - m*m);

			ImageCoordinates::LMToRaDec(l, m, _phaseCentreRA, _phaseCentreDec, ra, dec);

			vector3r_t itrfDirection;

			itrfDirection[0] = l*_l_vector_itrf[0] + m*_m_vector_itrf[0] + n*_n_vector_itrf[0];
			itrfDirection[1] = l*_l_vector_itrf[1] + m*_m_vector_itrf[1] + n*_n_vector_itrf[1];
			itrfDirection[2] = l*_l_vector_itrf[2] + m*_m_vector_itrf[2] + n*_n_vector_itrf[2];

			std::complex<float>* baseBuffer = buffer + (x + y*_height) * 4;

			std::complex<float>* antBufferPtr = baseBuffer + antennaIndex*valuesPerAntenna;
			double sbFrequency = _useChannelFrequency ? frequency : _subbandFrequency;
			matrix22c_t gainMatrix = _stations[antennaIndex]->response(time, frequency, itrfDirection, sbFrequency, _station0, _tile0);
			if(_useDifferentialBeam)
			{
				MC2x2F stationGains;
				stationGains[0] = gainMatrix[0][0];
				stationGains[1] = gainMatrix[0][1];
				stationGains[2] = gainMatrix[1][0];
				stationGains[3] = gainMatrix[1][1];
				MC2x2F::ATimesB(antBufferPtr, _inverseCentralGain[antennaIndex], stationGains);
			}
			else {
				antBufferPtr[0] = gainMatrix[0][0];
				antBufferPtr[1] = gainMatrix[0][1];
				antBufferPtr[2] = gainMatrix[1][0];
				antBufferPtr[3] = gainMatrix[1][1];
			}
		}
	}
}
