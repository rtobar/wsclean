#include "lofarbeamterm.h"

#include "../banddata.h"
#include "../matrix2x2.h"

#include "../units/imagecoordinates.h"

#include "../wsclean/logger.h"

#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>
#include <casacore/measures/Measures/MEpoch.h>

LofarBeamTerm::LofarBeamTerm(casa::MeasurementSet& ms, size_t width, size_t height, double dl, double dm, double phaseCentreDL, double phaseCentreDM, bool useDifferentialBeam) :
	_width(width),
	_height(height),
	_dl(dl), _dm(dm),
	_phaseCentreDL(phaseCentreDL), _phaseCentreDM(phaseCentreDM),
	_useDifferentialBeam(useDifferentialBeam)
{
	casacore::MSAntenna aTable(ms.antenna());
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

void LofarBeamTerm::Calculate(std::complex<float>* buffer, double time, double frequency)
{
	casacore::MEpoch timeEpoch(casacore::Quantity(time, "s"));
	casacore::MeasFrame frame(_arrayPos, timeEpoch);
	casacore::MDirection::Ref
		j2000Ref(casacore::MDirection::J2000, frame),
		itrfRef(casacore::MDirection::ITRF, frame);
	casacore::MDirection::Convert
		j2000ToITRFRef(j2000Ref, itrfRef);
	
	vector3r_t station0, tile0;
	setITRF(j2000ToITRFRef(_delayDir), station0);
	setITRF(j2000ToITRFRef(_tileBeamDir), tile0);
	
	std::vector<MC2x2F> inverseCentralGain;
	if(_useDifferentialBeam)
	{
		vector3r_t diffBeamCentre;
		setITRF(j2000ToITRFRef(_referenceDir), diffBeamCentre);
		inverseCentralGain.resize(_stations.size());
		for(size_t a=0; a!=_stations.size(); ++a)
		{
			matrix22c_t gainMatrix = _stations[a]->response(time, frequency, diffBeamCentre, _subbandFrequency, station0, tile0);
			inverseCentralGain[a][0] = gainMatrix[0][0];
			inverseCentralGain[a][1] = gainMatrix[0][1];
			inverseCentralGain[a][2] = gainMatrix[1][0];
			inverseCentralGain[a][3] = gainMatrix[1][1];
			if(!inverseCentralGain[a].Invert())
			{
				inverseCentralGain[a] = MC2x2F::NaN();
			}
		}
	}
	
	const size_t valuesPerAntenna = _width * _height * 4;
	
	size_t pxIndex = 0;
	for(size_t y=0; y!=_height; ++y)
	{
		for(size_t x=0; x!=_width; ++x)
		{
			double l, m, ra, dec;
			ImageCoordinates::XYToLM(x, y, _dl, _dm, _width, _height, l, m);
			l += _phaseCentreDL; m += _phaseCentreDM;
			ImageCoordinates::LMToRaDec(l, m, _phaseCentreRA, _phaseCentreDec, ra, dec);
			
			static const casacore::Unit radUnit("rad");
			casacore::MDirection imageDir(casacore::MVDirection(
				casacore::Quantity(ra, radUnit),
				casacore::Quantity(dec,radUnit)),
				j2000Ref);

			vector3r_t itrfDirection;
			setITRF(j2000ToITRFRef(imageDir), itrfDirection);
			
			std::complex<float>* baseBuffer = buffer + pxIndex * 4;
			
			for(size_t antennaIndex=0; antennaIndex!=_stations.size(); ++antennaIndex)
			{
				std::complex<float>* antBufferPtr = baseBuffer + antennaIndex*valuesPerAntenna;
				matrix22c_t gainMatrix = _stations[antennaIndex]->response(time, frequency, itrfDirection, _subbandFrequency, station0, tile0);
				if(_useDifferentialBeam)
				{
					MC2x2F stationGains;
					stationGains[0] = gainMatrix[0][0];
					stationGains[1] = gainMatrix[0][1];
					stationGains[2] = gainMatrix[1][0];
					stationGains[3] = gainMatrix[1][1];
					MC2x2F::ATimesB(antBufferPtr, inverseCentralGain[antennaIndex], stationGains);
				}
				else {
					antBufferPtr[0] = gainMatrix[0][0];
					antBufferPtr[1] = gainMatrix[0][1];
					antBufferPtr[2] = gainMatrix[1][0];
					antBufferPtr[3] = gainMatrix[1][1];
				}
				if(!Matrix2x2::Invert(antBufferPtr))
				{
					antBufferPtr[0] = 0.0; antBufferPtr[1] = 0.0;
					antBufferPtr[2] = 0.0; antBufferPtr[3] = 0.0;
				}
			}
			
			++pxIndex;
		}
	}
}
