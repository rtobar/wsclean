#include "lofarbeamterm.h"

#include "../banddata.h"
#include "../matrix2x2.h"

#include "../units/imagecoordinates.h"

#include "../wsclean/logger.h"

#include "../system.h"

#include "../lane.h"
#include "../uvector.h"

#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>
#include <casacore/measures/Measures/MEpoch.h>

#include <thread>

LofarBeamTerm::LofarBeamTerm(casacore::MeasurementSet& ms, size_t width, size_t height, double dl, double dm, double phaseCentreDL, double phaseCentreDM, bool useDifferentialBeam) :
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

struct LofarBeamTermThreadData
{
	std::complex<float>* buffer;
	ao::lane<size_t>* lane;
	casacore::MDirection::Ref j2000Ref;
	casacore::MDirection::Convert j2000ToITRFRef;
	std::vector<MC2x2F> inverseCentralGain;
	double time, frequency;
	vector3r_t station0, tile0;
};

void LofarBeamTerm::Calculate(std::complex<float>* buffer, double time, double frequency)
{
	size_t nCPUs = System::ProcessorCount();
	ao::lane<size_t> lane(nCPUs);
	
	size_t nThreads = 1; // nCPUs TODO
	std::vector<std::thread> threads(nThreads);
	std::vector<LofarBeamTermThreadData> threadData(nThreads);
	for(size_t i=0; i!=nThreads; ++i)
	{
		// Make a private copy of the data so that each thread has its local copy
		// (in particular to make sure casacore objects do not cause sync bugs)
		LofarBeamTermThreadData& data = threadData[i];
		data.buffer = buffer;
		data.lane = &lane;
		data.time = time;
		data.frequency = frequency;
		
		casacore::MEpoch timeEpoch(casacore::Quantity(time, "s"));
		casacore::MeasFrame frame(_arrayPos, timeEpoch);
		data.j2000Ref = casacore::MDirection::Ref(casacore::MDirection::J2000, frame);
		casacore::MDirection::Ref itrfRef(casacore::MDirection::ITRF, frame);
		data.j2000ToITRFRef = casacore::MDirection::Convert(data.j2000Ref, itrfRef);
		
		setITRF(data.j2000ToITRFRef(_delayDir), data.station0);
		setITRF(data.j2000ToITRFRef(_tileBeamDir), data.tile0);
		
		if(_useDifferentialBeam)
		{
			vector3r_t diffBeamCentre;
			setITRF(data.j2000ToITRFRef(_referenceDir), diffBeamCentre);
			data.inverseCentralGain.resize(_stations.size());
			for(size_t a=0; a!=_stations.size(); ++a)
			{
				matrix22c_t gainMatrix = _stations[a]->response(time, frequency, diffBeamCentre, _subbandFrequency, data.station0, data.tile0);
				data.inverseCentralGain[a][0] = gainMatrix[0][0];
				data.inverseCentralGain[a][1] = gainMatrix[0][1];
				data.inverseCentralGain[a][2] = gainMatrix[1][0];
				data.inverseCentralGain[a][3] = gainMatrix[1][1];
				if(!data.inverseCentralGain[a].Invert())
				{
					data.inverseCentralGain[a] = MC2x2F::Zero();
				}
			}
		}
		// It is necessary to use each converter once in the global thread, during which
		// it initializes itself. This initializes is not thread safe, apparently.
		threadData[i].j2000ToITRFRef(_delayDir);
	}
	for(size_t i=0; i!=nThreads; ++i)
	{
		threads[i] = std::thread(&LofarBeamTerm::calcThread, this, &threadData[i]);
	}
	for(size_t y=0; y!=_height; ++y)
	{
		lane.write(y);
	}
	lane.write_end();
	for(size_t i=0; i!=1; ++i)
		threads[i].join();
	
	static int index = 0;
	std::ostringstream f;
	f << "aterm" << index << ".fits";
	StoreATerms(f.str(), buffer);
	++index;
}

void LofarBeamTerm::calcThread(struct LofarBeamTermThreadData* data)
{
	const size_t valuesPerAntenna = _width * _height * 4;
	const casacore::Unit radUnit("rad");
	
	size_t y;
	while(data->lane->read(y))
	{
		for(size_t x=0; x!=_width; ++x)
		{
			double l, m, ra, dec;
            double l0, m0, l1, m1;
			ImageCoordinates::XYToLM(x, y, _dl, _dm, _width, _height, l0, m0);
			ImageCoordinates::XYToLM(x+1, y+1, _dl, _dm, _width, _height, l1, m1);
			l = (l0+l1)/2 + _phaseCentreDL;
            m = (m0+m1)/2 + _phaseCentreDM;
			ImageCoordinates::LMToRaDec(l, m, _phaseCentreRA, _phaseCentreDec, ra, dec);
			
			casacore::MDirection imageDir(casacore::MVDirection(
				casacore::Quantity(ra, radUnit),
				casacore::Quantity(dec,radUnit)),
				data->j2000Ref);

			vector3r_t itrfDirection;
			setITRF(data->j2000ToITRFRef(imageDir), itrfDirection);
			
			std::complex<float>* baseBuffer = data->buffer + (x + y*_height) * 4;
			
			for(size_t antennaIndex=0; antennaIndex!=_stations.size(); ++antennaIndex)
			{
				std::complex<float>* antBufferPtr = baseBuffer + antennaIndex*valuesPerAntenna;
				matrix22c_t gainMatrix = _stations[antennaIndex]->response(data->time, data->frequency, itrfDirection, _subbandFrequency, data->station0, data->tile0);
				if(_useDifferentialBeam)
				{
					MC2x2F stationGains;
					stationGains[0] = gainMatrix[0][0];
					stationGains[1] = gainMatrix[0][1];
					stationGains[2] = gainMatrix[1][0];
					stationGains[3] = gainMatrix[1][1];
					MC2x2F::ATimesB(antBufferPtr, data->inverseCentralGain[antennaIndex], stationGains);
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
}

#include "../fitswriter.h"
void LofarBeamTerm::StoreATerms(const std::string& filename, std::complex<float>* buffer)
{
	size_t ny = floor(sqrt(_stations.size())), nx = (_stations.size()+ny-1) / ny;
	Logger::Info << "Storing " << filename << " (" << _stations.size() << " ant, " << nx << " x " << ny << ")\n";
	ao::uvector<double> img(nx*ny * _width*_height, 0.0);
	for(size_t ant=0; ant!=_stations.size(); ++ant)
	{
		size_t xCorner = (ant%nx)*_width, yCorner = (ant/nx)*_height;
		for(size_t y=0; y!=_height; ++y)
		{
			for(size_t x=0; x!=_width; ++x)
			{
				std::complex<float> e1, e2;
				Matrix2x2::EigenValues(buffer+(_width*(ant*_height + y) + x)*4, e1, e2);
				double val = std::max(std::abs(e1), std::abs(e2));
				img[(yCorner + y)*_width*nx + x + xCorner] = val;
			}
		}
	}
	FitsWriter writer;
	writer.SetImageDimensions(nx*_width, ny*_height);
	writer.Write(filename, img.data());
}
