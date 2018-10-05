#include "mwabeam.h"

#include "../units/angle.h"
#include "../units/radeccoord.h"

#include "../banddata.h"
#include "../fftresampler.h"
#include "../fitsreader.h"
#include "../fitswriter.h"
#include "../units/imagecoordinates.h"
#include "../imageweights.h"
#include "../matrix2x2.h"
#include "../progressbar.h"
#include "../uvector.h"

#include "../wsclean/imageweightcache.h"
#include "../wsclean/logger.h"
#include "../msproviders/msprovider.h"
#include "../multibanddata.h"

#include "tilebeambase.h"
#include "tilebeam2016.h"

#include <casacore/ms/MeasurementSets/MSField.h>

#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/TableRecord.h>
#include <casacore/tables/Tables/ScalarColumn.h>

#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MCPosition.h>

#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>
#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>

#include <stdexcept>

void MWABeam::Make(PrimaryBeamImageSet& beamImages)
{
	_sampledWidth = _width / _undersample;
	_sampledHeight = _width / _undersample;
	_sPixelSizeX = _pixelSizeX * _undersample;
	_sPixelSizeY = _pixelSizeX * _undersample;

	beamImages.SetToZero();
	
	_totalWeightSum = 0.0;
	
	Logger::Debug << "Making beam for " << _msProviders.size() << " parts (" << _tableEntry->msData.size() << " ms)\n";
	for(const MSProviderInfo& msProviderInfo : _msProviders)
	{
		const ImagingTableEntry::MSInfo& msInfo = _tableEntry->msData[msProviderInfo.msIndex];
		const MSSelection& selection = *msProviderInfo.selection;
		casacore::MeasurementSet& ms = msProviderInfo.provider->MS();
		MultiBandData band(ms.spectralWindow(), ms.dataDescription());
		double centralFrequency = 0.0;
		for(size_t dataDescId=0; dataDescId!=band.DataDescCount(); ++dataDescId)
		{
			BandData subBand(band[dataDescId], selection.ChannelRangeStart(), selection.ChannelRangeEnd());
			centralFrequency += subBand.CentreFrequency();
		}
		centralFrequency /= msInfo.bands.size();
		makeBeamForMS(beamImages, *msProviderInfo.provider, centralFrequency);
	}

	for(size_t i=0; i!=8; ++i)
	{
		for(size_t j=0; j!=_sampledWidth*_sampledHeight; ++j)
		{
			beamImages[i][j] /= _totalWeightSum;
		}
	}
	
	if(_width!=_sampledWidth || _height!=_sampledHeight)
	{
		FFTResampler resampler(_sampledWidth, _sampledHeight, _width, _height, 1);
		ImageBufferAllocator::Ptr scratch;
		_allocator->Allocate(_width*_height, scratch);
		for(size_t p=0; p!=8; ++p)
		{
			resampler.RunSingle(&beamImages[p][0], scratch.data());
			memcpy(&beamImages[p][0], scratch.data(), sizeof(double)*_width*_height);
		}
	}
}

void MWABeam::makeBeamForMS(PrimaryBeamImageSet& beamImages, MSProvider& msProvider, double centralFrequency)
{
	/**
		* Read meta data from the measurement set
		*/
	casacore::MeasurementSet& ms = msProvider.MS();
	casacore::MSAntenna aTable = ms.antenna();
	if(aTable.nrow() == 0) throw std::runtime_error("No antennae in set");
	
	casacore::MPosition::ROScalarColumn antPosColumn(aTable, aTable.columnName(casacore::MSAntennaEnums::POSITION));
	casacore::MPosition arrayPos = antPosColumn(0);
		
	casacore::Table mwaTilePointing = ms.keywordSet().asTable("MWA_TILE_POINTING");
	casacore::ROArrayColumn<int> delaysCol(mwaTilePointing, "DELAYS");
	casacore::Array<int> delaysArr = delaysCol(0);
	casacore::Array<int>::contiter delaysArrPtr = delaysArr.cbegin();
	for(int i=0; i!=16; ++i)
		_delays[i] = delaysArrPtr[i];
		
	Logger::Debug << "Making MWA beam for frequency " << centralFrequency * 1e-6 << " MHz.\n";
	
	Logger::Debug << "Delays: [";
	for(int i=0; i!=16; ++i)
	{
		Logger::Debug << _delays[i];
		if(i != 15) Logger::Debug << ',';
	}
	Logger::Debug << "]\n";
	
	Logger::Debug << "Counting timesteps...\n";
	msProvider.Reset();
	size_t timestepCount = 0;
	double startTime = 0.0, endTime = 0.0;
	if(msProvider.CurrentRowAvailable())
	{
		MSProvider::MetaData meta;
		msProvider.ReadMeta(meta);
		startTime = meta.time;
		endTime = meta.time;
		++timestepCount;
		msProvider.NextRow();
		while(msProvider.CurrentRowAvailable())
		{
			msProvider.ReadMeta(meta);
			if(endTime != meta.time)
			{
				++timestepCount;
				endTime = meta.time;
			}
			msProvider.NextRow();
		}
	}
	const double totalSeconds = endTime - startTime;
	size_t intervalCount = (totalSeconds + _secondsBeforeBeamUpdate - 1) / _secondsBeforeBeamUpdate;
	if(intervalCount > timestepCount)
		intervalCount = timestepCount;
	Logger::Debug << "MS spans " << totalSeconds << " seconds, dividing in " << intervalCount << " intervals.\n";
	
	casacore::MEpoch::ROScalarColumn timeColumn(ms, ms.columnName(casacore::MSMainEnums::TIME));
	casacore::MEpoch midTime(casacore::MVEpoch((0.5/86400.0) * (startTime + endTime)), timeColumn(0).getRef());
	Logger::Debug << "Mid time for full selection: " << midTime << '\n';
	casacore::MeasFrame midFrame(arrayPos, midTime);
	const casacore::MDirection::Ref hadecRef(casacore::MDirection::HADEC, midFrame);
	const casacore::MDirection::Ref azelgeoRef(casacore::MDirection::AZELGEO, midFrame);
	const casacore::MDirection::Ref midJ2000Ref(casacore::MDirection::J2000, midFrame);
	
	msProvider.Reset();
	for(size_t intervalIndex=0; intervalIndex!=intervalCount; ++intervalIndex)
	{
		// Find the mid time step
		double firstTime = startTime + (endTime - startTime) * intervalIndex / intervalCount;
		double lastTime = startTime + (endTime - startTime) * (intervalIndex+1) / intervalCount;
		casacore::MEpoch timeEpoch = casacore::MEpoch(casacore::MVEpoch((0.5/86400.0)*(firstTime + lastTime)), timeColumn(0).getRef());
		Logger::Debug << "Mid time for this interval: " << timeEpoch << '\n';
		
		casacore::MeasFrame frame(arrayPos, timeEpoch);
		ao::uvector<double> singleImages[8];
		double *imgPtr[8];
		for(size_t i=0; i!=8; ++i)
		{
			singleImages[i].assign(_sampledWidth*_sampledHeight, 0.0);
			imgPtr[i] = &singleImages[i][0];
		}
	
		makeBeamSnapshot(imgPtr, centralFrequency, frame, arrayPos);
	
		_totalWeightSum += 1.0;
		for(size_t i=0; i!=8; ++i)
		{
			for(size_t j=0; j!=_sampledWidth*_sampledHeight; ++j)
			{
				beamImages[i][j] += singleImages[i][j];
			}
		}
	}
}

void MWABeam::makeBeamSnapshot(double** imgPtr, double frequency, casacore::MeasFrame frame, casacore::MPosition arrayPos)
{
	const casacore::MDirection::Ref hadecRef(casacore::MDirection::HADEC, frame);
	const casacore::MDirection::Ref azelgeoRef(casacore::MDirection::AZELGEO, frame);
	const casacore::MDirection::Ref j2000Ref(casacore::MDirection::J2000, frame);
	casacore::MDirection::Convert
		j2000ToHaDecRef(j2000Ref, hadecRef),
		j2000ToAzelGeoRef(j2000Ref, azelgeoRef);
	casacore::MPosition wgs = casacore::MPosition::Convert(arrayPos, casacore::MPosition::WGS84)();
	double arrLatitude = wgs.getValue().getLat();
	
	casacore::MDirection zenith(casacore::MVDirection(0.0, 0.0, 1.0), azelgeoRef);
	casacore::MDirection zenithHaDec = casacore::MDirection::Convert(zenith, hadecRef)();
	double zenithHa = zenithHaDec.getAngle().getValue()[0];
	double zenithDec = zenithHaDec.getAngle().getValue()[1];

	TileBeamBase<TileBeam2016> tilebeam(_delays, _frequencyInterpolation, _searchPath);
	ProgressBar progressBar("Constructing beam");
	for(size_t y=0;y!=_sampledHeight;++y)
	{
		for(size_t x=0;x!=_sampledWidth;++x)
		{
			double l, m, ra, dec;
			ImageCoordinates::XYToLM(x, y, _sPixelSizeX, _sPixelSizeY, _sampledWidth, _sampledHeight, l, m);
			l += _phaseCentreDL; m += _phaseCentreDM;
			ImageCoordinates::LMToRaDec(l, m, _phaseCentreRA, _phaseCentreDec, ra, dec);
			
			std::complex<double> gain[4];
			tilebeam.ArrayResponse(ra, dec, j2000Ref, j2000ToHaDecRef, j2000ToAzelGeoRef, arrLatitude, zenithHa, zenithDec, frequency, gain);
			
			for(size_t i=0; i!=4; ++i)
			{
				*imgPtr[i*2] = gain[i].real();
				*imgPtr[i*2 + 1] = gain[i].imag();
				++imgPtr[i*2];
				++imgPtr[i*2 + 1];
			}
		}
		progressBar.SetProgress(y, _sampledHeight);
	}
	progressBar.SetProgress(_sampledHeight, _sampledHeight);
}
