#include "lbeamimagemaker.h"

#ifndef HAVE_LOFAR_BEAM
#include <stdexcept>
void LBeamImageMaker::Make(std::vector<ImageBufferAllocator::Ptr>&)
{
	throw std::runtime_error("LOFAR beam imager maker called, but the software has been compiled without the LOFAR beam. Recompile your software and make sure that cmake finds the LOFAR station response library.");
}
#else

#include "../angle.h"
#include "../banddata.h"
#include "../fftresampler.h"
#include "../fitsreader.h"
#include "../fitswriter.h"
#include "../imagecoordinates.h"
#include "../imageweights.h"
#include "../matrix2x2.h"
#include "../progressbar.h"
#include "../uvector.h"

#include "../wsclean/imageweightcache.h"
#include "../wsclean/logger.h"
#include "../msproviders/msprovider.h"
#include "../multibanddata.h"

#include <StationResponse/LofarMetaDataUtil.h>

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

using namespace LOFAR::StationResponse;

static void dirToITRF(const casacore::MDirection& dir, casacore::MDirection::Convert& convert, vector3r_t& itrf);

void LBeamImageMaker::Make(std::vector<ImageBufferAllocator::Ptr>& beamImages)
{
	_sampledWidth = _width / _undersample;
	_sampledHeight = _width / _undersample;
	_sPixelSizeX = _pixelSizeX * _undersample;
	_sPixelSizeY = _pixelSizeX * _undersample;

	for(size_t i=0; i!=8; ++i)
	{
		for(size_t j=0; j!=_width*_height; ++j)
			beamImages[i][j] = 0.0;
	}
	
	_totalWeightSum = 0.0;
	
	for(size_t i=0; i!=_msProviders.size(); ++i)
	{
		makeBeamForMS(beamImages, *_msProviders[i].provider, _tableEntry->msData[i], *_msProviders[i].selection);
	}

	for(size_t i=0; i!=8; ++i)
	{
		for(size_t j=0; j!=_sampledWidth*_sampledHeight; ++j)
		{
			beamImages[i][j] /= _totalWeightSum;
		}
	}
	
	FFTResampler resampler(_sampledWidth, _sampledHeight, _width, _height, 1);
	ImageBufferAllocator::Ptr scratch;
	_allocator->Allocate(_width*_height, scratch);
	for(size_t p=0; p!=8; ++p)
	{
		resampler.RunSingle(&beamImages[p][0], scratch.data());
		memcpy(&beamImages[p][0], scratch.data(), sizeof(double)*_width*_height);
	}
}

void LBeamImageMaker::makeBeamForMS(std::vector<ImageBufferAllocator::Ptr>& beamImages, MSProvider& msProvider, const ImagingTableEntry::MSInfo& msInfo, const MSSelection& selection)
{
	/**
		* Read some meta data from the measurement set
		*/
	casacore::MeasurementSet& ms = msProvider.MS();
	if(ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
	
	casacore::MSAntenna aTable = ms.antenna();
	if(aTable.nrow() == 0) throw std::runtime_error("No antennae in set");
	casacore::MPosition::ROScalarColumn antPosColumn(aTable, aTable.columnName(casacore::MSAntennaEnums::POSITION));
	casacore::MPosition arrayPos = antPosColumn(0);
	
	MultiBandData band(ms.spectralWindow(), ms.dataDescription());
	double centralFrequency = 0.0;
	for(size_t b=0; b!=msInfo.bands.size(); ++b)
	{
		BandData subBand(band[b], selection.ChannelRangeStart(), selection.ChannelRangeEnd());
		centralFrequency += subBand.CentreFrequency();
	}
	centralFrequency /= msInfo.bands.size();
	Logger::Debug << "Making beam for frequency " << centralFrequency * 1e-6 << " MHz.\n";
	
	casacore::MSField fieldTable(ms.field());
	casacore::ROScalarMeasColumn<casacore::MDirection> delayDirColumn(fieldTable, casacore::MSField::columnName(casacore::MSFieldEnums::DELAY_DIR));
	if(fieldTable.nrow() != 1)
		throw std::runtime_error("Set has multiple fields");
	_delayDir = delayDirColumn(0);
	
	casacore::ROScalarMeasColumn<casa::MDirection> referenceDirColumn(fieldTable, casacore::MSField::columnName(casacore::MSFieldEnums::REFERENCE_DIR));
	_referenceDir = referenceDirColumn(0);
	
	if(fieldTable.tableDesc().isColumn("LOFAR_TILE_BEAM_DIR")) {
		casacore::ROArrayMeasColumn<casacore::MDirection> tileBeamDirColumn(fieldTable, "LOFAR_TILE_BEAM_DIR");
		_tileBeamDir = *(tileBeamDirColumn(0).data());
	} else {
		throw std::runtime_error("LOFAR_TILE_BEAM_DIR column not found");
	}
	
	std::vector<Station::Ptr> stations(aTable.nrow());
	readStations(ms, stations.begin());
	
	Logger::Debug << "Counting timesteps...\n";
	std::vector<size_t> idToMSRow;
	msProvider.MakeIdToMSRowMapping(idToMSRow);
	
	casacore::MEpoch::ROScalarColumn timeColumn(ms, ms.columnName(casacore::MSMainEnums::TIME));
	casacore::MEpoch time = timeColumn(idToMSRow[0]);
	casacore::MEpoch endTime = timeColumn(idToMSRow.back());
	double totalSeconds = endTime.get("s").getValue() - time.get("s").getValue();
	size_t intervalCount = (totalSeconds + _secondsBeforeBeamUpdate - 1) / _secondsBeforeBeamUpdate;
	std::vector<size_t> timestepIds(1, 0);
	for(size_t id=0;id!=idToMSRow.size();++id)
	{
		if(timeColumn(idToMSRow[id]).getValue() != time.getValue())
		{
			timestepIds.push_back(id);
			time = timeColumn(idToMSRow[id]);
		}
	}
	size_t timestepCount = timestepIds.size();
	timestepIds.push_back(idToMSRow.size());
	
	if(intervalCount > timestepCount)
		intervalCount = timestepCount;
	Logger::Debug << "MS spans " << totalSeconds << " seconds, dividing in " << intervalCount << " intervals.\n";
	
	casacore::MEpoch midTime(casacore::MVEpoch(0.5 * (timeColumn(idToMSRow[0]).getValue().get() + timeColumn(idToMSRow.back()).getValue().get())), timeColumn(idToMSRow[0]).getRef());
	Logger::Debug << "Mid time for full selection: " << midTime << '\n';
	casacore::MeasFrame midFrame(arrayPos, midTime);
	const casacore::MDirection::Ref hadecRef(casacore::MDirection::HADEC, midFrame);
	const casacore::MDirection::Ref azelgeoRef(casacore::MDirection::AZELGEO, midFrame);
	const casacore::MDirection::Ref midJ2000Ref(casacore::MDirection::J2000, midFrame);
	
	double refIntervalWeight = 0.0;
	msProvider.Reset();
	for(size_t intervalIndex=0; intervalIndex!=intervalCount; ++intervalIndex)
	{
		size_t timestepStart = intervalIndex*timestepCount/intervalCount;
		size_t timestepEnd = (intervalIndex+1)*timestepCount/intervalCount;
		size_t intervalStartId = timestepIds[timestepStart];
		size_t intervalEndId = timestepIds[timestepEnd];
	
		// Find the mid time step
		casacore::MEpoch firstTime = timeColumn(idToMSRow[intervalStartId]);
		casacore::MEpoch lastTime = timeColumn(idToMSRow[intervalEndId-1]);
		time = casacore::MEpoch(casacore::MVEpoch(0.5 * (firstTime.getValue().get() + lastTime.getValue().get())), firstTime.getRef());
		Logger::Debug << "Mid time for this interval: " << time << '\n';
		
		casacore::MeasFrame frame(arrayPos, time);
		const casacore::MDirection::Ref j2000Ref(casacore::MDirection::J2000, frame);
		
		if(_useDifferentialBeam)
			Logger::Debug << "Making differential snapshot beam for timesteps " << timestepStart << " - " << timestepEnd << "\n";
		else
			Logger::Debug << "Making snapshot beam for timesteps " << timestepStart << " - " << timestepEnd << "\n";
		
		double intervalWeight = 0.0;
		ao::uvector<double> stationWeights(stations.size(), 0.0);
		WeightMatrix baselineWeights(stations.size());
		calculateStationWeights(_imageWeightCache->Weights(), intervalWeight, stationWeights, baselineWeights, msProvider, selection, intervalStartId, intervalEndId);
		
		if(refIntervalWeight == 0.0)
			refIntervalWeight = intervalWeight;
		Logger::Debug << "Relative weight of this interval: " << ((intervalWeight==0.0) ? 0.0 : (intervalWeight * 100.0 / refIntervalWeight)) << " %\n";
		
		if(intervalWeight != 0.0)
		{
			ao::uvector<double> singleImages[8];
			double *imgPtr[8];
			for(size_t i=0; i!=8; ++i)
			{
				singleImages[i].assign(_sampledWidth*_sampledHeight, 0.0);
				imgPtr[i] = &singleImages[i][0];
			}
		
			makeBeamSnapshot(stations, stationWeights, baselineWeights, imgPtr, time.getValue().get()*86400.0, centralFrequency, centralFrequency, frame);
		
			_totalWeightSum += intervalWeight;
			for(size_t i=0; i!=8; ++i)
			{
				for(size_t j=0; j!=_sampledWidth*_sampledHeight; ++j)
				{
					double val = singleImages[i][j];
					beamImages[i][j] += val * intervalWeight;
				}
			}
		}
	}
}

static void dirToITRF(const casacore::MDirection& dir, casacore::MDirection::Convert& convert, vector3r_t& itrf)
{
	casacore::MDirection itrfDir = convert(dir);
	casacore::Vector<double> itrfVal = itrfDir.getValue().getValue();
	itrf[0] = itrfVal[0];
	itrf[1] = itrfVal[1];
	itrf[2] = itrfVal[2];
}

void LBeamImageMaker::makeBeamSnapshot(const std::vector<Station::Ptr>& stations, const ao::uvector<double>& weights, const WeightMatrix& baselineWeights, double** imgPtr, double time, double frequency, double subbandFrequency, const casacore::MeasFrame& frame)
{
	static const casacore::Unit radUnit("rad");
	const casacore::MDirection::Ref itrfRef(casacore::MDirection::ITRF, frame);
	const casacore::MDirection::Ref j2000Ref(casacore::MDirection::J2000, frame);
	casacore::MDirection::Convert
		j2000ToITRFRef(j2000Ref, itrfRef);
		
	vector3r_t station0, tile0, diffBeamCentre;
	dirToITRF(_delayDir, j2000ToITRFRef, station0);
	dirToITRF(_tileBeamDir, j2000ToITRFRef, tile0);
	
	std::vector<MC2x2> inverseCentralGain;
	if(_useDifferentialBeam)
	{
		dirToITRF(_referenceDir, j2000ToITRFRef, diffBeamCentre);
		inverseCentralGain.resize(stations.size());
		for(size_t a=0; a!=stations.size(); ++a)
		{
			Station::Ptr s=stations[a];
			matrix22c_t gainMatrix = s->response(time, frequency, diffBeamCentre, subbandFrequency, station0, tile0);
			inverseCentralGain[a][0] = gainMatrix[0][0];
			inverseCentralGain[a][1] = gainMatrix[0][1];
			inverseCentralGain[a][2] = gainMatrix[1][0];
			inverseCentralGain[a][3] = gainMatrix[1][1];
			if(!inverseCentralGain[a].Invert())
			{
				inverseCentralGain[a] = MC2x2::NaN();
			}
		}
	}

	ProgressBar progressBar("Constructing beam");
	for(size_t y=0;y!=_sampledHeight;++y)
	{
		for(size_t x=0;x!=_sampledWidth;++x)
		{
			double l, m, ra, dec;
			ImageCoordinates::XYToLM(x, y, _sPixelSizeX, _sPixelSizeY, _sampledWidth, _sampledHeight, l, m);
			l += _phaseCentreDL; m += _phaseCentreDM;
			ImageCoordinates::LMToRaDec(l, m, _phaseCentreRA, _phaseCentreDec, ra, dec);
			
			casacore::MDirection imageDir(casacore::MVDirection(
							casacore::Quantity(ra, radUnit),
							casacore::Quantity(dec,radUnit)),
							j2000Ref);
			
			vector3r_t itrfDirection;
			dirToITRF(imageDir, j2000ToITRFRef, itrfDirection);
			
			std::vector<MC2x2> stationGains(stations.size());
			for(size_t a=0; a!=stations.size(); ++a)
			{
				Station::Ptr s=stations[a];
				matrix22c_t gainMatrix = s->response(time, frequency, itrfDirection, subbandFrequency, station0, tile0);
				stationGains[a][0] = gainMatrix[0][0];
				stationGains[a][1] = gainMatrix[0][1];
				stationGains[a][2] = gainMatrix[1][0];
				stationGains[a][3] = gainMatrix[1][1];
				
				if(_useDifferentialBeam)
				{
					// The data have been premultiplied with the central beam C, and we want to
					// return a matrix that corrects the data for the full beam B. Given our data:
					//    data = Ci^-1 vis Cj^-*
					// we want to multiple data with a differential beam matrix D such that
					//    Di^-1 data Dj^-* = Bi^-1 vis Bj^-*
					// With B the full beam matrix. We can solve for D:
					//    Di^-1 Ci^-1 = Bi^-1  (and Cj^-* Dj^-* = Bj^-* )
					//          Di^-1 = Bi^-1 Ci
					//          Di    = Ci^-1 Bi
					// which means: diffBeam = inverseCentralGain * stationGain
					MC2x2 diffBeam;
					MC2x2::ATimesB(diffBeam, inverseCentralGain[a], stationGains[a]);
					stationGains[a] = diffBeam;
				}
			}
			
			MC2x2 gain = MC2x2::Zero();
			double totalWeight = 0.0;
			for(size_t a=0; a!=stations.size(); ++a)
			{
				double w = sqrt(weights[a]);
				gain.AddWithFactorAndAssign(stationGains[a], w);
				totalWeight += w;
			}
			gain *= (1.0/totalWeight);
		
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

void LBeamImageMaker::calculateStationWeights(const ImageWeights& imageWeights, double& totalWeight, ao::uvector<double>& weights, WeightMatrix& baselineWeights, MSProvider& msProvider, const MSSelection& selection, size_t intervalStartIdIndex, size_t intervalEndIdIndex)
{
	casacore::MeasurementSet& ms = msProvider.MS();
	casacore::MSAntenna antTable(ms.antenna());
	totalWeight = 0.0;
	weights.assign(antTable.nrow(), 0.0);
	
	MultiBandData multiband(ms.spectralWindow(), ms.dataDescription());
	size_t channelCount = selection.ChannelRangeEnd() - selection.ChannelRangeStart();
	ao::uvector<float> weightArr(channelCount);
	
	size_t id = intervalStartIdIndex;
	while(msProvider.CurrentRowAvailable() && id < intervalEndIdIndex)
	{
		size_t a1, a2, dataDescId;
		double uInM, vInM, wInM;
		msProvider.ReadMeta(uInM, vInM, wInM, dataDescId, a1, a2);
		const BandData &band(multiband[dataDescId]);
		msProvider.ReadWeights(weightArr.data());
		
		for(size_t ch=0; ch!=channelCount; ++ch)
		{
			double
				u = uInM / band.ChannelWavelength(ch),
				v = vInM / band.ChannelWavelength(ch);
			double iw = imageWeights.GetWeight(u, v);
			double w = weightArr[ch] * iw;
			totalWeight += w;
			weights[a1] += w;
			weights[a2] += w;
			baselineWeights.Value(a1, a2) += w;
		}
		msProvider.NextRow();
		id = msProvider.RowId();
	}
	
	if(Logger::IsVerbose())
		logWeights(ms, weights);
}

void LBeamImageMaker::logWeights(casacore::MeasurementSet& ms, const ao::uvector<double>& weights)
{
	casacore::MSAntenna antTable(ms.antenna());
	Logger::Debug << "Weights:";
	casacore::ROScalarColumn<casacore::String> antNameCol(antTable, casacore::MSAntenna::columnName(casacore::MSAntenna::NAME));
	double maxWeight = 0.0;
	for(size_t a=0; a!=antTable.nrow(); ++a)
		maxWeight=std::max(weights[a], maxWeight);
	
	for(size_t a=0; a!=antTable.nrow(); ++a)
	{
		std::string name = antNameCol(a);
		double w = (maxWeight==0.0) ? 0.0 : round(weights[a]*1000.0/maxWeight)*0.1;
		Logger::Debug << ' ' << name << '=' << w;
	}
	Logger::Debug << '\n';
}

#endif
