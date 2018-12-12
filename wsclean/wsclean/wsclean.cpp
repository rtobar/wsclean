#include "wsclean.h"

#include "directmsgridder.h"
#include "imageweightcache.h"
#include "measurementsetgridder.h"
#include "logger.h"
#include "wscfitswriter.h"
#include "wsmsgridder.h"
#include "primarybeam.h"
#include "imagefilename.h"

#include "../units/angle.h"

#include "../application.h"
#include "../areaset.h"
#include "../dftpredictionalgorithm.h"
#include "../fftresampler.h"
#include "../fitswriter.h"
#include "../gaussianfitter.h"
#include "../image.h"
#include "../imageweights.h"
#include "../modelrenderer.h"
#include "../msselection.h"
#include "../msproviders/contiguousms.h"
#include "../ndppp.h"
#include "../progressbar.h"
#include "../uvector.h"

#include "../deconvolution/deconvolutionalgorithm.h"
#include "../deconvolution/imageset.h"

#include "../idg/idgmsgridder.h"

#include "../lofar/lmspredicter.h"

#include "../model/areaparser.h"
#include "../model/model.h"

#include <iostream>
#include <memory>

std::string commandLine;

WSClean::WSClean() :
	_globalSelection(),
	_commandLine(),
	_inversionWatch(false), _predictingWatch(false), _deconvolutionWatch(false),
	_isFirstInversion(true), _doReorder(false),
	_majorIterationNr(0),
	_deconvolution(_settings)
{ }

WSClean::~WSClean()
{ }

void WSClean::multiplyImage(double factor, double* image) const
{
	size_t nPix = _settings.trimmedImageWidth * _settings.trimmedImageHeight;
	for(size_t i=0; i!=nPix; ++i)
		image[i] *= factor;
}

void WSClean::imagePSF(ImagingTableEntry& entry)
{
	size_t channelIndex = entry.outputChannelIndex;
	Logger::Info.Flush();
	Logger::Info << " == Constructing PSF ==\n";
	_inversionWatch.Start();
	_gridder->SetDoImagePSF(true);
	_gridder->SetDoSubtractModel(false);
	_gridder->SetVerbose(_isFirstInversion);
	_gridder->SetMetaDataCache(&_msGridderMetaCache[entry.index]);
	_gridder->SetStoreImagingWeights(_settings.writeImagingWeightSpectrumColumn);
	_gridder->Invert();
	
	size_t centralIndex = _settings.trimmedImageWidth/2 + (_settings.trimmedImageHeight/2) * _settings.trimmedImageWidth;
	
	double normFactor;
	if(_gridder->ImageRealResult()[centralIndex] != 0.0)
		normFactor = 1.0/_gridder->ImageRealResult()[centralIndex];
	else
		normFactor = 0.0;
	_infoPerChannel[channelIndex].psfNormalizationFactor = normFactor;
	multiplyImage(normFactor, _gridder->ImageRealResult());
	Logger::Debug << "Normalized PSF by factor of " << normFactor << ".\n";
		
	DeconvolutionAlgorithm::RemoveNaNsInPSF(_gridder->ImageRealResult(), _settings.trimmedImageWidth, _settings.trimmedImageHeight);
	_psfImages.SetFitsWriter(createWSCFitsWriter(entry, false, false).Writer());
	_psfImages.Store(_gridder->ImageRealResult(), *_settings.polarizations.begin(), channelIndex, false);
	_inversionWatch.Pause();
	
	_isFirstInversion = false;
	
	double bMaj, bMin, bPA;
	determineBeamSize(bMaj, bMin, bPA, _gridder->ImageRealResult(), _gridder->BeamSize());
	entry.imageWeight = _gridder->ImageWeight();
	entry.normalizationFactor = _gridder->NormalizationFactor();
	_infoPerChannel[channelIndex].theoreticBeamSize = _gridder->BeamSize();
	_infoPerChannel[channelIndex].beamMaj = bMaj;
	_infoPerChannel[channelIndex].beamMin = bMin;
	_infoPerChannel[channelIndex].beamPA = bPA;
	_infoPerChannel[channelIndex].weight = entry.imageWeight;
	_infoPerChannel[channelIndex].normalizationFactor = entry.normalizationFactor;
	_infoPerChannel[channelIndex].wGridSize = _gridder->WGridSize();
	_infoPerChannel[channelIndex].visibilityCount = _gridder->GriddedVisibilityCount();
	_infoPerChannel[channelIndex].effectiveVisibilityCount = _gridder->EffectiveGriddedVisibilityCount();
	_infoPerChannel[channelIndex].visibilityWeightSum = _gridder->VisibilityWeightSum();
		
	if(_settings.isUVImageSaved)
	{
		saveUVImage(_gridder->ImageRealResult(), *_settings.polarizations.begin(), entry, false, "uvpsf");
	}
	
	Logger::Info << "Writing psf image... ";
	Logger::Info.Flush();
	const std::string name(ImageFilename::GetPSFPrefix(_settings, channelIndex, entry.outputIntervalIndex) + "-psf.fits");
	WSCFitsWriter fitsFile = createWSCFitsWriter(entry, false, false);
	fitsFile.WritePSF(name, _gridder->ImageRealResult());
	Logger::Info << "DONE\n";
}

void WSClean::imageGridding()
{
	Logger::Info << "Writing gridding correction image... ";
	Logger::Info.Flush();
	double* gridding = _imageAllocator.Allocate(_gridder->ActualInversionWidth() * _gridder->ActualInversionHeight());
	_gridder->GetGriddingCorrectionImage(&gridding[0]);
	FitsWriter fitsWriter;
	fitsWriter.SetImageDimensions(_gridder->ActualInversionWidth(), _gridder->ActualInversionHeight());
	fitsWriter.Write(_settings.prefixName + "-gridding.fits", &gridding[0]);
	_imageAllocator.Free(gridding);
	Logger::Info << "DONE\n";
}

void WSClean::imageMainFirst(const ImagingTableEntry& entry)
{
	size_t joinedChannelIndex = entry.outputChannelIndex;
	
	Logger::Info.Flush();
	Logger::Info << " == Constructing image ==\n";
	_inversionWatch.Start();
	_gridder->SetDoImagePSF(false);
	_gridder->SetDoSubtractModel(_settings.subtractModel || _settings.continuedRun);
	_gridder->SetVerbose(_isFirstInversion);
	_gridder->SetMetaDataCache(&_msGridderMetaCache[entry.index]);
	_gridder->SetStoreImagingWeights(_settings.writeImagingWeightSpectrumColumn);
	_gridder->Invert();
	_inversionWatch.Pause();
	_gridder->SetVerbose(false);
	
	multiplyImage(_infoPerChannel[joinedChannelIndex].psfNormalizationFactor, _gridder->ImageRealResult());
	storeAndCombineXYandYX(_residualImages, entry.polarization, joinedChannelIndex, false, _gridder->ImageRealResult());
	if(Polarization::IsComplex(entry.polarization))
	{
		multiplyImage(_infoPerChannel[joinedChannelIndex].psfNormalizationFactor, _gridder->ImageImaginaryResult());
		storeAndCombineXYandYX(_residualImages, entry.polarization, joinedChannelIndex, true, _gridder->ImageImaginaryResult());
	}
}

void WSClean::imageMainNonFirst(const ImagingTableEntry& entry)
{
	size_t joinedChannelIndex = entry.outputChannelIndex;
	
	Logger::Info.Flush();
	Logger::Info << " == Constructing image ==\n";
	_inversionWatch.Start();
	_gridder->SetDoSubtractModel(true);
	_gridder->SetMetaDataCache(&_msGridderMetaCache[entry.index]);
	_gridder->Invert();
	_inversionWatch.Pause();
	
	multiplyImage(_infoPerChannel[joinedChannelIndex].psfNormalizationFactor, _gridder->ImageRealResult());
	storeAndCombineXYandYX(_residualImages, entry.polarization, joinedChannelIndex, false, _gridder->ImageRealResult());
	if(Polarization::IsComplex(entry.polarization))
	{
		multiplyImage(_infoPerChannel[joinedChannelIndex].psfNormalizationFactor, _gridder->ImageImaginaryResult());
		storeAndCombineXYandYX(_residualImages, entry.polarization, joinedChannelIndex, true, _gridder->ImageImaginaryResult());
	}
}

void WSClean::storeAndCombineXYandYX(CachedImageSet& dest, PolarizationEnum polarization, size_t joinedChannelIndex, bool isImaginary, const double* image)
{
	if(polarization == Polarization::YX && _settings.polarizations.count(Polarization::XY)!=0)
	{
		Logger::Info << "Adding XY and YX together...\n";
		double
			*xyImage = _imageAllocator.Allocate(_settings.trimmedImageWidth*_settings.trimmedImageHeight);
		dest.Load(xyImage, Polarization::XY, joinedChannelIndex, isImaginary);
		size_t count = _settings.trimmedImageWidth*_settings.trimmedImageHeight;
		if(isImaginary)
		{
			for(size_t i=0; i!=count; ++i)
				xyImage[i] = (xyImage[i]-image[i])*0.5;
		}
		else {
			for(size_t i=0; i!=count; ++i)
				xyImage[i] = (xyImage[i]+image[i])*0.5;
		}
		dest.Store(xyImage, Polarization::XY, joinedChannelIndex, isImaginary);
		_imageAllocator.Free(xyImage);
	}
	else {
		dest.Store(image, polarization, joinedChannelIndex, isImaginary);
	}
}

void WSClean::predict(const ImagingTableEntry& entry)
{
	Logger::Info.Flush();
	Logger::Info << " == Converting model image to visibilities ==\n";
	const size_t size = _settings.trimmedImageWidth*_settings.trimmedImageHeight;
	double
		*modelImageReal = _imageAllocator.Allocate(size),
		*modelImageImaginary = 0;
		
	if(entry.polarization == Polarization::YX)
	{
		_modelImages.Load(modelImageReal, Polarization::XY, entry.outputChannelIndex, false);
		modelImageImaginary = _imageAllocator.Allocate(size);
		_modelImages.Load(modelImageImaginary, Polarization::XY, entry.outputChannelIndex, true);
		for(size_t i=0; i!=size; ++i)
			modelImageImaginary[i] = -modelImageImaginary[i];
	}
	else {
		_modelImages.Load(modelImageReal, entry.polarization, entry.outputChannelIndex, false);
		if(Polarization::IsComplex(entry.polarization))
		{
			modelImageImaginary = _imageAllocator.Allocate(size);
			_modelImages.Load(modelImageImaginary, entry.polarization, entry.outputChannelIndex, true);
		}
	}
	
	_predictingWatch.Start();
	_gridder->SetAddToModel(false);
	_gridder->SetMetaDataCache(&_msGridderMetaCache[entry.index]);
	if(Polarization::IsComplex(entry.polarization))
		_gridder->Predict(modelImageReal, modelImageImaginary);
	else
		_gridder->Predict(modelImageReal);
	_predictingWatch.Pause();
	_imageAllocator.Free(modelImageReal);
	_imageAllocator.Free(modelImageImaginary);
}

void WSClean::initializeImageWeights(const ImagingTableEntry& entry)
{
	if(!_settings.mfsWeighting)
	{
		_imageWeightCache->Update(*_gridder, entry.outputChannelIndex, entry.outputIntervalIndex);
		if(_settings.isWeightImageSaved)
			_imageWeightCache->Weights().Save(_settings.prefixName+"-weights.fits");
	}
	_gridder->SetPrecalculatedWeightInfo(&_imageWeightCache->Weights());
}

void WSClean::initializeMFSImageWeights()
{
	Logger::Info << "Precalculating MFS weights for " << _settings.weightMode.ToString() << " weighting...\n";
	_imageWeightCache->ResetWeights();
	if(_doReorder)
	{
		for(size_t sg=0; sg!=_imagingTable.SquaredGroupCount(); ++sg)
		{
			const ImagingTable subTable = _imagingTable.GetSquaredGroup(sg);
			const ImagingTableEntry& entry = subTable.Front();
			for(size_t msIndex=0; msIndex!=_settings.filenames.size(); ++msIndex)
			{
				const ImagingTableEntry::MSInfo& ms = entry.msData[msIndex];
				for(size_t dataDescId=0; dataDescId!=_msBands[msIndex].DataDescCount(); ++dataDescId)
				{
					MSSelection partSelection(_globalSelection);
					partSelection.SetBandId(dataDescId);
					bool hasSelection = selectChannels(partSelection, msIndex, dataDescId, subTable.Front());
					if(hasSelection)
					{
						PolarizationEnum pol = _settings.useIDG ? Polarization::Instrumental : entry.polarization;
						PartitionedMS msProvider(_partitionedMSHandles[msIndex], ms.bands[dataDescId].partIndex, pol, dataDescId);
						_imageWeightCache->Weights().Grid(msProvider, partSelection);
					}
				}
			}
		}
	}
	else {
		for(size_t i=0; i!=_settings.filenames.size(); ++i)
		{
			for(size_t d=0; d!=_msBands[i].DataDescCount(); ++d)
			{
				PolarizationEnum pol = _settings.useIDG ? Polarization::Instrumental : *_settings.polarizations.begin();
				ContiguousMS msProvider(_settings.filenames[i], _settings.dataColumnName, _globalSelection, pol, d);
				_imageWeightCache->Weights().Grid(msProvider,  _globalSelection);
				Logger::Info << '.';
				Logger::Info.Flush();
			}
		}
	}
	_imageWeightCache->Weights().FinishGridding();
	_imageWeightCache->InitializeWeightTapers();
	if(_settings.isWeightImageSaved)
		_imageWeightCache->Weights().Save(_settings.prefixName+"-weights.fits");
}

void WSClean::prepareInversionAlgorithm(PolarizationEnum polarization)
{
	_gridder->SetGridMode(_settings.gridMode);
	_gridder->SetImageWidth(_settings.paddedImageWidth);
	_gridder->SetImageHeight(_settings.paddedImageHeight);
	_gridder->SetTrimSize(_settings.trimmedImageWidth, _settings.trimmedImageHeight);
	_gridder->SetNWSize(_settings.widthForNWCalculation, _settings.heightForNWCalculation);
	_gridder->SetPixelSizeX(_settings.pixelScaleX);
	_gridder->SetPixelSizeY(_settings.pixelScaleY);
	if(_settings.nWLayers != 0)
		_gridder->SetWGridSize(_settings.nWLayers);
	else
		_gridder->SetNoWGridSize();
	_gridder->SetAntialiasingKernelSize(_settings.antialiasingKernelSize);
	_gridder->SetOverSamplingFactor(_settings.overSamplingFactor);
	_gridder->SetPolarization(polarization);
	_gridder->SetIsComplex(polarization == Polarization::XY || polarization == Polarization::YX);
	_gridder->SetDataColumnName(_settings.dataColumnName);
	_gridder->SetWeighting(_settings.weightMode);
	_gridder->SetWLimit(_settings.wLimit/100.0);
	_gridder->SetSmallInversion(_settings.smallInversion);
	_gridder->SetVisibilityWeightingMode(_settings.visibilityWeightingMode);
}

std::unique_ptr<MSGridderBase> WSClean::createGridder() const
{
	if(_settings.useIDG)
		return std::unique_ptr<MSGridderBase>(new IdgMsGridder(_settings));
	else if(_settings.directFT)
	{
		switch(_settings.directFTPrecision) {
		case DirectFTPrecision::Half:
			throw std::runtime_error("Half precision is not implemented");
			//return std::unique_ptr<MSGridderBase>(new DirectMSGridder<half_float::half>(&_imageAllocator, _settings.threadCount));
			break;
		case DirectFTPrecision::Float:
			return std::unique_ptr<MSGridderBase>(new DirectMSGridder<float>(&_imageAllocator, _settings.threadCount));
			break;
		default:
		case DirectFTPrecision::Double:
			return std::unique_ptr<MSGridderBase>(new DirectMSGridder<double>(&_imageAllocator, _settings.threadCount));
			break;
		case DirectFTPrecision::LongDouble:
			return std::unique_ptr<MSGridderBase>(new DirectMSGridder<long double>(&_imageAllocator, _settings.threadCount));
			break;
		}
	}
	else
		return std::unique_ptr<MSGridderBase>(new WSMSGridder(&_imageAllocator, _settings.threadCount, _settings.memFraction, _settings.absMemLimit));
}

void WSClean::performReordering(bool isPredictMode)
{
	_partitionedMSHandles.clear();
	for(size_t i=0; i != _settings.filenames.size(); ++i)
	{
		std::vector<PartitionedMS::ChannelRange> channels;
		std::map<PolarizationEnum, size_t> nextIndex;
		for(size_t j=0; j!=_imagingTable.SquaredGroupCount(); ++j)
		{
			ImagingTable squaredGroup = _imagingTable.GetSquaredGroup(j);
			for(size_t s=0; s!=squaredGroup.EntryCount(); ++s)
			{
				ImagingTableEntry& entry =
					_imagingTable[squaredGroup[s].index];
				for(size_t d=0; d!=_msBands[i].DataDescCount(); ++d)
				{
					MSSelection selection(_globalSelection);
					if(selectChannels(selection, i, d, entry))
					{
						if(entry.polarization == *_settings.polarizations.begin())
						{
							PartitionedMS::ChannelRange r;
							r.dataDescId = d;
							r.start = selection.ChannelRangeStart();
							r.end = selection.ChannelRangeEnd();
							channels.push_back(r);
						}
						entry.msData[i].bands[d].partIndex = nextIndex[entry.polarization];
						++nextIndex[entry.polarization];
					}
				}
			}
		}
		
		bool useModel = _settings.deconvolutionMGain != 1.0 || isPredictMode || _settings.subtractModel || _settings.continuedRun;
		bool initialModelRequired = _settings.subtractModel || _settings.continuedRun;
		_partitionedMSHandles.push_back(PartitionedMS::Partition(_settings.filenames[i], channels, _globalSelection, _settings.dataColumnName, useModel, initialModelRequired, _settings));
	}
}

void WSClean::RunClean()
{
	// If no column specified, determine column to use
	if(_settings.dataColumnName.empty())
	{
		casacore::MeasurementSet ms(_settings.filenames.front());
		bool hasCorrected = ms.tableDesc().isColumn("CORRECTED_DATA");
		if(hasCorrected) {
			Logger::Info << "First measurement set has corrected data: tasks will be applied on the corrected data column.\n";
			_settings.dataColumnName = "CORRECTED_DATA";
		} else {
			Logger::Info << "No corrected data in first measurement set: tasks will be applied on the data column.\n";
			_settings.dataColumnName= "DATA";
		}
	}

	_settings.Propogate();
	
	_globalSelection = _settings.GetMSSelection();
	MSSelection fullSelection = _globalSelection;
	
	for(size_t intervalIndex=0; intervalIndex!=_settings.intervalsOut; ++intervalIndex)
	{
		makeImagingTable(intervalIndex);
		
		_globalSelection = selectInterval(fullSelection, intervalIndex);
		
		_doReorder = preferReordering();
		
		if(_doReorder) performReordering(false);
		
		_infoPerChannel.assign(_settings.channelsOut, OutputChannelInfo());
		
		_imageWeightCache = createWeightCache();
		
		if(_settings.mfsWeighting)
			initializeMFSImageWeights();
		
		_gridder = createGridder();
		
		for(size_t groupIndex=0; groupIndex!=_imagingTable.IndependentGroupCount(); ++groupIndex)
		{
			ImagingTable group = _imagingTable.GetIndependentGroup(groupIndex);
			runIndependentGroup(group);
		}

		// Needs to be destructed before image allocator, or image allocator will report error caused by leaked memory
		_gridder.reset();
	
		if(_settings.channelsOut > 1)
		{
			for(std::set<PolarizationEnum>::const_iterator pol=_settings.polarizations.begin(); pol!=_settings.polarizations.end(); ++pol)
			{
				bool psfWasMade = (_settings.deconvolutionIterationCount > 0 || _settings.makePSF || _settings.makePSFOnly) && pol == _settings.polarizations.begin();
				
				if(psfWasMade)
				{
					makeMFSImage("psf.fits", intervalIndex, *pol, false, true);
					if(_settings.savePsfPb)
						makeMFSImage("psf-pb.fits", intervalIndex, *pol, false, true);
				}
				
				if(!(*pol == Polarization::YX && _settings.polarizations.count(Polarization::XY)!=0) && !_settings.makePSFOnly)
				{
					if(_settings.isDirtySaved)
						makeMFSImage("dirty.fits", intervalIndex, *pol, false);
					if(_settings.deconvolutionIterationCount == 0)
					{
						makeMFSImage("image.fits", intervalIndex, *pol, false);
						if(_settings.applyPrimaryBeam || (_settings.gridWithBeam || !_settings.atermConfigFilename.empty()))
							makeMFSImage("image-pb.fits", intervalIndex, *pol, false);
					}
					else 
					{
						makeMFSImage("residual.fits", intervalIndex, *pol, false);
						makeMFSImage("model.fits", intervalIndex, *pol, false);
						renderMFSImage(intervalIndex, *pol, false, false);
						if(_settings.applyPrimaryBeam || (_settings.gridWithBeam || !_settings.atermConfigFilename.empty()))
						{
							makeMFSImage("residual-pb.fits", intervalIndex, *pol, false);
							makeMFSImage("model-pb.fits", intervalIndex, *pol, false);
							renderMFSImage(intervalIndex, *pol, false, true);
						}
					}
					if(Polarization::IsComplex(*pol))
					{
						if(_settings.isDirtySaved)
							makeMFSImage("dirty.fits", intervalIndex, *pol, true);
						if(_settings.deconvolutionIterationCount == 0)
						{
								makeMFSImage("image.fits", intervalIndex, *pol, true);
							if(_settings.applyPrimaryBeam || (_settings.gridWithBeam || !_settings.atermConfigFilename.empty()))
								makeMFSImage("image-pb.fits", intervalIndex, *pol, true);
						}
						else {
							makeMFSImage("residual.fits", intervalIndex, *pol, true);
							makeMFSImage("model.fits", intervalIndex, *pol, true);
							renderMFSImage(intervalIndex, *pol, true, false);
							if(_settings.applyPrimaryBeam || (_settings.gridWithBeam || !_settings.atermConfigFilename.empty()))
							{
								makeMFSImage("residual-pb.fits", intervalIndex, *pol, true);
								makeMFSImage("model-pb.fits", intervalIndex, *pol, true);
								renderMFSImage(intervalIndex, *pol, true, true);
							}
						}
					}
				}
			}
		}
	}
}

std::unique_ptr<ImageWeightCache> WSClean::createWeightCache()
{
	std::unique_ptr<ImageWeightCache> cache(new ImageWeightCache(
		_settings.weightMode,
		_settings.paddedImageWidth, _settings.paddedImageHeight,
		_settings.pixelScaleX, _settings.pixelScaleY,
		_settings.minUVInLambda, _settings.maxUVInLambda,
		_settings.rankFilterLevel, _settings.rankFilterSize,
		_settings.useWeightsAsTaper));
	cache->SetTaperInfo(
		_settings.gaussianTaperBeamSize,
		_settings.tukeyTaperInLambda, _settings.tukeyInnerTaperInLambda,
		_settings.edgeTaperInLambda, _settings.edgeTukeyTaperInLambda);
	return std::move(cache);
}

void WSClean::RunPredict()
{
	if(_settings.joinedFrequencyCleaning)
		throw std::runtime_error("Joined frequency cleaning specified for prediction: prediction doesn't clean, parameter invalid");
	if(_settings.joinedPolarizationCleaning)
		throw std::runtime_error("Joined polarization cleaning specified for prediction: prediction doesn't clean, parameter invalid");
	
	_settings.dataColumnName = "DATA";
	
	_settings.Propogate();
	
	_globalSelection = _settings.GetMSSelection();
	MSSelection fullSelection = _globalSelection;
	
	for(size_t intervalIndex=0; intervalIndex!=_settings.intervalsOut; ++intervalIndex)
	{
		makeImagingTable(intervalIndex);
		
		if(_settings.predictionChannels != 0)
		{
			// TODO
		}
		
		_infoPerChannel.assign(_settings.channelsOut, OutputChannelInfo());
		
		_globalSelection = selectInterval(fullSelection, intervalIndex);
		
		_doReorder = preferReordering();
		
		if(_doReorder) performReordering(true);
		
		_gridder = createGridder();
	
		for(size_t groupIndex=0; groupIndex!=_imagingTable.SquaredGroupCount(); ++groupIndex)
		{
			predictGroup(_imagingTable.GetSquaredGroup(groupIndex));
		}
	
		// Needs to be destructed before image allocator, or image allocator will report error caused by leaked memory
		_gridder.reset();
	}
}

bool WSClean::selectChannels(MSSelection& selection, size_t msIndex, size_t dataDescId, const ImagingTableEntry& entry)
{
	const BandData& band = _msBands[msIndex][dataDescId];
	double firstCh = band.ChannelFrequency(0), lastCh = band.ChannelFrequency(band.ChannelCount()-1);
	// Some mses have decreasing (i.e. reversed) channel frequencies in them
	bool isReversed = false;
	if(firstCh > lastCh) {
		std::swap(firstCh, lastCh);
		isReversed = true;
		Logger::Debug << "Warning: MS has reversed channel frequencies.\n";
	}
	if(band.ChannelCount()!=0 && entry.lowestFrequency <= lastCh && entry.highestFrequency >= firstCh)
	{
		size_t newStart, newEnd;
		if(isReversed)
		{
			BandData::const_reverse_iterator lowPtr, highPtr;
			lowPtr = std::lower_bound(band.rbegin(), band.rend(), entry.lowestFrequency);
			highPtr = std::lower_bound(lowPtr, band.rend(), entry.highestFrequency);
			
			if(highPtr == band.rend())
				--highPtr;
			newStart = band.ChannelCount() - 1 - (highPtr - band.rbegin());
			newEnd = band.ChannelCount() - (lowPtr - band.rbegin());
		}
		else {
			const double *lowPtr, *highPtr;
			lowPtr = std::lower_bound(band.begin(), band.end(), entry.lowestFrequency);
			highPtr = std::lower_bound(lowPtr, band.end(), entry.highestFrequency);
			
			if(highPtr == band.end())
				--highPtr;
			newStart = lowPtr - band.begin();
			newEnd = highPtr - band.begin() + 1;
		}
			
		selection.SetChannelRange(newStart, newEnd);
		return true;
	}
	else {
		return false;
	}
}

double WSClean::minTheoreticalBeamSize(const ImagingTable& table) const
{
	double beam = 0.0;
	for(size_t i=0; i!=table.EntryCount(); ++i)
	{
		const ImagingTableEntry& e = table[i];
		const OutputChannelInfo& info = _infoPerChannel[e.outputChannelIndex];
		if(std::isfinite(info.theoreticBeamSize) && (info.theoreticBeamSize < beam || beam == 0.0))
			beam = info.theoreticBeamSize;
	}
	return beam;
}

void WSClean::runIndependentGroup(ImagingTable& groupTable)
{
	WSCFitsWriter modelWriter(createWSCFitsWriter(groupTable.Front(), false, true));
	_modelImages.Initialize(modelWriter.Writer(), _settings.polarizations.size(), _settings.channelsOut, _settings.prefixName + "-model", _imageAllocator);
	WSCFitsWriter writer(createWSCFitsWriter(groupTable.Front(), false, false));
	_residualImages.Initialize(writer.Writer(), _settings.polarizations.size(), _settings.channelsOut, _settings.prefixName + "-residual", _imageAllocator);
	if(groupTable.Front().polarization == *_settings.polarizations.begin())
		_psfImages.Initialize(writer.Writer(), 1, groupTable.SquaredGroupCount(), _settings.prefixName + "-psf", _imageAllocator);
	
	const std::string rootPrefix = _settings.prefixName;
		
	for(size_t joinedIndex=0; joinedIndex!=groupTable.EntryCount(); ++joinedIndex)
	{
		ImagingTableEntry& entry = groupTable[joinedIndex];
		runFirstInversion(entry);
	}
	
	_deconvolution.InitializeDeconvolutionAlgorithm(groupTable, *_settings.polarizations.begin(), &_imageAllocator, minTheoreticalBeamSize(groupTable), _settings.threadCount);

	if(!_settings.makePSFOnly)
	{
		if(_settings.deconvolutionIterationCount > 0)
		{
			// Start major cleaning loop
			_majorIterationNr = 1;
			bool reachedMajorThreshold = false;
			do {
				_deconvolution.InitializeImages(_residualImages, _modelImages, _psfImages);
				_deconvolutionWatch.Start();
				_deconvolution.Perform(groupTable, reachedMajorThreshold, _majorIterationNr);
				_deconvolutionWatch.Pause();
				
				if(_majorIterationNr == 1 && _settings.deconvolutionMGain != 1.0 && _settings.isFirstResidualSaved)
					writeFirstResidualImages(groupTable);
		
				if(!reachedMajorThreshold)
					writeModelImages(groupTable);
		
				if(_settings.deconvolutionMGain != 1.0)
				{
					for(size_t sGroupIndex=0; sGroupIndex!=groupTable.SquaredGroupCount(); ++sGroupIndex)
					{
						const ImagingTable sGroupTable = groupTable.GetSquaredGroup(sGroupIndex);
						for(size_t e=0; e!=sGroupTable.EntryCount(); ++e)
						{
							prepareInversionAlgorithm(sGroupTable[e].polarization);
							initializeCurMSProviders(sGroupTable[e]);
							initializeImageWeights(sGroupTable[e]);
		
							predict(sGroupTable[e]);
							_currentPolMSes.clear();
						}
						for(size_t e=0; e!=sGroupTable.EntryCount(); ++e)
						{
							prepareInversionAlgorithm(sGroupTable[e].polarization);
							initializeCurMSProviders(sGroupTable[e]);
							initializeImageWeights(sGroupTable[e]);
							
							imageMainNonFirst(sGroupTable[e]);
							_currentPolMSes.clear();
						} // end of polarization loop
					} // end of joined channels loop
				}
				
				++_majorIterationNr;
			} while(reachedMajorThreshold);
			
			--_majorIterationNr;
			Logger::Info << _majorIterationNr << " major iterations were performed.\n";
		}
		
		_gridder->FreeImagingData();
		
		for(size_t joinedIndex=0; joinedIndex!=groupTable.EntryCount(); ++joinedIndex)
			saveRestoredImagesForGroup(groupTable[joinedIndex]);
		
		if(_settings.saveSourceList)
    {
			_deconvolution.SaveSourceList(groupTable, _gridder->PhaseCentreRA(), _gridder->PhaseCentreDec());
			if(_settings.applyPrimaryBeam)
			{
				_deconvolution.SavePBSourceList(groupTable, _gridder->PhaseCentreRA(), _gridder->PhaseCentreDec());
			}
    }
	}
	
	_deconvolution.FreeDeconvolutionAlgorithms();
	
	_imageAllocator.ReportStatistics();
	Logger::Info << "Inversion: " << _inversionWatch.ToString() << ", prediction: " << _predictingWatch.ToString() << ", deconvolution: " << _deconvolutionWatch.ToString() << '\n';
	
	_settings.prefixName = rootPrefix;
}

void WSClean::saveRestoredImagesForGroup(const ImagingTableEntry& tableEntry) const
{
	// Restore model to residual and save image
	size_t currentChannelIndex =
		tableEntry.outputChannelIndex;
	
	PolarizationEnum curPol = tableEntry.polarization;
	for(size_t imageIter=0; imageIter!=tableEntry.imageCount; ++imageIter)
	{
		bool isImaginary = (imageIter == 1);
		WSCFitsWriter writer(createWSCFitsWriter(tableEntry, isImaginary, false));
		double* restoredImage = _imageAllocator.Allocate(_settings.trimmedImageWidth*_settings.trimmedImageHeight);
		_residualImages.Load(restoredImage, curPol, currentChannelIndex, isImaginary);
		
		if(_settings.deconvolutionIterationCount != 0)
			writer.WriteImage("residual.fits", restoredImage);
		
		if(_settings.isUVImageSaved)
			saveUVImage(restoredImage, curPol, tableEntry, isImaginary, "uv");
		
		double* modelImage = _imageAllocator.Allocate(_settings.trimmedImageWidth*_settings.trimmedImageHeight);
		_modelImages.Load(modelImage, curPol, currentChannelIndex, isImaginary);
		double beamMaj = _infoPerChannel[currentChannelIndex].beamMaj;
		double beamMin, beamPA;
		std::string beamStr;
		if(std::isfinite(beamMaj))
		{
			beamMin = _infoPerChannel[currentChannelIndex].beamMin;
			beamPA = _infoPerChannel[currentChannelIndex].beamPA;
			beamStr = "(beam=" + Angle::ToNiceString(beamMin) + "-" +
			Angle::ToNiceString(beamMaj) + ", PA=" +
			Angle::ToNiceString(beamPA) + ")";
		}
		else {
			beamStr = "(beam is neither fitted nor estimated -- using delta scales!)";
			beamMaj = 0.0; beamMin = 0.0; beamPA = 0.0;
		}
		Logger::Info << "Rendering sources to restored image " + beamStr + "... ";
		Logger::Info.Flush();
		ModelRenderer::Restore(restoredImage, modelImage, _settings.trimmedImageWidth, _settings.trimmedImageHeight, beamMaj, beamMin, beamPA, _settings.pixelScaleX, _settings.pixelScaleY);
		Logger::Info << "DONE\n";
		_imageAllocator.Free(modelImage);
		
		Logger::Info << "Writing restored image... ";
		Logger::Info.Flush();
		writer.WriteImage("image.fits", restoredImage);
		Logger::Info << "DONE\n";
		_imageAllocator.Free(restoredImage);
		
		if(curPol == *_settings.polarizations.rbegin())
		{
			ImageFilename imageName = ImageFilename(currentChannelIndex, tableEntry.outputIntervalIndex);
			if(_settings.applyPrimaryBeam)
			{
				_primaryBeam->CorrectImages(writer.Writer(), imageName, "image", _imageAllocator);
				if(_settings.savePsfPb)
					_primaryBeam->CorrectImages(writer.Writer(), imageName, "psf", _imageAllocator);
				if(_settings.deconvolutionIterationCount != 0)
				{
					_primaryBeam->CorrectImages(writer.Writer(), imageName, "residual", _imageAllocator);
					_primaryBeam->CorrectImages(writer.Writer(), imageName, "model", _imageAllocator);
				}
			}
			else if(_settings.gridWithBeam || !_settings.atermConfigFilename.empty())
			{
				IdgMsGridder& idg = static_cast<IdgMsGridder&>(*_gridder);
				idg.SavePBCorrectedImages(writer.Writer(), imageName, "image", _imageAllocator);
				if(_settings.savePsfPb)
					idg.SavePBCorrectedImages(writer.Writer(), imageName, "psf", _imageAllocator);
				if(_settings.deconvolutionIterationCount != 0)
				{
					idg.SavePBCorrectedImages(writer.Writer(), imageName, "residual", _imageAllocator);
					idg.SavePBCorrectedImages(writer.Writer(), imageName, "model", _imageAllocator);
				}
			}
		}
	}
}

void WSClean::writeFirstResidualImages(const ImagingTable& groupTable) const
{
	Logger::Info << "Writing first iteration image(s)...\n";
	ImageBufferAllocator::Ptr ptr;
	_imageAllocator.Allocate(_settings.trimmedImageWidth*_settings.trimmedImageHeight, ptr);
	for(size_t e=0; e!=groupTable.EntryCount(); ++e)
	{
		const ImagingTableEntry& entry = groupTable[e];
		size_t ch = entry.outputChannelIndex;
		if(entry.polarization == Polarization::YX) {
			_residualImages.Load(ptr.data(), Polarization::XY, ch, true);
			WSCFitsWriter writer(createWSCFitsWriter(entry, Polarization::XY, true, false));
			writer.WriteImage("first-residual.fits", ptr.data());
		}
		else {
			_residualImages.Load(ptr.data(), entry.polarization, ch, false);
			WSCFitsWriter writer(createWSCFitsWriter(entry, false, false));
			writer.WriteImage("first-residual.fits", ptr.data());
		}
	}
}

void WSClean::writeModelImages(const ImagingTable& groupTable) const
{
	Logger::Info << "Writing model image...\n";
	ImageBufferAllocator::Ptr ptr;
	_imageAllocator.Allocate(_settings.trimmedImageWidth*_settings.trimmedImageHeight, ptr);
	for(size_t e=0; e!=groupTable.EntryCount(); ++e)
	{
		const ImagingTableEntry& entry = groupTable[e];
		size_t ch = entry.outputChannelIndex;
		if(entry.polarization == Polarization::YX) {
			_modelImages.Load(ptr.data(), Polarization::XY, ch, true);
			WSCFitsWriter writer(createWSCFitsWriter(entry, Polarization::XY, true, true));
			writer.WriteImage("model.fits", ptr.data());
		}
		else {
			_modelImages.Load(ptr.data(), entry.polarization, ch, false);
			WSCFitsWriter writer(createWSCFitsWriter(entry, false, true));
			writer.WriteImage("model.fits", ptr.data());
		}
	}
}

void WSClean::readEarlierModelImages(const ImagingTableEntry& entry)
{
	// load image(s) from disk and store them in the model-image cache.
	for(size_t i=0; i!=entry.imageCount; ++i)
	{
		std::string prefix = ImageFilename::GetPrefix(_settings, entry.polarization, entry.outputChannelIndex, entry.outputIntervalIndex, i==1);
		FitsReader reader(prefix + "-model.fits");
		Logger::Info << "Reading " << reader.Filename() << "...\n";
		if(_settings.trimmedImageWidth == 0 && _settings.trimmedImageHeight == 0)
		{
			_settings.trimmedImageWidth = reader.ImageWidth();
			_settings.trimmedImageHeight = reader.ImageHeight();
			_settings.RecalculatePaddedDimensions();
		}
		else if(reader.ImageWidth()!=_settings.trimmedImageWidth || reader.ImageHeight()!=_settings.trimmedImageHeight)
		{
			std::ostringstream msg;
			msg << "Inconsistent image size: dimensions of input image did not match, input: " << reader.ImageWidth() << " x " << reader.ImageHeight() << ", specified: " << _settings.trimmedImageWidth << " x " << _settings.trimmedImageHeight;
			throw std::runtime_error(msg.str());
		}
		
		if(reader.PixelSizeX()==0.0 || reader.PixelSizeY()==0.0)
			Logger::Warn << "Warning: input fits file misses the pixel size keywords.\n";
		else if(_settings.pixelScaleX == 0 && _settings.pixelScaleY == 0)
		{
			_settings.pixelScaleX = reader.PixelSizeX();
			_settings.pixelScaleY = reader.PixelSizeY();
			Logger::Debug << "Using pixel size of " << Angle::ToNiceString(_settings.pixelScaleX) << " x " << Angle::ToNiceString(_settings.pixelScaleY) << ".\n";
		}
		// Check if image corresponds with image dimensions of the settings
		// Here I require the pixel scale to be accurate enough so that the image is at most 1/10th pixel larger/smaller.
		else if(std::fabs(reader.PixelSizeX() - _settings.pixelScaleX) * _settings.trimmedImageWidth > 0.1 * _settings.pixelScaleX ||
			std::fabs(reader.PixelSizeY() - _settings.pixelScaleY) * _settings.trimmedImageHeight > 0.1 * _settings.pixelScaleY)
		{
			std::ostringstream msg;
			msg << "Inconsistent pixel size: pixel size of input image did not match. Input: " << reader.PixelSizeX() << " x " << reader.PixelSizeY() << ", specified: " << _settings.pixelScaleX << " x " << _settings.pixelScaleY;
			throw std::runtime_error(msg.str());
		}
		if(_settings.pixelScaleX == 0.0 || _settings.pixelScaleY == 0.0)
		{
			throw std::runtime_error("Could not determine proper pixel size. The input image did not provide proper pixel size values, and no or an invalid -scale was provided to WSClean");
		}
		
		// TODO check phase centre
		
		if(!_imageWeightCache)
		{
			// The construction of the weight cache is delayed in prediction mode, because only now
			// the image size and scale is known.
			_imageWeightCache = createWeightCache();
			if(_settings.mfsWeighting)
				initializeMFSImageWeights();
		}
		
		FitsWriter writer(reader);
		_modelImages.SetFitsWriter(writer);
		
		double* buffer = _imageAllocator.Allocate(_settings.trimmedImageWidth*_settings.trimmedImageHeight);
		reader.Read(buffer);
		for(size_t j=0; j!=_settings.trimmedImageWidth*_settings.trimmedImageHeight; ++j)
		{
			if(!std::isfinite(buffer[j]))
				throw std::runtime_error("The input image contains non-finite values -- can't predict from an image with non-finite values");
		}
		_modelImages.Store(buffer, entry.polarization, entry.outputChannelIndex, i==1);
		_imageAllocator.Free(buffer);
	}
}

void WSClean::predictGroup(const ImagingTable& imagingGroup)
{
	_modelImages.Initialize(
		createWSCFitsWriter(imagingGroup.Front(), false, true).Writer(),
		_settings.polarizations.size(), 1, _settings.prefixName + "-model", _imageAllocator
	);
	
	const std::string rootPrefix = _settings.prefixName;
		
	for(size_t e=0; e!=imagingGroup.EntryCount(); ++e)
	{
		const ImagingTableEntry& entry = imagingGroup[e];
		
		readEarlierModelImages(entry);
		
		prepareInversionAlgorithm(entry.polarization);
		initializeCurMSProviders(entry);
		initializeImageWeights(entry);

		predict(entry);
		
		_currentPolMSes.clear();
	} // end of polarization loop
	
	_imageAllocator.ReportStatistics();
	Logger::Info << "Inversion: " << _inversionWatch.ToString() << ", prediction: " << _predictingWatch.ToString() << ", cleaning: " << _deconvolutionWatch.ToString() << '\n';
	
	_settings.prefixName = rootPrefix;
}

std::unique_ptr<MSProvider> WSClean::initializeMSProvider(const ImagingTableEntry& entry, const MSSelection& selection, size_t filenameIndex, size_t dataDescId)
{
	PolarizationEnum pol = _settings.useIDG ? Polarization::Instrumental : entry.polarization;
	if(_doReorder)
		return std::unique_ptr<MSProvider>(new PartitionedMS(_partitionedMSHandles[filenameIndex], entry.msData[filenameIndex].bands[dataDescId].partIndex, pol, dataDescId));
	else
		return std::unique_ptr<MSProvider>(new ContiguousMS(_settings.filenames[filenameIndex], _settings.dataColumnName, selection, pol, dataDescId));
}

void WSClean::initializeCurMSProviders(const ImagingTableEntry& entry)
{
	_gridder->ClearMeasurementSetList();
	for(size_t i=0; i != _settings.filenames.size(); ++i)
	{
		for(size_t d=0; d!=_msBands[i].DataDescCount(); ++d)
		{
			MSSelection selection(_globalSelection);
			if(selectChannels(selection, i, d, entry))
			{
				std::unique_ptr<MSProvider> msProvider(
					initializeMSProvider(entry, selection, i, d));
				_gridder->AddMeasurementSet(msProvider.get(), selection);
				_currentPolMSes.emplace_back(std::move(msProvider));
			}
		}
	}
}

void WSClean::initializeMSProvidersForPB(const ImagingTableEntry& entry, PrimaryBeam& pb)
{
	for(size_t i=0; i != _settings.filenames.size(); ++i)
	{
		for(size_t d=0; d!=_msBands[i].DataDescCount(); ++d)
		{
			MSSelection selection(_globalSelection);
			if(selectChannels(selection, i, d, entry))
			{
				std::unique_ptr<MSProvider> msProvider = initializeMSProvider(entry, selection, i, d);
				pb.AddMS(msProvider.get(), selection);
				_currentPolMSes.emplace_back(std::move(msProvider));
			}
		}
	}
}

void WSClean::runFirstInversion(ImagingTableEntry& entry)
{
	initializeCurMSProviders(entry);
	initializeImageWeights(entry);
	
	prepareInversionAlgorithm(entry.polarization);
	
	const bool firstBeforePSF = _isFirstInversion;

	bool isFirstPol = entry.polarization == *_settings.polarizations.begin();
	bool isLastPol = entry.polarization == *_settings.polarizations.rbegin();
	bool doMakePSF = _settings.deconvolutionIterationCount > 0 || _settings.makePSF || _settings.makePSFOnly;
	if(doMakePSF && isFirstPol)
		imagePSF(entry);
	
	if(isLastPol)
	{
		ImageFilename imageName = ImageFilename(entry.outputChannelIndex, entry.outputIntervalIndex);
		if(_settings.applyPrimaryBeam)
		{
			_primaryBeam.reset(new PrimaryBeam(_settings));
			initializeMSProvidersForPB(entry, *_primaryBeam);
			// we don't have to call initializeImageWeights(entry), because they're still set ok.
			double ra, dec, dl, dm;
			MSGridderBase::GetPhaseCentreInfo(_currentPolMSes.front()->MS(), _settings.fieldId, ra, dec, dl, dm);
			_primaryBeam->SetPhaseCentre(ra, dec, dl, dm);
			_primaryBeam->MakeBeamImages(imageName, entry, _imageWeightCache.get(), _imageAllocator);
			_currentPolMSes.clear();
			initializeCurMSProviders(entry);
		}
	}
		
	if(!_settings.makePSFOnly)
	{
		FitsWriter mWriter(createWSCFitsWriter(entry, false, true).Writer());
		_modelImages.SetFitsWriter(mWriter);
		FitsWriter writer(createWSCFitsWriter(entry, false, false).Writer());
		_residualImages.SetFitsWriter(writer);
		
		imageMainFirst(entry);
		
		entry.imageWeight = _gridder->ImageWeight();
		entry.normalizationFactor = _gridder->NormalizationFactor();
		
		// If this was the first polarization of this channel, we need to set
		// the info for this channel
		if(isFirstPol)
		{
			_infoPerChannel[entry.outputChannelIndex].weight = entry.imageWeight;
			_infoPerChannel[entry.outputChannelIndex].normalizationFactor = entry.normalizationFactor;
			// If no PSF is made, also set the beam size. If the PSF was made, these would already be set
			// after imaging the PSF.
			if(!doMakePSF)
			{
				if(_settings.theoreticBeam) {
					_infoPerChannel[entry.outputChannelIndex].beamMaj = _gridder->BeamSize();
					_infoPerChannel[entry.outputChannelIndex].beamMin = _gridder->BeamSize();
					_infoPerChannel[entry.outputChannelIndex].beamPA = 0.0;
				}
				else if(_settings.manualBeamMajorSize != 0.0) {
					_infoPerChannel[entry.outputChannelIndex].beamMaj = _settings.manualBeamMajorSize;
					_infoPerChannel[entry.outputChannelIndex].beamMin = _settings.manualBeamMinorSize;
					_infoPerChannel[entry.outputChannelIndex].beamPA = _settings.manualBeamPA;
				}
				else {
					_infoPerChannel[entry.outputChannelIndex].beamMaj = std::numeric_limits<double>::quiet_NaN();
					_infoPerChannel[entry.outputChannelIndex].beamMin = std::numeric_limits<double>::quiet_NaN();
					_infoPerChannel[entry.outputChannelIndex].beamPA = std::numeric_limits<double>::quiet_NaN();
				}
			}
		}
		
		if(_settings.isGriddingImageSaved && firstBeforePSF && _gridder->HasGriddingCorrectionImage())
			imageGridding();
		
		_isFirstInversion = false;
		
		if(_settings.continuedRun)
		{
			readEarlierModelImages(entry);
		}
		else {
			// Set model to zero: already done if this is YX of XY/YX imaging combi
			if(!(entry.polarization == Polarization::YX && _settings.polarizations.count(Polarization::XY)!=0))
			{
				double* modelImage = _imageAllocator.Allocate(_settings.trimmedImageWidth * _settings.trimmedImageHeight);
				memset(modelImage, 0, _settings.trimmedImageWidth * _settings.trimmedImageHeight * sizeof(double));
				_modelImages.Store(modelImage, entry.polarization, entry.outputChannelIndex, false);
				if(Polarization::IsComplex(entry.polarization))
					_modelImages.Store(modelImage, entry.polarization, entry.outputChannelIndex, true);
				_imageAllocator.Free(modelImage);
			}
		}
		
		if(_settings.isDirtySaved)
		{
			for(size_t imageIndex=0; imageIndex!=entry.imageCount; ++imageIndex)
			{
				bool isImaginary = (imageIndex==1);
				WSCFitsWriter writer(createWSCFitsWriter(entry, isImaginary, false));
				double* dirtyImage = _imageAllocator.Allocate(_settings.trimmedImageWidth * _settings.trimmedImageHeight);
				_residualImages.Load(dirtyImage, entry.polarization, entry.outputChannelIndex, isImaginary);
				Logger::Info << "Writing dirty image...\n";
				writer.WriteImage("dirty.fits", dirtyImage);
				_imageAllocator.Free(dirtyImage);
			}
		}
	}

	if (isLastPol && (_settings.gridWithBeam || !_settings.atermConfigFilename.empty()))
	{
		ImageFilename imageName = ImageFilename(entry.outputChannelIndex, entry.outputIntervalIndex);
		static_cast<IdgMsGridder&>(*_gridder).SaveBeamImage(entry, imageName);
	}
	_currentPolMSes.clear();
}

void WSClean::makeMFSImage(const string& suffix, size_t intervalIndex, PolarizationEnum pol, bool isImaginary, bool isPSF)
{
	double lowestFreq = 0.0, highestFreq = 0.0;
	const size_t size = _settings.trimmedImageWidth * _settings.trimmedImageHeight;
	ao::uvector<double> mfsImage(size, 0.0), addedImage(size), weightImage(size, 0.0);
	double weightSum = 0.0;
	FitsWriter writer;
	for(size_t ch=0; ch!=_settings.channelsOut; ++ch)
	{
		std::string prefixStr = isPSF ?
			ImageFilename::GetPSFPrefix(_settings, ch, intervalIndex) :
			ImageFilename::GetPrefix(_settings, pol, ch, intervalIndex, isImaginary);
		const std::string name(prefixStr + '-' + suffix);
		FitsReader reader(name);
		if(ch == 0)
		{
			WSCFitsWriter wscWriter(reader);
			writer = wscWriter.Writer();
			lowestFreq = reader.Frequency() - reader.Bandwidth()*0.5;
			highestFreq = reader.Frequency() + reader.Bandwidth()*0.5;
		}
		else {
			lowestFreq = std::min(lowestFreq, reader.Frequency() - reader.Bandwidth()*0.5);
			highestFreq = std::max(highestFreq, reader.Frequency() + reader.Bandwidth()*0.5);
		}
		double weight;
		if(!reader.ReadDoubleKeyIfExists("WSCIMGWG", weight))
		{
			Logger::Error << "Error: image " << name << " did not have the WSCIMGWG keyword.\n";
			weight = 0.0;
		}
		weightSum += weight;
		reader.Read(addedImage.data());
		for(size_t i=0; i!=size; ++i)
		{
			if(std::isfinite(addedImage[i])) {
				mfsImage[i] += addedImage[i] * weight;
				weightImage[i] += weight;
			}
		}
	}
	for(size_t i=0; i!=size; ++i)
		mfsImage[i] /= weightImage[i];
	
	if(isPSF)
	{
		double smallestTheoreticBeamSize = 0.0;
		for(std::vector<OutputChannelInfo>::const_reverse_iterator r = _infoPerChannel.rbegin();
				r != _infoPerChannel.rend(); ++r)
		{
			if(std::isfinite(r->theoreticBeamSize)) {
				smallestTheoreticBeamSize = r->theoreticBeamSize;
				break;
			}
		}
		
		double bMaj, bMin, bPA;
		determineBeamSize(bMaj, bMin, bPA, mfsImage.data(), smallestTheoreticBeamSize);
		_infoForMFS.beamMaj = bMaj;
		_infoForMFS.beamMin = bMin;
		_infoForMFS.beamPA = bPA;
	}
	if(std::isfinite(_infoForMFS.beamMaj))
		writer.SetBeamInfo(_infoForMFS.beamMaj, _infoForMFS.beamMin, _infoForMFS.beamPA);
	else
		writer.SetNoBeamInfo();
	
	std::string mfsName(ImageFilename::GetMFSPrefix(_settings, pol, intervalIndex, isImaginary, isPSF) + '-' + suffix);
	Logger::Info << "Writing " << mfsName << "...\n";
	writer.SetFrequency((lowestFreq+highestFreq)*0.5, highestFreq-lowestFreq);
	writer.SetExtraKeyword("WSCIMGWG", weightSum);
	writer.RemoveExtraKeyword("WSCCHANS");
	writer.RemoveExtraKeyword("WSCCHANE");
	writer.Write(mfsName, mfsImage.data());
}

void WSClean::renderMFSImage(size_t intervalIndex, PolarizationEnum pol, bool isImaginary, bool isPBCorrected) const
{
	const size_t size = _settings.trimmedImageWidth * _settings.trimmedImageHeight;
	
	std::string mfsPrefix(ImageFilename::GetMFSPrefix(_settings, pol, intervalIndex, isImaginary, false));
	std::string postfix = isPBCorrected ? "-pb.fits" : ".fits";
	FitsReader residualReader(mfsPrefix + "-residual" + postfix);
	FitsReader modelReader(mfsPrefix + "-model" + postfix);
	ao::uvector<double> image(size), modelImage(size);
	residualReader.Read(image.data());
	modelReader.Read(modelImage.data());
	
	double beamMaj = _infoForMFS.beamMaj;
	double beamMin, beamPA;
	std::string beamStr;
	if(std::isfinite(beamMaj))
	{
		beamMin = _infoForMFS.beamMin;
		beamPA = _infoForMFS.beamPA;
		beamStr = "(beam=" + Angle::ToNiceString(beamMin) + "-" +
		Angle::ToNiceString(beamMaj) + ", PA=" +
		Angle::ToNiceString(beamPA) + ")";
	}
	else {
		beamStr = "(beam is neither fitted nor estimated -- using delta scales!)";
		beamMaj = 0.0; beamMin = 0.0; beamPA = 0.0;
	}
	Logger::Info << "Rendering sources to restored image " + beamStr + "... ";
	Logger::Info.Flush();
	ModelRenderer::Restore(image.data(), modelImage.data(), _settings.trimmedImageWidth, _settings.trimmedImageHeight, beamMaj, beamMin, beamPA, _settings.pixelScaleX, _settings.pixelScaleY);
	Logger::Info << "DONE\n";
	
	Logger::Info << "Writing " << mfsPrefix << "-image" << postfix << "...\n";
	FitsWriter imageWriter(residualReader);
	imageWriter.Write(mfsPrefix + "-image" + postfix, image.data());
}

MSSelection WSClean::selectInterval(MSSelection& fullSelection, size_t intervalIndex)
{
	if(_settings.intervalsOut == 1)
		return fullSelection;
	else {
		size_t tS, tE;
		if(fullSelection.HasInterval())
		{
			tS = fullSelection.IntervalStart();
			tE = fullSelection.IntervalEnd();
		}
		else {
			casacore::MeasurementSet ms(_settings.filenames[0]);
			Logger::Info << "Counting number of scans... ";
			Logger::Info.Flush();
			casacore::ROScalarColumn<double> timeColumn(ms, casacore::MS::columnName(casacore::MSMainEnums::TIME));
			double time = timeColumn(0);
			size_t timestepIndex = 1;
			for(size_t row = 0; row!=ms.nrow(); ++row)
			{
				if(time != timeColumn(row))
				{
					++timestepIndex;
					time = timeColumn(row);
				}
			}
			Logger::Info << "DONE (" << timestepIndex << ")\n";
			tS = 0;
			tE = timestepIndex;
			// Store the full interval in the selection, so that it doesn't need to be determined again.
			fullSelection.SetInterval(tS, tE);
		}
		if(_settings.intervalsOut > tE-tS)
		{
			std::ostringstream str;
			str << "Invalid interval selection: " << _settings.intervalsOut << " intervals requested, but measurement set has only " << tE-tS << " intervals.";
			throw std::runtime_error(str.str());
		}
		MSSelection newSelection(fullSelection);
		newSelection.SetInterval(
			tS + (tE-tS) * intervalIndex / _settings.intervalsOut,
			tS + (tE-tS) * (intervalIndex+1) / _settings.intervalsOut
		);
		return newSelection;
	}
}

void WSClean::saveUVImage(const double* image, PolarizationEnum pol, const ImagingTableEntry& entry, bool isImaginary, const std::string& prefix) const
{
	Image
		realUV(_settings.trimmedImageWidth, _settings.trimmedImageHeight, _imageAllocator),
		imagUV(_settings.trimmedImageWidth, _settings.trimmedImageHeight, _imageAllocator);
	FFTResampler fft(_settings.trimmedImageWidth, _settings.trimmedImageHeight, _settings.trimmedImageWidth, _settings.trimmedImageHeight, 1, true);
	fft.SingleFT(image, realUV.data(), imagUV.data());
	// Factors of 2 involved: because of SingleFT()
	// (also one from the fact that normF excludes a factor of two?)
	realUV *= _infoPerChannel[entry.outputChannelIndex].normalizationFactor / sqrt(0.5*_settings.trimmedImageWidth * _settings.trimmedImageHeight);
	imagUV *= _infoPerChannel[entry.outputChannelIndex].normalizationFactor / sqrt(0.5*_settings.trimmedImageWidth * _settings.trimmedImageHeight);
	WSCFitsWriter writer(createWSCFitsWriter(entry, isImaginary, false));
	writer.WriteUV(prefix+"-real.fits", realUV.data());
	writer.WriteUV(prefix+"-imag.fits", imagUV.data());
}

void WSClean::makeImagingTable(size_t outputIntervalIndex)
{
	std::set<ChannelInfo> channelSet;
	_msBands.assign(_settings.filenames.size(), MultiBandData());
	for(size_t i=0; i!=_settings.filenames.size(); ++i)
	{
		casacore::MeasurementSet ms(_settings.filenames[i]);
		_msBands[i] = MultiBandData(ms.spectralWindow(), ms.dataDescription());
		std::set<size_t> dataDescIds = _msBands[i].GetUsedDataDescIds(ms);
		if(dataDescIds.size() != _msBands[i].DataDescCount())
		{
			Logger::Debug << dataDescIds.size() << "/" << _msBands[i].DataDescCount() << " spws are used of " << _settings.filenames[i] << '\n';
		}
		
		// Apply user selection: remove unselected spws
		if(!_settings.spectralWindows.empty())
		{
			for(std::set<size_t>::iterator d = dataDescIds.begin(); d!=dataDescIds.end(); )
			{
				if(_settings.spectralWindows.find(_msBands[i].GetBandIndex(*d)) == _settings.spectralWindows.end())
					d = dataDescIds.erase(d);
				else
					++d;
			}
		}
		// accumulate channel info
		for(const size_t dataDescId : dataDescIds)
		{
			bool increasing = true;
			if(_msBands[i][dataDescId].ChannelCount() >= 2)
			{
				increasing = _msBands[i][dataDescId].Channel(1) > _msBands[i][dataDescId].Channel(0);
			}
			channelSet.insert(_msBands[i][dataDescId].Channel(0));
			for(size_t ch=1; ch!=_msBands[i][dataDescId].ChannelCount(); ++ch)
			{
				bool chanIncreasing = _msBands[i][dataDescId].Channel(ch) > _msBands[i][dataDescId].Channel(ch-1);
				if(chanIncreasing != increasing)
					throw std::runtime_error("Your measurement set has an incorrect frequency axis: the channels do neither only increase nor only decrease in frequency");
				if(_msBands[i][dataDescId].Channel(ch) == _msBands[i][dataDescId].Channel(ch-1))
					throw std::runtime_error("Your measurement set has an incorrect frequency axis: two adjacent channels had the same frequency. Channels should either strictly increase or strictly decrease in frequency.");
				channelSet.insert(_msBands[i][dataDescId].Channel(ch));
			}
		}
	}
	if(channelSet.size() < _settings.channelsOut)
	{
		std::ostringstream str;
		str << "Parameter '-channels-out' was set to an invalid value: " << _settings.channelsOut << " output channels requested, but combined in all specified measurement sets, there are only " << channelSet.size() << " unique channels.";
		throw std::runtime_error(str.str());
	}
	_inputChannelFrequencies.assign(channelSet.begin(), channelSet.end());
	Logger::Debug << "Total nr of channels found in measurement sets: " << _inputChannelFrequencies.size() << '\n';
	
	size_t joinedGroupIndex = 0, squaredGroupIndex = 0;
	_imagingTable.Clear();
	
	//for(size_t interval=0; interval!=_settings.intervalsOut; ++interval)
	//{
		if(_settings.joinedFrequencyCleaning)
		{
			size_t maxLocalJGI = joinedGroupIndex;
			for(size_t outChannelIndex=0; outChannelIndex!=_settings.channelsOut; ++outChannelIndex)
			{
				ImagingTableEntry freqTemplate;
				makeImagingTableEntry(_inputChannelFrequencies, outputIntervalIndex, outChannelIndex, freqTemplate);
				
				size_t localJGI = joinedGroupIndex;
				addPolarizationsToImagingTable(localJGI, squaredGroupIndex, outChannelIndex, freqTemplate);
				if(localJGI > maxLocalJGI)
					maxLocalJGI = localJGI;
			}
			joinedGroupIndex = maxLocalJGI;
		}
		else {
			for(size_t outChannelIndex=0; outChannelIndex!=_settings.channelsOut; ++outChannelIndex)
			{
				ImagingTableEntry freqTemplate;
				makeImagingTableEntry(_inputChannelFrequencies, outputIntervalIndex, outChannelIndex, freqTemplate);
				
				addPolarizationsToImagingTable(joinedGroupIndex, squaredGroupIndex, outChannelIndex, freqTemplate);
			}
		}
	//}
	_imagingTable.Update();
	_imagingTable.Print();
}

void WSClean::makeImagingTableEntry(const std::vector<ChannelInfo>& channels, size_t outIntervalIndex, size_t outChannelIndex, ImagingTableEntry& entry)
{
	size_t startCh, width;
	if(_settings.endChannel != 0)
	{
		if(_settings.endChannel > channels.size())
			throw std::runtime_error("Bad channel selection -- more channels selected than available");
		startCh = _settings.startChannel;
		width = _settings.endChannel-startCh;
	}
	else {
		startCh = 0;
		width = channels.size();
	}
	
	size_t chLowIndex, chHighIndex;
	if(_settings.divideChannelsByGaps)
	{
		std::multimap<double, size_t> gaps;
		for(size_t i = 1; i!=channels.size(); ++i)
		{
			double left = channels[i-1].Frequency();
			double right = channels[i].Frequency();
			gaps.insert(std::make_pair(right-left, i));
		}
		std::vector<size_t> orderedGaps;
		auto iter = gaps.rbegin();
		for(size_t i=0; i!=_settings.channelsOut-1; ++i)
		{
			orderedGaps.push_back(iter->second);
			++iter;
		}
		std::sort(orderedGaps.begin(), orderedGaps.end());
		if(outChannelIndex == 0)
			chLowIndex = 0;
		else
			chLowIndex = orderedGaps[outChannelIndex-1];
		if(outChannelIndex+1 == _settings.channelsOut)
			chHighIndex = width-1;
		else
			chHighIndex = orderedGaps[outChannelIndex]-1;
	}
	else {
		chLowIndex = startCh + outChannelIndex*width/_settings.channelsOut;
		chHighIndex = startCh + (outChannelIndex+1)*width/_settings.channelsOut - 1;
	}
	if(channels[chLowIndex].Frequency() > channels[chHighIndex].Frequency())
		std::swap(chLowIndex, chHighIndex);
	entry.inputChannelCount = chHighIndex+1 - chLowIndex;
	entry.lowestFrequency = channels[chLowIndex].Frequency();
	entry.highestFrequency = channels[chHighIndex].Frequency();
	entry.bandStartFrequency = entry.lowestFrequency - channels[chLowIndex].Width()*0.5;
	entry.bandEndFrequency = entry.highestFrequency + channels[chHighIndex].Width()*0.5;
	entry.outputIntervalIndex = outIntervalIndex;
	
	entry.msData.resize(_settings.filenames.size());
	for(size_t msIndex=0; msIndex!=_settings.filenames.size(); ++msIndex)
	{
		entry.msData[msIndex].bands.resize(_msBands[msIndex].DataDescCount());
	}
}

void WSClean::addPolarizationsToImagingTable(size_t& joinedGroupIndex, size_t& squaredGroupIndex, size_t outChannelIndex, const ImagingTableEntry& templateEntry)
{
	for(std::set<PolarizationEnum>::const_iterator p=_settings.polarizations.begin();
			p!=_settings.polarizations.end(); ++p)
	{
		ImagingTableEntry& entry = _imagingTable.AddEntry();
		entry = templateEntry;
		entry.index = _imagingTable.EntryCount()-1;
		entry.outputChannelIndex = outChannelIndex;
		entry.joinedGroupIndex = joinedGroupIndex;
		entry.squaredDeconvolutionIndex = squaredGroupIndex;
		entry.polarization = *p;
		if(*p == Polarization::XY)
			entry.imageCount = 2;
		else if(*p == Polarization::YX)
			entry.imageCount = 0;
		else
			entry.imageCount = 1;
		
		if(!_settings.joinedPolarizationCleaning)
		{
			++joinedGroupIndex;
			++squaredGroupIndex;
		}
	}
	
	if(_settings.joinedPolarizationCleaning)
	{
		++joinedGroupIndex;
		++squaredGroupIndex;
	}
}

void WSClean::fitBeamSize(double& bMaj, double& bMin, double& bPA, const double* image, double beamEstimate) const
{
	GaussianFitter beamFitter;
	Logger::Info << "Fitting beam... ";
	Logger::Info.Flush();
	if(_settings.circularBeam)
	{
		bMaj = beamEstimate;
		beamFitter.Fit2DCircularGaussianCentred(
			image,
			_settings.trimmedImageWidth, _settings.trimmedImageHeight,
			bMaj,
			_settings.beamFittingBoxSize);
		bMin = bMaj;
		bPA = 0.0;
	}
	else {
		beamFitter.Fit2DGaussianCentred(
			image,
			_settings.trimmedImageWidth, _settings.trimmedImageHeight,
			beamEstimate,
			bMaj, bMin, bPA,
			_settings.beamFittingBoxSize);
	}
	bMaj = bMaj*0.5*(_settings.pixelScaleX+_settings.pixelScaleY);
	bMin = bMin*0.5*(_settings.pixelScaleX+_settings.pixelScaleY);
}

void WSClean::determineBeamSize(double& bMaj, double& bMin, double& bPA, const double* image, double theoreticBeam) const
{
	double theoreticBeamWithTaper = theoreticBeam;
	if(_settings.gaussianTaperBeamSize != 0.0)
	{
		if(_settings.gaussianTaperBeamSize > theoreticBeamWithTaper)
		{
			theoreticBeamWithTaper = _settings.gaussianTaperBeamSize;
			Logger::Debug << "Beam is tapered; using " << Angle::ToNiceString(theoreticBeamWithTaper) << " as initial value in PSF fitting.\n";
		}
	}
	if(_settings.manualBeamMajorSize != 0.0)
	{
		bMaj = _settings.manualBeamMajorSize;
		bMin = _settings.manualBeamMinorSize;
		bPA = _settings.manualBeamPA;
	} else if(_settings.fittedBeam)
	{
		fitBeamSize(bMaj, bMin, bPA, image, theoreticBeamWithTaper*2.0/(_settings.pixelScaleX+_settings.pixelScaleY));
		Logger::Info << "major=" << Angle::ToNiceString(bMaj) << ", minor=" <<
		Angle::ToNiceString(bMin) << ", PA=" << Angle::ToNiceString(bPA) << ", theoretical=" <<
		Angle::ToNiceString(theoreticBeamWithTaper)<< ".\n";
	}
	else if(_settings.theoreticBeam) {
		bMaj = theoreticBeamWithTaper;
		bMin = theoreticBeamWithTaper;
		bPA = 0.0;
		Logger::Info << "Beam size is " << Angle::ToNiceString(theoreticBeamWithTaper) << '\n';
	} else {
		bMaj = std::numeric_limits<double>::quiet_NaN();
		bMin = std::numeric_limits<double>::quiet_NaN();
		bPA = std::numeric_limits<double>::quiet_NaN();
	}
}

WSCFitsWriter WSClean::createWSCFitsWriter(const ImagingTableEntry& entry, bool isImaginary, bool isModel) const
{
	return WSCFitsWriter(entry, isImaginary, _settings, _deconvolution, _majorIterationNr, *_gridder, _commandLine, _infoPerChannel[entry.outputChannelIndex], isModel);
}

WSCFitsWriter WSClean::createWSCFitsWriter(const ImagingTableEntry& entry, PolarizationEnum polarization, bool isImaginary, bool isModel) const
{
	return WSCFitsWriter(entry, polarization, isImaginary, _settings, _deconvolution, _majorIterationNr, *_gridder, _commandLine, _infoPerChannel[entry.outputChannelIndex], isModel);
}
