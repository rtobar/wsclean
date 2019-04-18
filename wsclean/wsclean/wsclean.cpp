#include "wsclean.h"

#include "directmsgridder.h"
#include "imagefilename.h"
#include "imageoperations.h"
#include "imageweightcache.h"
#include "griddingtaskmanager.h"
#include "logger.h"
#include "measurementsetgridder.h"
#include "primarybeam.h"
#include "wscfitswriter.h"

#include "../units/angle.h"

#include "../application.h"
#include "../areaset.h"
#include "../dftpredictionalgorithm.h"
#include "../fftresampler.h"
#include "../fitswriter.h"
#include "../image.h"
#include "../imageweights.h"
#include "../modelrenderer.h"
#include "../msselection.h"
#include "../msproviders/contiguousms.h"
#include "../progressbar.h"
#include "../uvector.h"

#include "../aocommon/parallelfor.h"

#include "../deconvolution/deconvolutionalgorithm.h"
#include "../deconvolution/imageset.h"

#include "../idg/idgmsgridder.h"

#include "../lofar/lmspredicter.h"

#include "../model/model.h"

#include <iostream>
#include <functional>
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
	if(factor != 1.0)
	{
		size_t nPix = _settings.trimmedImageWidth * _settings.trimmedImageHeight;
		for(size_t i=0; i!=nPix; ++i)
			image[i] *= factor;
	}
}

void WSClean::imagePSF(ImagingTableEntry& entry)
{
	Logger::Info.Flush();
	Logger::Info << " == Constructing PSF ==\n";
	_inversionWatch.Start();
	GriddingTask task;
	task.operation = GriddingTask::Invert;
	task.imagePSF = true;
	task.polarization = entry.polarization;
	task.subtractModel = false;
	task.verbose = _isFirstInversion;
	task.cache = &_msGridderMetaCache[entry.index];
	task.storeImagingWeights = _settings.writeImagingWeightSpectrumColumn;
	initializeCurMSProviders(entry, task);
	task.precalculatedWeightInfo = &initializeImageWeights(entry, task.msList);
	
	_griddingTaskManager->Run(task, std::bind(&WSClean::imagePSFCallback, this, std::ref(entry), std::placeholders::_1));
}

void WSClean::imagePSFCallback(ImagingTableEntry& entry, GriddingResult& result)
{
	size_t centralIndex = _settings.trimmedImageWidth/2 + (_settings.trimmedImageHeight/2) * _settings.trimmedImageWidth;
	double normFactor;
	if(result.imageRealResult[centralIndex] != 0.0)
		normFactor = 1.0/result.imageRealResult[centralIndex];
	else
		normFactor = 0.0;
	
	size_t channelIndex = entry.outputChannelIndex;
	_infoPerChannel[channelIndex].psfNormalizationFactor = normFactor;
	multiplyImage(normFactor, result.imageRealResult);
	Logger::Debug << "Normalized PSF by factor of " << normFactor << ".\n";
		
	DeconvolutionAlgorithm::RemoveNaNsInPSF(result.imageRealResult.data(), _settings.trimmedImageWidth, _settings.trimmedImageHeight);
	_psfImages.SetFitsWriter(createWSCFitsWriter(entry, false, false).Writer());
	_psfImages.Store(result.imageRealResult.data(), *_settings.polarizations.begin(), channelIndex, false);
	_inversionWatch.Pause();
	
	_phaseCentreRA = result.phaseCentreRA;
	_phaseCentreDec = result.phaseCentreDec;
	
	double bMaj, bMin, bPA;
	ImageOperations::DetermineBeamSize(_settings, bMaj, bMin, bPA, result.imageRealResult.data(), result.beamSize);
	entry.imageWeight = result.imageWeight;
	entry.normalizationFactor = result.normalizationFactor;
	_infoPerChannel[channelIndex].theoreticBeamSize = result.beamSize;
	_infoPerChannel[channelIndex].beamMaj = bMaj;
	_infoPerChannel[channelIndex].beamMin = bMin;
	_infoPerChannel[channelIndex].beamPA = bPA;
	_infoPerChannel[channelIndex].weight = entry.imageWeight;
	_infoPerChannel[channelIndex].normalizationFactor = entry.normalizationFactor;
	_infoPerChannel[channelIndex].wGridSize = result.actualWGridSize;
	_infoPerChannel[channelIndex].visibilityCount = result.griddedVisibilityCount;
	_infoPerChannel[channelIndex].effectiveVisibilityCount = result.effectiveGriddedVisibilityCount;
	_infoPerChannel[channelIndex].visibilityWeightSum = result.visibilityWeightSum;
	
	if(_settings.isUVImageSaved)
	{
		saveUVImage(result.imageRealResult.data(), entry, false, "uvpsf");
	}
	
	Logger::Info << "Writing psf image... ";
	Logger::Info.Flush();
	const std::string name(ImageFilename::GetPSFPrefix(_settings, channelIndex, entry.outputIntervalIndex) + "-psf.fits");
	WSCFitsWriter fitsFile = createWSCFitsWriter(entry, false, false);
	fitsFile.WritePSF(name, result.imageRealResult.data());
	Logger::Info << "DONE\n";
	
	if(_settings.isGriddingImageSaved && _isFirstInversion && result.hasGriddingCorrectionImage)
	{
		Logger::Info << "Writing gridding correction image... ";
		Logger::Info.Flush();
		double* gridding = _imageAllocator.Allocate(result.actualInversionWidth * result.actualInversionHeight);
		_griddingTaskManager->Gridder()->GetGriddingCorrectionImage(&gridding[0]);
		FitsWriter fitsWriter;
		fitsWriter.SetImageDimensions(result.actualInversionWidth, result.actualInversionHeight);
		fitsWriter.Write(_settings.prefixName + "-gridding.fits", &gridding[0]);
		_imageAllocator.Free(gridding);
		Logger::Info << "DONE\n";
	}
	
	_isFirstInversion = false;
}

void WSClean::imageMain(ImagingTableEntry& entry, bool isFirst, bool updateBeamInfo, bool isInitialInversion)
{
	Logger::Info.Flush();
	Logger::Info << " == Constructing image ==\n";
	_inversionWatch.Start();
	
	GriddingTask task;
	task.operation = GriddingTask::Invert;
	task.imagePSF = false;
	task.polarization = entry.polarization;
	task.subtractModel = !isFirst || _settings.subtractModel || _settings.continuedRun;
	task.verbose = isFirst && _isFirstInversion;
	task.cache = &_msGridderMetaCache[entry.index];
	task.storeImagingWeights = !isFirst && _settings.writeImagingWeightSpectrumColumn;
	initializeCurMSProviders(entry, task);
	task.precalculatedWeightInfo = &initializeImageWeights(entry, task.msList);
	
	_griddingTaskManager->Run(task, std::bind(&WSClean::imageMainCallback, this, std::ref(entry), std::placeholders::_1, updateBeamInfo, isInitialInversion));
	
	_inversionWatch.Pause();
}

void WSClean::imageMainCallback(ImagingTableEntry& entry, GriddingResult& result, bool updateBeamInfo, bool isInitialInversion)
{
	size_t joinedChannelIndex = entry.outputChannelIndex;
	
	multiplyImage(_infoPerChannel[joinedChannelIndex].psfNormalizationFactor, result.imageRealResult);
	storeAndCombineXYandYX(_residualImages, entry.polarization, joinedChannelIndex, false, result.imageRealResult.data());
	if(Polarization::IsComplex(entry.polarization))
	{
		multiplyImage(_infoPerChannel[joinedChannelIndex].psfNormalizationFactor, result.imageImaginaryResult);
		storeAndCombineXYandYX(_residualImages, entry.polarization, joinedChannelIndex, true, result.imageImaginaryResult.data());
	}
	entry.imageWeight = result.imageWeight;
	entry.normalizationFactor = result.normalizationFactor;
	_infoPerChannel[entry.outputChannelIndex].weight = result.imageWeight;
	_infoPerChannel[entry.outputChannelIndex].normalizationFactor = result.normalizationFactor;
	
	// If no PSF is made, also set the beam size. If the PSF was made, these would already be set
	// after imaging the PSF.
	if(updateBeamInfo)
	{
		if(_settings.theoreticBeam) {
			_infoPerChannel[entry.outputChannelIndex].beamMaj = result.beamSize;
			_infoPerChannel[entry.outputChannelIndex].beamMin = result.beamSize;
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
	
	if(isInitialInversion)
	{
		_modelImages.SetFitsWriter( createWSCFitsWriter(entry, false, true).Writer() );
		_residualImages.SetFitsWriter( createWSCFitsWriter(entry, false, false).Writer() );
		
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
	ImageBufferAllocator::Ptr
		modelImageReal(_imageAllocator.AllocatePtr(size)),
		modelImageImaginary;
		
	if(entry.polarization == Polarization::YX)
	{
		_modelImages.Load(modelImageReal.data(), Polarization::XY, entry.outputChannelIndex, false);
		modelImageImaginary = _imageAllocator.AllocatePtr(size);
		_modelImages.Load(modelImageImaginary.data(), Polarization::XY, entry.outputChannelIndex, true);
		for(size_t i=0; i!=size; ++i)
			modelImageImaginary[i] = -modelImageImaginary[i];
	}
	else {
		_modelImages.Load(modelImageReal.data(), entry.polarization, entry.outputChannelIndex, false);
		if(Polarization::IsComplex(entry.polarization))
		{
			modelImageImaginary = _imageAllocator.AllocatePtr(size);
			_modelImages.Load(modelImageImaginary.data(), entry.polarization, entry.outputChannelIndex, true);
		}
	}
	
	_predictingWatch.Start();
	GriddingTask task;
	task.operation = GriddingTask::Predict;
	task.polarization = entry.polarization;
	task.addToModel = false;
	task.cache = &_msGridderMetaCache[entry.index];
	task.verbose = false;
	task.storeImagingWeights = false;
	task.modelImageReal = std::move(modelImageReal);
	task.modelImageImaginary = std::move(modelImageImaginary);
	initializeCurMSProviders(entry, task);
	task.precalculatedWeightInfo = &initializeImageWeights(entry, task.msList);
	_griddingTaskManager->Run(task, std::bind(&WSClean::predictCallback, this, std::ref(entry), std::placeholders::_1));
	_predictingWatch.Pause();
}

void WSClean::predictCallback(const ImagingTableEntry&, GriddingResult&)
{ }

ImageWeights& WSClean::initializeImageWeights(const ImagingTableEntry& entry, std::vector<std::pair<std::unique_ptr<MSProvider>, MSSelection>>& msList)
{
	if(!_settings.mfWeighting)
	{
		_imageWeightCache->Update(msList, entry.outputChannelIndex, entry.outputIntervalIndex);
		if(_settings.isWeightImageSaved)
			_imageWeightCache->Weights().Save(_settings.prefixName+"-weights.fits");
	}
	return _imageWeightCache->Weights(); // task.precalculatedWeightInfo =
}

void WSClean::initializeMFSImageWeights()
{
	Logger::Info << "Precalculating MF weights for " << _settings.weightMode.ToString() << " weighting...\n";
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

void WSClean::performReordering(bool isPredictMode)
{
	std::mutex mutex;
	_partitionedMSHandles.resize(_settings.filenames.size());
	bool useModel = _settings.deconvolutionMGain != 1.0 || isPredictMode || _settings.subtractModel || _settings.continuedRun;
	bool initialModelRequired = _settings.subtractModel || _settings.continuedRun;
	
	if(_settings.parallelReordering!=1)
		Logger::Info << "Reordering...\n";
	
	ao::ParallelFor<size_t> loop(_settings.parallelReordering);
	loop.Run(0, _settings.filenames.size(), [&](size_t i, size_t)
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
		
		PartitionedMS::Handle partMS = PartitionedMS::Partition(_settings.filenames[i], channels, _globalSelection, _settings.dataColumnName, useModel, initialModelRequired, _settings);
		std::lock_guard<std::mutex> lock(mutex);
		_partitionedMSHandles[i] = std::move(partMS);
		if(_settings.parallelReordering!=1)
			Logger::Info << "Finished reordering " << _settings.filenames[i] << " [" << i << "]\n";
	});
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
		
		_msGridderMetaCache.clear();
		_imageWeightCache = createWeightCache();
		
		if(_settings.mfWeighting)
			initializeMFSImageWeights();
		
		_griddingTaskManager.reset(new GriddingTaskManager(_settings, _imageAllocator));
		
		std::unique_ptr<PrimaryBeam> primaryBeam;
		for(size_t groupIndex=0; groupIndex!=_imagingTable.IndependentGroupCount(); ++groupIndex)
		{
			ImagingTable group = _imagingTable.GetIndependentGroup(groupIndex);
			runIndependentGroup(group, primaryBeam);
		}

		// Needs to be destructed before image allocator, or image allocator will report error caused by leaked memory
		_griddingTaskManager.reset();
	
		if(_settings.channelsOut > 1)
		{
			for(std::set<PolarizationEnum>::const_iterator pol=_settings.polarizations.begin(); pol!=_settings.polarizations.end(); ++pol)
			{
				bool psfWasMade = (_settings.deconvolutionIterationCount > 0 || _settings.makePSF || _settings.makePSFOnly) && pol == _settings.polarizations.begin();
				
				if(psfWasMade)
				{
					ImageOperations::MakeMFSImage(_settings, _infoPerChannel, _infoForMFS, "psf.fits", intervalIndex, *pol, false, true);
					if(_settings.savePsfPb)
						ImageOperations::MakeMFSImage(_settings, _infoPerChannel, _infoForMFS, "psf-pb.fits", intervalIndex, *pol, false, true);
				}
				
				if(!(*pol == Polarization::YX && _settings.polarizations.count(Polarization::XY)!=0) && !_settings.makePSFOnly)
				{
					if(_settings.isDirtySaved)
						ImageOperations::MakeMFSImage(_settings, _infoPerChannel, _infoForMFS, "dirty.fits", intervalIndex, *pol, false);
					if(_settings.deconvolutionIterationCount == 0)
					{
						ImageOperations::MakeMFSImage(_settings, _infoPerChannel, _infoForMFS, "image.fits", intervalIndex, *pol, false);
						if(_settings.applyPrimaryBeam || (_settings.gridWithBeam || !_settings.atermConfigFilename.empty()))
							ImageOperations::MakeMFSImage(_settings, _infoPerChannel, _infoForMFS, "image-pb.fits", intervalIndex, *pol, false);
					}
					else 
					{
						ImageOperations::MakeMFSImage(_settings, _infoPerChannel, _infoForMFS, "residual.fits", intervalIndex, *pol, false);
						ImageOperations::MakeMFSImage(_settings, _infoPerChannel, _infoForMFS, "model.fits", intervalIndex, *pol, false);
						ImageOperations::RenderMFSImage(_settings, _infoForMFS, intervalIndex, *pol, false, false);
						if(_settings.applyPrimaryBeam || (_settings.gridWithBeam || !_settings.atermConfigFilename.empty()))
						{
							ImageOperations::MakeMFSImage(_settings, _infoPerChannel, _infoForMFS, "residual-pb.fits", intervalIndex, *pol, false);
							ImageOperations::MakeMFSImage(_settings, _infoPerChannel, _infoForMFS, "model-pb.fits", intervalIndex, *pol, false);
							ImageOperations::RenderMFSImage(_settings, _infoForMFS, intervalIndex, *pol, false, true);
						}
					}
					if(Polarization::IsComplex(*pol))
					{
						if(_settings.isDirtySaved)
							ImageOperations::MakeMFSImage(_settings, _infoPerChannel, _infoForMFS, "dirty.fits", intervalIndex, *pol, true);
						if(_settings.deconvolutionIterationCount == 0)
						{
								ImageOperations::MakeMFSImage(_settings, _infoPerChannel, _infoForMFS, "image.fits", intervalIndex, *pol, true);
							if(_settings.applyPrimaryBeam || (_settings.gridWithBeam || !_settings.atermConfigFilename.empty()))
								ImageOperations::MakeMFSImage(_settings, _infoPerChannel, _infoForMFS, "image-pb.fits", intervalIndex, *pol, true);
						}
						else {
							ImageOperations::MakeMFSImage(_settings, _infoPerChannel, _infoForMFS, "residual.fits", intervalIndex, *pol, true);
							ImageOperations::MakeMFSImage(_settings, _infoPerChannel, _infoForMFS, "model.fits", intervalIndex, *pol, true);
							ImageOperations::RenderMFSImage(_settings, _infoForMFS, intervalIndex, *pol, true, false);
							if(_settings.applyPrimaryBeam || (_settings.gridWithBeam || !_settings.atermConfigFilename.empty()))
							{
								ImageOperations::MakeMFSImage(_settings, _infoPerChannel, _infoForMFS, "residual-pb.fits", intervalIndex, *pol, true);
								ImageOperations::MakeMFSImage(_settings, _infoPerChannel, _infoForMFS, "model-pb.fits", intervalIndex, *pol, true);
								ImageOperations::RenderMFSImage(_settings, _infoForMFS, intervalIndex, *pol, true, true);
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
		_msGridderMetaCache.clear();
		
		_globalSelection = selectInterval(fullSelection, intervalIndex);
		
		_doReorder = preferReordering();
		
		if(_doReorder) performReordering(true);
		
		_griddingTaskManager.reset(new GriddingTaskManager(_settings, _imageAllocator));
	
		for(size_t groupIndex=0; groupIndex!=_imagingTable.SquaredGroupCount(); ++groupIndex)
		{
			predictGroup(_imagingTable.GetSquaredGroup(groupIndex));
		}
	
		// Needs to be destructed before image allocator, or image allocator will report error caused by leaked memory
		_griddingTaskManager.reset();
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

void WSClean::runIndependentGroup(ImagingTable& groupTable, std::unique_ptr<PrimaryBeam>& primaryBeam)
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
		bool isFirstPol = entry.polarization == *_settings.polarizations.begin();
		bool doMakePSF = _settings.deconvolutionIterationCount > 0 || _settings.makePSF || _settings.makePSFOnly;
		if(doMakePSF && isFirstPol)
			imagePSF(entry);
	}
	_griddingTaskManager->Finish();
	
	for(size_t joinedIndex=0; joinedIndex!=groupTable.EntryCount(); ++joinedIndex)
	{
		runFirstInversion(groupTable[joinedIndex], primaryBeam);
	}
	_griddingTaskManager->Finish();

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
					if(_settings.useIDG)
					{
						// For now in the case of IDG we have to directly ask for all Four polarizations. This can't 
						// be parallelized yet in the current structure.
						for(size_t sGroupIndex=0; sGroupIndex!=groupTable.SquaredGroupCount(); ++sGroupIndex)
						{
							ImagingTable sGroupTable = groupTable.GetSquaredGroup(sGroupIndex);
							for(size_t e=0; e!=sGroupTable.EntryCount(); ++e)
							{
								predict(sGroupTable[e]);
							}
							_griddingTaskManager->Finish();
							for(size_t e=0; e!=sGroupTable.EntryCount(); ++e)
							{
								imageMain(sGroupTable[e], false, false, false);
							} // end of polarization loop
							_griddingTaskManager->Finish();
						} // end of joined channels loop
					}
					
					else { // not using IDG
						for(size_t sGroupIndex=0; sGroupIndex!=groupTable.SquaredGroupCount(); ++sGroupIndex)
						{
							ImagingTable sGroupTable = groupTable.GetSquaredGroup(sGroupIndex);
							for(size_t e=0; e!=sGroupTable.EntryCount(); ++e)
							{
								predict(sGroupTable[e]);
							}
						}
						_griddingTaskManager->Finish();
						
						for(size_t sGroupIndex=0; sGroupIndex!=groupTable.SquaredGroupCount(); ++sGroupIndex)
						{
							ImagingTable sGroupTable = groupTable.GetSquaredGroup(sGroupIndex);
							for(size_t e=0; e!=sGroupTable.EntryCount(); ++e)
							{
								imageMain(sGroupTable[e], false, false, false);
							} // end of polarization loop
						} // end of joined channels loop
						_griddingTaskManager->Finish();
					}
				}
				
				++_majorIterationNr;
			} while(reachedMajorThreshold);
			
			--_majorIterationNr;
			Logger::Info << _majorIterationNr << " major iterations were performed.\n";
		}
		
		//_gridder->FreeImagingData();
		
		for(size_t joinedIndex=0; joinedIndex!=groupTable.EntryCount(); ++joinedIndex)
			saveRestoredImagesForGroup(groupTable[joinedIndex], primaryBeam);
		
		if(_settings.saveSourceList)
    {
			_deconvolution.SaveSourceList(groupTable, _phaseCentreRA, _phaseCentreDec);
			if(_settings.applyPrimaryBeam)
			{
				_deconvolution.SavePBSourceList(groupTable, _phaseCentreRA, _phaseCentreDec);
			}
    }
	}
	
	_deconvolution.FreeDeconvolutionAlgorithms();
	
	_imageAllocator.ReportStatistics();
	Logger::Info << "Inversion: " << _inversionWatch.ToString() << ", prediction: " << _predictingWatch.ToString() << ", deconvolution: " << _deconvolutionWatch.ToString() << '\n';
	
	_settings.prefixName = rootPrefix;
}

void WSClean::saveRestoredImagesForGroup(const ImagingTableEntry& tableEntry, std::unique_ptr<PrimaryBeam>& primaryBeam) const
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
			saveUVImage(restoredImage, tableEntry, isImaginary, "uv");
		
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
				primaryBeam->CorrectImages(writer.Writer(), imageName, "image", _imageAllocator);
				if(_settings.savePsfPb)
					primaryBeam->CorrectImages(writer.Writer(), imageName, "psf", _imageAllocator);
				if(_settings.deconvolutionIterationCount != 0)
				{
					primaryBeam->CorrectImages(writer.Writer(), imageName, "residual", _imageAllocator);
					primaryBeam->CorrectImages(writer.Writer(), imageName, "model", _imageAllocator);
				}
			}
			else if(_settings.gridWithBeam || !_settings.atermConfigFilename.empty())
			{
				IdgMsGridder& idg = static_cast<IdgMsGridder&>(*_griddingTaskManager->Gridder());
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
			if(_settings.mfWeighting)
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
		
		predict(entry);
	} // end of polarization loop
	
	_griddingTaskManager->Finish();
	
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

void WSClean::initializeCurMSProviders(const ImagingTableEntry& entry, GriddingTask& task)
{
	task.msList.clear();
	for(size_t i=0; i != _settings.filenames.size(); ++i)
	{
		for(size_t d=0; d!=_msBands[i].DataDescCount(); ++d)
		{
			MSSelection selection(_globalSelection);
			if(selectChannels(selection, i, d, entry))
			{
				std::unique_ptr<MSProvider> msProvider(
					initializeMSProvider(entry, selection, i, d));
				task.msList.emplace_back(std::move(msProvider), selection);
			}
		}
	}
}

void WSClean::initializeMSProvidersForPB(const ImagingTableEntry& entry, std::vector<std::pair<std::unique_ptr<MSProvider>, MSSelection>>& msList, PrimaryBeam& pb)
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
				msList.emplace_back(std::move(msProvider), selection);
			}
		}
	}
}

void WSClean::runFirstInversion(ImagingTableEntry& entry, std::unique_ptr<class PrimaryBeam>& primaryBeam)
{
	bool isLastPol = entry.polarization == *_settings.polarizations.rbegin();
	bool doMakePSF = _settings.deconvolutionIterationCount > 0 || _settings.makePSF || _settings.makePSFOnly;
	
	if(isLastPol)
	{
		ImageFilename imageName = ImageFilename(entry.outputChannelIndex, entry.outputIntervalIndex);
		if(_settings.applyPrimaryBeam)
		{
			primaryBeam.reset(new PrimaryBeam(_settings));
			std::vector<std::pair<std::unique_ptr<MSProvider>, MSSelection>> pbmsList;
			initializeMSProvidersForPB(entry, pbmsList, *primaryBeam);
			initializeImageWeights(entry, pbmsList);
			double ra, dec, dl, dm;
			casacore::MeasurementSet ms(pbmsList.front().first->MS());
			MSGridderBase::GetPhaseCentreInfo(ms, _settings.fieldId, ra, dec, dl, dm);
			primaryBeam->SetPhaseCentre(ra, dec, dl, dm);
			primaryBeam->MakeBeamImages(imageName, entry, _imageWeightCache.get(), _imageAllocator);
		}
	}
	
	if(!_settings.makePSFOnly)
	{
		imageMain(entry, true, !doMakePSF, true);
		
		_isFirstInversion = false;
	}

	if (isLastPol && (_settings.gridWithBeam || !_settings.atermConfigFilename.empty()))
	{
		ImageFilename imageName(entry.outputChannelIndex, entry.outputIntervalIndex);
		static_cast<IdgMsGridder&>(*_griddingTaskManager->Gridder()).SaveBeamImage(entry, imageName);
	}
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

void WSClean::saveUVImage(const double* image, const ImagingTableEntry& entry, bool isImaginary, const std::string& prefix) const
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
	std::vector<ChannelInfo> inputChannelFrequencies(channelSet.begin(), channelSet.end());
	Logger::Debug << "Total nr of channels found in measurement sets: " << inputChannelFrequencies.size() << '\n';
	
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
				makeImagingTableEntry(inputChannelFrequencies, outputIntervalIndex, outChannelIndex, freqTemplate);
				
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
				makeImagingTableEntry(inputChannelFrequencies, outputIntervalIndex, outChannelIndex, freqTemplate);
				
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

WSCFitsWriter WSClean::createWSCFitsWriter(const ImagingTableEntry& entry, bool isImaginary, bool isModel) const
{
	return WSCFitsWriter(entry, isImaginary, _settings, _deconvolution, _majorIterationNr, *_griddingTaskManager->Gridder(), _commandLine, _infoPerChannel[entry.outputChannelIndex], isModel);
}

WSCFitsWriter WSClean::createWSCFitsWriter(const ImagingTableEntry& entry, PolarizationEnum polarization, bool isImaginary, bool isModel) const
{
	return WSCFitsWriter(entry, polarization, isImaginary, _settings, _deconvolution, _majorIterationNr, *_griddingTaskManager->Gridder(), _commandLine, _infoPerChannel[entry.outputChannelIndex], isModel);
}
