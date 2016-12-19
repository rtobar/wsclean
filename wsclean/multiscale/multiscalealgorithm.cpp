#include "multiscalealgorithm.h"

#include "clarkloop.h"
#include "multiscaletransforms.h"

#include "../deconvolution/simpleclean.h"
#include "../fftwmultithreadenabler.h"
#include "../wsclean/logger.h"

MultiScaleAlgorithm::MultiScaleAlgorithm(ImageBufferAllocator& allocator, double beamSize, double pixelScaleX, double pixelScaleY) :
	_allocator(allocator),
	_width(0),
	_height(0),
	_beamSizeInPixels(beamSize / std::max(pixelScaleX, pixelScaleY)),
	_trackPerScaleMasks(false), _usePerScaleMasks(false),
	_fastSubMinorLoop(true)
{
}

MultiScaleAlgorithm::~MultiScaleAlgorithm()
{
	Logger::Info << "Multi-scale cleaning summary:\n";
	size_t sumComponents = 0;
	double sumFlux = 0.0;
	for(size_t scaleIndex=0; scaleIndex!=_scaleInfos.size(); ++scaleIndex)
	{
		const ScaleInfo& scaleEntry = _scaleInfos[scaleIndex];
		Logger::Info << "- Scale " << round(scaleEntry.scale) << " px, nr of components cleaned: " << scaleEntry.nComponentsCleaned << " (" << scaleEntry.totalFluxCleaned << " Jy)\n";
		sumComponents += scaleEntry.nComponentsCleaned;
		sumFlux += scaleEntry.totalFluxCleaned;
	}
	Logger::Info << "Total: " << sumComponents << " components (" << sumFlux << " Jy)\n";
}


void MultiScaleAlgorithm::ExecuteMajorIteration(DynamicSet& dirtySet, DynamicSet& modelSet, const ao::uvector<const double*>& psfs, size_t width, size_t height, bool& reachedMajorThreshold)
{
	// Rough overview of the procedure:
	// Convolve integrated image (all scales)
	// Find integrated peak & scale
	// Minor loop:
	// - Convolve individual images at fixed scale
	// - Subminor loop:
	//   - Measure individual peaks per individually convolved image
	//   - Subtract convolved PSF from individual images
	//   - Subtract double convolved PSF from individually convolved images
	//   - Find integrated peak at fixed scale
	// - Convolve integrated image (all scales)
	// - Find integrated peak & scale
	//
	// (This excludes creating the convolved PSFs and double-convolved PSFs
	//  at the appropriate moments).
	
	_width = width;
	_height = height;
		
	// The threads always need to be stopped at the end of this function, so we use a scoped
	// unique ptr.
	std::unique_ptr<ThreadedDeconvolutionTools> tools(new ThreadedDeconvolutionTools(_threadCount));
	_tools = tools.get();
	
	initializeScaleInfo();
	
	if(_trackPerScaleMasks && _scaleMasks.empty())
	{
		_scaleMasks.resize(_scaleInfos.size());
		for(size_t s=0; s!=_scaleInfos.size(); ++s)
			_scaleMasks[s].assign(_width*_height, false);
	}
	
	ImageBufferAllocator::Ptr scratch, integratedScratch;
	_allocator.Allocate(_width*_height, scratch);
	_allocator.Allocate(_width*_height, integratedScratch);
	std::unique_ptr<std::unique_ptr<ImageBufferAllocator::Ptr[]>[]> convolvedPSFs(
		new std::unique_ptr<ImageBufferAllocator::Ptr[]>[dirtySet.PSFCount()]);
	dirtySet.GetIntegratedPSF(integratedScratch.data(), psfs);
	convolvePSFs(convolvedPSFs[0], integratedScratch.data(), scratch.data(), true);

	// If there's only one, the integrated equals the first, so we can skip this
	if(dirtySet.PSFCount() > 1)
	{
		for(size_t i=0; i!=dirtySet.PSFCount(); ++i)
		{
			convolvePSFs(convolvedPSFs[i], psfs[i], scratch.data(), false);
		}
	}
	
	MultiScaleTransforms msTransforms(_width, _height);
	
	size_t scaleWithPeak;
	findActiveScaleConvolvedMaxima(dirtySet, integratedScratch.data(), true);
	sortScalesOnMaxima(scaleWithPeak);
	
	double mGainThreshold = std::fabs(_scaleInfos[scaleWithPeak].maxImageValue * _scaleInfos[scaleWithPeak].biasFactor) * (1.0 - _mGain);
	const double firstThreshold = std::max(_threshold, mGainThreshold);
	
	Logger::Info << "Starting multi-scale cleaning. Start peak="
		<< _scaleInfos[scaleWithPeak].maxImageValue * _scaleInfos[scaleWithPeak].biasFactor
		<< " Jy, major iteration threshold=" << firstThreshold << "\n";
	
	std::unique_ptr<ImageBufferAllocator::Ptr[]> doubleConvolvedPSFs(
		new ImageBufferAllocator::Ptr[dirtySet.PSFCount()]);
	for(size_t i=0; i!=dirtySet.PSFCount(); ++i)
	{
		_allocator.Allocate(_width*_height, doubleConvolvedPSFs[i]);
	}
	
	DynamicSet individualConvolvedImages(&dirtySet.Table(), dirtySet.Allocator(), dirtySet.ChannelsInDeconvolution(), dirtySet.SquareJoinedChannels(), _width, _height);
	
	//
	// The minor iteration loop
	//
	while(_iterationNumber < MaxNIter() && std::fabs(_scaleInfos[scaleWithPeak].maxImageValue * _scaleInfos[scaleWithPeak].biasFactor) > firstThreshold)
	{
		// Create double-convolved PSFs & individually convolved images for this scale
		ao::uvector<double*> transformList;
		for(size_t i=0; i!=dirtySet.PSFCount(); ++i)
		{
			double* psf = getConvolvedPSF(i, scaleWithPeak, psfs, scratch.data(), convolvedPSFs);
			memcpy(doubleConvolvedPSFs[i].data(), psf, _width*_height*sizeof(double));
			transformList.push_back(doubleConvolvedPSFs[i].data());
		}
		for(size_t i=0; i!=dirtySet.size(); ++i)
		{
			memcpy(individualConvolvedImages[i], dirtySet[i], _width*_height*sizeof(double));
			transformList.push_back(individualConvolvedImages[i]);
		}
		if(scaleWithPeak != 0)
		{
			_tools->MultiScaleTransform(&msTransforms, transformList, scratch.data(), _scaleInfos[scaleWithPeak].scale);
			//msTransforms.Transform(transformList, scratch.data(), _scaleInfos[scaleWithPeak].scale);
		}
		
		//
		// The sub-minor iteration loop for this scale
		//
		double firstSubIterationThreshold = std::max(
			std::fabs(_scaleInfos[scaleWithPeak].maxImageValue * _scaleInfos[scaleWithPeak].biasFactor) * (1.0 - _gain),
			firstThreshold);
		// TODO we could chose to run the non-fast loop until we hit e.g. 10 iterations in a scale,
		// because the fast loop takes more constant time and is only efficient when doing
		// many iterations.
		if(_fastSubMinorLoop)
		{
			FFTWMultiThreadEnabler fftwMT(false);
			size_t subMinorStartIteration = _iterationNumber;
			ClarkLoop clarkLoop(_width, _height);
			clarkLoop.SetIterationInfo(_iterationNumber, MaxNIter());
			clarkLoop.SetThreshold(firstSubIterationThreshold / _scaleInfos[scaleWithPeak].biasFactor);
			clarkLoop.SetGain(_scaleInfos[scaleWithPeak].gain);
			clarkLoop.SetAllowNegativeComponents(AllowNegativeComponents());
			if(_usePerScaleMasks)
				clarkLoop.SetMask(_scaleMasks[scaleWithPeak].data());
			else if(_cleanMask)
				clarkLoop.SetMask(_cleanMask);
			clarkLoop.SetSpectralFitter(&Fitter());
			
			std::vector<const double*> clarkPSFs(dirtySet.PSFCount());
			for(size_t psfIndex=0; psfIndex!=clarkPSFs.size(); ++psfIndex)
				clarkPSFs[psfIndex] = doubleConvolvedPSFs[psfIndex].data();
			
			clarkLoop.Run(individualConvolvedImages, clarkPSFs);
			
			_iterationNumber = clarkLoop.CurrentIteration();
			_scaleInfos[scaleWithPeak].nComponentsCleaned += (_iterationNumber - subMinorStartIteration);
			_scaleInfos[scaleWithPeak].totalFluxCleaned += clarkLoop.FluxCleaned();
			
			for(size_t imageIndex=0; imageIndex!=dirtySet.size(); ++imageIndex)
			{
				// TODO this can be multi-threaded if each thread has its own temporaries
				double *psf = getConvolvedPSF(dirtySet.PSFIndex(imageIndex), scaleWithPeak, psfs, scratch.data(), convolvedPSFs);
				clarkLoop.CorrectResidualDirty(integratedScratch.data(), scratch.data(), imageIndex, dirtySet[imageIndex],  psf);
				
				clarkLoop.GetFullIndividualModel(imageIndex, scratch.data());
				if(_trackPerScaleMasks && imageIndex==0)
				{
					// TODO this should be set from the sparse model to prevent traversing over all values
					clarkLoop.UpdateAutoMask(_scaleMasks[scaleWithPeak].data());
				}
				if(_scaleInfos[scaleWithPeak].scale != 0.0)
					_tools->MultiScaleTransform(&msTransforms, ao::uvector<double*>{scratch.data()}, integratedScratch.data(), _scaleInfos[scaleWithPeak].scale);
				double* model = modelSet[imageIndex];
				for(size_t i=0; i!=_width*_height; ++i)
					model[i] += scratch.data()[i];
			}
			
		}
		else { // don't use the Clark optimization
			while(_iterationNumber < MaxNIter() && std::fabs(_scaleInfos[scaleWithPeak].maxImageValue * _scaleInfos[scaleWithPeak].biasFactor) > firstSubIterationThreshold)
			{
				ao::uvector<double> componentValues;
				measureComponentValues(componentValues, scaleWithPeak, individualConvolvedImages);
				PerformSpectralFit(componentValues.data());
				
				for(size_t imgIndex=0; imgIndex!=dirtySet.size(); ++imgIndex)
				{
					// Subtract component from individual, non-deconvolved images
					double componentGain = componentValues[imgIndex] * _scaleInfos[scaleWithPeak].gain;
					
					double* psf = getConvolvedPSF(dirtySet.PSFIndex(imgIndex), scaleWithPeak, psfs, scratch.data(), convolvedPSFs);
					tools->SubtractImage(dirtySet[imgIndex], psf, _width, _height, _scaleInfos[scaleWithPeak].maxImageValueX, _scaleInfos[scaleWithPeak].maxImageValueY, componentGain);
					
					// Subtract double convolved PSFs from convolved images
					tools->SubtractImage(individualConvolvedImages[imgIndex], doubleConvolvedPSFs[dirtySet.PSFIndex(imgIndex)].data(), _width, _height, _scaleInfos[scaleWithPeak].maxImageValueX, _scaleInfos[scaleWithPeak].maxImageValueY, componentGain);
					// TODO this is incorrect, but why is the residual without Cotton-Schwab still OK ?
					// Should test
					//tools->SubtractImage(individualConvolvedImages[imgIndex], psf, _width, _height, _scaleInfos[scaleWithPeak].maxImageValueX, _scaleInfos[scaleWithPeak].maxImageValueY, componentGain);
					
					// Adjust model
					addComponentToModel(modelSet[imgIndex], scaleWithPeak, componentValues[imgIndex]);
				}
				
				// Find maximum for this scale
				individualConvolvedImages.GetLinearIntegrated(integratedScratch.data());
				findSingleScaleMaximum(integratedScratch.data(), scaleWithPeak);
				Logger::Debug << "Scale now " << std::fabs(_scaleInfos[scaleWithPeak].maxImageValue * _scaleInfos[scaleWithPeak].biasFactor) << '\n';
				
				++_iterationNumber;
			}
		}
		
		activateScales(scaleWithPeak);
		
		findActiveScaleConvolvedMaxima(dirtySet, integratedScratch.data(), false);
		sortScalesOnMaxima(scaleWithPeak);
		
		Logger::Info << "Iteration " << _iterationNumber << ", scale " << round(_scaleInfos[scaleWithPeak].scale) << " px : " << _scaleInfos[scaleWithPeak].maxImageValue*_scaleInfos[scaleWithPeak].biasFactor << " Jy at " << _scaleInfos[scaleWithPeak].maxImageValueX << ',' << _scaleInfos[scaleWithPeak].maxImageValueY << '\n';
	}
	
	reachedMajorThreshold =
		_iterationNumber < MaxNIter() && std::fabs(_scaleInfos[scaleWithPeak].maxImageValue * _scaleInfos[scaleWithPeak].biasFactor) > _threshold;
}

void MultiScaleAlgorithm::initializeScaleInfo()
{
	if(_scaleInfos.empty())
	{
		if(_manualScaleList.empty())
		{
			size_t scaleIndex = 0;
			double scale = _beamSizeInPixels * 2.0;
			while(scale < std::min(_width, _height)*0.5)
			{
				_scaleInfos.push_back(ScaleInfo());
				ScaleInfo& newEntry = _scaleInfos.back();
				if(scaleIndex == 0)
					newEntry.scale = 0.0;
				else
					newEntry.scale = scale;
				newEntry.kernelPeak = MultiScaleTransforms::KernelPeakValue(scale);
				
				scale *= 2.0;
				++scaleIndex;
			}
		}
		else {
			std::sort(_manualScaleList.begin(), _manualScaleList.end());
			for(size_t scaleIndex = 0; scaleIndex != _manualScaleList.size(); ++scaleIndex)
			{
				_scaleInfos.push_back(ScaleInfo());
				ScaleInfo& newEntry = _scaleInfos.back();
				newEntry.scale = _manualScaleList[scaleIndex];
				newEntry.kernelPeak = MultiScaleTransforms::KernelPeakValue(newEntry.scale);
			}
		}
	}
}

void MultiScaleAlgorithm::convolvePSFs(std::unique_ptr<ImageBufferAllocator::Ptr[]>& convolvedPSFs, const double* psf, double* tmp, bool isIntegrated)
{
	MultiScaleTransforms msTransforms(_width, _height);
	convolvedPSFs.reset(new ImageBufferAllocator::Ptr[_scaleInfos.size()]);
	if(isIntegrated)
		Logger::Info << "Scale info:\n";
	const double firstAutoScaleSize = _beamSizeInPixels * 2.0;
	for(size_t scaleIndex=0; scaleIndex!=_scaleInfos.size(); ++scaleIndex)
	{
		ScaleInfo& scaleEntry = _scaleInfos[scaleIndex];
		
		_allocator.Allocate(_width*_height, convolvedPSFs[scaleIndex]);
		memcpy(convolvedPSFs[scaleIndex].data(), psf, _width*_height*sizeof(double));
		
		if(isIntegrated)
		{
			msTransforms.Transform(convolvedPSFs[scaleIndex].data(), tmp, scaleEntry.scale);
			
			scaleEntry.psfPeak = convolvedPSFs[scaleIndex][_width/2 + (_height/2)*_width];
			// We normalize this factor to 1 for scale 0, so:
			// factor = (psf / kernel) / (psf0 / kernel0) = psf * kernel0 / (kernel * psf0)
			//scaleEntry.biasFactor = std::max(1.0,
			//	scaleEntry.psfPeak * scaleInfos[0].kernelPeak /
			//	(scaleEntry.kernelPeak * scaleInfos[0].psfPeak));
			double responseNormalization = _multiscaleNormalizeResponse ? scaleEntry.psfPeak : 1.0;
			double expTerm;
			if(scaleEntry.scale == 0.0 || _scaleInfos.size() < 2)
				expTerm = 0.0;
			else
				expTerm = log2(scaleEntry.scale / firstAutoScaleSize);
			scaleEntry.biasFactor = pow(_multiscaleScaleBias, -double(expTerm)) * 1.0 / responseNormalization;
			
			// I tried this, but wasn't perfect:
			// _gain * _scaleInfos[0].kernelPeak / scaleEntry.kernelPeak;
			scaleEntry.gain = _gain * _scaleInfos[0].psfPeak / scaleEntry.psfPeak;
			
			scaleEntry.isActive = true;
			
			if(scaleEntry.scale == 0.0)
				memcpy(convolvedPSFs[scaleIndex].data(), psf, _width*_height*sizeof(double));
			
			Logger::Info << "- Scale " << round(scaleEntry.scale) << ", bias factor=" << round(scaleEntry.biasFactor*10.0)/10.0 << ", psfpeak=" << scaleEntry.psfPeak << ", gain=" << scaleEntry.gain << ", kernel peak=" << scaleEntry.kernelPeak << '\n';
		}
		else {
			if(scaleEntry.scale != 0.0)
				msTransforms.Transform(convolvedPSFs[scaleIndex].data(), tmp, scaleEntry.scale);
		}
	}
}

void MultiScaleAlgorithm::findActiveScaleConvolvedMaxima(const DynamicSet& imageSet, double* integratedScratch, bool reportRMS)
{
	MultiScaleTransforms msTransforms(_width, _height);
	//ImageBufferAllocator::Ptr convolvedImage;
	//_allocator.Allocate(_width*_height, convolvedImage);
	imageSet.GetLinearIntegrated(integratedScratch);
	ao::uvector<double> transformScales;
	ao::uvector<size_t> transformIndices;
	std::vector<ao::uvector<bool>> transformScaleMasks;
	for(size_t scaleIndex=0; scaleIndex!=_scaleInfos.size(); ++scaleIndex)
	{
		ScaleInfo& scaleEntry = _scaleInfos[scaleIndex];
		if(scaleEntry.isActive)
		{
			if(scaleEntry.scale == 0)
			{
				// Don't convolve scale 0: this is the delta function scale
				findSingleScaleMaximum(integratedScratch, scaleIndex);
				scaleEntry.maxImageValue = findPeak(integratedScratch, scaleEntry.maxImageValueX, scaleEntry.maxImageValueY, scaleIndex);
				if(reportRMS)
					scaleEntry.rms = ThreadedDeconvolutionTools::RMS(integratedScratch, _width*_height);
			} else {
				transformScales.push_back(scaleEntry.scale);
				transformIndices.push_back(scaleIndex);
				if(_usePerScaleMasks)
					transformScaleMasks.push_back(_scaleMasks[scaleIndex]);
			}
		}
	}
	std::vector<ThreadedDeconvolutionTools::PeakData> results;
	
	_tools->FindMultiScalePeak(&msTransforms, &_allocator, integratedScratch, transformScales, results, _allowNegativeComponents, _cleanMask, transformScaleMasks, _cleanBorderRatio, reportRMS);
	
	for(size_t i=0; i!=results.size(); ++i)
	{
		ScaleInfo& scaleEntry = _scaleInfos[transformIndices[i]];
		scaleEntry.maxImageValue = results[i].value;
		scaleEntry.maxImageValueX = results[i].x;
		scaleEntry.maxImageValueY = results[i].y;
		if(reportRMS)
			scaleEntry.rms = results[i].rms;
	}
	if(reportRMS)
	{
		Logger::Info << "RMS per scale: {";
		for(size_t scaleIndex=0; scaleIndex!=_scaleInfos.size(); ++scaleIndex)
		{
			ScaleInfo& scaleEntry = _scaleInfos[scaleIndex];
			//double rmsBias = _scaleInfos[0].rms / scaleEntry.rms;
			if(scaleIndex != 0)
				Logger::Info << ", ";
			Logger::Info << round(scaleEntry.scale) << ": " << scaleEntry.rms << " Jy";
			// This can be made an option later:
			// scaleEntry.biasFactor = rmsBias;
			// However, at large scales the RMS is not a good estimator of the significance, because
			// at large scales also the signal looks noise like, and increases thereby the RMS.
		}
		Logger::Info << "}\n";
	}
}

void MultiScaleAlgorithm::findSingleScaleMaximum(const double* convolvedImage, size_t scaleIndex)
{
	ScaleInfo& scaleEntry = _scaleInfos[scaleIndex];
	scaleEntry.maxImageValue = findPeak(convolvedImage, scaleEntry.maxImageValueX, scaleEntry.maxImageValueY, scaleIndex);
}

void MultiScaleAlgorithm::sortScalesOnMaxima(size_t& scaleWithPeak)
{
	// Find max component
	std::map<double,size_t> peakToScaleMap;
	for(size_t i=0; i!=_scaleInfos.size(); ++i)
	{
		if(_scaleInfos[i].isActive)
			peakToScaleMap.insert(std::make_pair(std::fabs(_scaleInfos[i].maxImageValue * _scaleInfos[i].biasFactor), i));
	}
	std::map<double,size_t>::const_reverse_iterator mapIter = peakToScaleMap.rbegin();
	scaleWithPeak = mapIter->second;
	//++mapIter;
	//size_t runnerUp;
	//if(mapIter != peakToScaleMap.rend())
	//	runnerUp = mapIter->second;
	//else
	//	runnerUp = scaleWithPeak;
}

void MultiScaleAlgorithm::activateScales(size_t scaleWithLastPeak)
{
	for(size_t i=0; i!=_scaleInfos.size(); ++i)
	{
		bool doActivate = i == scaleWithLastPeak || /*i == runnerUp ||*/ std::fabs(_scaleInfos[i].maxImageValue) * _scaleInfos[i].biasFactor > std::fabs(_scaleInfos[scaleWithLastPeak].maxImageValue) * (1.0-_gain) * _scaleInfos[scaleWithLastPeak].biasFactor;
		if(!_scaleInfos[i].isActive && doActivate)
		{
			Logger::Debug << "Scale " << _scaleInfos[i].scale << " is now significant and is activated.\n";
			_scaleInfos[i].isActive = true;
		}
		else if(_scaleInfos[i].isActive && !doActivate) {
			Logger::Debug << "Scale " << _scaleInfos[i].scale << " is insignificant and is deactivated.\n";
			_scaleInfos[i].isActive = false;
		}
	}
}

void MultiScaleAlgorithm::measureComponentValues(ao::uvector<double>& componentValues, size_t scaleIndex, DynamicSet& imageSet)
{
	const ScaleInfo& scale = _scaleInfos[scaleIndex];
	componentValues.resize(imageSet.size());
	Logger::Debug << "Measuring " << scale.maxImageValueX << ',' << scale.maxImageValueY << ", scale " << scale.scale << ", integrated=" << scale.maxImageValue << ":";
	for(size_t i=0; i!=imageSet.size(); ++i)
	{
		componentValues[i] = imageSet[i][scale.maxImageValueX + scale.maxImageValueY*_width];
		Logger::Debug << ' ' << componentValues[i];
	}
	Logger::Debug << '\n';
}

void MultiScaleAlgorithm::addComponentToModel(double* model, size_t scaleWithPeak, double componentValue)
{
	double componentGain = componentValue * _scaleInfos[scaleWithPeak].gain;
	const size_t
		x = _scaleInfos[scaleWithPeak].maxImageValueX,
		y = _scaleInfos[scaleWithPeak].maxImageValueY;
	if(_scaleInfos[scaleWithPeak].scale == 0.0)
		model[x + _width*y] += componentGain;
	else
		MultiScaleTransforms::AddShapeComponent(model, _width, _height, _scaleInfos[scaleWithPeak].scale, x, y, componentGain);
	
	_scaleInfos[scaleWithPeak].nComponentsCleaned++;
	_scaleInfos[scaleWithPeak].totalFluxCleaned += componentGain;
	
	if(_trackPerScaleMasks)
	{
		_scaleMasks[scaleWithPeak][x + _width*y] = true;
	}
}

double* MultiScaleAlgorithm::getConvolvedPSF(size_t psfIndex, size_t scaleIndex, const ao::uvector<const double*>& psfs, double* scratch,const std::unique_ptr<std::unique_ptr<ImageBufferAllocator::Ptr[]>[]>& convolvedPSFs)
{
	return convolvedPSFs[psfIndex][scaleIndex].data();
}

double MultiScaleAlgorithm::findPeak(const double* image, size_t& x, size_t& y, size_t scaleIndex)
{
	if(_usePerScaleMasks)
		return SimpleClean::FindPeak(image, _width, _height, x, y, _allowNegativeComponents, 0, _height, _scaleMasks[scaleIndex].data(), _cleanBorderRatio);
	else if(_cleanMask == 0)
		return SimpleClean::FindPeak(image, _width, _height, x, y, _allowNegativeComponents, 0, _height, _cleanBorderRatio);
	else
		return SimpleClean::FindPeak(image, _width, _height, x, y, _allowNegativeComponents, 0, _height, _cleanMask, _cleanBorderRatio);
}
