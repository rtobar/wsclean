#include "multiscalealgorithm.h"

#include "multiscaletransforms.h"

#include "../deconvolution/simpleclean.h"

MultiScaleAlgorithm::MultiScaleAlgorithm(ImageBufferAllocator& allocator, double beamSize, double pixelScaleX, double pixelScaleY) :
	_allocator(allocator),
	_width(0),
	_height(0),
	_beamSizeInPixels(beamSize / std::max(pixelScaleX, pixelScaleY)),
	_verbose(false)
{
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
	//   - Subtract PSF from individual & individually convolved images
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
	findActiveScaleConvolvedMaxima(dirtySet, integratedScratch.data());
	sortScalesOnMaxima(scaleWithPeak);
	
	double mGainThreshold = std::fabs(_scaleInfos[scaleWithPeak].maxImageValue * _scaleInfos[scaleWithPeak].factor) * (1.0 - _mGain);
	const double firstThreshold = std::max(_threshold, mGainThreshold);
	
	std::cout << "Starting multi-scale cleaning. Start peak="
		<< _scaleInfos[scaleWithPeak].maxImageValue * _scaleInfos[scaleWithPeak].factor
		<< " Jy, major iteration threshold=" << firstThreshold << "\n";
	
	std::unique_ptr<ImageBufferAllocator::Ptr[]> doubleConvolvedPSFs(
		new ImageBufferAllocator::Ptr[dirtySet.PSFCount()]);
	for(size_t i=0; i!=dirtySet.PSFCount(); ++i)
	{
		_allocator.Allocate(_width*_height, doubleConvolvedPSFs[i]);
	}
	
	DynamicSet individualConvolvedImages(&dirtySet.Table(), dirtySet.Allocator(), dirtySet.ChannelsInDeconvolution(), _width, _height);
	
	//
	// The minor iteration loop
	//
	while(_iterationNumber < MaxNIter() && std::fabs(_scaleInfos[scaleWithPeak].maxImageValue * _scaleInfos[scaleWithPeak].factor) > firstThreshold)
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
			std::fabs(_scaleInfos[scaleWithPeak].maxImageValue * _scaleInfos[scaleWithPeak].factor) * (1.0 - _gain),
			firstThreshold);
		while(_iterationNumber < MaxNIter() && std::fabs(_scaleInfos[scaleWithPeak].maxImageValue * _scaleInfos[scaleWithPeak].factor) > firstSubIterationThreshold)
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
				
				// Adjust model
				addComponentToModel(modelSet[imgIndex], scaleWithPeak, componentValues[imgIndex]);
			}
			
			// Find maximum for this scale
			individualConvolvedImages.GetLinearIntegrated(integratedScratch.data());
			findSingleScaleMaximum(integratedScratch.data(), scaleWithPeak);
			
			++_iterationNumber;
		}
		
		activateScales(scaleWithPeak);
		
		findActiveScaleConvolvedMaxima(dirtySet, integratedScratch.data());
		sortScalesOnMaxima(scaleWithPeak);
		
		std::cout << "Iteration " << _iterationNumber << ", scale " << round(_scaleInfos[scaleWithPeak].scale) << " px : " << _scaleInfos[scaleWithPeak].maxImageValue*_scaleInfos[scaleWithPeak].factor << " Jy at " << _scaleInfos[scaleWithPeak].maxImageValueX << ',' << _scaleInfos[scaleWithPeak].maxImageValueY << '\n';
	}
	
	reachedMajorThreshold =
		_iterationNumber < MaxNIter() && std::fabs(_scaleInfos[scaleWithPeak].maxImageValue * _scaleInfos[scaleWithPeak].factor) > _threshold;
}

void MultiScaleAlgorithm::initializeScaleInfo()
{
	if(_scaleInfos.empty())
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
}

void MultiScaleAlgorithm::convolvePSFs(std::unique_ptr<ImageBufferAllocator::Ptr[]>& convolvedPSFs, const double* psf, double* tmp, bool isIntegrated)
{
	MultiScaleTransforms msTransforms(_width, _height);
	convolvedPSFs.reset(new ImageBufferAllocator::Ptr[_scaleInfos.size()]);
	if(isIntegrated)
		std::cout << "Scale info:\n";
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
			//scaleEntry.factor = std::max(1.0,
			//	scaleEntry.psfPeak * scaleInfos[0].kernelPeak /
			//	(scaleEntry.kernelPeak * scaleInfos[0].psfPeak));
			scaleEntry.factor = pow(_multiscaleScaleBias, -double(scaleIndex));
			
			// I tried this, but wasn't perfect:
			// _gain * _scaleInfos[0].kernelPeak / scaleEntry.kernelPeak;
			scaleEntry.gain = _gain * _scaleInfos[0].psfPeak / scaleEntry.psfPeak;
			
			scaleEntry.isActive = true;
			
			if(scaleIndex == 0)
				memcpy(convolvedPSFs[scaleIndex].data(), psf, _width*_height*sizeof(double));
			
			std::cout << "- Scale " << round(scaleEntry.scale) << ", bias factor=" << round(scaleEntry.factor*10.0)/10.0 << ", response=" << scaleEntry.psfPeak << ", gain=" << scaleEntry.gain << ", kernel peak=" << scaleEntry.kernelPeak << '\n';
		}
		else {
			if(scaleIndex != 0)
				msTransforms.Transform(convolvedPSFs[scaleIndex].data(), tmp, scaleEntry.scale);
		}
	}
}

void MultiScaleAlgorithm::findActiveScaleConvolvedMaxima(const DynamicSet& imageSet, double* integratedScratch)
{
	MultiScaleTransforms msTransforms(_width, _height);
	//ImageBufferAllocator::Ptr convolvedImage;
	//_allocator.Allocate(_width*_height, convolvedImage);
	imageSet.GetLinearIntegrated(integratedScratch);
	ao::uvector<double> transformScales;
	ao::uvector<size_t> transformIndices;
	for(size_t scaleIndex=0; scaleIndex!=_scaleInfos.size(); ++scaleIndex)
	{
		ScaleInfo& scaleEntry = _scaleInfos[scaleIndex];
		if(scaleEntry.isActive)
		{
			//memcpy(convolvedImage.data(), integratedScratch, _width*_height*sizeof(double));
				
			if(scaleIndex == 0)
			{
				// Don't convolve scale 0: this is the delta function scale
				findSingleScaleMaximum(integratedScratch, scaleIndex);
				scaleEntry.maxImageValue = findPeak(integratedScratch, scaleEntry.maxImageValueX, scaleEntry.maxImageValueY);
			} else {
				transformScales.push_back(scaleEntry.scale);
				transformIndices.push_back(scaleIndex);
				//msTransforms.Transform(convolvedImage.data(), scratch, scaleEntry.scale);
			}
		}
	}
	std::vector<ThreadedDeconvolutionTools::PeakData> results;
	_tools->FindMultiScalePeak(&msTransforms, &_allocator, integratedScratch, transformScales, results, _allowNegativeComponents, _cleanMask, _cleanBorderRatio);
	
	for(size_t i=0; i!=results.size(); ++i)
	{
		ScaleInfo& scaleEntry = _scaleInfos[transformIndices[i]];
		scaleEntry.maxImageValue = results[i].value;
		scaleEntry.maxImageValueX = results[i].x;
		scaleEntry.maxImageValueY = results[i].y;
	}
}

void MultiScaleAlgorithm::findSingleScaleMaximum(const double* convolvedImage, size_t scaleIndex)
{
	ScaleInfo& scaleEntry = _scaleInfos[scaleIndex];
	scaleEntry.maxImageValue = findPeak(convolvedImage, scaleEntry.maxImageValueX, scaleEntry.maxImageValueY);
}

void MultiScaleAlgorithm::sortScalesOnMaxima(size_t& scaleWithPeak)
{
	// Find max component
	std::map<double,size_t> peakToScaleMap;
	for(size_t i=0; i!=_scaleInfos.size(); ++i)
	{
		if(_scaleInfos[i].isActive)
			peakToScaleMap.insert(std::make_pair(std::fabs(_scaleInfos[i].maxImageValue * _scaleInfos[i].factor), i));
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
		bool doActivate = i == scaleWithLastPeak || /*i == runnerUp ||*/ std::fabs(_scaleInfos[i].maxImageValue) * _scaleInfos[i].factor > std::fabs(_scaleInfos[scaleWithLastPeak].maxImageValue) * (1.0-_gain) * _scaleInfos[scaleWithLastPeak].factor;
		if(!_scaleInfos[i].isActive && doActivate)
		{
			if(_verbose)
				std::cout << "Scale " << _scaleInfos[i].scale << " is now significant and is activated.\n";
			_scaleInfos[i].isActive = true;
		}
		else if(_scaleInfos[i].isActive && !doActivate) {
			if(_verbose)
				std::cout << "Scale " << _scaleInfos[i].scale << " is insignificant and is deactivated.\n";
			_scaleInfos[i].isActive = false;
		}
	}
}

void MultiScaleAlgorithm::measureComponentValues(ao::uvector<double>& componentValues, size_t scaleIndex, DynamicSet& imageSet)
{
	const ScaleInfo& scale = _scaleInfos[scaleIndex];
	componentValues.resize(imageSet.size());
	if(_verbose)
		std::cout << "Measuring " << scale.maxImageValueX << ',' << scale.maxImageValueY << ", scale " << scale.scale << ", integrated=" << scale.maxImageValue << ":";
	for(size_t i=0; i!=imageSet.size(); ++i)
	{
		componentValues[i] = imageSet[i][scale.maxImageValueX + scale.maxImageValueY*_width];
		if(_verbose)
			std::cout << ' ' << componentValues[i];
	}
	if(_verbose)
		std::cout << '\n';
}

void MultiScaleAlgorithm::addComponentToModel(double* model, size_t scaleWithPeak, double componentValue)
{
	double componentGain = componentValue * _scaleInfos[scaleWithPeak].gain;
	if(scaleWithPeak == 0)
		model[_scaleInfos[scaleWithPeak].maxImageValueX + _width*_scaleInfos[scaleWithPeak].maxImageValueY]
			+= componentGain;
	else
		MultiScaleTransforms::AddShapeComponent(model, _width, _height, _scaleInfos[scaleWithPeak].scale, _scaleInfos[scaleWithPeak].maxImageValueX, _scaleInfos[scaleWithPeak].maxImageValueY, componentGain);
}

double* MultiScaleAlgorithm::getConvolvedPSF(size_t psfIndex, size_t scaleIndex, const ao::uvector<const double*>& psfs, double* scratch,const std::unique_ptr<std::unique_ptr<ImageBufferAllocator::Ptr[]>[]>& convolvedPSFs)
{
	return convolvedPSFs[psfIndex][scaleIndex].data();
}

double MultiScaleAlgorithm::findPeak(const double* image, size_t& x, size_t& y)
{
	if(_cleanMask == 0)
		return SimpleClean::FindPeak(image, _width, _height, x, y, _allowNegativeComponents, 0, _height, _cleanBorderRatio);
	else
		return SimpleClean::FindPeak(image, _width, _height, x, y, _allowNegativeComponents, 0, _height, _cleanMask, _cleanBorderRatio);
}
