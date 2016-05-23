#include "dynamicjoinedclean.h"

#include "../lane.h"

#include "../multiscale/threadeddeconvolutiontools.h"

#include <boost/thread/thread.hpp>

DynamicJoinedClean::DynamicJoinedClean(ImageBufferAllocator& allocator) :
	_allocator(allocator)
{ }

void DynamicJoinedClean::ExecuteMajorIteration(DynamicSet& dirtySet, DynamicSet& modelSet, const ao::uvector<const double*>& psfs, size_t width, size_t height, bool& reachedMajorThreshold)
{
	Logger::Debug << "Running dynamically-joined clean algorithm.\n";
	
	ThreadedDeconvolutionTools tools(_threadCount);
	
	if(_stopOnNegativeComponent)
		_allowNegativeComponents = true;
	_width = width;
	_height = height;
	
	ImageBufferAllocator::Ptr integrated, scratch;
	_allocator.Allocate(width*height, integrated);
	_allocator.Allocate(width*height, scratch);
	dirtySet.GetSquareIntegrated(integrated.data(), scratch.data());
	size_t componentX=0, componentY=0;
	double peakNormalized = findPeak(integrated.data(), componentX, componentY);
	Logger::Info << "Initial peak: " << peakDescription(integrated.data(), componentX, componentY) << '\n';
	
	size_t peakIndex = componentX + componentY*_width;
	double firstThreshold = this->_threshold;
	double stopGainThreshold = std::fabs(peakNormalized)*(1.0-this->_mGain);
	if(stopGainThreshold > firstThreshold)
	{
		firstThreshold = stopGainThreshold;
		Logger::Info << "Next major iteration at: " << stopGainThreshold << '\n';
	}
	else if(this->_mGain != 1.0) {
		Logger::Info << "Major iteration threshold reached global threshold of " << this->_threshold << ": final major iteration.\n";
	}

	ao::uvector<double> peakValues(dirtySet.size());
	
	while(fabs(peakNormalized) > firstThreshold && this->_iterationNumber < this->_maxIter && !(peakNormalized<0.0 && this->_stopOnNegativeComponent))
	{
		if(this->_iterationNumber <= 10 ||
			(this->_iterationNumber <= 100 && this->_iterationNumber % 10 == 0) ||
			(this->_iterationNumber <= 1000 && this->_iterationNumber % 100 == 0) ||
			this->_iterationNumber % 1000 == 0)
			Logger::Info << "Iteration " << this->_iterationNumber << ": " << peakDescription(integrated.data(), componentX, componentY) << '\n';
		
		for(size_t i=0; i!=dirtySet.size(); ++i)
			peakValues[i] = dirtySet[i][peakIndex];
		
		PerformSpectralFit(peakValues.data());
		
		for(size_t i=0; i!=dirtySet.size(); ++i)
		{
			peakValues[i] *= this->_gain;
			modelSet[i][peakIndex] += peakValues[i];
			
			size_t psfIndex = dirtySet.PSFIndex(i);
			
			tools.SubtractImage(dirtySet[i], psfs[psfIndex], width, height, componentX, componentY, peakValues[i]);
		}
		
		dirtySet.GetSquareIntegrated(integrated.data(), scratch.data());
		peakNormalized = findPeak(integrated.data(), componentX, componentY);
		
		peakIndex = componentX + componentY*_width;
		
		++this->_iterationNumber;
	}
	Logger::Info << "Stopped on peak " << peakNormalized << '\n';
	reachedMajorThreshold = std::fabs(peakNormalized) <= stopGainThreshold && (peakNormalized != 0.0);
}

std::string DynamicJoinedClean::peakDescription(const double* image, size_t& x, size_t& y)
{
	std::ostringstream str;
	size_t index = x + y*_width;
	double peak = image[index];
	str << peak << " Jy at " << x << "," << y;
	return str.str();
}

double DynamicJoinedClean::findPeak(const double* image, size_t& x, size_t& y)
{
	if(_cleanMask == 0)
		return SimpleClean::FindPeak(image, _width, _height, x, y, _allowNegativeComponents, 0, _height, _cleanBorderRatio);
	else
		return SimpleClean::FindPeak(image, _width, _height, x, y, _allowNegativeComponents, 0, _height, _cleanMask, _cleanBorderRatio);
}
