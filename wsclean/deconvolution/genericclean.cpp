#include "genericclean.h"
#include "../multiscale/clarkloop.h"

#include "../lane.h"

#include "../multiscale/threadeddeconvolutiontools.h"

#include <boost/thread/thread.hpp>

GenericClean::GenericClean(ImageBufferAllocator& allocator, bool clarkOptimization) :
	_convolutionPadding(1.1),
	_useClarkOptimization(clarkOptimization),
	_allocator(allocator)
{
}

void GenericClean::ExecuteMajorIteration(ImageSet& dirtySet, ImageSet& modelSet, const ao::uvector<const double*>& psfs, size_t width, size_t height, bool& reachedMajorThreshold)
{
	Logger::Debug << "Running dynamically-joined clean algorithm.\n";
	
	const size_t iterationCounterAtStart = _iterationNumber;
	if(_stopOnNegativeComponent)
		_allowNegativeComponents = true;
	_width = width;
	_height = height;
	_convolutionWidth = ceil(_convolutionPadding * _width);
	_convolutionHeight = ceil(_convolutionPadding * _height);
	if(_convolutionWidth%2 != 0)
		++_convolutionWidth;
	if(_convolutionHeight%2 != 0)
		++_convolutionHeight;
	
	ImageBufferAllocator::Ptr integrated, scratchA, scratchB;
	_allocator.Allocate(width*height, integrated);
	_allocator.Allocate(_convolutionWidth*_convolutionHeight, scratchA);
	_allocator.Allocate(_convolutionWidth*_convolutionHeight, scratchB);
	dirtySet.GetSquareIntegrated(integrated.data(), scratchA.data());
	size_t componentX=0, componentY=0;
	double maxValue = findPeak(integrated.data(), scratchA.data(), componentX, componentY);
	Logger::Info << "Initial peak: " << peakDescription(integrated.data(), componentX, componentY) << '\n';
	double firstThreshold = this->_threshold;
	double stopGainThreshold = std::fabs(maxValue)*(1.0-this->_mGain);
	if(stopGainThreshold > firstThreshold)
	{
		firstThreshold = stopGainThreshold;
		Logger::Info << "Next major iteration at: " << stopGainThreshold << '\n';
	}
	else if(this->_mGain != 1.0) {
		Logger::Info << "Major iteration threshold reached global threshold of " << this->_threshold << ": final major iteration.\n";
	}
		
	if(_useClarkOptimization)
	{
		size_t startIteration = _iterationNumber;
		ClarkLoop clarkLoop(_width, _height, _convolutionWidth, _convolutionHeight);
		clarkLoop.SetIterationInfo(_iterationNumber, MaxNIter());
		clarkLoop.SetThreshold(firstThreshold);
		clarkLoop.SetGain(Gain());
		clarkLoop.SetAllowNegativeComponents(AllowNegativeComponents());
		clarkLoop.SetSpectralFitter(&Fitter());
		if(!_rmsFactorImage.empty())
			clarkLoop.SetRMSFactorImage(_rmsFactorImage);
		if(_cleanMask)
			clarkLoop.SetMask(_cleanMask);
		const size_t
			horBorderSize = round(_width * CleanBorderRatio()),
			vertBorderSize = round(_height * CleanBorderRatio());
		clarkLoop.SetCleanBorders(horBorderSize, vertBorderSize);
		
		maxValue = clarkLoop.Run(dirtySet, psfs);
		
		_iterationNumber = clarkLoop.CurrentIteration();
		
		Logger::Info << "Performed " << (_iterationNumber - startIteration) << " / " << _iterationNumber << " iterations with Clark optimization in this major iteration.\n";
		
		for(size_t imageIndex=0; imageIndex!=dirtySet.size(); ++imageIndex)
		{
			// TODO this can be multi-threaded if each thread has its own temporaries
			const double *psf = psfs[dirtySet.PSFIndex(imageIndex)];
			clarkLoop.CorrectResidualDirty(scratchA.data(), scratchB.data(), integrated.data(), imageIndex, dirtySet[imageIndex],  psf);
			
			clarkLoop.GetFullIndividualModel(imageIndex, scratchA.data());
			double* model = modelSet[imageIndex];
			for(size_t i=0; i!=_width*_height; ++i)
				model[i] += scratchA.data()[i];
		}
	}
	else {
		ThreadedDeconvolutionTools tools(_threadCount);
		size_t peakIndex = componentX + componentY*_width;

		ao::uvector<double> peakValues(dirtySet.size());
		
		while(fabs(maxValue) > firstThreshold && this->_iterationNumber < this->_maxIter && !(maxValue<0.0 && this->_stopOnNegativeComponent))
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
			
			dirtySet.GetSquareIntegrated(integrated.data(), scratchA.data());
			maxValue = findPeak(integrated.data(), scratchA.data(), componentX, componentY);
			
			peakIndex = componentX + componentY*_width;
			
			++this->_iterationNumber;
		}
	}
	Logger::Info << "Stopped on peak " << maxValue << '\n';
	reachedMajorThreshold = std::fabs(maxValue) <= stopGainThreshold && (maxValue != 0.0) && (_iterationNumber-iterationCounterAtStart)!=0;
}

std::string GenericClean::peakDescription(const double* image, size_t& x, size_t& y)
{
	std::ostringstream str;
	size_t index = x + y*_width;
	double peak = image[index];
	str << peak << " Jy at " << x << "," << y;
	return str.str();
}

double GenericClean::findPeak(const double* image, double* scratch, size_t& x, size_t& y)
{
	if(_rmsFactorImage.empty())
	{
		if(_cleanMask == 0)
			return SimpleClean::FindPeak(image, _width, _height, x, y, _allowNegativeComponents, 0, _height, _cleanBorderRatio);
		else
			return SimpleClean::FindPeakWithMask(image, _width, _height, x, y, _allowNegativeComponents, 0, _height, _cleanMask, _cleanBorderRatio);
	}
	else {
		for(size_t i=0; i!=_width*_height; ++i)
		{
			scratch[i] = image[i] * _rmsFactorImage[i];
		}
		double maxValue;
		if(_cleanMask == 0)
			maxValue = SimpleClean::FindPeak(scratch, _width, _height, x, y, _allowNegativeComponents, 0, _height, _cleanBorderRatio);
		else
			maxValue = SimpleClean::FindPeakWithMask(scratch, _width, _height, x, y, _allowNegativeComponents, 0, _height, _cleanMask, _cleanBorderRatio);
		return maxValue / _rmsFactorImage[x + y*_width];
	}
}
