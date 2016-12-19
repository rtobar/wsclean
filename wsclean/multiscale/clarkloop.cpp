#include "clarkloop.h"

#include "../deconvolution/spectralfitter.h"

#include "../fftconvolver.h"

#include "../wsclean/logger.h"

template<bool AllowNegatives>
size_t ClarkModel::GetMaxComponent(double* scratch, double& maxValue) const
{
	_residual->GetLinearIntegrated(scratch);
	size_t maxComponent = 0;
	maxValue = scratch[0];
	for(size_t i=0; i!=size(); ++i)
	{
		double value;
		if(AllowNegatives)
			value = std::fabs(scratch[i]);
		else
			value = scratch[i];
		if(value > maxValue)
		{
			maxComponent = i;
			maxValue = value;
		}
	}
	return maxComponent;
}

void ClarkLoop::Run(DynamicSet& convolvedResidual, const std::vector<const double*>& doubleConvolvedPsfs)
{
	_clarkModel = ClarkModel(_width, _height);
	
	findPeakPositions(convolvedResidual);
	
	_clarkModel.MakeSets(convolvedResidual);
	Logger::Debug << "Number of components selected > " << _threshold << ": " << _clarkModel.size() << '\n';
	
	ao::uvector<double> scratch(_clarkModel.size());
	double maxValue;
	size_t maxComponent = _clarkModel.GetMaxComponent(scratch.data(), maxValue, _allowNegativeComponents);
		
	while(std::fabs(maxValue) > _threshold && _currentIteration < _maxIterations)
	{
		ao::uvector<double> componentValues(_clarkModel.Residual().size());
		for(size_t imgIndex=0; imgIndex!=_clarkModel.Residual().size(); ++imgIndex)
			componentValues[imgIndex] = _clarkModel.Residual()[imgIndex][maxComponent] * _gain;
		_fluxCleaned += maxValue * _gain;
		
		if(_fitter)
			_fitter->FitAndEvaluate(componentValues.data());
			
		for(size_t imgIndex=0; imgIndex!=_clarkModel.Model().size(); ++imgIndex)
			_clarkModel.Model()[imgIndex][maxComponent] += componentValues[imgIndex];
		
		size_t
			x = _clarkModel.X(maxComponent),
			y = _clarkModel.Y(maxComponent);
		for(size_t imgIndex=0; imgIndex!=_clarkModel.Residual().size(); ++imgIndex)
		{
			double* image = _clarkModel.Residual()[imgIndex];
			const double* psf = doubleConvolvedPsfs[_clarkModel.Residual().PSFIndex(imgIndex)];
			double psfFactor = componentValues[imgIndex];
			for(size_t px=0; px!=_clarkModel.size(); ++px)
			{
				int psfX = _clarkModel.X(px) - x + _width/2;
				int psfY = _clarkModel.Y(px) - y + _height/2;
				if(psfX >= 0 && psfX < int(_width) && psfY >= 0 && psfY < int(_height))
					image[px] -= psf[psfX + psfY*_height] * psfFactor;
			}
		}
		
		maxComponent = _clarkModel.GetMaxComponent(scratch.data(), maxValue, _allowNegativeComponents);
		++_currentIteration;
	}
}

void ClarkModel::MakeSets(const DynamicSet& residualSet)
{
	_residual.reset(new DynamicSet(&residualSet.Table(), residualSet.Allocator(), residualSet.ChannelsInDeconvolution(), residualSet.SquareJoinedChannels(), size(), 1));
	_model.reset(new DynamicSet(&residualSet.Table(), residualSet.Allocator(), residualSet.ChannelsInDeconvolution(), residualSet.SquareJoinedChannels(), size(), 1));
	for(size_t imgIndex=0; imgIndex!=_model->size(); ++imgIndex)
	{
		std::fill((*_model)[imgIndex], (*_model)[imgIndex]+size(), 0.0);
		
		const double* sourceResidual = residualSet[imgIndex];
		double* destResidual = (*_residual)[imgIndex];
		for(size_t pxIndex=0; pxIndex!=size(); ++pxIndex)
		{
			size_t srcIndex = _positions[pxIndex].second*_width + _positions[pxIndex].first;
			destResidual[pxIndex] = sourceResidual[srcIndex];
		}
	}
}

void ClarkLoop::findPeakPositions(DynamicSet& convolvedResidual)
{
	ImageBufferAllocator::Ptr integratedScratch;
	convolvedResidual.Allocator().Allocate(_width * _height, integratedScratch);
	convolvedResidual.GetLinearIntegrated(integratedScratch.data());
	
	double* imagePtr = integratedScratch.data();
	if(_mask)
	{
		const bool* maskPtr = _mask;
		for(size_t y=0; y!=_height; ++y)
		{
			for(size_t x=0; x!=_width; ++x)
			{
				double value;
				if(_allowNegativeComponents)
					value = fabs(*imagePtr);
				else
					value = *imagePtr;
				if(value >= _threshold && *maskPtr)
					_clarkModel.AddPosition(x, y);
				++imagePtr;
				++maskPtr;
			}
		}
	}
	else {
		for(size_t y=0; y!=_height; ++y)
		{
			for(size_t x=0; x!=_width; ++x)
			{
				double value;
				if(_allowNegativeComponents)
					value = fabs(*imagePtr);
				else
					value = *imagePtr;
				if(value >= _threshold)
					_clarkModel.AddPosition(x, y);
				++imagePtr;
			}
		}
	}
}

void ClarkLoop::GetFullIndividualModel(size_t imageIndex, double* individualModelImg) const
{
	std::fill(individualModelImg, individualModelImg + _width*_height, 0.0);
	const double* data = _clarkModel.Model()[imageIndex];
	for(size_t px=0; px!=_clarkModel.size(); ++px)
	{
		individualModelImg[_clarkModel.FullIndex(px)] = data[px];
	}
}

void ClarkLoop::CorrectResidualDirty(double* scratchA, double* scratchB, size_t imageIndex, double* residual, const double* singleConvolvedPsf) const
{
	GetFullIndividualModel(imageIndex, scratchA);
	
	FFTConvolver::PrepareKernel(scratchB, singleConvolvedPsf, _width, _height);
	FFTConvolver::ConvolveSameSize(scratchA, scratchB, _width, _height);
	
	for(size_t i=0; i!=_width*_height; ++i)
		residual[i] -= scratchA[i];
}

void ClarkLoop::UpdateAutoMask(bool* mask) const
{
	for(size_t px=0; px!=_clarkModel.size(); ++px)
	{
		if(_clarkModel.Model()[0][px] != 0.0)
			mask[_clarkModel.FullIndex(px)] = true;
	}
}
