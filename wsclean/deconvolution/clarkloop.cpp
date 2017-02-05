#include "clarkloop.h"

#include "../deconvolution/spectralfitter.h"
#include "../deconvolution/componentlist.h"

#include "../fftconvolver.h"
#include "../image.h"

#include "../wsclean/logger.h"

template<bool AllowNegatives>
size_t ClarkModel::GetMaxComponent(double* scratch, double& maxValue) const
{
	_residual->GetLinearIntegrated(scratch);
	if(!_rmsFactorImage.empty())
	{
		for(size_t i=0; i!=size(); ++i)
			scratch[i] *= _rmsFactorImage[i];
	}
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

double ClarkLoop::Run(ImageSet& convolvedResidual, const ao::uvector<const double*>& doubleConvolvedPsfs)
{
	_clarkModel = ClarkModel(_width, _height);
	
	findPeakPositions(convolvedResidual);
	
	_clarkModel.MakeSets(convolvedResidual);
	if(!_rmsFactorImage.empty())
		_clarkModel.MakeRMSFactorImage(_rmsFactorImage);
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
					image[px] -= psf[psfX + psfY*_width] * psfFactor;
			}
		}
		
		maxComponent = _clarkModel.GetMaxComponent(scratch.data(), maxValue, _allowNegativeComponents);
		++_currentIteration;
	}
	return maxValue;
}

void ClarkModel::MakeSets(const ImageSet& residualSet)
{
	_residual.reset(new ImageSet(&residualSet.Table(), residualSet.Allocator(), residualSet.ChannelsInDeconvolution(), residualSet.SquareJoinedChannels(), size(), 1));
	_model.reset(new ImageSet(&residualSet.Table(), residualSet.Allocator(), residualSet.ChannelsInDeconvolution(), residualSet.SquareJoinedChannels(), size(), 1));
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

void ClarkModel::MakeRMSFactorImage(Image& rmsFactorImage)
{
	_rmsFactorImage = Image(size(), 1, _residual->Allocator());
	for(size_t pxIndex=0; pxIndex!=size(); ++pxIndex)
	{
		size_t srcIndex = _positions[pxIndex].second*_width + _positions[pxIndex].first;
		_rmsFactorImage[pxIndex] = rmsFactorImage[srcIndex];
	}
}

void ClarkLoop::findPeakPositions(ImageSet& convolvedResidual)
{
	Image integratedScratch(_width, _height, convolvedResidual.Allocator());
	convolvedResidual.GetLinearIntegrated(integratedScratch.data());
	
	if(!_rmsFactorImage.empty())
	{
		integratedScratch *= _rmsFactorImage;
	}
	
	const size_t
		xiStart = _horizontalBorder, xiEnd = std::max<long>(xiStart, _width - _horizontalBorder),
		yiStart = _verticalBorder, yiEnd = std::max<long>(yiStart, _height - _verticalBorder);
	
	if(_mask)
	{
		for(size_t y=yiStart; y!=yiEnd; ++y)
		{
			const bool* maskPtr = _mask + y*_width;
			double* imagePtr = integratedScratch.data() + y*_width;
			for(size_t x=xiStart; x!=xiEnd; ++x)
			{
				double value;
				if(_allowNegativeComponents)
					value = fabs(imagePtr[x]);
				else
					value = imagePtr[x];
				if(value >= _threshold && maskPtr[x])
					_clarkModel.AddPosition(x, y);
			}
		}
	}
	else {
		for(size_t y=yiStart; y!=yiEnd; ++y)
		{
			double* imagePtr = integratedScratch.data() + y*_width;
			for(size_t x=xiStart; x!=xiEnd; ++x)
			{
				double value;
				if(_allowNegativeComponents)
					value = fabs(imagePtr[x]);
				else
					value = imagePtr[x];
				if(value >= _threshold)
					_clarkModel.AddPosition(x, y);
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

void ClarkLoop::CorrectResidualDirty(double* scratchA, double* scratchB, double* scratchC, size_t imageIndex, double* residual, const double* singleConvolvedPsf) const
{
	// Get padded kernel in scratchB
	Image::Untrim(scratchA, _untrimmedWidth, _untrimmedHeight, singleConvolvedPsf, _width, _height);
	FFTConvolver::PrepareKernel(scratchB, scratchA, _untrimmedWidth, _untrimmedHeight);
	
	// Get padded model image in scratchA
	GetFullIndividualModel(imageIndex, scratchC);
	Image::Untrim(scratchA, _untrimmedWidth, _untrimmedHeight, scratchC, _width, _height);
	
	// Convolve and store in scratchA
	FFTConvolver::ConvolveSameSize(scratchA, scratchB, _untrimmedWidth, _untrimmedHeight);
	
	//Trim the result into scratchC
	Image::Trim(scratchC, _width, _height, scratchA, _untrimmedWidth, _untrimmedHeight);
	
	for(size_t i=0; i!=_width*_height; ++i)
		residual[i] -= scratchC[i];
}

void ClarkLoop::UpdateAutoMask(bool* mask) const
{
	for(size_t imageIndex=0; imageIndex!=_clarkModel.Model().size(); ++imageIndex)
	{
		const double* image = _clarkModel.Model()[imageIndex];
		for(size_t px=0; px!=_clarkModel.size(); ++px)
		{
			if(image[px] != 0.0)
				mask[_clarkModel.FullIndex(px)] = true;
		}
	}
}

void ClarkLoop::UpdateComponentList(class ComponentList& list, size_t scaleIndex) const
{
	ao::uvector<double> values(_clarkModel.Model().size());
	for(size_t px=0; px!=_clarkModel.size(); ++px)
	{
		bool isNonZero = false;
		for(size_t imageIndex=0; imageIndex!=_clarkModel.Model().size(); ++imageIndex)
		{
			values[imageIndex] = _clarkModel.Model()[imageIndex][px];
			if(values[imageIndex] != 0.0)
				isNonZero = true;
		}
		if(isNonZero)
		{
			size_t posIndex = _clarkModel.FullIndex(px);
			size_t x = posIndex % _width, y = posIndex / _width;
			list.Add(x, y, scaleIndex, values.data());
		}
	}
}
