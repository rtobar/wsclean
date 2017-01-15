#ifndef BINNED_UV_GRID_H
#define BINNED_UV_GRID_H

#include "wstackinggridder.h"

#include "../imageweights.h"

class BinnedUVGrid
{
public:
	BinnedUVGrid(size_t width, size_t height, double pixelScaleX, double pixelScaleY, size_t fftThreadCount, ImageBufferAllocator* allocator, size_t kernelSize=7, size_t overSamplingFactor=63) :
		_fullGridder(new WStackingGridder(width, height, pixelScaleX, pixelScaleY, fftThreadCount, allocator, kernelSize, overSamplingFactor)),
		_imageWeights(WeightMode::UniformWeighted, width, height, pixelScaleX, pixelScaleY)
	{
		_fullGridder->SetGridMode(WStackingGridder::NearestNeighbour);
	}

	/**
	 * Prepare the w-layers of the uv grids. This method should be called before calling
	 * @ref GridSample(). If the w-layers do not fit in the specified amount of memory
	 * in one pass, an exception is thrown.
	 * @param nWLayers Number of uv grids at different w-values, should be >= 1.
	 * @param maxMem Allowed memory in bytes. The gridder will try to set the number
	 * of passes such that this value is not exceeded. Note that this is approximate.
	 * @param minW The smallest w-value to be inverted/predicted.
	 * @param maxW The largest w-value to be inverted/predicted.
	 */
	void PrepareWLayers(size_t nWLayers, double maxMem, double minW, double maxW)
	{
		_maxMem = maxMem;
		_minW = minW;
		_maxW = maxW;
		_fullGridder->PrepareWLayers(nWLayers, maxMem*double(nWLayers)/double(nWLayers+1), minW, maxW);
		_fullGridder->StartInversionPass(0);
		
		if(_fullGridder->NPasses()!=1)
		{
			throw std::runtime_error("The specified amount of memory in BinnedUVGrid::PrepareWLayers() was not enough to hold all w-layers.");
		}
		_totalWeight = 0.0;
	}

	/**
	 * Grid a single sample on the UV grid.
	 */
	void GridSample(const std::complex<float>& value, double weight, double uInLambda, double vInLambda, double wInLambda)
	{
		_fullGridder->AddDataSample(value * float(weight), uInLambda, vInLambda, wInLambda);
		_imageWeights.Grid(uInLambda, vInLambda, weight);
		_totalWeight += weight;
	}
	
	void FinishGridding()
	{
		size_t width = _fullGridder->Width(), height = _fullGridder->Height();
		
		_fullGridder->FinishInversionPass();
		double predictFactor = 1.0 / width / height;
		_fullGridder->FinalizeImage(predictFactor, false);
		
		WStackingGridder *predictingGridder = new WStackingGridder(
			width, height,
			_fullGridder->PixelSizeX(), _fullGridder->PixelSizeY(),
			_fullGridder->NFFTThreads(), _fullGridder->Allocator());
		predictingGridder->SetGridMode(WStackingGridder::NearestNeighbour);
		predictingGridder->PrepareWLayers(1, _maxMem, _minW, _maxW);
		predictingGridder->InitializePrediction(_fullGridder->RealImage());
		std::cout << "Image centre: " << _fullGridder->RealImage()[width/2 + height*width/2]/predictFactor << '\n';
		predictingGridder->StartPredictionPass(0);
		_fullGridder.reset(predictingGridder);
		
		_imageWeights.FinishGridding();
		_imageWeightGrid.resize(_imageWeights.Width() * _imageWeights.Height());
		_imageWeights.GetGrid(_imageWeightGrid.data());
		
		_nonZeroBinCount = countNonZeroBins();
		_currentX = 0;
		_currentY = 0;
	}
	
	/**
	 * Count the number of bins that contain at least one sample. This is different from
	 * grid points that are non-zero, since the anti-alias kernel sets several uv-pixels
	 * around each sample to non-zero.
	 */
	size_t NonZeroBinCount()
	{
		return _nonZeroBinCount;
	}
	
	bool GetNextNonZeroBin(double& uInLambda, double& vInLambda, double& wInLambda,
												 std::complex<double>& value, double& weight)
	{
		size_t
			width = _imageWeights.Width(),
			height = _imageWeights.Height();
		ao::uvector<double> image(width * height);
		const double* ptr = _imageWeightGrid.data() +
			_currentY * width + _currentX;
		
		while(_currentY != height/2)
		{
			while(_currentX != width)
			{
				if((*ptr) != 0.0)
				{
					uInLambda = double(int(width)/2-int(_currentX)-1) / (_fullGridder->PixelSizeX() * width);
					vInLambda = double(int(height)/2-int(_currentY)-1) / (_fullGridder->PixelSizeY() * height);
					wInLambda = 0.0;
					_fullGridder->SampleDataSample(value, uInLambda, vInLambda, wInLambda);
					value *= (*ptr);
					weight = 1.0 / (*ptr);
					++_currentX;
					return true;
				}
				
				++_currentX;
				++ptr;
			}
			_currentX = 0;
			++_currentY;
		}
		return false;
	}
	
private:
	size_t _nonZeroBinCount;
	std::unique_ptr<WStackingGridder> _fullGridder;
	ImageWeights _imageWeights;
	size_t _currentX, _currentY;
	std::complex<double>* _grid;
	ao::uvector<double> _imageWeightGrid;
	double _maxMem, _minW, _maxW, _totalWeight;
	
	size_t countNonZeroBins() const
	{
		size_t count = 0;
		const double* ptr = _imageWeightGrid.data();
		for(size_t y=0; y!=_imageWeights.Height()/2; ++y)
		{
			for(size_t x=0; x!=_imageWeights.Width(); ++x)
			{
				if((*ptr) != 0.0) ++count;
				++ptr;
			}
		}
		return count;
	}
};

#endif
