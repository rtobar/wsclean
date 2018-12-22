#include "wstackinggridder.h"
#include "imagebufferallocator.h"
#include "logger.h"

#include <fftw3.h>

#include <iostream>
#include <fstream>

WStackingGridder::WStackingGridder(size_t width, size_t height, double pixelSizeX, double pixelSizeY, size_t fftThreadCount, ImageBufferAllocator* allocator, size_t kernelSize, size_t overSamplingFactor) :
	_width(width),
	_height(height),
	_pixelSizeX(pixelSizeX),
	_pixelSizeY(pixelSizeY),
	_nWLayers(0),
	_nPasses(0),
	_curLayerRangeIndex(0),
	_minW(0.0),
	_maxW(0.0),
	_phaseCentreDL(0.0),
	_phaseCentreDM(0.0),
	_isComplex(false),
	_imageConjugatePart(false),
	_gridMode(KaiserBesselKernel),
	_overSamplingFactor(overSamplingFactor),
	_kernelSize(kernelSize),
	_imageData(fftThreadCount, 0),
	_imageDataImaginary(fftThreadCount, 0),
	_nFFTThreads(fftThreadCount),
	_imageBufferAllocator(allocator)
{
	makeKernels();
}

WStackingGridder::~WStackingGridder()
{
	try {
		for(size_t i=0; i!=_nFFTThreads; ++i)
		{
			_imageBufferAllocator->Free(_imageData[i]);
			_imageBufferAllocator->Free(_imageDataImaginary[i]);
		}
		freeLayeredUVData();
		fftw_cleanup();
	} catch(std::exception& e) { }
}

void WStackingGridder::PrepareWLayers(size_t nWLayers, double maxMem, double minW, double maxW)
{
	_minW = minW;
	_maxW = maxW;
	_nWLayers = nWLayers;
	
	if(_minW == _maxW)
	{
		// All values have the same w-value. Some computations divide by _maxW-_minW, so prevent
		// division by zero. By changing only maxW, one makes sure that layer 0 is still at the
		// exact w value of all visibilies.
		_maxW += 1.0;
	}
	
	size_t nrCopies = _nFFTThreads;
	if(nrCopies > _nWLayers) nrCopies = _nWLayers;
	double memPerImage = _width * _height * sizeof(double);
	double memPerCore = memPerImage * 5.0; // two complex ones for FFT, one for projecting on
	double remainingMem = maxMem - nrCopies * memPerCore;
	if(remainingMem <= memPerImage * _nFFTThreads)
	{
		_nFFTThreads = size_t(maxMem*3.0/(5.0*memPerCore)); // times 3/5 to use 3/5 of mem for FFTing at most
		if(_nFFTThreads==0) _nFFTThreads = 1;
		remainingMem = maxMem - _nFFTThreads * memPerCore;
		
		Logger::Warn <<
			"WARNING: the amount of available memory is too low for the image size,\n"
			"       : not all cores might be used.\n"
			"       : nr buffers avail for FFT: " << _nFFTThreads << " remaining mem: " << round(remainingMem/1.0e8)/10.0 << " GB \n";
	}
	
	// Allocate FFT buffers
	size_t imgSize = _height * _width;
	for(size_t i=0; i!=_nFFTThreads; ++i)
	{
		_imageData[i] = _imageBufferAllocator->Allocate(imgSize);
		memset(_imageData[i], 0, imgSize * sizeof(double));
		if(_isComplex)
		{
			_imageDataImaginary[i] = _imageBufferAllocator->Allocate(imgSize);
			memset(_imageDataImaginary[i], 0, imgSize * sizeof(double));
		}
	}
	
	// Calculate nr wlayers per pass from remaining memory
	int maxNWLayersPerPass = int((double) remainingMem / (2.0*memPerImage));
	if(maxNWLayersPerPass < 1)
		maxNWLayersPerPass=1;
	_nPasses = (nWLayers+maxNWLayersPerPass-1)/maxNWLayersPerPass;
	if(_nPasses == 0) _nPasses = 1;
	Logger::Info << "Will process " << (_nWLayers / _nPasses) << "/" << _nWLayers << " w-layers per pass.\n";
	
	_curLayerRangeIndex = 0;
}

void WStackingGridder::initializeLayeredUVData(size_t n)
{
	while(_layeredUVData.size() > n)
	{
		_imageBufferAllocator->Free(_layeredUVData.back());
		_layeredUVData.pop_back();
	}
	while(_layeredUVData.size() < n)
		_layeredUVData.push_back(_imageBufferAllocator->AllocateComplex(_width * _height));
}

void WStackingGridder::StartInversionPass(size_t passIndex)
{
	initializeSqrtLMLookupTable();
	
	_curLayerRangeIndex = passIndex;
	size_t nLayersInPass = layerRangeStart(passIndex+1) - layerRangeStart(passIndex);
	initializeLayeredUVData(nLayersInPass);
	for(size_t i=0; i!=nLayersInPass; ++i)
		std::uninitialized_fill_n(_layeredUVData[i], _width*_height, 0.0);
}

void WStackingGridder::StartPredictionPass(size_t passIndex)
{
	initializeSqrtLMLookupTableForSampling();
	
	_curLayerRangeIndex = passIndex;
	size_t layerOffset = layerRangeStart(passIndex);
	size_t nLayersInPass = layerRangeStart(passIndex+1) - layerOffset;
	initializeLayeredUVData(nLayersInPass);
	
	std::stack<size_t> layers;
	for(size_t layer=0; layer!=nLayersInPass; ++layer)
		layers.push(layer);
	
	std::mutex mutex;
	std::vector<std::thread> threadGroup;
	for(size_t i=0; i!=std::min(_nFFTThreads, nLayersInPass); ++i)
		threadGroup.emplace_back(&WStackingGridder::fftToUVThreadFunction, this, &mutex, &layers);
	for(std::thread& thr : threadGroup)
		thr.join();
}

void WStackingGridder::fftToImageThreadFunction(std::mutex *mutex, std::stack<size_t> *tasks, size_t threadIndex)
{
	const size_t imgSize = _width * _height;
	std::complex<double> *fftwIn = _imageBufferAllocator->AllocateComplex(imgSize);
	// reinterpret_cast<std::complex<double>*>(fftw_malloc(imgSize * sizeof(double) * 2));
	std::complex<double> *fftwOut = _imageBufferAllocator->AllocateComplex(imgSize);
	// reinterpret_cast<std::complex<double>*>(fftw_malloc(imgSize * sizeof(double) * 2));
	
	std::unique_lock<std::mutex> lock(*mutex);
	fftw_plan plan =
		fftw_plan_dft_2d(_height, _width,
			reinterpret_cast<fftw_complex*>(fftwIn), reinterpret_cast<fftw_complex*>(fftwOut),
			FFTW_BACKWARD, FFTW_ESTIMATE);
		
	const size_t layerOffset = layerRangeStart(_curLayerRangeIndex);

	while(!tasks->empty())
	{
		size_t layer = tasks->top();
		tasks->pop();
		lock.unlock();
		
		// Fourier transform the layer
		std::complex<double> *uvData = _layeredUVData[layer];
		memcpy(fftwIn, uvData, imgSize * sizeof(double) * 2);
		fftw_execute(plan);
		
		// Add layer to full image
		if(_isComplex)
			projectOnImageAndCorrect<true>(fftwOut, LayerToW(layer + layerOffset), threadIndex);
		else
			projectOnImageAndCorrect<false>(fftwOut, LayerToW(layer + layerOffset), threadIndex);
		
		// lock for accessing tasks in guard
		lock.lock();
	}
	// Lock is still required for destroying plan
	fftw_destroy_plan(plan);
	lock.unlock();
	_imageBufferAllocator->Free(fftwIn);
	_imageBufferAllocator->Free(fftwOut);
}

void WStackingGridder::fftToUVThreadFunction(std::mutex *mutex, std::stack<size_t> *tasks)
{
	const size_t imgSize = _width * _height;
	std::complex<double> *fftwIn = _imageBufferAllocator->AllocateComplex(imgSize);
		// reinterpret_cast<std::complex<double>*>(fftw_malloc(imgSize * sizeof(double) * 2));
	std::complex<double> *fftwOut = _imageBufferAllocator->AllocateComplex(imgSize);
		// reinterpret_cast<std::complex<double>*>(fftw_malloc(imgSize * sizeof(double) * 2));
	
	std::unique_lock<std::mutex> lock(*mutex);
	fftw_plan plan =
		fftw_plan_dft_2d(_height, _width,
			reinterpret_cast<fftw_complex*>(fftwIn), reinterpret_cast<fftw_complex*>(fftwOut),
			FFTW_FORWARD, FFTW_ESTIMATE);
		
	const size_t layerOffset = layerRangeStart(_curLayerRangeIndex);

	while(!tasks->empty())
	{
		size_t layer = tasks->top();
		tasks->pop();
		lock.unlock();
		
		// Make copy of input and w-correct it
		if(_isComplex)
			copyImageToLayerAndInverseCorrect<true>(fftwIn, LayerToW(layer + layerOffset));
		else
			copyImageToLayerAndInverseCorrect<false>(fftwIn, LayerToW(layer + layerOffset));
		
		// Fourier transform the layer
		fftw_execute(plan);
		std::complex<double> *uvData = _layeredUVData[layer];
		memcpy(uvData, fftwOut, imgSize * sizeof(double) * 2);
		
		// lock for accessing tasks in guard
		lock.lock();
	}
	// Lock is still required for destroying plan
	fftw_destroy_plan(plan);
	lock.unlock();
	
	_imageBufferAllocator->Free(fftwIn);
	_imageBufferAllocator->Free(fftwOut);
}

void WStackingGridder::FinishInversionPass()
{
	size_t layerOffset = layerRangeStart(_curLayerRangeIndex);
	size_t nLayersInPass = layerRangeStart(_curLayerRangeIndex+1) - layerOffset;
	std::stack<size_t> layers;
	for(size_t layer=0; layer!=nLayersInPass; ++layer)
		layers.push(layer);
	
	std::mutex mutex;
	std::vector<std::thread> threadGroup;
	for(size_t i=0; i!=std::min(_nFFTThreads, nLayersInPass); ++i)
		threadGroup.emplace_back(&WStackingGridder::fftToImageThreadFunction, this, &mutex, &layers, i);
	for(std::thread& thr : threadGroup)
		thr.join();
}

void WStackingGridder::makeKernels()
{
	_griddingKernels.resize(_overSamplingFactor);
	_1dKernel.resize(_kernelSize*_overSamplingFactor);
	const double alpha = 8.6;
	
	switch(_gridMode)
	{
		case NearestNeighbourGridding:
		case KaiserBesselKernel:
			makeKaiserBesselKernel(_1dKernel, alpha, _overSamplingFactor, true);
			break;
		case KaiserBesselWithoutSinc:
			makeKaiserBesselKernel(_1dKernel, alpha, _overSamplingFactor, false);
			break;
		case RectangularKernel:
			makeRectangularKernel(_1dKernel, _overSamplingFactor);
			break;
		case GaussianKernel:
			makeGaussianKernel(_1dKernel, _overSamplingFactor, true);
			break;
		case GaussianKernelWithoutSinc:
			makeGaussianKernel(_1dKernel, _overSamplingFactor, false);
			break;
		case BlackmanNuttallKernel:
			makeBlackmanNutallKernel(_1dKernel, _overSamplingFactor, true);
			break;
		case BlackmanNuttallKernelWithoutSinc:
			makeBlackmanNutallKernel(_1dKernel, _overSamplingFactor, false);
			break;
	}
	
	std::vector<std::vector<double>>::iterator gridKernelIter = _griddingKernels.begin();
	for(size_t i=0; i!=_overSamplingFactor; ++i)
	{
		std::vector<double> &kernel = _griddingKernels[_overSamplingFactor - i - 1];
		kernel.resize(_kernelSize);
		std::vector<double>::iterator kernelValueIter = kernel.begin();
		for(size_t x=0; x!=_kernelSize; ++x)
		{
			size_t xIndex = x*_overSamplingFactor + i;
			*kernelValueIter = _1dKernel[xIndex];
			++kernelValueIter;
		}
		++gridKernelIter;
	}
}

void WStackingGridder::GetKernel(enum GridModeEnum gridMode, double* kernel, size_t oversampling, size_t size)
{
	double alpha = 8.6;
	std::vector<double> v(oversampling * size);
	switch(gridMode)
	{
		case KaiserBesselKernel:
			makeKaiserBesselKernel(v, alpha, oversampling, true);
			break;
		case KaiserBesselWithoutSinc:
			makeKaiserBesselKernel(v, alpha, oversampling, false);
			break;
		case GaussianKernel:
			makeGaussianKernel(v, oversampling, true);
			break;
		case GaussianKernelWithoutSinc:
			makeGaussianKernel(v, oversampling, false);
			break;
		case RectangularKernel:
			makeRectangularKernel(v, oversampling);
			break;
		case NearestNeighbourGridding:
			v[oversampling * size/2] = 1.0;
			break;
		case BlackmanNuttallKernel:
			makeBlackmanNutallKernel(v, oversampling, true);
			break;
		case BlackmanNuttallKernelWithoutSinc:
			makeBlackmanNutallKernel(v, oversampling, false);
			break;
	}
	for(size_t i=0; i!=oversampling*size; ++i)
		kernel[i] = v[i];
}

void WStackingGridder::makeKaiserBesselKernel(std::vector<double> &kernel, double alpha, size_t overSamplingFactor, bool withSinc)
{
	size_t
		n = kernel.size(),
		mid = n/2;
	std::vector<double> sincKernel(mid+1);
	const double filterRatio = 1.0 / double(overSamplingFactor); // FILTER POINT / TOTAL BANDWIDTH
	sincKernel[0] = filterRatio;
	for(size_t i=1; i!=mid+1; i++)
	{
		double x = i;
		sincKernel[i] = withSinc ? (sin(M_PI*filterRatio*x)/(M_PI*x)) : filterRatio;
	}
	const double normFactor = double(overSamplingFactor) / bessel0(alpha, 1e-8);
	for(size_t i=0; i!=mid+1; i++)
	{
		double term = double(i)/mid;
		kernel[mid+i] = sincKernel[i] * bessel0(alpha * sqrt(1.0-(term*term)), 1e-8) * normFactor;
	}
	for(size_t i=0; i!=mid; i++)
		kernel[i] = kernel[n-1-i];
}

void WStackingGridder::makeRectangularKernel(std::vector<double> &kernel, size_t overSamplingFactor)
{
	size_t
		n = kernel.size(),
		mid = n/2;
	const double filterRatio = 1.0 / double(overSamplingFactor); // FILTER POINT / TOTAL BANDWIDTH
	kernel[mid] = 1.0;
	const double normFactor = double(overSamplingFactor);
	for(size_t i=1; i!=mid+1; i++)
	{
		double x = i;
		kernel[mid+i] = normFactor * sin(M_PI*filterRatio*x)/(M_PI*x);
	}
	for(size_t i=0; i!=mid; i++)
		kernel[i] = kernel[n-1-i];
}

void WStackingGridder::makeGaussianKernel(std::vector<double> &kernel, size_t overSamplingFactor, bool withSinc)
{
	size_t
		n = kernel.size(),
		mid = n/2;
	const double filterRatio = 1.0 / double(overSamplingFactor);
	kernel[mid] = 1.0;
	const double normFactor = double(overSamplingFactor);
	for(size_t i=1; i!=mid+1; i++)
	{
		double x = i, y = x*x*3.0*3.0/(mid*mid);
		kernel[mid+i] = withSinc ? normFactor * sin(M_PI*filterRatio*x)/(M_PI*x) * exp(y/-2.0) : exp(y/-2.0);
	}
	for(size_t i=0; i!=mid; i++)
		kernel[i] = kernel[n-1-i];
}

void WStackingGridder::makeBlackmanNutallKernel(std::vector<double> &kernel, size_t overSamplingFactor, bool withSinc)
{
	size_t
		n = kernel.size(),
		mid = n/2;
	const double filterRatio = 1.0 / double(overSamplingFactor);
	kernel[mid] = 1.0;
	const double normFactor = double(overSamplingFactor);
	const double
		a0 = 0.3635819,
		a1 = 0.4891775,
		a2 = 0.1365995,
		a3 = 0.0106411;
	for(size_t i=1; i!=mid+1; i++)
	{
		double x = i, y = mid + i;
		double bn = a0 - a1 * cos(2.0*M_PI*y/(n-1)) + a2 * cos(4.0*M_PI*y/(n-1)) - a3 * cos(6.0*M_PI*y/(n-1));
		if(withSinc)
			kernel[mid+i] = normFactor * sin(M_PI*filterRatio*x)/(M_PI*x) * bn;
		else
			kernel[mid+i] = bn;
	}
	for(size_t i=0; i!=mid; i++)
		kernel[i] = kernel[n-1-i];
}

double WStackingGridder::bessel0(double x, double precision)
{
	// Calculate I_0 = SUM of m 0 -> inf [ (x/2)^(2m) ]
	// This is the unnormalized bessel function of order 0.
	double
		d = 0.0,
		ds = 1.0,
		sum = 1.0;
	do
	{
		d += 2.0;
		ds *= x*x/(d*d);
		sum += ds;
	} while (ds > sum*precision);
	return sum;
}

void WStackingGridder::AddDataSample(std::complex<float> sample, double uInLambda, double vInLambda, double wInLambda)
{
 	const size_t
		layerOffset = layerRangeStart(_curLayerRangeIndex),
		layerRangeEnd = layerRangeStart(_curLayerRangeIndex+1);
	if(_imageConjugatePart)
	{
		uInLambda = -uInLambda;
		vInLambda = -vInLambda;
		sample = std::conj(sample);
	}
	if(wInLambda < 0.0 && !_isComplex)
	{
		uInLambda = -uInLambda;
		vInLambda = -vInLambda;
		wInLambda = -wInLambda;
		sample = std::conj(sample);
	}
	size_t
		wLayer = WToLayer(wInLambda);
	if(wLayer >= layerOffset && wLayer < layerRangeEnd)
	{
		size_t layerIndex = wLayer - layerOffset;
		std::complex<double>* uvData = _layeredUVData[layerIndex];
		if(_gridMode == NearestNeighbourGridding)
		{
			int
				x = int(round(uInLambda * _pixelSizeX * _width)),
				y = int(round(vInLambda * _pixelSizeY * _height));
			if(x > -int(_width)/2 && y > -int(_height)/2 && x <= int(_width)/2 && y <= int(_height)/2)
			{
				if(x < 0) x += _width;
				if(y < 0) y += _height;
				uvData[x + y*_width] += sample;
			}
		}
		else {
			double
				xExact = uInLambda * _pixelSizeX * _width,
				yExact = vInLambda * _pixelSizeY * _height;
			int
				x = round(xExact),
				y = round(yExact),
				xKernelIndex = round((xExact - double(x)) * _overSamplingFactor),
				yKernelIndex = round((yExact - double(y)) * _overSamplingFactor);
			xKernelIndex = (xKernelIndex + (_overSamplingFactor*3)/2) % _overSamplingFactor;
			yKernelIndex = (yKernelIndex + (_overSamplingFactor*3)/2) % _overSamplingFactor;
			const std::vector<double>& xKernel = _griddingKernels[xKernelIndex];
			const std::vector<double>& yKernel = _griddingKernels[yKernelIndex];
			int mid = _kernelSize / 2;
			if(x > -int(_width)/2 && y > -int(_height)/2 && x <= int(_width)/2 && y <= int(_height)/2)
			{
				if(x < 0) x += _width;
				if(y < 0) y += _height;
				// Are we on the edge?
				if(x < mid || x+mid+1 >= int(_width) || y < mid || y+mid+1 >= int(_height))
				{
					for(size_t j=0; j!=_kernelSize; ++j)
					{
						const double yKernelValue = yKernel[j];
						size_t cy = ((y+j+_height-mid) % _height) * _width;
						for(size_t i=0; i!=_kernelSize; ++i)
						{
							size_t cx = (x+i+_width-mid) % _width;
							std::complex<double> *uvRowPtr = &uvData[cx + cy];
							const double kernelValue = yKernelValue * xKernel[i];
							*uvRowPtr += std::complex<double>(sample.real() * kernelValue, sample.imag() * kernelValue);
						}
					}
				}
				else {
					x -= mid;
					y -= mid;
					for(size_t j=0; j!=_kernelSize; ++j)
					{
						const double yKernelValue = yKernel[j];
						std::complex<double> *uvRowPtr = &uvData[x + y*_width];
						for(size_t i=0; i!=_kernelSize; ++i)
						{
							const double kernelValue = yKernelValue * xKernel[i];
							*uvRowPtr += std::complex<double>(sample.real() * kernelValue, sample.imag() * kernelValue);
							++uvRowPtr;
						}
						++y;
					}
				}
			}
		}
	}
}

void WStackingGridder::SampleDataSample(std::complex<double>& value, double uInLambda, double vInLambda, double wInLambda)
{
	const size_t
		layerOffset = layerRangeStart(_curLayerRangeIndex),
		layerRangeEnd = layerRangeStart(_curLayerRangeIndex+1);
		
	bool isConjugated = (wInLambda < 0.0 && !_isComplex);
	if(isConjugated)
	{
		uInLambda = -uInLambda;
		vInLambda = -vInLambda;
		wInLambda = -wInLambda;
	}
	size_t
		wLayer = WToLayer(wInLambda);
	if(wLayer >= layerOffset && wLayer < layerRangeEnd)
	{
		size_t layerIndex = wLayer - layerOffset;
		std::complex<double>* uvData = _layeredUVData[layerIndex];
		std::complex<double> sample;
		if(_gridMode == NearestNeighbourGridding)
		{
			int
				x = int(round(uInLambda * _pixelSizeX * _width)),
				y = int(round(vInLambda * _pixelSizeY * _height));
			if(x > -int(_width)/2 && y > -int(_height)/2 && x <= int(_width)/2 && y <= int(_height)/2)
			{
				if(x < 0) x += _width;
				if(y < 0) y += _height;
				sample = uvData[x + y*_width];
			} else {
				sample = std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
				//std::cout << "Sampling outside uv-plane (" << x << "," << y << ")\n";
			}
		}
		else {
			sample = 0.0;
			double
				xExact = uInLambda * _pixelSizeX * _width,
				yExact = vInLambda * _pixelSizeY * _height;
			int
				x = round(xExact),
				y = round(yExact),
				xKernelIndex = round((xExact - double(x)) * _overSamplingFactor),
				yKernelIndex = round((yExact - double(y)) * _overSamplingFactor);
			xKernelIndex = (xKernelIndex + (_overSamplingFactor*3)/2) % _overSamplingFactor;
			yKernelIndex = (yKernelIndex + (_overSamplingFactor*3)/2) % _overSamplingFactor;
			const std::vector<double> &xKernel = _griddingKernels[xKernelIndex];
			const std::vector<double> &yKernel = _griddingKernels[yKernelIndex];
			int mid = _kernelSize / 2;
			if(x > -int(_width)/2 && y > -int(_height)/2 && x <= int(_width)/2 && y <= int(_height)/2)
			{
				if(x < 0) x += _width;
				if(y < 0) y += _height;
				// Are we on the edge?
				if(x < mid || x+mid+1 >= int(_width) || y < mid || y+mid+1 >= int(_height))
				{
					for(size_t j=0; j!=_kernelSize; ++j)
					{
						const double yKernelValue = yKernel[j];
						size_t cy = ((y+j+_height-mid) % _height) * _width;
						for(size_t i=0; i!=_kernelSize; ++i)
						{
							const double kernelValue = xKernel[i] * yKernelValue;
							size_t cx = (x+i+_width-mid) % _width;
							std::complex<double> *uvRowPtr = &uvData[cx + cy];
							sample += std::complex<double>(uvRowPtr->real() * kernelValue, uvRowPtr->imag() * kernelValue);
						}
					}
				}
				else {
					x -= mid;
					y -= mid;
					for(size_t j=0; j!=_kernelSize; ++j)
					{
						const double yKernelValue = yKernel[j];
						std::complex<double> *uvRowPtr = &uvData[x + y*_width];
						for(size_t i=0; i!=_kernelSize; ++i)
						{
							const double kernelValue = xKernel[i] * yKernelValue;
							sample += std::complex<double>(uvRowPtr->real() * kernelValue, uvRowPtr->imag() * kernelValue);
							++uvRowPtr;
						}
						++y;
					}
				}
			}
			else {
				sample = std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
				//std::cout << "Sampling outside uv-plane (" << x << "," << y << ")\n";
			}
		}
		if(isConjugated)
			value = sample;
		else
			value = std::conj(sample);
	} else {
		value = std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
	}
}

void WStackingGridder::FinalizeImage(double multiplicationFactor, bool correctFFTFactor)
{
	freeLayeredUVData();
	if(correctFFTFactor)
	{
		multiplicationFactor /= sqrt(_width*_height);
	}
	finalizeImage(multiplicationFactor, _imageData);
	if(_isComplex)
		finalizeImage(multiplicationFactor, _imageDataImaginary);
}

void WStackingGridder::finalizeImage(double multiplicationFactor, std::vector<double*>& dataArray)
{
	for(size_t i=1;i!=_nFFTThreads;++i)
	{
		double *primaryData = dataArray[0];
		double *endPtr = dataArray[i] + (_width * _height);
		for(double *dataPtr = dataArray[i]; dataPtr!=endPtr; ++dataPtr)
		{
			*primaryData += *dataPtr;
			++primaryData;
		}
	}
	double *dataPtr = dataArray[0];
	for(size_t y=0;y!=_height;++y)
	{
		//double m = ((double) y-(_height/2)) * _pixelSizeY + _phaseCentreDM;
		for(size_t x=0;x!=_width;++x)
		{
			//double l = ((_width/2)-(double) x) * _pixelSizeX + _phaseCentreDL;
			*dataPtr *= multiplicationFactor;
			++dataPtr;
		}
	}
	
	if(_gridMode != NearestNeighbourGridding)
		correctImageForKernel<false>(dataArray[0]);
}

template<bool Inverse>
void WStackingGridder::correctImageForKernel(double *image) const
{
	const size_t nX = _width * _overSamplingFactor, nY = _height * _overSamplingFactor;
	
	double
		*fftwInX = reinterpret_cast<double*>(fftw_malloc(nX/2 * sizeof(double))),
		*fftwOutX = reinterpret_cast<double*>(fftw_malloc(nX/2 * sizeof(double)));
	double
		*fftwOutY;
	fftw_plan planX = fftw_plan_r2r_1d(nX/2, fftwInX, fftwOutX, FFTW_REDFT01, FFTW_ESTIMATE);
	memset(fftwInX, 0, nX/2 * sizeof(double));
	memcpy(fftwInX, &_1dKernel[_kernelSize*_overSamplingFactor/2], (_kernelSize*_overSamplingFactor/2+1) * sizeof(double));
	fftw_execute(planX);
	fftw_free(fftwInX);
	fftw_destroy_plan(planX);
	if(_width == _height)
	{
		fftwOutY = fftwOutX;
	}
	else {
		double *fftwInY = reinterpret_cast<double*>(fftw_malloc(nY/2 * sizeof(double)));
		fftwOutY = reinterpret_cast<double*>(fftw_malloc(nY/2 * sizeof(double)));
		fftw_plan planY = fftw_plan_r2r_1d(nY/2, fftwInY, fftwOutY, FFTW_REDFT01, FFTW_ESTIMATE);
		memset(fftwInY, 0, nY/2 * sizeof(double));
		memcpy(fftwInY, &_1dKernel[_kernelSize*_overSamplingFactor/2], (_kernelSize*_overSamplingFactor/2+1) * sizeof(double));
		fftw_execute(planY);
		fftw_free(fftwInY);
		fftw_destroy_plan(planY);
	}
	
	double normFactor = 1.0 / (_overSamplingFactor * _overSamplingFactor);
	for(size_t y=0; y!=_height; ++y)
	{
		for(size_t x=0; x!=_width; ++x)
		{
			double xVal = (x>=_width/2) ? fftwOutX[x-_width/2] : fftwOutX[_width/2-x];
			double yVal = (y>=_height/2) ? fftwOutY[y-_height/2] : fftwOutY[_height/2-y];
			if(Inverse)
				*image *= xVal * yVal * normFactor;
			else
				*image /= xVal * yVal * normFactor;
			++image;
		}
	}
	
	fftw_free(fftwOutX);
	if(_width != _height)
	{
		fftw_free(fftwOutY);
	}
}

void WStackingGridder::GetGriddingCorrectionImage(double *image) const
{
	for(size_t i=0; i!=_width*_height; ++i)
		image[i] = 1.0;
	correctImageForKernel<true>(image);
}

void WStackingGridder::initializePrediction(const double* image, std::vector<double*>& dataArray)
{
	double *dataPtr = dataArray[0];
	const double *inPtr = image;
	for(size_t y=0;y!=_height;++y)
	{
		double m = ((double) y-(_height/2)) * _pixelSizeY + _phaseCentreDM;
		for(size_t x=0;x!=_width;++x)
		{
			double l = ((_width/2)-(double) x) * _pixelSizeX + _phaseCentreDL;
			if(std::isfinite(*dataPtr) && l*l + m*m < 1.0)
				*dataPtr = *inPtr;
			else
				*dataPtr = 0.0;
			++dataPtr;
			++inPtr;
		}
	}
	if(_gridMode != NearestNeighbourGridding)
	{
		correctImageForKernel<false>(dataArray[0]);
	}
}

void WStackingGridder::initializeSqrtLMLookupTable()
{
	_sqrtLMLookupTable.resize(_width * _height);
	std::vector<double>::iterator iter = _sqrtLMLookupTable.begin();
	for(size_t y=0;y!=_height;++y)
	{
		size_t ySrc = (_height - y) + _height / 2;
		if(ySrc >= _height) ySrc -= _height;
		double m = ((double) ySrc-(_height/2)) * _pixelSizeY + _phaseCentreDM;
		
		for(size_t x=0;x!=_width;++x)
		{
			size_t xSrc = x + _width / 2;
			if(xSrc >= _width) xSrc -= _width;
			double l = ((_width/2)-(double) xSrc) * _pixelSizeX + _phaseCentreDL;
			
			if(l*l + m*m < 1.0)
				*iter = sqrt(1.0 - l*l - m*m) - 1.0;
			else
				*iter = 0.0;
			++iter;
		}
	}
}

template<bool IsComplexImpl>
void WStackingGridder::projectOnImageAndCorrect(const std::complex<double> *source, double w, size_t threadIndex)
{
	double *dataReal = _imageData[threadIndex], *dataImaginary;
	if(IsComplexImpl)
		dataImaginary = _imageDataImaginary[threadIndex];
	
	const double twoPiW = -2.0 * M_PI * w;
	std::vector<double>::const_iterator sqrtLMIter = _sqrtLMLookupTable.begin();
	for(size_t y=0;y!=_height;++y)
	{
		size_t ySrc = (_height - y) + _height / 2;
		if(ySrc >= _height) ySrc -= _height;
		
		for(size_t x=0;x!=_width;++x)
		{
			size_t xSrc = x + _width / 2;
			if(xSrc >= _width) xSrc -= _width;
			
			double
				rad = twoPiW * *sqrtLMIter,
				s = sin(rad),
				c = cos(rad);
			/*std::complex<double> val = std::complex<double>(
				source->real() * c - source->imag() * s,
				source->real() * s + source->imag() * c
			);*/
			dataReal[xSrc + ySrc*_width] += source->real()*c - source->imag()*s;
			if(IsComplexImpl)
			{
				if(_imageConjugatePart)
					dataImaginary[xSrc + ySrc*_width] += -source->real()*s + source->imag()*c;
				else
					dataImaginary[xSrc + ySrc*_width] += source->real()*s + source->imag()*c;
			}
			
			++source;
			++sqrtLMIter;
		}
	}
}

void WStackingGridder::initializeSqrtLMLookupTableForSampling()
{
	_sqrtLMLookupTable.resize(_width * _height);
	std::vector<double>::iterator iter = _sqrtLMLookupTable.begin();
	for(size_t y=0;y!=_height;++y)
	{
		//size_t yDest = (_height - y) + _height / 2;
		size_t yDest = y + _height / 2;
		if(yDest >= _height) yDest -= _height;
		double m = ((double) yDest-(_height/2)) * _pixelSizeY + _phaseCentreDM;
		
		for(size_t x=0;x!=_width;++x)
		{
			//size_t xDest = x + _width / 2;
			size_t xDest = (_width - x) + _width / 2;
			if(xDest >= _width) xDest -= _width;
			double l = ((_width/2)-(double) xDest) * _pixelSizeX + _phaseCentreDL;
			
			if(l*l + m*m < 1.0)
				*iter = sqrt(1.0 - l*l - m*m) - 1.0;
			else
				*iter = 0.0;
			++iter;
		}
	}
}

template<bool IsComplexImpl>
void WStackingGridder::copyImageToLayerAndInverseCorrect(std::complex<double> *dest, double w)
{
	double *dataReal = _imageData[0], *dataImaginary;
	if(IsComplexImpl)
		dataImaginary = _imageDataImaginary[0];
	
	const double twoPiW = 2.0 * M_PI * w;
	std::vector<double>::const_iterator sqrtLMIter = _sqrtLMLookupTable.begin();
	for(size_t y=0;y!=_height;++y)
	{
		// The fact that yDest is different than ySrc as above, is because of the way
		// fftw expects the data to be ordered.
		//size_t ySrc = (_height - y) + _height / 2;
		//if(ySrc >= _height) ySrc -= _height;
		size_t yDest = y + _height / 2;
		if(yDest >= _height) yDest -= _height;
		
		for(size_t x=0;x!=_width;++x)
		{
			size_t xDest = (_width - x) + _width / 2;
			if(xDest >= _width) xDest -= _width;
			//size_t xSrc = x + _width / 2;
			//if(xSrc >= _width) xSrc -= _width;
			
			double
				rad = twoPiW * *sqrtLMIter,
				s = sin(rad),
				c = cos(rad);
			double realVal = dataReal[xDest + yDest*_width];
			if(IsComplexImpl)
			{
				double imagVal = -dataImaginary[xDest + yDest*_width];
				*dest = std::complex<double>(realVal*c + imagVal*s, imagVal*c - realVal*s);
			}
			else
				*dest = std::complex<double>(realVal*c, -realVal*s);
			
			++dest;
			++sqrtLMIter;
		}
	}
}

void WStackingGridder::ReplaceRealImageBuffer(double* newBuffer)
{
	_imageBufferAllocator->Free(_imageData[0]);
	_imageData[0] = newBuffer;
}

void WStackingGridder::ReplaceImaginaryImageBuffer(double* newBuffer)
{
	_imageBufferAllocator->Free(_imageDataImaginary[0]);
	_imageDataImaginary[0] = newBuffer;
}

#ifndef AVOID_CASACORE
void WStackingGridder::AddData(const std::complex<float>* data, size_t dataDescId, double uInM, double vInM, double wInM)
{
	const BandData& curBand = _bandData[dataDescId];
	for(size_t ch=0; ch!=curBand.ChannelCount(); ++ch)
	{
		double
			wavelength = curBand.ChannelWavelength(ch),
			u = uInM / wavelength,
			v = vInM / wavelength,
			w = wInM / wavelength;
		AddDataSample(data[ch], u, v, w);
	}
}

void WStackingGridder::SampleData(std::complex<float>* data, size_t dataDescId, double uInM, double vInM, double wInM)
{
	const BandData& curBand(_bandData[dataDescId]);
	for(size_t ch=0; ch!=curBand.ChannelCount(); ++ch)
	{
		double
			wavelength = curBand.ChannelWavelength(ch),
			u = uInM / wavelength,
			v = vInM / wavelength,
			w = wInM / wavelength;
		SampleDataSample(data[ch], u, v, w);
	}
}

#endif

