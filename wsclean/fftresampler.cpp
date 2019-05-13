#include "fftresampler.h"
#include "wsclean/logger.h"

#include <complex>
#include <iostream>

FFTResampler::FFTResampler(size_t inWidth, size_t inHeight, size_t outWidth, size_t outHeight, size_t cpuCount, bool verbose) :
	_inputWidth(inWidth), _inputHeight(inHeight),
	_outputWidth(outWidth), _outputHeight(outHeight),
	_fftWidth(std::max(inWidth, outWidth)), _fftHeight(std::max(inHeight, outHeight)),
	_windowFunction(WindowFunction::Rectangular),
	_tukeyInsetSize(0.0),
	_correctWindow(false),
	_tasks(cpuCount),
	_verbose(verbose)
{
	double* inputData = reinterpret_cast<double*>(fftw_malloc(_fftWidth*_fftHeight * sizeof(double)));
	fftw_complex* fftData = reinterpret_cast<fftw_complex*>(fftw_malloc(_fftWidth*_fftHeight * sizeof(fftw_complex)));
	_inToFPlan =
		fftw_plan_dft_r2c_2d(_inputHeight, _inputWidth,
			inputData, fftData, FFTW_ESTIMATE);
	_fToOutPlan =
		fftw_plan_dft_c2r_2d(_outputHeight, _outputWidth,
			fftData, inputData, FFTW_ESTIMATE);
	fftw_free(fftData);
	fftw_free(inputData);
}

FFTResampler::~FFTResampler()
{
	Finish();
	fftw_destroy_plan(_inToFPlan);
	fftw_destroy_plan(_fToOutPlan);
}

void FFTResampler::runThread()
{
	Task task;
	while(_tasks.read(task))
	{
		runSingle(task, false);
	}
}

void FFTResampler::runSingle(const Task& task, bool skipWindow) const
{
	double *endPtr = task.input + _inputWidth*_inputHeight;
	for(double* i=task.input; i!=endPtr; ++i)
	{
		if(!std::isfinite(*i))
			*i = 0.0;
	}
	
	if(_windowFunction != WindowFunction::Rectangular && !skipWindow)
		applyWindow(task.input);
	
	size_t fftInWidth = _inputWidth/2+1;
	std::complex<double>
		*fftData = reinterpret_cast<std::complex<double>*>(fftw_malloc(fftInWidth*_inputHeight*sizeof(std::complex<double>)));
	if(_verbose) Logger::Debug << "FFT " << _inputWidth << " x " << _inputHeight << " real -> complex...\n";
	fftw_execute_dft_r2c(_inToFPlan, task.input, reinterpret_cast<fftw_complex*>(fftData));
	
	size_t fftOutWidth = _outputWidth/2+1;
	// TODO this can be done without allocating more mem!
	std::complex<double>
		*newfftData = reinterpret_cast<std::complex<double>*>(fftw_malloc(fftOutWidth*_outputHeight*sizeof(std::complex<double>)));
	std::uninitialized_fill_n(newfftData, fftOutWidth*_outputHeight, std::complex<double>(0));
	
	size_t oldMidX = _inputWidth/2;
	size_t newMidX = _outputWidth/2;
	
	size_t minWidth = std::min(_inputWidth, _outputWidth);
	size_t minHeight = std::min(_inputHeight, _outputHeight);
	
	size_t minMidX = minWidth/2;
	size_t minMidY = minHeight/2;
	
	double factor = 1.0 / (minWidth*minHeight);
	
	for(size_t y=0; y!=minHeight; ++y)
	{
		size_t oldY = y-minMidY + _inputHeight;
		size_t newY = y-minMidY + _outputHeight;
		if(oldY >= _inputHeight) oldY -= _inputHeight;
		if(newY >= _outputHeight) newY -= _outputHeight;
		
		// The last dimension is stored half
		for(size_t x=0; x!=minMidX; ++x)
		{
			size_t oldX = x;
			size_t newX = x;
			size_t oldIndex = oldX + oldY * (oldMidX+1);
			size_t newIndex = newX + newY * (newMidX+1);
			
			newfftData[newIndex] = fftData[oldIndex] * factor;
			
			//if((x == 0 && newY == 0) || (x==minMidX-1 && y==minHeight-1))
			//	std::cout << newfftData[newIndex] << " (" << oldX << " , " << oldY << ") - (" << newX << " , " << newY << ")\n";
		}
		if(_inputWidth >= _outputWidth)
		{
			size_t oldIndex = _inputWidth/2 + oldY * (oldMidX+1);
			size_t newIndex = _outputWidth/2 + newY * (newMidX+1);
			newfftData[newIndex] = fftData[oldIndex] * factor;
		}
	}
	
	fftw_free(fftData);
	
	if(_verbose) Logger::Debug << "FFT " << _outputWidth << " x " << _outputHeight << " complex -> real...\n";
	fftw_execute_dft_c2r(_fToOutPlan, reinterpret_cast<fftw_complex*>(newfftData), task.output);
	
	fftw_free(newfftData);
	
	if(_correctWindow && _windowFunction != WindowFunction::Rectangular && !skipWindow)
		unapplyWindow(task.output);
}

void FFTResampler::SingleFT(const double* input, double* realOutput, double* imaginaryOutput)
{
	ao::uvector<double> data(_inputWidth*_inputHeight);
	size_t
		halfWidth = _inputWidth/2,
		halfHeight = _inputHeight/2;
	for(size_t y=0; y!=_inputHeight; ++y)
	{
		size_t yIn = y + halfHeight;
		if(yIn >= _inputHeight) yIn -= _inputHeight;
		double* rowOutPtr = &data[y*_inputWidth];
		const double* rowInPtr = &input[yIn*_inputWidth];
		for(size_t x=0; x!=_inputWidth; ++x)
		{
			size_t xIn = x + halfWidth;
			if(xIn >= _inputWidth) xIn -= _inputWidth;
			if(std::isfinite(rowInPtr[xIn]))
				rowOutPtr[x] = rowInPtr[xIn];
			else
				rowOutPtr[x] = 0.0;
		}
	}
	
	size_t fftInWidth = _inputWidth/2+1;
	std::complex<double>
		*fftData = reinterpret_cast<std::complex<double>*>(fftw_malloc(fftInWidth*_inputHeight*sizeof(std::complex<double>)));
	if(_verbose) Logger::Debug << "FFT " << _inputWidth << " x " << _inputHeight << " real -> complex...\n";
	fftw_execute_dft_r2c(_inToFPlan, data.data(), reinterpret_cast<fftw_complex*>(fftData));
	
	size_t midX = _inputWidth/2;
	size_t midY = _inputHeight/2;
	
	double factor = 1.0 / sqrt(_inputWidth*_inputHeight);
	
	for(size_t y=0; y!=_inputHeight; ++y)
	{
		size_t oldY = y+midY;
		if(oldY >= _inputHeight) oldY -= _inputHeight;
		
		// The last dimension is stored half
		for(size_t x=0; x!=midX+1; ++x)
		{
			size_t oldIndex = x + oldY * (midX+1);
			size_t newIndex1 = midX - x + y * _inputWidth;
			
			const std::complex<double>& val = fftData[oldIndex] * factor;
			
			realOutput[newIndex1] = val.real();
			imaginaryOutput[newIndex1] = val.imag();
			if(x != midX)
			{
				size_t yTo = _inputHeight-y;
				if(yTo == _inputHeight) yTo = 0;
				size_t newIndex2 = midX + x + yTo * _inputWidth;
				realOutput[newIndex2] = val.real();
				imaginaryOutput[newIndex2] = -val.imag();
			}
		}
	}
	
	fftw_free(fftData);
}

void FFTResampler::makeWindow(ao::uvector<double>& data, size_t width) const
{
	if(_windowFunction == WindowFunction::Tukey)
		makeTukeyWindow(data, width);
	else {
		data.resize(width);
		for(size_t x=0; x!=width; ++x)
			data[x] = WindowFunction::Evaluate(_windowFunction, width, x);
	}
}

void FFTResampler::makeTukeyWindow(ao::uvector<double>& data, size_t width) const
{
	// Make a Tukey window, which consists of
	// left: a cosine going from 0 to 1
	// mid: all 1
	// right: a cosine going from 1 to 0
	data.resize(width);
	for(size_t x=0; x!=width; ++x)
	{
		// left part of Tukey window
		double xSh = (0.5+x)*2;
		if(xSh < width - _tukeyInsetSize)
		{
			double pos = xSh / (width - _tukeyInsetSize);
			data[x] = (std::cos((pos + 1.0) * M_PI) + 1.0) * 0.5;
		}
		else if(xSh < width + _tukeyInsetSize)
		{
			data[x] = 1.0;
		}
		else {
			double pos = (xSh - (width + _tukeyInsetSize)) / (width - _tukeyInsetSize);
			data[x] = (std::cos(pos * M_PI) + 1.0) * 0.5;
		}
	}
}

void FFTResampler::applyWindow(double* data) const
{
	if(_windowRowIn.empty())
	{
		makeWindow(_windowRowIn, _inputWidth);
		makeWindow(_windowColIn, _inputHeight);
		if(_correctWindow)
		{
			ao::uvector<double> windowImgIn(_inputWidth * _inputHeight);
			_windowOut.resize(_outputWidth * _outputHeight);
			Task task;
			double *inPtr = windowImgIn.data();
			for(size_t y=0; y!=_inputHeight; ++y)
			{
				for(size_t x=0; x!=_inputWidth; ++x)
				{
					*inPtr = _windowRowIn[x] * _windowColIn[y];
					++inPtr;
				}
			}
			
			task.input = windowImgIn.data();
			task.output = _windowOut.data();
			runSingle(task, true);
		}
	}
	for(size_t y=0; y!=_inputHeight; ++y)
	{
		for(size_t x=0; x!=_inputWidth; ++x)
		{
			*data *= _windowRowIn[x] * _windowColIn[y];
			++data;
		}
	}
}

void FFTResampler::unapplyWindow(double* data) const
{
	size_t n = _outputWidth * _outputHeight;
	for(size_t i=0; i!=n; ++i)
	{
		data[i] /= _windowOut[i];
	}
}
