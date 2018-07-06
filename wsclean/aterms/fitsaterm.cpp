#include "fitsaterm.h"

#include "../units/imagecoordinates.h"

FitsATerm::FitsATerm(size_t nAntenna, size_t width, size_t height, double dl, double dm, double phaseCentreDL, double phaseCentreDM) :
	_nAntenna(nAntenna),
	_width(width),
	_height(height),
	_dl(dl), _dm(dm),
	_phaseCentreDL(phaseCentreDL),
	_phaseCentreDM(phaseCentreDM)
{
}

void FitsATerm::OpenTECFile(const std::string& filename)
{
	_reader.reset(new FitsReader(filename, true, true));
	if(_reader->NFrequencies() != 1)
		throw std::runtime_error("Fits file for TEC A-terms has multiple frequencies in it");
	double time = _reader->DateObs(); // TODO this has to be the time corresponding to the image
	_timesteps.emplace_back(time);
}

void FitsATerm::Calculate(std::complex<float>* buffer, double time, double frequency)
{
	bool newImage = false;
	// Have we arrived at the next timestep?
	if(_curTimeindex+1 < _timesteps.size())
	{
		double curTime = _timesteps[_curTimeindex];
		double nextTime = _timesteps[_curTimeindex+1];
		// If we are closer to the next timestep, use the next.
		if(std::fabs(nextTime - time) < std::fabs(curTime - time))
		{
			newImage = true;
			++_curTimeindex;
			_bufferCache.clear();
		}
	}
	// Do we have this frequency in the cache?
	auto iter = _bufferCache.find(frequency);
	if(iter == _bufferCache.end())
	{
		newImage = true;
	}
	
	if(newImage)
	{
		readImages(buffer, _curTimeindex, frequency);
		std::vector<std::complex<float>>& cacheEntry = _bufferCache[frequency];
		cacheEntry.assign(buffer, buffer + _nAntenna * 4 * _width * _height);
	}
	else {
		// Use a previously calculated buffer
		std::copy(iter->second.begin(), iter->second.end(), buffer);
	}
}

void FitsATerm::readImages(std::complex<float>* buffer, size_t timeIndex, double frequency)
{
	for(size_t antennaIndex = 0; antennaIndex != _nAntenna; ++antennaIndex)
	{
		// TODO if we are in the same timestep but at a different frequency, it would
		// be possible to script reading and resampling, and immediately call evaluateTEC()
		// with the "scratch" data still there.
		
		ao::uvector<double> image(_reader->ImageWidth() * _reader->ImageHeight());
		_reader->ReadIndex(image.data(), antennaIndex + timeIndex*_nAntenna);
		
		resample(image.data(), _scratch.data());
		
		std::complex<float>* antennaBuffer = buffer + antennaIndex * _width*_height*4;
		
		evaluateTEC(antennaBuffer, _scratch.data(), frequency);
	}
}

void FitsATerm::resample(double* dest, const double* source)
{
	size_t inWidth = _reader->ImageWidth(), inHeight = _reader->ImageHeight();
	double inPixelSizeX = _reader->PixelSizeX(), inPixelSizeY = _reader->PixelSizeY();
	double inPhaseCentreDL = _reader->PhaseCentreDL(), inPhaseCentreDM = _reader->PhaseCentreDM();
	/**
	 * We assume that the phase centra of input and output are the same, i.e. they have the same
	 * tangential plane. This saves a few calculations.
	 */
	//double inPhaseCentreRA = _reader->PhaseCentreRA(), inPhaseCentreDec = _reader->PhaseCentreDec();
	
	size_t index = 0;
	for(size_t y=0; y!=_height; ++y)
	{
		for(size_t x=0; x!=_width; ++x)
		{
			double l, m;
			ImageCoordinates::XYToLM(x, y, _dl, _dm, _width, _height, l, m);
			l += _phaseCentreDL - inPhaseCentreDL;
			m += _phaseCentreDM - inPhaseCentreDM;
			int inX, inY;
			ImageCoordinates::LMToXY(l, m, inPixelSizeX, inPixelSizeY, inWidth, inHeight, inX, inY);
			if(inX < 0 || inY < 0 || inX >= int(inWidth) || inY >= int(inHeight))
				dest[index] = 0;
			else {
				dest[index] = source[inX + inY * inWidth];
			}
			++index;
		}
	}
}

void FitsATerm::evaluateTEC(std::complex<float>* dest, const double* source, double frequency)
{
	for(size_t pixel = 0; pixel != _width*_height*4; ++pixel)
	{
		dest[pixel] = std::polar(1.0, source[pixel] * -8.44797245e9 / frequency);
	}
}
