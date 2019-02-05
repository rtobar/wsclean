#include "fitsaterm.h"

#include "../units/imagecoordinates.h"

#include "../wsclean/logger.h"

FitsATerm::FitsATerm(size_t nAntenna, size_t width, size_t height, double ra, double dec, double dl, double dm, double phaseCentreDL, double phaseCentreDM) :
	_nAntenna(nAntenna),
	_width(width),
	_height(height),
	_ra(ra), _dec(dec),
	_dl(dl), _dm(dm),
	_phaseCentreDL(phaseCentreDL),
	_phaseCentreDM(phaseCentreDM)
{ }

void FitsATerm::OpenTECFile(const std::string& filename)
{
	_mode = TECMode;
	_reader.reset(new FitsReader(filename, true, true));
	if(_reader->NFrequencies() != 1)
		throw std::runtime_error("FITS file for TEC A-terms has multiple frequencies in it");
	_nFrequencies = 1;
	if(_reader->NAntennas() != _nAntenna)
	{
		std::ostringstream str;
		str << "FITS file for TEC A-terms has incorrect number of antennas. Measurement set has "
		<< _nAntenna << " antennas, a-term FITS file has " << _reader->NAntennas() << " antennas.";
		throw std::runtime_error(str.str());
	}
	double time0 = _reader->TimeDimensionStart();
	for(size_t i=0; i!=_reader->NTimesteps(); ++i)
		_timesteps.emplace_back(time0 + i*_reader->TimeDimensionIncr());
	_curTimeindex = std::numeric_limits<size_t>::max();
}

void FitsATerm::OpenDiagGainFile(const std::string& filename)
{
	_mode = DiagonalMode;
	_reader.reset(new FitsReader(filename, true, true));
	_nFrequencies = _reader->NFrequencies();
	if(_reader->NMatrixElements() != 4)
		throw std::runtime_error("FITS file for diagonal gains did not have 4 matrix elements in it");
	if(_reader->NAntennas() != _nAntenna)
	{
		std::ostringstream str;
		str << "FITS file for diagonal gains A-terms has incorrect number of antennas. Measurement set has "
		<< _nAntenna << " antennas, a-term FITS file has " << _reader->NAntennas() << " antennas.";
		throw std::runtime_error(str.str());
	}
	double time0 = _reader->TimeDimensionStart();
	for(size_t i=0; i!=_reader->NTimesteps(); ++i)
		_timesteps.emplace_back(time0 + i*_reader->TimeDimensionIncr());
	_curTimeindex = std::numeric_limits<size_t>::max();
}

bool FitsATerm::Calculate(std::complex<float>* buffer, double time, double frequency)
{
	bool newImage = false;
	if(_curTimeindex == std::numeric_limits<size_t>::max())
	{
		newImage = true;
		_bufferCache.clear();
		_curTimeindex = 0;
	}
	
	bool finishedSearch = false;
	while(_curTimeindex+1 < _timesteps.size() && !finishedSearch)
	{
		// Do we need to calculate a next timestep?
		double curTime = _timesteps[_curTimeindex];
		double nextTime = _timesteps[_curTimeindex+1];
		/*std::cout.precision(15);
		std::cout << 
			"Searching index for time " << time << 
			", time[" << _curTimeindex << "]=" << curTime <<
			", next=" << nextTime << '\n';*/
		// If we are closer to the next timestep, use the next.
		if(std::fabs(nextTime - time) < std::fabs(curTime - time))
		{
			++_curTimeindex;
			newImage = true;
			_bufferCache.clear();
			finishedSearch = false;
		}
		else {
			finishedSearch = true;
		}
	}
	//std::cout << "Selected time index = " << _curTimeindex << '\n';
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
		return true;
	}
	return false;
}

void FitsATerm::readImages(std::complex<float>* buffer, size_t timeIndex, double frequency)
{
	_scratch.resize(_width*_height);
	ao::uvector<double> image(_reader->ImageWidth() * _reader->ImageHeight());
	for(size_t antennaIndex = 0; antennaIndex != _nAntenna; ++antennaIndex)
	{
		std::complex<float>* antennaBuffer = buffer + antennaIndex * _width*_height*4;
		
		switch(_mode)
		{
			case TECMode: { 
				// TODO When we are in the same timestep but at a different frequency, it would
				// be possible to skip reading and resampling, and immediately call evaluateTEC()
				// with the "scratch" data still there.
				
				_reader->ReadIndex(image.data(), antennaIndex + timeIndex*_nAntenna);
				
				resample(_scratch.data(), image.data());
				
				evaluateTEC(antennaBuffer, _scratch.data(), frequency);
			} break;
			
			case DiagonalMode: {
				for(size_t p=0; p!=2; ++p)
				{
					_reader->ReadIndex(image.data(), (antennaIndex + timeIndex*_nAntenna) * 4 + p*2);
					resample(_scratch.data(), image.data());
					copyToRealPolarization(antennaBuffer, _scratch.data(), p*3);
					
					_reader->ReadIndex(image.data(), (antennaIndex + timeIndex*_nAntenna) * 4 + p*2 + 1);
					resample(_scratch.data(), image.data());
					copyToImaginaryPolarization(antennaBuffer, _scratch.data(), p*3);
				}
				setPolarization(antennaBuffer, 1, std::complex<float>(0.0, 0.0));
				setPolarization(antennaBuffer, 2, std::complex<float>(0.0, 0.0));
			} break;
		}
	}
}

void FitsATerm::resample(double* dest, const double* source)
{
	size_t inWidth = _reader->ImageWidth(), inHeight = _reader->ImageHeight();
	double inPixelSizeX = _reader->PixelSizeX(), inPixelSizeY = _reader->PixelSizeY();
	double inPhaseCentreDL = _reader->PhaseCentreDL(), inPhaseCentreDM = _reader->PhaseCentreDM();
	
	double inPhaseCentreRA = _reader->PhaseCentreRA(), inPhaseCentreDec = _reader->PhaseCentreDec();
	/**
	 * If phase centra of input and output are the same, i.e. they have the same
	 * tangential plane, a few calculations can be saved.
	 */
	bool samePlane = inPhaseCentreRA == _ra && inPhaseCentreDec == _dec;
	
	size_t index = 0;
	for(size_t y=0; y!=_height; ++y)
	{
		for(size_t x=0; x!=_width; ++x)
		{
			double l, m;
			ImageCoordinates::XYToLM(x, y, _dl, _dm, _width, _height, l, m);
			l += _phaseCentreDL;
			m += _phaseCentreDM;
			if(!samePlane)
			{
				double pixra, pixdec;
				ImageCoordinates::LMToRaDec(l, m, _ra, _dec, pixra, pixdec);
				ImageCoordinates::RaDecToLM(pixra, pixdec, inPhaseCentreRA, inPhaseCentreDec, l, m);
			}
			l -= inPhaseCentreDL;
			m -= inPhaseCentreDM;
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
	for(size_t pixel = 0; pixel != _width*_height; ++pixel)
	{
		dest[pixel*4 ] = std::polar(1.0, source[pixel] * -8.44797245e9 / frequency);
		dest[pixel*4 + 1] = 0.0;
		dest[pixel*4 + 2] = 0.0;
		dest[pixel*4 + 3] = dest[pixel*4];
	}
}

void FitsATerm::copyToRealPolarization(std::complex<float>* dest, const double* source, size_t polIndex)
{
	dest += polIndex;
	for(size_t i=0; i!=_width*_height; ++i)
	{
		dest[i*4].real(source[i]);
	}
}


void FitsATerm::copyToImaginaryPolarization(std::complex<float>* dest, const double* source, size_t polIndex)
{
	dest += polIndex;
	for(size_t i=0; i!=_width*_height; ++i)
	{
		dest[i*4].imag(source[i]);
	}
}

void FitsATerm::setPolarization(std::complex<float>* dest, size_t polIndex, std::complex<float> value)
{
	dest += polIndex;
	for(size_t i=0; i!=_width*_height; ++i)
	{
		dest[i*4] = value;
	}
}
