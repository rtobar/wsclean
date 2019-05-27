#include "fitsaterm.h"

#include "../units/imagecoordinates.h"

#include "../wsclean/logger.h"

#include "../fftresampler.h"

FitsATerm::FitsATerm(size_t nAntenna, size_t width, size_t height, double ra, double dec, double dl, double dm, double phaseCentreDL, double phaseCentreDM, size_t atermSize) :
	_nAntenna(nAntenna),
	_width(width),
	_height(height),
	_ra(ra), _dec(dec),
	_dl(dl), _dm(dm),
	_phaseCentreDL(phaseCentreDL),
	_phaseCentreDM(phaseCentreDM),
	_atermWidth(atermSize),
	_atermHeight(atermSize),
	_window(WindowFunction::Rectangular),
	_padding(1.0),
	_cache(nAntenna * 4 * width * height)
{ }

FitsATerm::~FitsATerm()
{ }

void FitsATerm::OpenTECFiles(const std::vector<std::string>& filenames)
{
	_mode = TECMode;
	_readers.reserve(filenames.size());
	for(const std::string& filename : filenames)
	{
		_readers.emplace_back(filename, true, true);
		if(_readers.back().NFrequencies() != 1)
			throw std::runtime_error("FITS file for TEC A-terms has multiple frequencies in it");
	}
	initializeFromFile();
}

void FitsATerm::OpenDiagGainFiles(const std::vector<std::string>& filenames)
{
	_mode = DiagonalMode;
	_readers.reserve(filenames.size());
	for(const std::string& filename : filenames)
	{
		_readers.emplace_back(filename, true, true);
		if(_readers.back().NMatrixElements() != 4)
			throw std::runtime_error("FITS file for diagonal gains did not have 4 matrix elements in it");
	}
	initializeFromFile();
}

void FitsATerm::initializeFromFile()
{
	// Sort the readers on observation time
	std::sort(_readers.begin(), _readers.end(), [](const FitsReader& a, const FitsReader& b)->bool
		{ return a.TimeDimensionStart() < b.TimeDimensionStart(); }
	);
	_nFrequencies = _readers.front().NFrequencies();
	for(size_t readerIndex=0; readerIndex!=_readers.size(); ++readerIndex)
	{
		const FitsReader& reader = _readers[readerIndex];
		if(_nFrequencies != reader.NFrequencies())
			throw std::runtime_error("A-term FITS files have inconsistent number of frequencies");
		if(reader.NAntennas() != _nAntenna)
		{
			std::ostringstream str;
			str << "FITS file for A-terms has incorrect number of antennas. Measurement set has "
			<< _nAntenna << " antennas, a-term FITS file has " << reader.NAntennas() << " antennas.";
			throw std::runtime_error(str.str());
		}
		double time0 = reader.TimeDimensionStart();
		if(!_timesteps.empty() && time0 < _timesteps.back().time)
			throw std::runtime_error("Time axis of FITS files seem to overlap (start of fitsfile " + std::to_string(readerIndex) +
				" (t=" + std::to_string(time0) + " was before end of previous fitsfile)");
		for(size_t i=0; i!=reader.NTimesteps(); ++i)
			_timesteps.emplace_back(Timestep{time0 + i*reader.TimeDimensionIncr(), readerIndex, i});
	}
	_curTimeindex = std::numeric_limits<size_t>::max();
	_curFrequency = std::numeric_limits<double>::max();
}

bool FitsATerm::Calculate(std::complex<float>* buffer, double time, double frequency)
{
	bool requiresRead = false;
	if(_curTimeindex == std::numeric_limits<size_t>::max())
	{
		requiresRead = true;
		_cache.Reset();
		_curTimeindex = 0;
	}
	
	bool finishedSearch = false;
	while(_curTimeindex+1 < _timesteps.size() && !finishedSearch)
	{
		// Do we need to calculate a next timestep?
		double curTime = _timesteps[_curTimeindex].time;
		double nextTime = _timesteps[_curTimeindex+1].time;
		// If we are closer to the next timestep, use the next.
		if(std::fabs(nextTime - time) < std::fabs(curTime - time))
		{
			++_curTimeindex;
			requiresRead = true;
			_cache.Reset();
			finishedSearch = false;
		}
		else {
			finishedSearch = true;
		}
	}
	
	if(!requiresRead)
	{
		// If we are here, it means that the timestep didn't
		// change. So if the frequency also didn't change, we're done...
		if(_curFrequency == frequency)
			return false;
		// If it did change: do we have this frequency in the cache?
		size_t cacheIndex = _cache.Find(frequency);
		if(cacheIndex == Cache::NOT_FOUND)
			requiresRead = true;
		else {
			_cache.Get(cacheIndex, buffer);
			_curFrequency = frequency;
			return true;
		}
	}
	
	if(requiresRead)
	{
		readImages(buffer, _curTimeindex, frequency);
		_curFrequency = frequency;
		_cache.Store(frequency, buffer);
		return true;
	}
	
	return false;
}

void FitsATerm::readImages(std::complex<float>* buffer, size_t timeIndex, double frequency)
{
	_scratchA.resize(_atermWidth * _atermHeight);
	_scratchB.resize(_width * _height);
	const size_t freqIndex = round((frequency - _readers.front().FrequencyDimensionStart())
		/ _readers.front().FrequencyDimensionIncr());
	const size_t imgIndex = _timesteps[timeIndex].imgIndex * _nFrequencies + freqIndex;
	FitsReader& reader = _readers[_timesteps[timeIndex].readerIndex];
	ao::uvector<double> image(reader.ImageWidth() * reader.ImageHeight());
	if(!_resampler)
	{
		_resampler.reset(new FFTResampler(_atermWidth, _atermHeight, _width, _height, 1, false));
		if(_window == WindowFunction::Tukey)
			_resampler->SetTukeyWindow(double(_atermWidth) / _padding, false);
		else
			_resampler->SetWindowFunction(_window, true);
	}
	// TODO do this in parallel. Needs to fix Resampler too, as currently it can't run in
	// parallel when a window is used.
	for(size_t antennaIndex = 0; antennaIndex != _nAntenna; ++antennaIndex)
	{
		std::complex<float>* antennaBuffer = buffer + antennaIndex * _width*_height*4;
		
		switch(_mode)
		{
			case TECMode: { 
				// TODO When we are in the same timestep but at a different frequency, it would
				// be possible to skip reading and resampling, and immediately call evaluateTEC()
				// with the "scratch" data still there.
				reader.ReadIndex(image.data(), antennaIndex + imgIndex*_nAntenna);
				resample(reader, _scratchA.data(), image.data());
				_resampler->Resample(_scratchA.data(), _scratchB.data());
				evaluateTEC(antennaBuffer, _scratchB.data(), frequency);
			} break;
			
			case DiagonalMode: {
				for(size_t p=0; p!=2; ++p)
				{
					reader.ReadIndex(image.data(), (antennaIndex + imgIndex*_nAntenna) * 4 + p*2);
					resample(reader, _scratchA.data(), image.data());
					_resampler->Resample(_scratchA.data(), _scratchB.data());
					copyToRealPolarization(antennaBuffer, _scratchB.data(), p*3);
					
					reader.ReadIndex(image.data(), (antennaIndex + imgIndex*_nAntenna) * 4 + p*2 + 1);
					resample(reader, _scratchA.data(), image.data());
					_resampler->Resample(_scratchA.data(), _scratchB.data());
					copyToImaginaryPolarization(antennaBuffer, _scratchB.data(), p*3);
				}
				setPolarization(antennaBuffer, 1, std::complex<float>(0.0, 0.0));
				setPolarization(antennaBuffer, 2, std::complex<float>(0.0, 0.0));
			} break;
		}
	}
}

void FitsATerm::resample(const FitsReader& reader, double* dest, const double* source)
{
	size_t inWidth = reader.ImageWidth(), inHeight = reader.ImageHeight();
	double inPixelSizeX = reader.PixelSizeX(), inPixelSizeY = reader.PixelSizeY();
	double inPhaseCentreDL = reader.PhaseCentreDL(), inPhaseCentreDM = reader.PhaseCentreDM();
	double inPhaseCentreRA = reader.PhaseCentreRA(), inPhaseCentreDec = reader.PhaseCentreDec();
	double atermDL = _dl * _width / _atermWidth;
	double atermDM = _dm * _height / _atermHeight;
	/**
	 * If phase centra of input and output are the same, i.e. they have the same
	 * tangential plane, a few calculations can be saved.
	 */
	bool samePlane = inPhaseCentreRA == _ra && inPhaseCentreDec == _dec;
	
	size_t index = 0;
	for(size_t y=0; y!=_atermHeight; ++y)
	{
		for(size_t x=0; x!=_atermWidth; ++x)
		{
			double l, m;
			ImageCoordinates::XYToLM(x, y, atermDL, atermDM, _atermWidth, _atermHeight, l, m);
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
