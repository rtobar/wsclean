#ifndef PRIMARY_BEAM_IMAGE_SET
#define PRIMARY_BEAM_IMAGE_SET

#include <vector>

#include "../matrix2x2.h"

#include "../wsclean/imagebufferallocator.h"

class PrimaryBeamImageSet
{
public:
	PrimaryBeamImageSet(size_t width, size_t height, ImageBufferAllocator& allocator) :
		_beamImages(8),
		_width(width),
		_height(height)
	{
		for(size_t i=0; i!=8; ++i)
			allocator.Allocate(width*height, _beamImages[i]);
	}
	
	void SetToZero()
	{
		for(size_t i=0; i!=8; ++i)
			std::fill(_beamImages[i].data(), _beamImages[i].data() + _width*_height, 0.0);
	}
	
	double GetUnpolarizedCorrectionFactor(size_t x, size_t y)
	{
		size_t index = y * _width + x;
		MC2x2 val, squared;
		val[0] = std::complex<double>(_beamImages[0][index], _beamImages[1][index]);
		val[1] = std::complex<double>(_beamImages[2][index], _beamImages[3][index]);
		val[2] = std::complex<double>(_beamImages[4][index], _beamImages[5][index]);
		val[3] = std::complex<double>(_beamImages[6][index], _beamImages[7][index]);
		MC2x2::ATimesHermB(squared, val, val);
		double value;
		if(squared.Invert())
			value = 0.5 * (squared[0].real() + squared[3].real());
		else
			value = std::numeric_limits<double>::quiet_NaN();
		return value;
	}
	
	const ImageBufferAllocator::Ptr& operator[](size_t index) const
	{
		return _beamImages[index];
	}
	
	PrimaryBeamImageSet& operator+=(const PrimaryBeamImageSet& rhs)
	{
		for(size_t i=0; i!=8; ++i)
		{
			for(size_t j=0; j!=_width*_height; ++j)
				_beamImages[i][j] += rhs._beamImages[i][j];
		}
		return *this;
	}
	
	PrimaryBeamImageSet& operator*=(double factor)
	{
		for(size_t i=0; i!=8; ++i)
		{
			for(size_t j=0; j!=_width*_height; ++j)
				_beamImages[i][j] *= factor;
		}
		return *this;
	}
private:
	std::vector<ImageBufferAllocator::Ptr> _beamImages;
	size_t _width, _height;
};

#endif

