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
	
	void ApplyStokesI(double* stokesI) const
	{
		// If Iu is uncorrected and Ic is corrected:
		// Iu = B Ic B^*
		// when I is unpolarized (diagonal, constant)
		// Iu = Ic B B^*
		// Ic = Iu (B B^*)^-1
		// Since we have measured Iu_xx + Iu_yy, and want to know Ic_xx + Ic_yy, let B2 = (B B^*)^-1 and Iu_xx = Iu_yy :
		// Ic_xx + Ic_yy = Iu_xx B2_xx + Iu_yy B2_yy = (Iu_xx + Iu_yy) (B2_xx + B2_yy)
		
		size_t size = _width * _height;
		for(size_t j=0; j!=size; ++j)
		{
			MC2x2 val, squared;
			val[0] = std::complex<double>(_beamImages[0][j], _beamImages[1][j]);
			val[1] = std::complex<double>(_beamImages[2][j], _beamImages[3][j]);
			val[2] = std::complex<double>(_beamImages[4][j], _beamImages[5][j]);
			val[3] = std::complex<double>(_beamImages[6][j], _beamImages[7][j]);
			MC2x2::ATimesHermB(squared, val, val);
			if(squared.Invert())
				stokesI[j] = stokesI[j] * 0.5 * (squared[0].real() + squared[3].real());
			else
				stokesI[j] = std::numeric_limits<double>::quiet_NaN();
		}
	}
	
	void ApplyFullStokes(double* images[4]) const
	{
		size_t size = _width * _height;
		for(size_t j=0; j!=size; ++j)
		{
			MC2x2 beamVal;
			beamVal[0] = std::complex<double>(_beamImages[0][j], _beamImages[1][j]);
			beamVal[1] = std::complex<double>(_beamImages[2][j], _beamImages[3][j]);
			beamVal[2] = std::complex<double>(_beamImages[4][j], _beamImages[5][j]);
			beamVal[3] = std::complex<double>(_beamImages[6][j], _beamImages[7][j]);
			if(beamVal.Invert())
			{
				double stokesVal[4] = { images[0][j], images[1][j], images[2][j], images[3][j] };
				MC2x2 linearVal, scratch;
				Polarization::StokesToLinear(stokesVal, linearVal.Data());
				MC2x2::ATimesB(scratch, beamVal, linearVal);
				MC2x2::ATimesHermB(linearVal, scratch, beamVal);
				Polarization::LinearToStokes(linearVal.Data(), stokesVal);
				for(size_t p=0; p!=4; ++p)
					images[p][j] = stokesVal[p];
			}
			else {
				for(size_t p=0; p!=4; ++p)
					images[p][j] = std::numeric_limits<double>::quiet_NaN();
			}
		}
	}
private:
	std::vector<ImageBufferAllocator::Ptr> _beamImages;
	size_t _width, _height;
};

#endif

