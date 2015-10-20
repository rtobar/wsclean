#ifndef SPECTRAL_IMAGE_FITTER
#define SPECTRAL_IMAGE_FITTER

#include "../wsclean/imagebufferallocator.h"

#include "../uvector.h"
#include "spectralfitter.h"

class SpectralImageFitter
{
public:
	SpectralImageFitter(size_t imgWidth, size_t imgHeight, ImageBufferAllocator &allocator, SpectralFitter& fitter)
		: _width(imgWidth), _height(imgHeight), _allocator(allocator), _fitter(fitter)
	{
	}
	
	void AddImage(const double* image, double frequencyHz)
	{
		double* newImg = _allocator.Allocate(_width*_height);
		memcpy(image, newImg, _width*_height*sizeof(double));
		_images.push_back(newImg);
		_frequencies.push_back(frequencyHz);
	}
	
	void PerformFit()
	{
		for(size_t i=0; i!=_fitter.NTerms(); ++i)
			_terms.push_back(_allocator.Allocate(_width*_height));
		
		ao::uvector<double> values(_images.size());
		ao::uvector<double> terms(_fitter.NTerms());
		for(size_t px=0; px!=_width*_height; ++px)
		{
			double isZero = true;
			for(size_t s=0; s!=_images.size(); ++s)
			{
				double value = _images[s][px];
				values[s] = value;
				isZero = isZero && (value == 0.0);
			}
			
			if(isZero)
				terms.assign(_fitter.NTerms(), 0.0);
			else
				_fitter.Fit(terms, values.data());
			
			
		}
	}
	
	void Interpolate(double* destination, double frequency)
	{
	}
	
private:
	size_t _width, _height;
	ImageBufferAllocator& _allocator;
	SpectralFitter& _fitter;
	ao::uvector<double*> _images;
	ao::uvector<double*> _terms;
	ao::uvector<double> _frequencies;
};

#endif
