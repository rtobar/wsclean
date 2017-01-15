#ifndef GENERIC_CLEAN_H
#define GENERIC_CLEAN_H

#include "deconvolutionalgorithm.h"
#include "imageset.h"
#include "simpleclean.h"

#include "../uvector.h"

/**
 * This class implements a generalized version of Högbom clean. It performs a single-channel
 * or joined cleaning, depending on the number of images provided. It can use the Clark optimization
 * to speed up the cleaning. When multiple frequencies are provided, it can perform spectral fitting.
 */
class GenericClean : public DeconvolutionAlgorithm
{
public:
	explicit GenericClean(class ImageBufferAllocator& allocator, bool useClarkOptimization);
	
	virtual void ExecuteMajorIteration(ImageSet& dirtySet, ImageSet& modelSet, const ao::uvector<const double*>& psfs, size_t width, size_t height, bool& reachedMajorThreshold) final override;
	
private:
	size_t _width, _height, _convolutionWidth, _convolutionHeight;
	double _convolutionPadding;
	bool _useClarkOptimization;
	
	double findPeak(const double *image, size_t &x, size_t &y);
	
	std::string peakDescription(const double* image, size_t& x, size_t& y);
	
	void subtractImage(double *image, const double *psf, size_t x, size_t y, double factor, size_t startY, size_t endY) const
	{
		SimpleClean::PartialSubtractImage(image, _width, _height, psf, _width, _height, x, y, factor, startY, endY);
	}
	
	class ImageBufferAllocator& _allocator;
};

#endif
