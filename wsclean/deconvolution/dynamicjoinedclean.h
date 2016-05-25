#ifndef DYNAMIC_JOINED_CLEAN_H
#define DYNAMIC_JOINED_CLEAN_H

#include "deconvolutionalgorithm.h"
#include "dynamicset.h"
#include "simpleclean.h"

#include "../uvector.h"

namespace ao {
	template<typename T> class lane;
}

class DynamicJoinedClean : public UntypedDeconvolutionAlgorithm
{
public:
	DynamicJoinedClean(class ImageBufferAllocator& allocator);
	
	virtual void ExecuteMajorIteration(DynamicSet& dirtySet, DynamicSet& modelSet, const ao::uvector<const double*>& psfs, size_t width, size_t height, bool& reachedMajorThreshold);
	
private:
	size_t _width, _height;
	
	double findPeak(const double *image, size_t &x, size_t &y);
	
	std::string peakDescription(const double* image, size_t& x, size_t& y);
	
	void subtractImage(double *image, const double *psf, size_t x, size_t y, double factor, size_t startY, size_t endY) const
	{
		SimpleClean::PartialSubtractImage(image, _width, _height, psf, _width, _height, x, y, factor, startY, endY);
	}
	
	class ImageBufferAllocator& _allocator;
};

#endif
