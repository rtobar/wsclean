#include "simpleclean.h"

#include "../lane.h"

#include <boost/thread/thread.hpp>

#ifdef __SSE__
#define USE_INTRINSICS
#endif

#ifdef USE_INTRINSICS
#include <emmintrin.h>
#include <immintrin.h>
#endif

#include <iostream>
#include <limits>

boost::optional<double> SimpleClean::FindPeakSimple(const double *image, size_t width, size_t height, size_t &x, size_t &y, bool allowNegativeComponents, size_t startY, size_t endY, size_t horizontalBorder, size_t verticalBorder)
{
	double peakMax = std::numeric_limits<double>::min();
	size_t peakIndex = width * height;
	
	size_t xiStart = horizontalBorder, xiEnd = width - horizontalBorder;
	size_t yiStart = std::max(startY, verticalBorder), yiEnd = std::min(endY, height - verticalBorder);
	if(xiEnd < xiStart) xiEnd = xiStart;
	if(yiEnd < yiStart) yiEnd = yiStart;

	for(size_t yi=yiStart; yi!=yiEnd; ++yi)
	{
		size_t index = yi*width + xiStart;
		for(size_t xi=xiStart; xi!=xiEnd; ++xi)
		{
			double value = image[index];
			if(allowNegativeComponents) value = std::fabs(value);
			if(value > peakMax)
			{
				peakIndex = index;
				peakMax = std::fabs(value);
			}
			++value;
			++index;
		}
	}
	if(peakIndex == width * height)
	{
		x = width; y = height;
		return boost::optional<double>();
	}
	else {
		x = peakIndex % width;
		y = peakIndex / width;
		return image[x + y*width];
	}
}

boost::optional<double> SimpleClean::FindPeakWithMask(const double *image, size_t width, size_t height, size_t &x, size_t &y, bool allowNegativeComponents, size_t startY, size_t endY, const bool* cleanMask, size_t horizontalBorder, size_t verticalBorder)
{
	double peakMax = std::numeric_limits<double>::min();
	x = width; y = height;
	
	size_t xiStart = horizontalBorder, xiEnd = width - horizontalBorder;
	size_t yiStart = std::max(startY, verticalBorder), yiEnd = std::min(endY, height - verticalBorder);
	if(xiEnd < xiStart) xiEnd = xiStart;
	if(yiEnd < yiStart) yiEnd = yiStart;
	
	for(size_t yi=yiStart; yi!=yiEnd; ++yi)
	{
		const double *imgIter = &image[yi*width+xiStart];
		const bool* cleanMaskPtr = &cleanMask[yi*width+xiStart];
		for(size_t xi=xiStart; xi!=xiEnd; ++xi)
		{
			double value = *imgIter;
			if(allowNegativeComponents) value = std::fabs(value);
			if(value > peakMax && *cleanMaskPtr)
			{
				x = xi;
				y = yi;
				peakMax = std::fabs(value);
			}
			++imgIter;
			++cleanMaskPtr;
		}
	}
	if(y == height)
		return boost::optional<double>();
	else
		return image[x + y*width];
}

#if defined __AVX__ && defined USE_INTRINSICS && !defined FORCE_NON_AVX
template<bool AllowNegativeComponent>
boost::optional<double> SimpleClean::FindPeakAVX(const double *image, size_t width, size_t height, size_t& x, size_t& y, size_t startY, size_t endY, size_t horizontalBorder, size_t verticalBorder)
{
	double peakMax = std::numeric_limits<double>::min();
	size_t peakIndex = 0;
	
	__m256d mPeakMax = _mm256_set1_pd(peakMax);
	
	size_t xiStart = horizontalBorder, xiEnd = width - horizontalBorder;
	size_t yiStart = std::max(startY, verticalBorder), yiEnd = std::min(endY, height - verticalBorder);
	if(xiEnd < xiStart) xiEnd = xiStart;
	if(yiEnd < yiStart) yiEnd = yiStart;
	
	for(size_t yi=yiStart; yi!=yiEnd; ++yi)
	{
		size_t index = yi*width + xiStart;
		const double* const endPtr = image + yi*width + xiEnd - 4;
		const double *i=image + index;
		for(; i<endPtr; i+=4)
		{
			__m256d val = _mm256_loadu_pd(i);
			if(AllowNegativeComponent) {
				__m256d negVal = _mm256_sub_pd(_mm256_set1_pd(0.0), val);
				val = _mm256_max_pd(val, negVal);
			}
			int mask = _mm256_movemask_pd(_mm256_cmp_pd(val, mPeakMax, _CMP_GT_OQ));
			if(mask != 0)
			{
				for(size_t di=0; di!=4; ++di)
				{
					double value = i[di];
					if(AllowNegativeComponent) value = std::fabs(value);
					if(value > peakMax)
					{
						peakIndex = index+di;
						peakMax = std::fabs(i[di]);
						mPeakMax = _mm256_set1_pd(peakMax);
					}
				}
			}
			index+=4;
		}
		for(; i!=endPtr+4; ++i)
		{
			double value = *i;
			if(AllowNegativeComponent) value = std::fabs(value);
			if(value > peakMax)
			{
				peakIndex = index;
				peakMax = std::fabs(*i);
			}
			++index;
		}
	}
	x = peakIndex % width;
	y = peakIndex / width;
	return image[x + y*width];
}

template
boost::optional<double> SimpleClean::FindPeakAVX<false>(const double *image, size_t width, size_t height, size_t &x, size_t &y, size_t startY, size_t endY, size_t horizontalBorder, size_t verticalBorder);
template
boost::optional<double> SimpleClean::FindPeakAVX<true>(const double *image, size_t width, size_t height, size_t &x, size_t &y, size_t startY, size_t endY, size_t horizontalBorder, size_t verticalBorder);
#else
#warning "Not using AVX optimized version of FindPeak()!"
#endif // __AVX__

void SimpleClean::SubtractImage(double *image, const double *psf, size_t width, size_t height, size_t x, size_t y, double factor)
{
	size_t startX, startY, endX, endY;
	int offsetX = (int) x - width/2, offsetY = (int) y - height/2;
	
	if(offsetX > 0)
		startX = offsetX;
	else
		startX = 0;
	
	if(offsetY > 0)
		startY = offsetY;
	else
		startY = 0;
	
	endX = x + width/2;
	if(endX > width) endX = width;
	
	bool isAligned = ((endX - startX) % 2) == 0;
	if(!isAligned) --endX;
	
	endY = y + height/2;
	if(endY > height) endY = height;
	
	for(size_t ypos = startY; ypos != endY; ++ypos)
	{
		double *imageIter = image + ypos * width + startX;
		const double *psfIter = psf + (ypos - offsetY) * width + startX - offsetX;
		for(size_t xpos = startX; xpos != endX; xpos++)
		{
			// I've SSE-ified this, but it didn't improve speed at all :-/
			// (Compiler probably already did it)
			*imageIter -= (*psfIter * factor);
			//*(imageIter+1) = *(imageIter+1) - (*(psfIter+1) * factor);
			++imageIter;
			++psfIter;
		}
	}
}

void SimpleClean::PartialSubtractImage(double *image, const double *psf, size_t width, size_t height, size_t x, size_t y, double factor, size_t startY, size_t endY)
{
	size_t startX, endX;
	int offsetX = (int) x - width/2, offsetY = (int) y - height/2;
	
	if(offsetX > 0)
		startX = offsetX;
	else
		startX = 0;
	
	if(offsetY > (int) startY)
		startY = offsetY;
	
	endX = x + width/2;
	if(endX > width) endX = width;
	
	bool isAligned = ((endX - startX) % 2) == 0;
	if(!isAligned) --endX;
	
	endY = std::min(y + height/2, endY);
	
	for(size_t ypos = startY; ypos < endY; ++ypos)
	{
		double *imageIter = image + ypos * width + startX;
		const double *psfIter = psf + (ypos - offsetY) * width + startX - offsetX;
		for(size_t xpos = startX; xpos != endX; xpos+=2)
		{
			*imageIter = *imageIter - (*psfIter * factor);
			*(imageIter+1) = *(imageIter+1) - (*(psfIter+1) * factor);
			imageIter+=2;
			psfIter+=2;
		}
		if(!isAligned)
			*imageIter -= *psfIter * factor;
	}
}

void SimpleClean::PartialSubtractImage(double *image, size_t imgWidth, size_t imgHeight, const double *psf, size_t psfWidth, size_t psfHeight, size_t x, size_t y, double factor, size_t startY, size_t endY)
{
	size_t startX, endX;
	int offsetX = (int) x - psfWidth/2, offsetY = (int) y - psfHeight/2;
	
	if(offsetX > 0)
		startX = offsetX;
	else
		startX = 0;
	
	if(offsetY > (int) startY)
		startY = offsetY;
	
	endX = std::min(x + psfWidth/2, imgWidth);
	
	bool isAligned = ((endX - startX) % 2) == 0;
	if(!isAligned) --endX;
	
	endY = std::min(y + psfHeight/2, endY);
	
	for(size_t ypos = startY; ypos < endY; ++ypos)
	{
		double *imageIter = image + ypos * imgWidth + startX;
		const double *psfIter = psf + (ypos - offsetY) * psfWidth + startX - offsetX;
		for(size_t xpos = startX; xpos != endX; xpos+=2)
		{
			*imageIter = *imageIter - (*psfIter * factor);
			*(imageIter+1) = *(imageIter+1) - (*(psfIter+1) * factor);
			imageIter+=2;
			psfIter+=2;
		}
		if(!isAligned)
			*imageIter -= *psfIter * factor;
	}
}

#if defined __AVX__ && defined USE_INTRINSICS
void SimpleClean::PartialSubtractImageAVX(double *image, size_t imgWidth, size_t imgHeight, const double *psf, size_t psfWidth, size_t psfHeight, size_t x, size_t y, double factor, size_t startY, size_t endY)
{
	size_t startX, endX;
	int offsetX = (int) x - psfWidth/2, offsetY = (int) y - psfHeight/2;
	
	if(offsetX > 0)
		startX = offsetX;
	else
		startX = 0;
	
	if(offsetY > (int) startY)
		startY = offsetY;
	
	endX = std::min(x + psfWidth/2, imgWidth);
	
	size_t unAlignedCount = (endX - startX) % 4;
	endX -= unAlignedCount;
	
	endY = std::min(y + psfHeight/2, endY);
	
	const __m256d mFactor = _mm256_set1_pd(-factor);
	for(size_t ypos = startY; ypos < endY; ++ypos)
	{
		double *imageIter = image + ypos * imgWidth + startX;
		const double *psfIter = psf + (ypos - offsetY) * psfWidth + startX - offsetX;
		for(size_t xpos = startX; xpos != endX; xpos+=4)
		{
			__m256d
				imgVal = _mm256_loadu_pd(imageIter),
				psfVal = _mm256_loadu_pd(psfIter);
#ifdef __AVX2__
			_mm256_storeu_pd(imageIter, _mm256_fmadd_pd(psfVal, mFactor, imgVal));
#else
			_mm256_storeu_pd(imageIter, _mm256_add_pd(imgVal, _mm256_mul_pd(psfVal, mFactor)));
#endif
			imageIter+=4;
			psfIter+=4;
		}
		for(size_t xpos = endX; xpos!=endX + unAlignedCount; ++xpos)
		{
			*imageIter -= *psfIter * factor;
			++imageIter;
			++psfIter;
		}
	}
}

#endif
