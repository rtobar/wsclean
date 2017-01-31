#ifndef IMAGE_H
#define IMAGE_H

#include <cstring>

#include "wsclean/imagebufferallocator.h"

#include "uvector.h"

class Image
{
public:
	Image(size_t width, size_t height, class ImageBufferAllocator& allocator);
	
	~Image();
	
	Image(Image&& source);
	Image& operator=(Image&& source);
	
	Image(const Image&) = delete;
	Image& operator=(const Image&) = delete;
	
	double* data() { return _data; }
	const double* data() const { return _data; }
	
	size_t Width() const { return _width; }
	size_t Height() const { return _height; }
	
	Image& operator*=(double factor);
	
	/** Cut-off the borders of an image.
	 * @param outWidth Should be &lt;= inWidth.
	 * @param outHeight Should be &lt;= inHeight.
	 */
	static void Trim(double* output, size_t outWidth, size_t outHeight, const double* input, size_t inWidth, size_t inHeight);
	
	static void TrimBox(bool* output, size_t x1, size_t y1, size_t boxWidth, size_t boxHeight, const bool* input, size_t inWidth, size_t inHeight);
	
	/** Extend an image with zeros, complement of Trim.
	 * @param outWidth Should be &gt;= inWidth.
	 * @param outHeight Should be &gt;= inHeight.
	 */
	static void Untrim(double* output, size_t outWidth, size_t outHeight, const double* input, size_t inWidth, size_t inHeight);
	
	static double Median(const double* data, size_t size)
	{
		ao::uvector<double> copy;
		return median_with_copy(data, size, copy);
	}
	
	static double MAD(const double* data, size_t size);
	
	static double StdDevFromMAD(const double* data, size_t size)
	{
		// norminv(0.75) x MAD
		return 1.48260221850560 * MAD(data, size);
	}
	
	static double RMS(const double* data, size_t size)
	{
		double sum = 0.0;
		for(size_t i=0; i!=size; ++i)
			sum += data[i]*data[i];
		return sqrt(sum/size);
	}
private:
	double* _data;
	size_t _width, _height;
	ImageBufferAllocator* _allocator;
	
	static double median_with_copy(const double* data, size_t size, ao::uvector<double>& copy);
};

#endif
