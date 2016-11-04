#ifndef IMAGE_H
#define IMAGE_H

#include <cstring>

#include "uvector.h"

class Image
{
public:
	// Cut-off the borders of an image.
	// @param outWidth Should be <= inWidth.
	// @param outHeight Should be <= inHeight.
	static void Trim(double* output, size_t outWidth, size_t outHeight, const double* input, size_t inWidth, size_t inHeight);
	
	static void TrimBox(bool* output, size_t x1, size_t y1, size_t boxWidth, size_t boxHeight, const bool* input, size_t inWidth, size_t inHeight);
	
	// Extend an image with zeros, complement of Trim.
	// @param outWidth Should be >= inWidth.
	// @param outHeight Should be >= inHeight.
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
	
private:
	static double median_with_copy(const double* data, size_t size, ao::uvector<double>& copy);
};

#endif
