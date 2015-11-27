#ifndef IMAGE_OPERATIONS_H
#define IMAGE_OPERATIONS_H

class ImageOperations
{
public:
	// Cut-off the borders of an image.
	// @param outWidth Should be <= inWidth.
	// @param outHeight Should be <= inHeight.
	static void Trim(double* output, size_t outWidth, size_t outHeight, double* input, size_t inWidth, size_t inHeight)
	{
		size_t startX = (inWidth - outWidth) / 2;
		size_t startY = (inHeight - outHeight) / 2;
		size_t endY = (inHeight + outHeight) / 2;
		for(size_t y=startY; y!=endY; ++y)
		{
			memcpy(&output[(y-startY)*outWidth], &input[y*inWidth + startX], outWidth*sizeof(double));
		}
	}
	
	// Extend an image with zeros, complement of Trim.
	// @param outWidth Should be >= inWidth.
	// @param outHeight Should be >= inHeight.
	static void Untrim(double* output, size_t outWidth, size_t outHeight, double* input, size_t inWidth, size_t inHeight)
	{
		size_t startX = (outWidth - inWidth) / 2;
		size_t endX = (outWidth + inWidth) / 2;
		size_t startY = (outHeight - inHeight) / 2;
		size_t endY = (outHeight + inHeight) / 2;
		for(size_t y=0; y!=startY; ++y)
		{
			double* ptr = &output[y*outWidth];
			for(size_t x=0; x!=outWidth; ++x)
				ptr[x] = 0.0;
		}
		for(size_t y=startY; y!=endY; ++y)
		{
			double* ptr = &output[y*outWidth];
			for(size_t x=0; x!=startX; ++x)
				ptr[x] = 0.0;
			memcpy(&output[y*outWidth + startX], &input[(y-startY)*inWidth], inWidth*sizeof(double));
			for(size_t x=endX; x!=outWidth; ++x)
				ptr[x] = 0.0;
		}
		for(size_t y=endY; y!=outHeight; ++y)
		{
			double* ptr = &output[y*outWidth];
			for(size_t x=0; x!=outWidth; ++x)
				ptr[x] = 0.0;
		}
	}
};

#endif
