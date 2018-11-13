#ifndef IMAGING_INTERFACE_H
#define IMAGING_INTERFACE_H

/**
 * @file imaginginterface.h
 * This file is the generic part of an interface between imaging operators.
 * It is used in @ref wscleaninterface.h. A client program can use this file,
 * put it in their code repository and link to the WSClean library to avoid
 * needing to compile C++ code.
 */

#ifndef DCOMPLEX

#ifdef __cplusplus
#include <complex>
#define DCOMPLEX std::complex<double>
#else
/** In C++, this is defined as std::complex<double>, while in
 * plain C, this is defined as 'double complex'. */
#define DCOMPLEX double complex
#endif
	
#endif

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Class that defines the image properties.
 * Used by @ref wsclean_initialize().
 */
typedef struct
{
	/** Path to the measurement set */
	const char* msPath;
	/** Image width in pixels */
	int imageWidth;
	/** Image height in pixels */
	int imageHeight;
	/** Horizontal size of a pixel in radians */
	double pixelScaleX;
	/** Vertical size of a pixel in radians */
	double pixelScaleY;
	/** A string containing command line parameters that are passed to WSClean,
	 * for example "-weight natural" to use natural weighting. */
	const char* extraParameters;
	/** Whether the PSF should be normalized to one, as is the normal case.
	 * When set to true, the image is in Jy/by. When set to false, the image
	 * is unnormalized. */
	int doNormalize;
} imaging_parameters;

/**
 * Holds information about the visibility data.
 */
typedef struct
{
	/** Number of visibilities in the set (that are selected). */
	long dataSize;
	enum { DATA_TYPE_DOUBLE, DATA_TYPE_COMPLEX_DOUBLE }
		lhs_data_type,
		rhs_data_type;
		
	/*void (*deinitialize_function)(void* userData);
	void (*read_function)(void* userData, DCOMPLEX* data, double* weights);
	void (*write_function)(void* userData, const double* image);
	void (*operator_A_function)(void* userData, void* dataOut, void* dataIn);
	void (*operator_At_function)(void* userData, void* dataOut, void* dataIn);*/
} imaging_data;

#ifdef __cplusplus
}
#endif

#endif
