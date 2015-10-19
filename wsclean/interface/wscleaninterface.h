#ifndef WSCLEAN_INTERFACE_H
#define WSCLEAN_INTERFACE_H

#include <complex.h>

#include "imaginginterface.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Initialize WSClean for use as measurement operator in compressed sensing
 * application. This should be called before any other function.
 * 
 * @p userdata should be a pointer to a void pointer which will be set to
 * a structure that WSClean internally uses. The way to pass this pointer
 * is like this:
 * void* userdata;
 * wsclean_initialize(&userdata, ...)
 * ...
 * wsclean_deinitialize(userdata);
 * @p domain_info domain specific information, containing the measurement set.
 * @p data_info will be filled with info describing the data.
 */
void wsclean_initialize(
	void** userData,
	const imaging_parameters* domain_info,
	imaging_data* data_info
);

/**
 * Clean up.
 */
void wsclean_deinitialize(void* userData);

/**
 * Initializes the data array
 * @p data An already allocated array of data which will be set to the selected data
 * in the measurement set.
 * @p weights An already allocated array which will be set to the weights.
 */
void wsclean_read(void* userData, DCOMPLEX* data, double* weights);

/**
 * Write the final image out.
 * @p image The image data.
 */
void wsclean_write(void* userData, const double* image);

/**
 * Calculate the visibilities given the image data.
 */
void wsclean_operator_A(void* userData, DCOMPLEX* dataOut, const double* dataIn);

void wsclean_operator_At(void* userData, double* dataOut, const DCOMPLEX* dataIn);

double wsclean_parse_angle(const char* angle);

#ifdef __cplusplus
}
#endif

#endif
