#ifndef WSCLEAN_INTERFACE_H
#define WSCLEAN_INTERFACE_H

/**
 * @file wscleaninterface.h
 * This is a low-level interface useful for using WSClean as the
 * "operator" in functions such as compressed sensing. It allows using the
 * WSClean gridding / prediction operations. This header can be used both
 * for C and for C++ programs.
 * 
 * This interface tries to abstract the measurement set, such that the interfacing
 * program does not need to read/write to the measurement set at all.
 * 
 * The general order to call this is:
 * - Fill an @ref imaging_parameters structure.
 * - Call @ref wsclean_initialize()
 * - Call @ref wsclean_read() to read the visibilities from a measurement set.
 * - Call @ref wsclean_operator_A() and / or @ref wsclean_operator_At() as many
 * times as you like to predict / image data.
 * - Call @ref wsclean_write() once or more to save image results.
 * - Call @ref wsclean_deinitialize()
 * 
 * Methods are thread safe, in that one can e.g. call the operator methods at the
 * same time from different threads, and things will still work. However, this is
 * accomplished by using a global lock, such that this will actually not speed up
 * processing.
 * 
 * @todo Currently, these methods write the visibilities to disk before performing the
 * imaging or prediction operation. This is unnecessary overhead.
 */

/**
 * @example interfaceexample.c
 * This example demonstrates how to use this interface in a C program. This file
 * can be found in the directory wsclean/examples.
 */
#include <complex.h>

#include "imaginginterface.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Initialize WSClean for use as measurement operator in compressed sensing
 * application. This should be called before any other function.
 * 
 * The way to pass the userdata pointer is as in the following example:
 * \code{.cpp}
 * void* userdata;
 * wsclean_initialize(&userdata, ...)
 * // ...
 * wsclean_deinitialize(userdata);
 * \endcode
 * 
 * @param userData should be a pointer to a void pointer which will be set to
 * a structure that WSClean internally uses.
 * @param parameters domain specific information, containing the measurement set.
 * @param imgData will be filled with info describing the data.
 */
void wsclean_initialize(
	void** userData,
	const imaging_parameters* parameters,
	imaging_data* imgData
);

/**
 * Release all resources. After this call, the userData should no longer be used.
 * Every call to @ref wsclean_initialize() should be followed by a call to 
 * wsclean_deinitialize().
 * @param userData A wsclean userdata struct as returned by @ref wsclean_initialize().
 */
void wsclean_deinitialize(void* userData);

/**
 * Reads the visibility data array from the measurement set. The returned data are
 * unweighted. The weights can be used by the client program to converge optimally,
 * but should not be applied (or should be unapplied) to the data before calling the
 * operator functions.
 * @param userData A wsclean userdata struct as returned by @ref wsclean_initialize().
 * @param data An already allocated array of data which will be set to the selected data
 * in the measurement set. The array should be of the size that @ref wsclean_initialize()
 * returned in the @ref imaging_data struct.
 * @param weights An already allocated array which will be set to the weights, of equal size
 * as the data.
 */
void wsclean_read(void* userData, DCOMPLEX* data, double* weights);

/**
 * Write the final image out.
 * @param userData A wsclean userdata struct as returned by @ref wsclean_initialize().
 * @param filename Filename of fits output file.
 * @param image The image data of size width x height.
 */
void wsclean_write(void* userData, const char* filename, const double* image);

/**
 * Calculate the unweighted visibilities for the given image data.
 * @param userData A wsclean userdata struct as returned by @ref wsclean_initialize().
 * @param dataOut Array that will be filled with the predicted visibilities.
 * @param dataIn The image data: array of size width x height.
 */
void wsclean_operator_A(void* userData, DCOMPLEX* dataOut, const double* dataIn);

/**
 * Calculate the dirty image from the visibilities. The weights will be applied
 * while imaging, so the client should not already have multiplied the visibilities
 * with the weights (see @ref wsclean_read() ). If the client does not want to apply
 * the weights, the client can divide the visibilities by the weights beforehand, and
 * multiply the image by the sum of weights afterwards.
 * @param userData A wsclean userdata struct as returned by @ref wsclean_initialize().
 * @param dataOut Array of size width x height that will be filled with the dirty image.
 * @param dataIn The visibility data to image.
 */
void wsclean_operator_At(void* userData, double* dataOut, const DCOMPLEX* dataIn);

/**
 * Convert a string with units to an angle in radians. A client program can use
 * this to convert a string like "10asec" or "1deg" to a numeric angle that can
 * be passed to @ref wsclean_initialize().
 * @param angle A string specifying an angle.
 * @return the angle converted to double, in radians.
 */
double wsclean_parse_angle(const char* angle);

#ifdef __cplusplus
}
#endif

#endif
