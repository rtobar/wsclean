#ifndef IMAGING_INTERFACE_H
#define IMAGING_INTERFACE_H

#ifndef DCOMPLEX

#ifdef __cplusplus
#include <complex>
#define DCOMPLEX std::complex<double>
#else
#define DCOMPLEX double complex
#endif
	
#endif

#ifdef __cplusplus
extern "C" {
#endif
	
typedef struct
{
	const char* msPath;
	int imageWidth;
	int imageHeight;
	double pixelScaleX;
	double pixelScaleY;
	const char* extraParameters;
} imaging_parameters;

typedef struct
{
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
