#include "multiscaletransforms.h"

#include "../fftconvolver.h"

void MultiScaleTransforms::Transform(const ao::uvector<double*>& images, double* scratch, double scale)
{
	ao::uvector<double> shape;
	size_t kernelSize;
	MakeShapeFunction(scale, shape, kernelSize);
	
	std::fill_n(scratch, _width*_height, 0.0);
	
	FFTConvolver::PrepareSmallKernel(scratch, _width, _height, shape.data(), kernelSize);
	for(double*const* imageIter = images.begin(); imageIter!=images.end(); ++imageIter)
		FFTConvolver::ConvolveSameSize(_fftwManager, *imageIter, scratch, _width, _height);
}

void MultiScaleTransforms::PrepareTransform(double* kernel, double scale)
{
	ao::uvector<double> shape;
	size_t kernelSize;
	MakeShapeFunction(scale, shape, kernelSize);
	
	std::fill_n(kernel, _width*_height, 0.0);
	
	FFTConvolver::PrepareSmallKernel(kernel, _width, _height, shape.data(), kernelSize);
}

void MultiScaleTransforms::FinishTransform(double* image, const double* kernel)
{
	FFTConvolver::ConvolveSameSize(_fftwManager, image, kernel, _width, _height);
}
