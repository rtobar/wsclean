#include <iostream>

#include "../fitsreader.h"
#include "../fitswriter.h"

#include "multiscalealgorithm.h"

#include "../wsclean/imagebufferallocator.h"

double rms(const double* image, size_t size)
{
	double sum = 0.0;
	for(const double* i=image; i!=image+size; ++i)
		sum += (*i)*(*i);
	return sqrt(sum/size);
}

int main(int argc, char* argv[])
{
	if(argc != 3)
	{
		std::cout << "Syntax: mscaleexample <image> <psf>\n";
	}
	else {
		FitsReader imgReader(argv[1]), psfReader(argv[2]);
		double beamScale = imgReader.BeamMajorAxisRad() / imgReader.PixelSizeX();
		size_t width = imgReader.ImageWidth(), height = imgReader.ImageHeight();
		
		ImageBufferAllocator allocator;
		ImageBufferAllocator::Ptr image, psf;
		allocator.Allocate(width*height, image);
		allocator.Allocate(width*height, psf);
		imgReader.Read(image.data());
		psfReader.Read(psf.data());
	
		ImageBufferAllocator::Ptr model;
		allocator.Allocate(width*height, model);
		for(size_t i=0; i!=width*height; ++i)
			model[i] = 0.0;
	
		double gain = 0.1;
		bool allowNegativeComponents = true;
		double borderRatio = 0.05;
		
		MultiScaleAlgorithm algorithm(allocator, beamScale, gain, allowNegativeComponents, borderRatio);
		algorithm.ExecuteMajorIteration(image.data(), model.data(), psf.data(), width, height);
		
		FitsWriter writer(imgReader);
		writer.Write("result.fits", image.data());
		writer.Write("model.fits", model.data());
	}
}
