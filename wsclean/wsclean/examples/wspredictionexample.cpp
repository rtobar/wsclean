#include <iostream>

#include "../imagebufferallocator.h"
#include "../wstackinggridder.h"

#include "../../fftwmultithreadenabler.h"
#include "../../system.h"

int main(int argc, char* argv[])
{
	size_t width = 4000, height = 4000;
	double pixelScale = 1.0/60.0*(M_PI/180.0);          // one arcmin in radians
	size_t threadCount = System::ProcessorCount(); // number of CPUs in system
	
	// Calculate available memory
	long int pageCount = sysconf(_SC_PHYS_PAGES), pageSize = sysconf(_SC_PAGE_SIZE);
	int64_t memSize = (int64_t) pageCount * (int64_t) pageSize;
	
	// Initialize an image with all zero except one pixel
	std::vector<double> image(width*height, 0.0);
	size_t sourceX = width/2 - width/17, sourceY = height/2-height/21;
	image[sourceX + sourceY*width] = 100.0;
	
	ImageBufferAllocator allocator;
	
	// Not strictly required, but the following line will turn on FFTW multi-threading.
	// In this particular scenario, without w-term correction, this speeds up the
	// FFT a bit.
	FFTWMultiThreadEnabler fftwMT;
	
	// Initialize the gridder
	// Here, larger and more accurate values are used for the kernel size (here: 15) and
	// oversampling factor (here: 201). Often, the default values are good enough and provide a
	// better compromise between accuracy and performance. To get the default values, the last
	// two parameters can be left out.
	WStackingGridder gridder(width, height, pixelScale, pixelScale, threadCount, &allocator, 15, 201);
	
	// Prepare the w-layers. This example disables w-term correction by setting the
	// number of w-layers to 1, and therefore this call becomes quite simple.
	// For w-term correction, actually proper w-extremes should have been given.
	// Specifying a memsize will allow the gridder to warn when not enough memory
	// would be available. Alternatively, a very large number can be given to avoid this.
	gridder.PrepareWLayers(1, memSize, -1.0, 1.0);
	
	// To keep things simple, we do not handle situations in which multiple passes would
	// have been required. Since only one w-layer is requested, this test is not necessary,
	// but here for good style.
	if(gridder.NPasses() != 1)
	{
		std::cerr << "Error; w-layers do not all fit in memory at once.\n";
		return -1;
	}
	
	// Now we supply the image.
	gridder.InitializePrediction(&image[0]);
	
	// Start the first (and only) pass. This will perform the FFTs for this pass.
	// Since there is only one w-layer, only one FFT will be performed.
	gridder.StartPredictionPass(0);
	
	// We can now start to sample visibilities
	// We calculate some random uvw's that are likely within the
	// gridded area:
	double
		u = 100.0 , // These are in units of number of wavelengths
		v =  10.0 ,
		w =  10.0 ;
		
	for(size_t i=0; i!=20; ++i, u/=2.0, v/=2.0, w/=2.0)
	{
		std::complex<float> value;
		gridder.SampleDataSample(value, u, v, w);
		std::cout << "uvw = (" << u << ", " << v << ", " << w << "), visibility = " << value << '\n';
	}
}
