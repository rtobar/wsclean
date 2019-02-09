#include "atermbase.h"

#include "../fitswriter.h"
#include "../matrix2x2.h"
#include "../uvector.h"

#include "../wsclean/logger.h"

void ATermBase::StoreATermsEigenvalues(const std::string& filename, const std::complex<float>* buffer, size_t nStations, size_t width, size_t height)
{
	size_t ny = floor(sqrt(nStations)), nx = (nStations+ny-1) / ny;
	Logger::Info << "Storing " << filename << " (" << nStations << " ant, " << nx << " x " << ny << ")\n";
	ao::uvector<double> img(nx*ny * width*height, 0.0);
	for(size_t ant=0; ant!=nStations; ++ant)
	{
		size_t xCorner = (ant%nx)*width, yCorner = (ant/nx)*height;
		for(size_t y=0; y!=height; ++y)
		{
			for(size_t x=0; x!=width; ++x)
			{
				std::complex<float> e1, e2;
				Matrix2x2::EigenValues(buffer+(width*(ant*height + y) + x)*4, e1, e2);
				double val = std::max(std::abs(e1), std::abs(e2));
				img[(yCorner + y)*width*nx + x + xCorner] = val;
			}
		}
	}
	FitsWriter writer;
	writer.SetImageDimensions(nx*width, ny*height);
	writer.Write(filename, img.data());
}

void ATermBase::StoreATermsReal(const std::string& filename, const std::complex<float>* buffer, size_t nStations, size_t width, size_t height)
{
	size_t ny = floor(sqrt(nStations)), nx = (nStations+ny-1) / ny;
	Logger::Info << "Storing " << filename << " (" << nStations << " ant, " << nx << " x " << ny << ")\n";
	ao::uvector<double> img(nx*ny * width*height, 0.0);
	for(size_t ant=0; ant!=nStations; ++ant)
	{
		size_t xCorner = (ant%nx)*width, yCorner = (ant/nx)*height;
		for(size_t y=0; y!=height; ++y)
		{
			for(size_t x=0; x!=width; ++x)
			{
				std::complex<float> xx = (buffer+(width*(ant*height + y) + x)*4)[0];
				img[(yCorner + y)*width*nx + x + xCorner] = xx.real();
			}
		}
	}
	FitsWriter writer;
	writer.SetImageDimensions(nx*width, ny*height);
	writer.Write(filename, img.data());
}
