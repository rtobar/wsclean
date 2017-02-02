#include "rmsimage.h"
#include "modelrenderer.h"

void RMSImage::Make(Image& rmsOutput, const Image& inputImage, double windowSize, long double beamMaj, long double beamMin, long double beamPA, long double pixelScaleL, long double pixelScaleM)
{
	Image image(inputImage);
	rmsOutput = Image(image.Width(), image.Height(), 0.0, image.Allocator());
	
	for(double& val : image)
		val *= val;
	
	ModelRenderer::Restore(rmsOutput.data(), image.data(), image.Width(), image.Height(), beamMaj*windowSize, beamMin*windowSize, beamPA, pixelScaleL, pixelScaleM);
	
	double s = sqrt(2.0 * M_PI);
	const long double sigmaMaj = beamMaj / (2.0L * sqrtl(2.0L * logl(2.0L)));
	const long double sigmaMin = beamMin / (2.0L * sqrtl(2.0L * logl(2.0L)));
	const double norm = 1.0 / (s * sigmaMaj/pixelScaleL * windowSize * s * sigmaMin/pixelScaleL * windowSize);
	for(double& val : rmsOutput)
		val = sqrt(val * norm);
}

void RMSImage::SlidingMinimum(Image& output, const Image& input, size_t windowSize)
{
	const size_t width = input.Width();
	output = Image(width, input.Height(), input.Allocator());
	Image temp(output);
	for(size_t y=0; y!=input.Height(); ++y)
	{
		double *outRowptr = &temp[y*width];
		const double* inRowptr = &input[y*width];
		for(size_t x=0; x!=width; ++x)
		{
			size_t left = std::max(x, windowSize/2) - windowSize/2;
			size_t right = std::min(x, width-windowSize/2) + windowSize/2;
			outRowptr[x] = *std::min_element(inRowptr + left, inRowptr + right);
		}
	}
	for(size_t x=0; x!=width; ++x)
	{
		ao::uvector<double> vals;
		for(size_t y=0; y!=input.Height(); ++y)
		{
			size_t top = std::max(y, windowSize/2) - windowSize/2;
			size_t bottom = std::min(y, input.Height()-windowSize/2) + windowSize/2;
			vals.clear();
			for(size_t winY=top; winY!=bottom; ++winY)
				vals.push_back(temp[winY*width + x]);
			output[y*width + x] = *std::min_element(vals.begin(), vals.end());
		}
	}
}

void RMSImage::MakeWithNegativityLimit(Image& rmsOutput, const Image& inputImage, double windowSize, long double beamMaj, long double beamMin, long double beamPA, long double pixelScaleL, long double pixelScaleM)
{
	Make(rmsOutput, inputImage, windowSize, beamMaj, beamMin, beamPA, pixelScaleL, pixelScaleM);
	Image slidingMinimum(inputImage.Width(), inputImage.Height(), inputImage.Allocator());
	double beamInPixels = std::max(beamMaj / pixelScaleL, 1.0L);
	SlidingMinimum(slidingMinimum, inputImage, windowSize * beamInPixels);
	for(size_t i=0; i!=rmsOutput.size(); ++i)
	{
		rmsOutput[i] = std::max(rmsOutput[i], std::abs(slidingMinimum[i]) * (1.5/5.0) );
	}
}
