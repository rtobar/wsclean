#include "rmsimage.h"
#include "modelrenderer.h"

void RMSImage::Make(Image& rmsOutput, const Image& inputImage, long double beamMaj, long double beamMin, long double beamPA, long double pixelScaleL, long double pixelScaleM)
{
	const double windowSize = 25.0;
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
