#include <boost/test/unit_test.hpp>

#include "../gaussianfitter.h"
#include "../image.h"
#include "../modelrenderer.h"

BOOST_AUTO_TEST_SUITE(gaussian_fitter)

BOOST_AUTO_TEST_CASE( fit )
{
	const size_t
		width = 512,
		height = 512;
	ImageBufferAllocator allocator;
	Image
		model(width, height, 0.0, allocator),
		restored(width, height, 0.0, allocator);
	model[((height/2)*width) + (width/2)] = 1.0;
	long double pixelScale = 1 /*amin*/ * (M_PI/180.0/60.0);
	long double beamMaj = 20*pixelScale, beamMin = 5*pixelScale, beamPA = 1.0;
	ModelRenderer::Restore(restored.data(), model.data(), width, height, beamMaj, beamMin, beamPA, pixelScale, pixelScale);
	
	GaussianFitter fitter;
	double fitMaj, fitMin, fitPA;
	fitter.Fit2DGaussianCentred(restored.data(), 512, 512, 5.0, fitMaj, fitMin, fitPA, false);
	
	BOOST_CHECK_CLOSE_FRACTION(fitMaj, 20.0, 1e-4);
	BOOST_CHECK_CLOSE_FRACTION(fitMin, 5.0, 1e-4);
	BOOST_CHECK_CLOSE_FRACTION(fitPA, 1.0, 1e-4);
}

BOOST_AUTO_TEST_SUITE_END()
