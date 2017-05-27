#include <boost/test/unit_test.hpp>

#include "../gaussianfitter.h"
#include "../image.h"
#include "../modelrenderer.h"

#include "../model/model.h"
#include "../model/powerlawsed.h"

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
	fitter.Fit2DGaussianCentred(restored.data(), width, height, 5.0, fitMaj, fitMin, fitPA, false);
	
	BOOST_CHECK_CLOSE_FRACTION(fitMaj, 20.0, 1e-4);
	BOOST_CHECK_CLOSE_FRACTION(fitMin, 5.0, 1e-4);
	BOOST_CHECK_CLOSE_FRACTION(fitPA, 1.0, 1e-4);
}

BOOST_AUTO_TEST_CASE( fit_oversize )
{
	const size_t
		width = 64,
		height = 64;
	ImageBufferAllocator allocator;
	Image restored(width, height, 0.0, allocator);
	PowerLawSED sed(150.0e6, 1.0);
	ModelComponent component;
	component.SetPosDec(0.0);
	component.SetPosRA(0.0);
	component.SetSED(sed);
	ModelSource source;
	source.AddComponent(component);
	Model model;
	model.AddSource(source);
	long double pixelScale = 1 /*amin*/ * (M_PI/180.0/60.0);
	long double beamMaj = 4*pixelScale, beamMin = 4*pixelScale, beamPA = 0.0;
	long double estimatedBeamPx = 1.0; // this is on purpose way off
	ModelRenderer renderer(0.0, 0.0, pixelScale, pixelScale);
	renderer.Restore(restored.data(), width, height, model, beamMaj, beamMin, beamPA, 100e6, 200e6, Polarization::StokesI);
	
	GaussianFitter fitter;
	double fitMajor, fitMinor, fitPA;
	fitter.Fit2DGaussianCentred(restored.data(), width, height, estimatedBeamPx, fitMajor, fitMinor, fitPA, false);
	
	BOOST_CHECK_CLOSE_FRACTION(fitMajor, 4.0, 1e-4);
	BOOST_CHECK_CLOSE_FRACTION(fitMinor, 4.0, 1e-4);
}

BOOST_AUTO_TEST_SUITE_END()
