#include <boost/test/unit_test.hpp>

#include "../gaussianfitter.h"
#include "../image.h"
#include "../modelrenderer.h"

#include "../model/model.h"
#include "../model/powerlawsed.h"

BOOST_AUTO_TEST_SUITE(gaussian_fitter)

BOOST_AUTO_TEST_CASE( conversions )
{
	const long double sigmaToBeam = 2.0L * sqrtl(2.0L * logl(2.0L));
	double gMax, gMin, gAngle;
	double outMax, outMin, outAngle;
	double sxsx, sxsy, sysy;
	
	gMax = sigmaToBeam; // sx = 1
	gMin = sigmaToBeam;
	gAngle = 0.0;
	GaussianFitter::ToCovariance(gMax, gMin, gAngle, sxsx, sxsy, sysy);
	
	BOOST_CHECK_CLOSE_FRACTION(sxsx, 1.0, 1e-4);  // a = 1 / (2 sx ^2) = 0.5
	BOOST_CHECK_CLOSE_FRACTION(sxsy, 0.0, 1e-4);
	BOOST_CHECK_CLOSE_FRACTION(sysy, 1.0, 1e-4);
	
	GaussianFitter::FromCovariance(sxsx, sxsy, sysy, outMax, outMin, outAngle);
	BOOST_CHECK_CLOSE_FRACTION(outMax, sigmaToBeam, 1e-4);
	BOOST_CHECK_CLOSE_FRACTION(outMin, sigmaToBeam, 1e-4);
	
	gMax = sigmaToBeam;     // sx = 1
	gMin = sigmaToBeam*0.5; // sy = 0.5
	gAngle = 0.0;
	GaussianFitter::ToCovariance(gMax, gMin, gAngle, sxsx, sxsy, sysy);
	
	BOOST_CHECK_CLOSE_FRACTION(sxsx, 1.0, 1e-4);
	BOOST_CHECK_CLOSE_FRACTION(sxsy, 0.0, 1e-4);
	BOOST_CHECK_CLOSE_FRACTION(sysy, 0.25, 1e-4);
	
	GaussianFitter::ToCovariance(1.0, 0.5, gAngle, sxsx, sxsy, sysy);
	GaussianFitter::FromCovariance(sxsx, sxsy, sysy, outMax, outMin, outAngle);
	BOOST_CHECK_CLOSE_FRACTION(outMax, 1.0, 1e-4);
	BOOST_CHECK_CLOSE_FRACTION(outMin, 0.5, 1e-4);
	outAngle += 2.0*M_PI;
	while(outAngle > 0.5*M_PI) outAngle -= M_PI;
	BOOST_CHECK_CLOSE_FRACTION(outAngle, 0.0, 1e-4);
	
	for(double x=0.0; x<2.0; x+=0.1)
	{
		gAngle = x * M_PI;
		GaussianFitter::ToCovariance(1.0, 0.5, gAngle, sxsx, sxsy, sysy);
		GaussianFitter::FromCovariance(sxsx, sxsy, sysy, outMax, outMin, outAngle);
		BOOST_CHECK_CLOSE_FRACTION(outMax, 1.0, 1e-4);
		BOOST_CHECK_CLOSE_FRACTION(outMin, 0.5, 1e-4);
		outAngle += 2.0*M_PI;
		while(outAngle > gAngle + 0.5*M_PI) outAngle -= M_PI;
		BOOST_CHECK_CLOSE_FRACTION(outAngle*180.0/M_PI, gAngle*180.0/M_PI, 1e-4);
	}
}

BOOST_AUTO_TEST_CASE( fit )
{
	for(size_t beamPAindex=0;  beamPAindex!=10; ++beamPAindex)
	{
		const size_t
			width = 512,
			height = 512;
		ImageBufferAllocator allocator;
		Image
			model(width, height, 0.0, allocator),
			restored(width, height, 0.0, allocator);
		model[((height/2)*width) + (width/2)] = 1.0;
		long double
			pixelScale = 1 /*amin*/ * (M_PI/180.0/60.0),
			beamMaj = 20*pixelScale, beamMin = 5*pixelScale,
			beamPA = beamPAindex*M_PI / 10.0;
		ModelRenderer::Restore(restored.data(), model.data(), width, height, beamMaj, beamMin, beamPA, pixelScale, pixelScale);

		GaussianFitter fitter;
		double fitMaj, fitMin, fitPA;
		fitter.Fit2DGaussianCentred(restored.data(), width, height, 5.0, fitMaj, fitMin, fitPA, 10.0, false);
		fitPA = fmod((fitPA + 2.0 * M_PI), M_PI);

		BOOST_CHECK_CLOSE_FRACTION(fitMaj, 20.0, 1e-3);
		BOOST_CHECK_CLOSE_FRACTION(fitMin, 5.0, 1e-3);
		BOOST_CHECK_CLOSE_FRACTION(fitPA, beamPA, 1e-3);
	}
}

BOOST_AUTO_TEST_CASE( fit_with_bad_initial_value )
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
	fitter.Fit2DGaussianCentred(restored.data(), width, height, estimatedBeamPx, fitMajor, fitMinor, fitPA, 10.0, false);
	
	BOOST_CHECK_CLOSE_FRACTION(fitMajor, 4.0, 1e-4);
	BOOST_CHECK_CLOSE_FRACTION(fitMinor, 4.0, 1e-4);
}

BOOST_AUTO_TEST_CASE( fit_circular )
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
	double fitMajor = estimatedBeamPx;
	fitter.Fit2DCircularGaussianCentred(restored.data(), width, height, fitMajor);
	
	BOOST_CHECK_CLOSE_FRACTION(fitMajor, 4.0, 1e-4);
}

BOOST_AUTO_TEST_CASE( fit_small_beam )
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
	long double beamMaj = 4*pixelScale, beamMin = 0.5*pixelScale, beamPA = 0.0;
	long double estimatedBeamPx = 1.0; // this is on purpose way off
	ModelRenderer renderer(0.0, 0.0, pixelScale, pixelScale);
	renderer.Restore(restored.data(), width, height, model, beamMaj, beamMin, beamPA, 100e6, 200e6, Polarization::StokesI);
	
	GaussianFitter fitter;
	double fitMajor = estimatedBeamPx, fitMinor = estimatedBeamPx, fitPA = 0.0;
	fitter.Fit2DGaussianCentred(restored.data(), width, height, estimatedBeamPx, fitMajor, fitMinor, fitPA, 10.0, false);
	
	BOOST_CHECK_CLOSE_FRACTION(fitMinor, 0.5, 1e-4);
}

BOOST_AUTO_TEST_SUITE_END()
