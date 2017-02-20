#include "../polynomialchannelfitter.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

BOOST_AUTO_TEST_SUITE(polynomial_channel_fitter)

BOOST_AUTO_TEST_CASE( fit_line )
{
	PolynomialChannelFitter fitter;
	fitter.AddChannel(1.0, 2.0);
	fitter.AddChannel(2.0, 3.0);
	
	ao::uvector<double> terms;
	fitter.AddDataPoint(0, 0.0);
	fitter.AddDataPoint(0, 1.0);
	fitter.AddDataPoint(1, 1.0);
	fitter.AddDataPoint(1, 2.0);
	// This should fit a line through (1.5, 0.5) and (2.5, 1.5),
	// hence y = x - 1
	fitter.Fit(terms, 2);
	
	BOOST_CHECK_CLOSE_FRACTION(terms[0], -1.0, 1e-3);
	BOOST_CHECK_CLOSE_FRACTION(terms[1], 1.0, 1e-3);
}

BOOST_AUTO_TEST_CASE( fit_squared_function )
{
	PolynomialChannelFitter fitter;
	fitter.AddChannel(1.0, 2.0);
	fitter.AddChannel(2.0, 3.0);
	fitter.AddChannel(3.0, 4.0);
	
	// Fit to function x^2
	ao::uvector<double> terms;
	fitter.AddDataPoint(0, 7.0/3.0); // int_{x=1}^2 x^2 / 1 = 1/3 (2^3 - 1^3) = 7/3
	fitter.AddDataPoint(1, 19.0/3.0); // int_{x=2}^2 x^2 / 1 = 1/3 (3^3 - 2^3) = 19/3
	fitter.AddDataPoint(2, 37.0/3.0); // int_{x=3}^4 x^2 / 1 = 1/3 (4^3 - 3^3) = 37/3
	fitter.Fit(terms, 3);
	
	BOOST_CHECK_LT(terms[0], 1e-3);
	BOOST_CHECK_LT(terms[1], 1e-3);
	BOOST_CHECK_CLOSE_FRACTION(terms[2], 1.0, 1e-3);
}

BOOST_AUTO_TEST_CASE( fit_irregular_channels )
{
	PolynomialChannelFitter fitter;
	fitter.AddChannel(0.0, 2.0);
	fitter.AddChannel(2.0, 3.0);
	fitter.AddChannel(3.0, 6.0);
	
	// Fit to function x^2
	ao::uvector<double> terms;
	fitter.AddDataPoint(0, 8.0/6.0); // int_{x=0}^2 x^2 / 2 = 1/3 (2^3 - 0^3) = 8/6
	fitter.AddDataPoint(1, 19.0/3.0); // int_{x=2}^3 x^2 / 1 = 1/3 (3^3 - 2^3) = 19/3
	fitter.AddDataPoint(2, 189.0/9.0); // int_{x=3}^6 x^2 / 3 = 1/3 (6^3 - 3^3) = 189/9
	fitter.Fit(terms, 3);
	
	BOOST_CHECK_LT(terms[0], 1e-3);
	BOOST_CHECK_LT(terms[1], 1e-3);
	BOOST_CHECK_CLOSE_FRACTION(terms[2], 1.0, 1e-3);
}

BOOST_AUTO_TEST_CASE( fit_squared_missing_bandwidth )
{
	PolynomialChannelFitter fitter;
	fitter.AddChannel(0.0, 2.0);
	fitter.AddChannel(10.0, 12.0);
	fitter.AddChannel(20.0, 22.0);
	
	// Fit to function x^2
	ao::uvector<double> terms;
	fitter.AddDataPoint(0,    8.0/6.0); // int_{x=0}^2 x^2 / 2 = 1/2 (2^3 - 0^3) = 8/6
	fitter.AddDataPoint(1,  728.0/6.0); // int_{x=10}^12 x^2 / 1 = 1/2 (12^3 - 10^3) = 728/6
	fitter.AddDataPoint(2, 2648.0/6.0); // int_{x=3}^6 x^2 / 3 = 1/2 (22^3 - 20^3) = 2648/6
	fitter.Fit(terms, 3);
	
	BOOST_CHECK_LT(terms[0], 1e-3);
	BOOST_CHECK_LT(terms[1], 1e-3);
	BOOST_CHECK_CLOSE_FRACTION(terms[2], 1.0, 1e-3);
}

double integrate(double xa, double xb, double a, double b, double c, double d)
{
	return
		((a*xb + b*xb*xb/2.0 + c*xb*xb*xb/3.0 + d*xb*xb*xb*xb/4.0) -
		 (a*xa + b*xa*xa/2.0 + c*xa*xa*xa/3.0 + d*xa*xa*xa*xa/4.0)) /
		       (xb - xa);
}

BOOST_AUTO_TEST_CASE( fit_2nd_order )
{
	double a = 0.5, b = 2.0, c = 1.2, d = 0.0;
	PolynomialChannelFitter fitter;
	fitter.AddChannel(0.0, 1.0);
	fitter.AddChannel(3.0, 4.0);
	fitter.AddChannel(6.0, 7.0);
	
	// Fit to function x^2
	ao::uvector<double> terms;
	fitter.AddDataPoint(0, integrate(0.0, 1.0, a, b, c, d));
	fitter.AddDataPoint(1, integrate(3.0, 4.0, a, b, c, d));
	fitter.AddDataPoint(2, integrate(6.0, 7.0, a, b, c, d));
	fitter.Fit(terms, 3);
	
	BOOST_CHECK_CLOSE_FRACTION(terms[0], a, 1e-3);
	BOOST_CHECK_CLOSE_FRACTION(terms[1], b, 1e-3);
	BOOST_CHECK_CLOSE_FRACTION(terms[2], c, 1e-3);
}

BOOST_AUTO_TEST_CASE( fit_higher_order )
{
	double a = 0.5, b = 2.0, c = 1.2, d = 0.1;
	PolynomialChannelFitter fitter;
	fitter.AddChannel(1.0, 1.1);
	fitter.AddChannel(1.3, 1.4);
	fitter.AddChannel(1.4, 1.5);
	fitter.AddChannel(1.8, 1.9);
	fitter.AddChannel(2.0, 2.1);
	
	ao::uvector<double> terms;
	fitter.AddDataPoint(0, integrate(1.0, 1.1, a, b, c, d));
	fitter.AddDataPoint(1, integrate(1.3, 1.4, a, b, c, d));
	fitter.AddDataPoint(2, integrate(1.4, 1.5, a, b, c, d));
	fitter.AddDataPoint(3, integrate(1.8, 1.9, a, b, c, d));
	fitter.AddDataPoint(4, integrate(2.0, 2.1, a, b, c, d));
	fitter.Fit(terms, 4);
	
	BOOST_CHECK_CLOSE_FRACTION(terms[0], a, 1e-3);
	BOOST_CHECK_CLOSE_FRACTION(terms[1], b, 1e-3);
	BOOST_CHECK_CLOSE_FRACTION(terms[2], c, 1e-3);
	BOOST_CHECK_CLOSE_FRACTION(terms[3], d, 1e-3);
}

BOOST_AUTO_TEST_SUITE_END()

