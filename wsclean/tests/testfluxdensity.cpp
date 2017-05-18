
#include <boost/test/unit_test.hpp>

#include "../units/fluxdensity.h"

BOOST_AUTO_TEST_SUITE(flux_density)

BOOST_AUTO_TEST_CASE( adding_values )
{
	BOOST_CHECK_EQUAL(FluxDensity::Parse("0Jy", "") , 0.0);
	BOOST_CHECK_CLOSE_FRACTION(FluxDensity::Parse("1Jy", ""), 1.0, 1e-7);
	BOOST_CHECK_CLOSE_FRACTION(FluxDensity::Parse("1KJy", ""), 1e3, 1e-7);
	BOOST_CHECK_CLOSE_FRACTION(FluxDensity::Parse("1mJy", ""), 1e-3, 1e-7);
	BOOST_CHECK_CLOSE_FRACTION(FluxDensity::Parse("1µJy", ""), 1e-6, 1e-7);
	BOOST_CHECK_CLOSE_FRACTION(FluxDensity::Parse("1nJy", ""), 1e-9, 1e-7);
}

BOOST_AUTO_TEST_SUITE_END()
