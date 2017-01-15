#include <boost/test/unit_test.hpp>

#include "../deconvolution/simpleclean.h"

#include "../uvector.h"

#include <random>

BOOST_AUTO_TEST_SUITE(clean_algorithms)

const size_t nRepeats = 3; /* This should be set to 100 to assert the performance */

#if defined __AVX__ && !defined FORCE_NON_AVX
struct CleanTestFixture
{
	size_t x, y;
	ao::uvector<double> img;

	CleanTestFixture() :
		x(size_t(-1)), y(size_t(-1)),
		img(16, 0)
	{	}
	void findPeak(size_t width=4, size_t height=2, size_t ystart=0, size_t yend=2)
	{
		SimpleClean::FindPeakAVX(img.data(), width, height, x, y, true, ystart, yend, 0, 0);
	}
};
#endif

struct NoiseFixture
{
	NoiseFixture() :
		n(2048),
		psf(n * n, 0.0),
		img(n * n),
		normal_dist(0.0, 1.0)
	{
		mt.seed(42);
		for(size_t i=0; i!=n*n; ++i)
		{
			img[i] = normal_dist(mt);
		}
	}
	
	size_t n;
	ao::uvector<double> psf, img;
	std::mt19937 mt;
	std::normal_distribution<double> normal_dist;
};

#if defined __AVX__ && !defined FORCE_NON_AVX
BOOST_AUTO_TEST_CASE( findPeakAVX1 )
{
	CleanTestFixture f;
	f.img[0] = 1;
	f.findPeak();
	BOOST_CHECK_EQUAL(f.x, 0);
	BOOST_CHECK_EQUAL(f.y, 0);
}

BOOST_AUTO_TEST_CASE( findPeakAVX2 )
{
	CleanTestFixture f;
	f.img[0] = 1;
	f.img[1] = 2;
	f.findPeak();
	BOOST_CHECK_EQUAL(f.x, 1);
	BOOST_CHECK_EQUAL(f.y, 0);
}

BOOST_AUTO_TEST_CASE( findPeakAVX3 )
{
	CleanTestFixture f;
	f.img[0] = 1;
	f.img[1] = 2;
	f.img[4] = 3;
	f.findPeak();
	BOOST_CHECK_EQUAL(f.x, 0);
	BOOST_CHECK_EQUAL(f.y, 1);
}

BOOST_AUTO_TEST_CASE( findPeakAVX4 )
{
	CleanTestFixture f;
	f.img[0] = 1;
	f.img[1] = 2;
	f.img[4] = 3;
	f.img[7] = 4;
	f.findPeak();
	BOOST_CHECK_EQUAL(f.x, 3);
	BOOST_CHECK_EQUAL(f.y, 1);
}

BOOST_AUTO_TEST_CASE( findPeakAVX5 )
{
	CleanTestFixture f;
	f.img[0] = 1;
	f.img[1] = 2;
	f.img[4] = 3;
	f.img[7] = 4;
	f.img[15] = 6;
	f.findPeak(4,4,0,4);
	BOOST_CHECK_EQUAL(f.x, 3);
	BOOST_CHECK_EQUAL(f.y, 3);
}

BOOST_AUTO_TEST_CASE( findPeakAVX6 )
{
	CleanTestFixture f;
	f.img[0] = 1;
	f.img[1] = 2;
	f.img[4] = 3;
	f.img[7] = 4;
	f.img[15] = 6;
	f.img[14] = 5;
	f.findPeak(3,5,0,5);
	BOOST_CHECK_EQUAL(f.x, 2);
	BOOST_CHECK_EQUAL(f.y, 4);
}
#endif

BOOST_AUTO_TEST_CASE( findPeakPerformance )
{
	NoiseFixture f;
	for(size_t repeat=0; repeat!=nRepeats; ++repeat)
	{
		size_t x, y;
		SimpleClean::FindPeak(f.img.data(), f.n, f.n, x, y, true, 0, f.n/2, 0.0);
	}
	BOOST_CHECK(true);
}

BOOST_AUTO_TEST_CASE( findPeakSimplePerformance )
{
	NoiseFixture f;
	for(size_t repeat=0; repeat!=nRepeats; ++repeat)
	{
		size_t x, y;
		SimpleClean::FindPeakSimple(f.img.data(), f.n, f.n, x, y, true, 0, f.n/2, 0, 0);
	}
	BOOST_CHECK(true);
}
	
#if defined __AVX__ && !defined FORCE_NON_AVX
BOOST_AUTO_TEST_CASE( findPeakAVXPerformance )
{
	NoiseFixture f;
	for(size_t repeat=0; repeat!=nRepeats; ++repeat)
	{
		size_t x, y;
		SimpleClean::FindPeakAVX(f.img.data(), f.n, f.n, x, y, true, 0, f.n/2, 0, 0);
	}
	BOOST_CHECK(true);
}
#endif
	
BOOST_AUTO_TEST_CASE( partialSubtractImagePerformance )
{
	NoiseFixture f;
	size_t x=f.n/2, y=f.n/2;
	for(size_t repeat=0; repeat!=nRepeats; ++repeat)
		SimpleClean::PartialSubtractImage(f.img.data(), f.n, f.n, f.psf.data(), f.n, f.n, x, y, 0.5, 0, f.n/2);
	BOOST_CHECK(true);
}
	
#if defined __AVX__ && !defined FORCE_NON_AVX
BOOST_AUTO_TEST_CASE( partialSubtractImageAVXPerformance )
{
	NoiseFixture f;
	size_t x=f.n/2, y=f.n/2;
	for(size_t repeat=0; repeat!=nRepeats; ++repeat)
		SimpleClean::PartialSubtractImageAVX(f.img.data(), f.n, f.n, f.psf.data(), f.n, f.n, x, y, 0.5, 0, f.n/2);
	BOOST_CHECK(true);
}
#endif

BOOST_AUTO_TEST_SUITE_END()
