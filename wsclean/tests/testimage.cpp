#include <boost/test/unit_test.hpp>

#include "../image.h"

#include "../uvector.h"

#include <random>

BOOST_AUTO_TEST_SUITE(image_operations)

BOOST_AUTO_TEST_CASE( median_empty )
{
	BOOST_CHECK_EQUAL(Image::Median(nullptr, 0), 0.0);
}

BOOST_AUTO_TEST_CASE( median_single )
{
	ao::uvector<double> arr(1, 1.0);
	BOOST_CHECK_EQUAL(Image::Median(arr.data(), arr.size()), 1.0);
	
	arr[0] = std::numeric_limits<double>::quiet_NaN();
	Image::Median(arr.data(), arr.size()); // undefined -- just make sure it doesn't crash
}

BOOST_AUTO_TEST_CASE( median_two_elements )
{
	ao::uvector<double> arr1(2, 1.0);
	BOOST_CHECK_CLOSE_FRACTION(Image::Median(arr1.data(), arr1.size()), 1.0, 1e-5);
	
	ao::uvector<double> arr2(2, 0.0);
	arr2[1] = 2.0;
	BOOST_CHECK_CLOSE_FRACTION(Image::Median(arr2.data(), arr2.size()), 1.0, 1e-5);
	
	ao::uvector<double> arr3(2, 1.0);
	arr3[1] = -1.0;
	BOOST_CHECK_CLOSE_FRACTION(Image::Median(arr3.data(), arr3.size()), 0.0, 1e-5);

	ao::uvector<double> arr4(2, 13.0);
	arr3[1] = std::numeric_limits<double>::quiet_NaN();
	BOOST_CHECK_CLOSE_FRACTION(Image::Median(arr4.data(), arr4.size()), 13.0, 1e-5);
}

BOOST_AUTO_TEST_CASE( median_three_elements )
{
	ao::uvector<double> arr(3, 1.0);
	BOOST_CHECK_CLOSE_FRACTION(Image::Median(arr.data(), arr.size()), 1.0, 1e-5);
	
	arr[0] = 0.0; arr[1] = 1.0; arr[2] = 2.0;
	BOOST_CHECK_CLOSE_FRACTION(Image::Median(arr.data(), arr.size()), 1.0, 1e-5);

	arr[0] = 3.0; arr[1] = -3.0; arr[2] = 2.0;
	BOOST_CHECK_CLOSE_FRACTION(Image::Median(arr.data(), arr.size()), 2.0, 1e-5);
	
	arr[1] = std::numeric_limits<double>::quiet_NaN();
	BOOST_CHECK_CLOSE_FRACTION(Image::Median(arr.data(), arr.size()), 2.5, 1e-5);
	
	arr[0] = std::numeric_limits<double>::quiet_NaN();
	BOOST_CHECK_CLOSE_FRACTION(Image::Median(arr.data(), arr.size()), 2.0, 1e-5);
}

BOOST_AUTO_TEST_CASE( mad_empty )
{
	BOOST_CHECK_EQUAL(Image::MAD(nullptr, 0), 0.0);
}

BOOST_AUTO_TEST_CASE( mad_single )
{
	ao::uvector<double> arr(1, 1.0);
	BOOST_CHECK_EQUAL(Image::MAD(arr.data(), arr.size()), 0.0);
}

BOOST_AUTO_TEST_CASE( mad_two_elements )
{
	ao::uvector<double> arr(2, 1.0);
	BOOST_CHECK_EQUAL(Image::MAD(arr.data(), arr.size()), 0.0);
	
	arr[0] = 0.0; arr[1] = 2.0;
	BOOST_CHECK_CLOSE_FRACTION(Image::Median(arr.data(), arr.size()), 1.0, 1e-5);
	
	arr[0] = 1.0; arr[1] = -1.0;
	BOOST_CHECK_CLOSE_FRACTION(Image::Median(arr.data(), arr.size()), 0.0, 1e-5);

	arr[0] = 13.0; arr[1] = std::numeric_limits<double>::quiet_NaN();
	BOOST_CHECK_CLOSE_FRACTION(Image::Median(arr.data(), arr.size()), 13.0, 1e-5);
}

BOOST_AUTO_TEST_CASE( mad_three_elements )
{
	ao::uvector<double> arr(3, 1.0);
	BOOST_CHECK_CLOSE_FRACTION(Image::MAD(arr.data(), arr.size()), 0.0, 1e-5);
	
	arr[0] = 0.0; arr[1] = 1.0; arr[2] = 2.0;
	BOOST_CHECK_CLOSE_FRACTION(Image::MAD(arr.data(), arr.size()), 1.0, 1e-5);

	arr[0] = 3.0; arr[1] = -3.0; arr[2] = 2.0;
	BOOST_CHECK_CLOSE_FRACTION(Image::MAD(arr.data(), arr.size()), 1.0, 1e-5);
	
	arr[1] = std::numeric_limits<double>::quiet_NaN();
	BOOST_CHECK_CLOSE_FRACTION(Image::MAD(arr.data(), arr.size()), 0.5, 1e-5);
	
	arr[0] = std::numeric_limits<double>::quiet_NaN();
	BOOST_CHECK_CLOSE_FRACTION(Image::MAD(arr.data(), arr.size()), 0.0, 1e-5);
}

BOOST_AUTO_TEST_CASE( stddev_from_mad )
{
	std::mt19937 rnd;
	std::normal_distribution<double> dist(1.0, 5.0);
	ao::uvector<double> data(10000);
	for(size_t i=0; i!=data.size(); ++i)
		data[i] = dist(rnd);
	BOOST_CHECK_CLOSE_FRACTION(Image::StdDevFromMAD(data.data(), data.size()), 5.0, 0.05);
}

BOOST_AUTO_TEST_SUITE_END()
