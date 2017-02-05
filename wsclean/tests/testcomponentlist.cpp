#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include "../wsclean/imagebufferallocator.h"
#include "../deconvolution/componentlist.h"

BOOST_AUTO_TEST_SUITE(component_list)

BOOST_AUTO_TEST_CASE( adding_values )
{
	ImageBufferAllocator allocator;
	ComponentList list(512, 512, 4, 3, allocator);
	ao::uvector<double> values;
	values = { 1.0, 2.0, 3.0 };
	list.Add(256, 256, 1, values.data() );
	values = { 5.0, 6.0, 7.0 };
	list.Add(256, 256, 1, values.data() );
	values = { 8.0, 9.0, 10.0 };
	list.Add(511, 511, 0, values.data() );
	values = { 11.0, 12.0, 13.0 };
	list.Add(13, 42, 3, values.data() );
	
	list.MergeDuplicates();
	
	BOOST_CHECK_EQUAL(list.ComponentCount(0) , 1);
	BOOST_CHECK_EQUAL(list.ComponentCount(1) , 1);
	BOOST_CHECK_EQUAL(list.ComponentCount(2) , 0);
	BOOST_CHECK_EQUAL(list.ComponentCount(3) , 1);
	
	size_t x, y;
	
	list.GetComponent(0, 0, x, y, values.data());
	BOOST_CHECK_EQUAL(x, 511);
	BOOST_CHECK_EQUAL(y, 511);
	BOOST_CHECK_CLOSE_FRACTION(values[0], 8.0, 1e-5);
	BOOST_CHECK_CLOSE_FRACTION(values[1], 9.0, 1e-5);
	BOOST_CHECK_CLOSE_FRACTION(values[2], 10.0, 1e-5);
	
	list.GetComponent(1, 0, x, y, values.data());
	BOOST_CHECK_EQUAL(x, 256);
	BOOST_CHECK_EQUAL(y, 256);
	BOOST_CHECK_CLOSE_FRACTION(values[0], 6.0, 1e-5);
	BOOST_CHECK_CLOSE_FRACTION(values[1], 8.0, 1e-5);
	BOOST_CHECK_CLOSE_FRACTION(values[2], 10.0, 1e-5);

	list.GetComponent(3, 0, x, y, values.data());
	BOOST_CHECK_EQUAL(x, 13);
	BOOST_CHECK_EQUAL(y, 42);
	BOOST_CHECK_CLOSE_FRACTION(values[0], 11.0, 1e-5);
	BOOST_CHECK_CLOSE_FRACTION(values[1], 12.0, 1e-5);
	BOOST_CHECK_CLOSE_FRACTION(values[2], 13.0, 1e-5);
}

BOOST_AUTO_TEST_SUITE_END()
