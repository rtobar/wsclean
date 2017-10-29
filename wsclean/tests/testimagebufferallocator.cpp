#include <boost/test/unit_test.hpp>

#include "../wsclean/imagebufferallocator.h"

BOOST_AUTO_TEST_SUITE(imagebufferallocator)

BOOST_AUTO_TEST_CASE( allocation )
{
	ImageBufferAllocator a;
	double* x = a.Allocate(1000);
	double* y = a.Allocate(2000);
	double* z = a.Allocate(1000);
	BOOST_CHECK_NE(x, y);
	BOOST_CHECK_NE(y, z);
	a.Free(z);
	a.Free(x);
	a.Free(y);
}

BOOST_AUTO_TEST_CASE( autoptr )
{
	ImageBufferAllocator a;
	ImageBufferAllocator::Ptr x, y, z;
	a.Allocate(1000, x);
	a.Allocate(2000, y);
	a.Allocate(1000, z);
	BOOST_CHECK_NE(x.data(), y.data());
	BOOST_CHECK_NE(y.data(), z.data());
}

BOOST_AUTO_TEST_CASE( zerosized )
{
	ImageBufferAllocator a;
	ImageBufferAllocator::Ptr w, x, y, z;
	a.Allocate(1000, w);
	a.Allocate(0, x);
	a.Allocate(0, y);
	a.Allocate(0, z);
	BOOST_CHECK_NE(w.data(), x.data());
	BOOST_CHECK_NE(x.data(), y.data());
	BOOST_CHECK_NE(y.data(), z.data());
}

BOOST_AUTO_TEST_CASE( combination )
{
	ImageBufferAllocator a;
	std::complex<double>* w = a.AllocateComplex(1000);
	std::complex<double>* x = a.AllocateComplex(0);
	a.Free(w);
	double* f = a.Allocate(1000);
	BOOST_CHECK_NE(f, (double*) x);
	a.Free(f);
	a.Free(x);
	double* g = a.Allocate(1000);
	std::complex<double>* y = a.AllocateComplex(0);
	std::complex<double>* z = a.AllocateComplex(0);
	BOOST_CHECK_NE(g, (double*)(y));
	BOOST_CHECK_NE(y, z);
	a.Free(g);
	a.Free(y);
	a.Free(z);
}

BOOST_AUTO_TEST_SUITE_END()

