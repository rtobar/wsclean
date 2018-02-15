#include <boost/test/unit_test.hpp>

#include "../matrix2x2.h"

BOOST_AUTO_TEST_SUITE(matrix2x2)

BOOST_AUTO_TEST_CASE( eigenvalue1 )
{
	double unit[4] = { 1.0, 0.0, 0.0, 1.0 };
	double e1, e2;
	Matrix2x2::EigenValues(unit, e1, e2);
	BOOST_CHECK_CLOSE(e1, 1.0, 1e-6);
	BOOST_CHECK_CLOSE(e2, 1.0, 1e-6);
}

BOOST_AUTO_TEST_CASE( eigenvalue2 )
{
	double unit[4] = { 0.0, 1.0, -2.0, -3.0 };
	double e1, e2;
	Matrix2x2::EigenValues(unit, e1, e2);
	if(e1 < e2) std::swap(e1, e2);
	BOOST_CHECK_CLOSE(e1, -1.0, 1e-6);
	BOOST_CHECK_CLOSE(e2, -2.0, 1e-6);
}

BOOST_AUTO_TEST_CASE( eigenvalue3 )
{
	double unit[4] = { 0.0, -2.0, 1.0, -3.0 };
	double e1, e2;
	Matrix2x2::EigenValues(unit, e1, e2);
	if(e1 < e2) std::swap(e1, e2);
	BOOST_CHECK_CLOSE(e1, -1.0, 1e-6);
	BOOST_CHECK_CLOSE(e2, -2.0, 1e-6);
}

BOOST_AUTO_TEST_CASE( eigenvalue4 )
{
	double unit[4] = { 0.0, 1.0, -1.0, 0.0 };
	double e1, e2;
	Matrix2x2::EigenValues(unit, e1, e2);
	if(e1 < e2) std::swap(e1, e2);
	BOOST_CHECK(!std::isfinite(e1));
	BOOST_CHECK(!std::isfinite(e2));
}

BOOST_AUTO_TEST_CASE( eigenvector2 )
{
	double unit[4] = { 0.0, 1.0, -2.0, -3.0 };
	double e1, e2, vec1[2], vec2[2];
	Matrix2x2::EigenValuesAndVectors(unit, e1, e2, vec1, vec2);
	if(e1 < e2) {
		std::swap(e1, e2);
		std::swap(vec1, vec2);
	}
	BOOST_CHECK_CLOSE(e1, -1.0, 1e-6);
	BOOST_CHECK_CLOSE(vec1[0]/vec1[1], -1.0, 1e-6); // vec1 = c [-1, 1]
	BOOST_CHECK_CLOSE(e2, -2.0, 1e-6);
	BOOST_CHECK_CLOSE(vec2[0]/vec2[1], -0.5, 1e-6); // vec2 = c [-1, 2]
}

BOOST_AUTO_TEST_CASE( eigenvector3 )
{
	double unit[4] = { 0.0, -2.0, 1.0, -3.0 };
	double e1, e2, vec1[2], vec2[2];
	Matrix2x2::EigenValuesAndVectors(unit, e1, e2, vec1, vec2);
	if(e1 < e2) {
		std::swap(e1, e2);
		std::swap(vec1, vec2);
	}
	BOOST_CHECK_CLOSE(e1, -1.0, 1e-6);
	BOOST_CHECK_CLOSE(vec1[0]/vec1[1], 2.0, 1e-6); // vec1 = c [2, 1]
	BOOST_CHECK_CLOSE(e2, -2.0, 1e-6);
	BOOST_CHECK_CLOSE(vec2[0]/vec2[1], 1.0, 1e-6); // vec2 = c [1, 1]
}

BOOST_AUTO_TEST_CASE( eigenvector4 )
{
	double unit[4] = { 1.0, 2.0, 3.0, -4.0 };
	double e1, e2, vec1[2], vec2[2];
	Matrix2x2::EigenValuesAndVectors(unit, e1, e2, vec1, vec2);
	if(e1 < e2) {
		std::swap(e1, e2);
		std::swap(vec1, vec2);
	}
	BOOST_CHECK_CLOSE(e1, 2.0, 1e-6);
	BOOST_CHECK_CLOSE(vec1[0]/vec1[1], 2.0, 1e-6); // vec1 = c [2, 1]
	BOOST_CHECK_CLOSE(e2, -5.0, 1e-6);
	BOOST_CHECK_CLOSE(vec2[1]/vec2[0], -3.0, 1e-6); // vec2 = c [-2, 6]
}

BOOST_AUTO_TEST_CASE( cholesky_real )
{
	std::complex<double> matrix[4] = {1., 2., 2., 13.};
	BOOST_CHECK(Matrix2x2::Cholesky(matrix));
	std::complex<double> answer[4] = {1., 0., 2., 3.};
	for(size_t i=0; i!=4; ++i) {
		BOOST_CHECK_CLOSE(matrix[i].real(), answer[i].real(), 1e-6);
		BOOST_CHECK_CLOSE(matrix[i].imag(), answer[i].imag(), 1e-6);
	}
}

BOOST_AUTO_TEST_CASE( cholesky_complex )
{
	std::complex<double> matrix[4] = {{1., 0.}, {2., -5.}, {2., 5.}, {38., 0.}};
	BOOST_CHECK(Matrix2x2::Cholesky(matrix));
	std::complex<double> answer[4] = {{1., 0.}, {0., 0.}, {2., 5.}, {3., 0.}};

	for(size_t i=0; i!=4; ++i) {
		BOOST_CHECK_CLOSE(matrix[i].real(), answer[i].real(), 1e-6);
		BOOST_CHECK_CLOSE(matrix[i].imag(), answer[i].imag(), 1e-6);
	}
}

BOOST_AUTO_TEST_CASE( cholesky_not_positive )
{
	std::complex<double> diag_not_positive[4] = {{0., 0.}, {0., 0.}, {0., 0.}, {1., 0.}}; // diagonal not positive
	BOOST_CHECK(!Matrix2x2::CheckedCholesky(diag_not_positive));
	std::complex<double> diag_not_real[4] = {{1., 0.}, {0., 0.}, {0., 0.}, {1., 1.}}; // diagonal not real
	BOOST_CHECK(!Matrix2x2::CheckedCholesky(diag_not_real));
	std::complex<double> not_hermitian[4] = {{1., 0.}, {1., 0.}, {2., 0.}, {1., 0.}}; // not hermitian
	BOOST_CHECK(!Matrix2x2::CheckedCholesky(not_hermitian));
}

BOOST_AUTO_TEST_SUITE_END()
