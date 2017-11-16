#ifndef ATERM_STUB
#define ATERM_STUB

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

class ATermStub
{
public:
	ATermStub(casacore::MeasurementSet&, size_t, size_t, double, double, double, double, bool) { }
	void Calculate(std::complex<float>*, double, double) { };
};

#endif

