#ifndef ATERM_STUB
#define ATERM_STUB

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

class ATermStub
{
public:
	ATermStub(casacore::MeasurementSet&, size_t, size_t, double, double, double, double, bool)
	{
		throw std::runtime_error("ATerm not implemented -- did you forget to turn specific beam options on during the compilation?");
	}
	void Calculate(std::complex<float>*, double, double) { };
};

#endif

