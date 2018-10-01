#ifndef ATERM_STUB
#define ATERM_STUB

#include "atermbase.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

class ATermStub : public ATermBase
{
public:
	ATermStub(casacore::MeasurementSet&, size_t, size_t, double, double, double, double, bool)
	{
		throw std::runtime_error("ATerm not implemented -- did you forget to turn specific beam options on during the compilation?");
	}
	virtual bool Calculate(std::complex<float>*, double, double) final override { return false; };
};

#endif

