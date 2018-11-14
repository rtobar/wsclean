#ifndef ATERM_STUB
#define ATERM_STUB

#include "atermbeam.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

class ATermStub : public ATermBeam
{
public:
	ATermStub(casacore::MeasurementSet&, size_t /*width*/, size_t /*height*/, double /*dl*/, double /*dm*/, double /*phaseCentreDL*/, double /*phaseCentreDM*/, const std::string& /*dataColumnName*/)
	{
		throw std::runtime_error("ATerm not implemented -- did you forget to turn specific beam options on during the compilation?");
	}
	virtual bool calculateBeam(std::complex<float>* /*buffer*/, double /*time*/, double /*frequency*/)
	{
		throw std::runtime_error("ATerm not implemented -- did you forget to turn specific beam options on during the compilation?");
	}
	
	void SetUseDifferentialBeam(bool) { }
	void SetUseChannelFrequency(bool) { }
};

#endif

