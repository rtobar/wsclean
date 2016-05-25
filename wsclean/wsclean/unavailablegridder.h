#ifndef NOT_IMPLEMENTED_GRIDDER_H
#define NOT_IMPLEMENTED_GRIDDER_H

#include "msgridderbase.h"

#include <stdexcept>

class UnavailableGridder : public MSGridderBase
{
public:
	UnavailableGridder() { doThrow(); }
	
	virtual ~UnavailableGridder() { doThrow(); }
	
	virtual void Invert() { doThrow(); }
	
	virtual void Predict(double* real) { doThrow(); }
	
	virtual void Predict(double* real, double* imaginary) { doThrow(); }
	
	virtual double* ImageRealResult() { doThrow(); }
	
	virtual double* ImageImaginaryResult() { doThrow(); }
	
	virtual double BeamSize() const { doThrow(); }
	
	virtual void GetGriddingCorrectionImage(double* image) const { doThrow(); }
	
	virtual bool HasGriddingCorrectionImage() const { doThrow(); }
	
private:
	virtual size_t getSuggestedWGridSize() const { doThrow(); }
	
	void doThrow() const
	{
		throw std::runtime_error("This gridder is not available, because WSClean was not compiled to have this gridder. Use a different gridder or recompile WSClean and make sure the necessary prerequisites are present.");
	}
};

#endif
