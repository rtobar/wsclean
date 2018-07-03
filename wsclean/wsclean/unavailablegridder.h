#ifndef NOT_IMPLEMENTED_GRIDDER_H
#define NOT_IMPLEMENTED_GRIDDER_H

#include "msgridderbase.h"

#include <stdexcept>
#include <string>

class UnavailableGridder : public MSGridderBase
{
public:
	UnavailableGridder(const WSCleanSettings&) { doThrow(); }
	
	virtual ~UnavailableGridder() final override { doThrow(); }
	
	virtual void Invert() final override { doThrow(); }
	
	virtual void Predict(double* real) final override { doThrow(); }
	
	virtual void Predict(double* real, double* imaginary) final override { doThrow(); }
	
	virtual double* ImageRealResult() final override { doThrow(); return 0; }
	
	virtual double* ImageImaginaryResult() final override { doThrow(); return 0; }
	
	virtual void GetGriddingCorrectionImage(double* image) const final override { doThrow(); }
	
	virtual bool HasGriddingCorrectionImage() const final override { doThrow(); return false; }
	
	void SavePBCorrectedImages(class FitsWriter& /*writer*/, class ImageFilename& /*filename*/, const std::string& /*filenameKind*/, class ImageBufferAllocator& /*allocator*/) const
	{ }
	
	void SaveBeamImage(const class ImagingTableEntry& /*entry*/, class ImageFilename& /*filename*/) const
	{ }

private:
	virtual size_t getSuggestedWGridSize() const final override { doThrow(); return 0; }
	
	void doThrow() const
	{
		throw std::runtime_error("This gridder is not available, because WSClean was not compiled to have this gridder. Use a different gridder or recompile WSClean and make sure the necessary prerequisites are satisfied.");
	}
};

#endif
