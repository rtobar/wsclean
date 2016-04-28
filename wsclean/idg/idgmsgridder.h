#ifndef IDG_MS_GRIDDER_H
#define IDG_MS_GRIDDER_H

#ifdef HAVE_IDG

#include <common/Proxy.h>

#include "../wsclean/msgridderbase.h"

class IdgMsGridder : public MSGridderBase
{
public:
	IdgMsGridder();
	
	virtual ~IdgMsGridder();
	
	virtual void Invert() override;
	
	virtual void Predict(double* real) final override;
	
	virtual void Predict(double* real, double* imaginary) final override;
	
	virtual double* ImageRealResult() const final override;
	
	virtual double* ImageImaginaryResult() const final override;
	
	virtual double BeamSize() const final override { return 0.0; }
	
	virtual double ImageWeight() const final override;
	
	virtual void GetGriddingCorrectionImage(double* image) const final override;
	
	virtual bool HasGriddingCorrectionImage() const final override;
	
private:
	idg::proxy::Proxy *_proxy;
	size_t _gridSize;
	size_t _imageSize;
};

#endif // HAVE IDG

#endif
