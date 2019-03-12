#ifndef DIRECT_MS_GRIDDER_H
#define DIRECT_MS_GRIDDER_H

#include "../lane.h"

#include "imagebufferallocator.h"
#include "msgridderbase.h"

template<typename num_t>
class DirectMSGridder : public MSGridderBase
{
public:
	const static size_t num_t_factor = (sizeof(num_t) + sizeof(double) - 1) / sizeof(double);
	
	DirectMSGridder(class ImageBufferAllocator* imageAllocator, size_t nThreads);

	virtual void Invert() final override;
	
	virtual void Predict(ImageBufferAllocator::Ptr image) final override;
	virtual void Predict(ImageBufferAllocator::Ptr /*real*/, ImageBufferAllocator::Ptr /*imaginary*/) final override
	{
		throw std::runtime_error("Direct FT imager can not predict complex images");
	}
	
	virtual ImageBufferAllocator::Ptr ImageRealResult() final override { return std::move(_image); }
	virtual ImageBufferAllocator::Ptr ImageImaginaryResult() final override {
		throw std::runtime_error("Direct FT imager can not make complex images");
	}
	virtual bool HasGriddingCorrectionImage() const final override { return false; }
	virtual void GetGriddingCorrectionImage(double *) const final override { }
	virtual size_t getSuggestedWGridSize() const override final { return 1; }
	
private:
	struct InversionSample {
		num_t uInLambda, vInLambda, wInLambda;
		std::complex<float> sample;
	};
	size_t _nThreads;
	ImageBufferAllocator::Ptr _image;
	num_t* _sqrtLMTable;
	std::vector<num_t*> _layers;
	ao::lane<InversionSample> _inversionLane;
	ImageBufferAllocator* _imageAllocator;
	
	void invertMeasurementSet(const MSData& msData, class ProgressBar& progress, size_t msIndex);
	void inversionWorker(size_t layer);
	void gridSample(const InversionSample& sample, size_t layer);
	void initializeSqrtLMLookupTable();
	
	num_t* allocate()
	{
		return new(_imageAllocator->Allocate(ImageWidth() * ImageHeight() * num_t_factor)) num_t[ImageWidth() * ImageHeight()];
	}
	void freeImg(num_t* ptr)
	{
		_imageAllocator->Free(reinterpret_cast<double*>(ptr));
	}
};

#endif
