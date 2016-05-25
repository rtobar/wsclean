#ifndef IUWT_DECONVOLUTION_H
#define IUWT_DECONVOLUTION_H

#include <memory>
#include <string>

#include "../uvector.h"

#include "deconvolutionalgorithm.h"
#include "dynamicset.h"

#include "../iuwt/iuwtdeconvolutionalgorithm.h"

#include "../wsclean/imagingtable.h"

class IUWTDeconvolution : public UntypedDeconvolutionAlgorithm
{
public:
	virtual void ExecuteMajorIteration(DynamicSet& dataImage, DynamicSet& modelImage, const ao::uvector<const double*>& psfImages, size_t width, size_t height, bool& reachedMajorThreshold) override final
	{
		IUWTDeconvolutionAlgorithm algorithm(width, height, _gain, _mGain, _cleanBorderRatio, _allowNegativeComponents, _cleanMask, _threshold);
		algorithm.PerformMajorIteration(_iterationNumber, MaxNIter(), modelImage, dataImage, psfImages, reachedMajorThreshold);
		if(_iterationNumber >= MaxNIter())
			reachedMajorThreshold = false;
	}
	
private:
};

#endif
