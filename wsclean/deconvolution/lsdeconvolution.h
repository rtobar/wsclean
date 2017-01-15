#ifndef LSDECONVOLUTION_H
#define LSDECONVOLUTION_H

#include <memory>
#include <string>

#include "../uvector.h"

#include "deconvolutionalgorithm.h"
#include "dynamicset.h"

struct LSDeconvolutionData;

class LSDeconvolution : public UntypedDeconvolutionAlgorithm
{
	public:
		LSDeconvolution();
		~LSDeconvolution();
		
    virtual void ExecuteMajorIteration(DynamicSet& dataImage, DynamicSet& modelImage, const ao::uvector<const double*>& psfImages, size_t width, size_t height, bool& reachedMajorThreshold)
		{
      _allocator = &dataImage.Allocator();
			ExecuteMajorIteration(dataImage[0], modelImage[0], psfImages[0], width, height, reachedMajorThreshold);
		}
		
		void ExecuteMajorIteration(double* dataImage, double* modelImage, const double* psfImage, size_t width, size_t height, bool& reachedMajorThreshold)
		{
			nonLinearFit(dataImage, modelImage, psfImage, width, height, reachedMajorThreshold);
		}
	private:
		void getMaskPositions(ao::uvector<std::pair<size_t, size_t>>& maskPositions, const bool* mask, size_t width, size_t height);
		
		void linearFit(double* dataImage, double* modelImage, const double* psfImage, size_t width, size_t height, bool& reachedMajorThreshold);
		
		void nonLinearFit(double* dataImage, double* modelImage, const double* psfImage, size_t width, size_t height, bool& reachedMajorThreshold);
		
		ImageBufferAllocator* _allocator;
		std::unique_ptr<LSDeconvolutionData> _data;
};

#endif
