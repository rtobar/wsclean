#ifndef MULTISCALE_ALGORITHM_H
#define MULTISCALE_ALGORITHM_H

#include <cstring>
#include <vector>

#include "threadeddeconvolutiontools.h"

#include "../uvector.h"

#include "../deconvolution/dynamicset.h"
#include "../deconvolution/deconvolutionalgorithm.h"

#include "../wsclean/imagebufferallocator.h"

class MultiScaleAlgorithm : public UntypedDeconvolutionAlgorithm
{
public:
	MultiScaleAlgorithm(class ImageBufferAllocator& allocator, double beamSize, double pixelScaleX, double pixelScaleY);
	
	void SetCleanMask(const bool* cleanMask) { _cleanMask = cleanMask; }
	
	//void PerformMajorIteration(size_t& iterCounter, size_t nIter, DynamicSet& modelSet, DynamicSet& dirtySet, const ao::uvector<const double*>& psfs, bool& reachedMajorThreshold);
	
	virtual void ExecuteMajorIteration(DynamicSet& dataImage, DynamicSet& modelImage, const ao::uvector<const double*>& psfImages, size_t width, size_t height, bool& reachedMajorThreshold);	
private:
	class ImageBufferAllocator& _allocator;
	size_t _width, _height;
	double _beamSizeInPixels;
	bool _verbose;
	ThreadedDeconvolutionTools* _tools;
	
	struct ScaleInfo
	{
		double scale;
		double psfPeak, kernelPeak, factor, gain;
		
		double maxImageValue;
		size_t maxImageValueX, maxImageValueY;
		bool isActive;
	};
	std::vector<MultiScaleAlgorithm::ScaleInfo> _scaleInfos;

	void initializeScaleInfo();
	void convolvePSFs(std::unique_ptr<ImageBufferAllocator::Ptr[]>& convolvedPSFs, const double* psf, double* tmp, bool isIntegrated);
	void findActiveScaleConvolvedMaxima(const DynamicSet& imageSet, double* integratedScratch);
	void findSingleScaleMaximum(const double* convolvedImage, size_t scaleIndex);
	void sortScalesOnMaxima(size_t& scaleWithPeak);
	void activateScales(size_t scaleWithLastPeak);
	void measureComponentValues(ao::uvector<double>& componentValues, size_t scaleIndex, DynamicSet& imageSet);
	void addComponentToModel(double* model, size_t scaleWithPeak, double componentValue);
	double findPeak(const double *image, size_t &x, size_t &y);
	
	double* getConvolvedPSF(size_t psfIndex, size_t scaleIndex, const ao::uvector<const double*>& psfs, double* scratch, const std::unique_ptr<std::unique_ptr<ImageBufferAllocator::Ptr[]>[]>& convolvedPSFs);
	
};

#endif
