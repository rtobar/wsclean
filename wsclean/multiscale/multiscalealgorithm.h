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
	~MultiScaleAlgorithm();
	
	void SetManualScaleList(const ao::uvector<double>& scaleList) { _manualScaleList = scaleList; }
	
	//void PerformMajorIteration(size_t& iterCounter, size_t nIter, DynamicSet& modelSet, DynamicSet& dirtySet, const ao::uvector<const double*>& psfs, bool& reachedMajorThreshold);
	
	virtual void ExecuteMajorIteration(DynamicSet& dataImage, DynamicSet& modelImage, const ao::uvector<const double*>& psfImages, size_t width, size_t height, bool& reachedMajorThreshold);
	
	void SetAutoMaskMode(bool trackPerScaleMasks, bool usePerScaleMasks) {
		_trackPerScaleMasks = trackPerScaleMasks;
		_usePerScaleMasks = usePerScaleMasks; 
	}
	void SetUseFastSubMinorLoop(bool fastSubMinorLoop) {
		_fastSubMinorLoop = fastSubMinorLoop;
	}
private:
	class ImageBufferAllocator& _allocator;
	size_t _width, _height;
	double _beamSizeInPixels;
	ThreadedDeconvolutionTools* _tools;
	
	struct ScaleInfo
	{
		ScaleInfo() :
			scale(0.0), psfPeak(0.0),
			kernelPeak(0.0), biasFactor(0.0),
			gain(0.0), maxImageValue(0.0),
			rms(0.0),
			maxImageValueX(0), maxImageValueY(0),
			isActive(false),
			nComponentsCleaned(0),
			totalFluxCleaned(0.0)
		{ }
		
		double scale;
		double psfPeak, kernelPeak, biasFactor, gain;
		
		double maxImageValue, rms;
		size_t maxImageValueX, maxImageValueY;
		bool isActive;
		size_t nComponentsCleaned;
		double totalFluxCleaned;
	};
	std::vector<MultiScaleAlgorithm::ScaleInfo> _scaleInfos;
	ao::uvector<double> _manualScaleList;
	
	bool _trackPerScaleMasks, _usePerScaleMasks, _fastSubMinorLoop;
	std::vector<ao::uvector<bool>> _scaleMasks;

	void initializeScaleInfo();
	void convolvePSFs(std::unique_ptr<ImageBufferAllocator::Ptr[]>& convolvedPSFs, const double* psf, double* tmp, bool isIntegrated);
	void findActiveScaleConvolvedMaxima(const DynamicSet& imageSet, double* integratedScratch, bool reportRMS);
	void findSingleScaleMaximum(const double* convolvedImage, size_t scaleIndex);
	void sortScalesOnMaxima(size_t& scaleWithPeak);
	void activateScales(size_t scaleWithLastPeak);
	void measureComponentValues(ao::uvector<double>& componentValues, size_t scaleIndex, DynamicSet& imageSet);
	void addComponentToModel(double* model, size_t scaleWithPeak, double componentValue);
	double findPeak(const double *image, size_t &x, size_t &y, size_t scaleIndex);
	
	double* getConvolvedPSF(size_t psfIndex, size_t scaleIndex, const ao::uvector<const double*>& psfs, double* scratch, const std::unique_ptr<std::unique_ptr<ImageBufferAllocator::Ptr[]>[]>& convolvedPSFs);
	
};

#endif
