#ifndef DECONVOLUTION_H
#define DECONVOLUTION_H

#include "../uvector.h"
#include "../wsclean/imagebufferallocator.h"
#include "../polarizationenum.h"

#include <cstring>

class Deconvolution
{
public:
	Deconvolution(const class WSCleanSettings& settings);
	~Deconvolution();
	
	void Perform(const class ImagingTable& groupTable, bool& reachedMajorThreshold, size_t majorIterationNr);
	
	void InitializeDeconvolutionAlgorithm(const ImagingTable& groupTable, PolarizationEnum psfPolarization, ImageBufferAllocator* imageAllocator, size_t imgWidth, size_t imgHeight, double pixelScaleX, double pixelScaleY, size_t outputChannels, double beamSize, size_t threadCount);
	
	void InitializeImages(class CachedImageSet& residuals, CachedImageSet& models, CachedImageSet& psfs)
	{
		_residualImages = &residuals;
		_modelImages = &models;
		_psfImages = &psfs;
	}
	
	void FreeDeconvolutionAlgorithms();
	
	class DeconvolutionAlgorithm& GetAlgorithm()
	{
		return *_cleanAlgorithm;
	}
	
	bool IsInitialized() const { return _cleanAlgorithm != 0; }
private:
	void performSimpleClean(class DynamicSet& residual, DynamicSet& model, const ao::uvector<const double*>& psfs, bool& reachedMajorThreshold, size_t majorIterationNr);
	void performSimpleClean(size_t currentChannelIndex, bool& reachedMajorThreshold, size_t majorIterationNr, PolarizationEnum polarization);
	
	template<size_t PolCount>
	void performJoinedPolClean(DynamicSet& residual, DynamicSet& model, const ao::uvector<const double*>& psfs, bool& reachedMajorThreshold, size_t majorIterationNr);
	template<size_t PolCount>
	void performJoinedPolClean(size_t currentChannelIndex, bool& reachedMajorThreshold, size_t majorIterationNr);
	
	template<size_t PolCount>
	void performJoinedPolFreqClean(bool& reachedMajorThreshold, size_t majorIterationNr);
	template<size_t PolCount>
	void performJoinedPolFreqClean(DynamicSet& residual, DynamicSet& model, const ao::uvector<const double*>& psfs, bool& reachedMajorThreshold, size_t majorIterationNr);

	void calculateDeconvolutionFrequencies(const ImagingTable& groupTable, ao::uvector<double>& frequencies);
	
	const class WSCleanSettings& _settings;
	
	std::unique_ptr<class DeconvolutionAlgorithm> _cleanAlgorithm;
	
	ao::uvector<bool> _cleanMask;
	
	size_t _summedCount, _squaredCount;
	std::set<PolarizationEnum> _polarizations;
	PolarizationEnum _psfPolarization;
	size_t _imgWidth, _imgHeight;
	ImageBufferAllocator* _imageAllocator;
	CachedImageSet *_psfImages, *_modelImages, *_residualImages;
};

#endif
