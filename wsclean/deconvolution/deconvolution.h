#ifndef DECONVOLUTION_H
#define DECONVOLUTION_H

#include "../uvector.h"
#include "../wsclean/imagebufferallocator.h"
#include "../polarization.h"

#include <cstring>

class Deconvolution
{
public:
	explicit Deconvolution(const class WSCleanSettings& settings);
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
	
	class DeconvolutionAlgorithm& GetAlgorithm() { return *_cleanAlgorithm; }
	const DeconvolutionAlgorithm& GetAlgorithm() const { return *_cleanAlgorithm; }
	
	bool IsInitialized() const { return _cleanAlgorithm != 0; }
	
	void SaveComponentList(const class ImagingTable& table, long double phaseCentreRA, long double phaseCentreDec) const;
private:
	void calculateDeconvolutionFrequencies(const ImagingTable& groupTable, ao::uvector<double>& frequencies);
	
	const class WSCleanSettings& _settings;
	
	std::unique_ptr<class DeconvolutionAlgorithm> _cleanAlgorithm;
	
	ao::uvector<bool> _cleanMask;
	
	bool _autoMaskIsFinished;
	size_t _summedCount, _squaredCount;
	std::set<PolarizationEnum> _polarizations;
	PolarizationEnum _psfPolarization;
	size_t _imgWidth, _imgHeight;
	ImageBufferAllocator* _imageAllocator;
	CachedImageSet *_psfImages, *_modelImages, *_residualImages;
	ao::uvector<bool> _autoMask;
	double _beamSize, _pixelScaleX, _pixelScaleY;
};

#endif
