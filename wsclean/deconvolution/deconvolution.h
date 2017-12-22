#ifndef DECONVOLUTION_H
#define DECONVOLUTION_H

#include "paralleldeconvolution.h"

#include "../wsclean/imagebufferallocator.h"
#include "../polarization.h"
#include "../uvector.h"

#include <cstring>

class Deconvolution
{
public:
	explicit Deconvolution(const class WSCleanSettings& settings);
	~Deconvolution();
	
	void Perform(const class ImagingTable& groupTable, bool& reachedMajorThreshold, size_t majorIterationNr);
	
	void InitializeDeconvolutionAlgorithm(const ImagingTable& groupTable, PolarizationEnum psfPolarization, ImageBufferAllocator* imageAllocator, double beamSize, size_t threadCount);
	
	void InitializeImages(class CachedImageSet& residuals, CachedImageSet& models, CachedImageSet& psfs)
	{
		_residualImages = &residuals;
		_modelImages = &models;
		_psfImages = &psfs;
	}
	
	void FreeDeconvolutionAlgorithms()
	{
		_parallelDeconvolution.FreeDeconvolutionAlgorithms();
	}
	
	class DeconvolutionAlgorithm& GetAlgorithm()
	{
		return _parallelDeconvolution.FirstAlgorithm();
	}
	const DeconvolutionAlgorithm& GetAlgorithm() const
	{
		return _parallelDeconvolution.FirstAlgorithm();
	}
	
	bool IsInitialized() const { return _parallelDeconvolution.IsInitialized(); }
	
	void SaveSourceList(const class ImagingTable& table, long double phaseCentreRA, long double phaseCentreDec)
	{
		_parallelDeconvolution.SaveSourceList(*_modelImages, table, phaseCentreRA, phaseCentreDec);
	}
	
	void SavePBSourceList(const class ImagingTable& table, long double phaseCentreRA, long double phaseCentreDec)
	{
		_parallelDeconvolution.SavePBSourceList(*_modelImages, table, phaseCentreRA, phaseCentreDec);
	}

private:
	void calculateDeconvolutionFrequencies(const ImagingTable& groupTable, ao::uvector<double>& frequencies, ao::uvector<double>& weights);
	
	void correctChannelForPB(class ComponentList& list, const class ImagingTableEntry& entry) const;
	
	void readMask(const ImagingTable& groupTable);
	
	const class WSCleanSettings& _settings;
	
	ParallelDeconvolution _parallelDeconvolution;
	
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
