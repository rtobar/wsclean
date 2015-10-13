#ifndef DECONVOLUTION_H
#define DECONVOLUTION_H

#include "../uvector.h"
#include "../wsclean/imagebufferallocator.h"

#include <cstring>

#include "deconvolutionalgorithm.h"

class Deconvolution
{
public:
	Deconvolution();
	~Deconvolution();
	
	void Perform(const class ImagingTable& groupTable, bool& reachedMajorThreshold, size_t majorIterationNr);
	
	void SetGain(double gain) { _gain = gain; }
	double Gain() const { return _gain; }
	
	void SetMGain(double mGain) { _mGain = mGain; }
	double MGain() const { return _mGain; }
	
	void SetNIter(size_t nIter) { _nIter = nIter; }
	size_t NIter() const { return _nIter; }
	
	void SetThreshold(double threshold) { _threshold = threshold; }
	double Threshold() const { return _threshold; }
	
	void SetCleanBorderRatio(double borderRatio) { _cleanBorderRatio = borderRatio; }
	
	void SetFitsMask(const std::string& fitsMask) { _fitsMask = fitsMask; }
	
	void SetCASAMask(const std::string& casaMask) { _casaMask = casaMask; }
	
	void SetAllowNegativeComponents(bool allowNegative) { _allowNegative = allowNegative; }
	bool AllowNegativeComponents() const { return _allowNegative; }
	
	void SetStopOnNegativeComponents(bool stopOnNegative) { _stopOnNegative = stopOnNegative; }
	bool StopOnNegativeComponents() const { return _stopOnNegative; }
	
	void SetMultiscale(bool multiscale) { _multiscale = multiscale; }
	void SetFastMultiscale(bool fastMultiscale) { _fastMultiscale = fastMultiscale; }
	void SetMultiscaleThresholdBias(double thresholdBias)
	{ _multiscaleThresholdBias = thresholdBias; }
	void SetMultiscaleScaleBias(double scaleBias)
	{ _multiscaleScaleBias = scaleBias; }
	void SetUseMoreSane(bool useMoreSane) { _useMoreSane = useMoreSane; }
	void SetUseIUWT(bool useIUWT) { _useIUWT = useIUWT; }
	void SetMoreSaneLocation(const std::string& location) { _moreSaneLocation = location; }
	void SetMoreSaneArgs(const std::string& arguments) { _moreSaneArgs = arguments; }
	void SetMoreSaneSigmaLevels(const std::vector<std::string> &slevels) { _moreSaneSigmaLevels = slevels; }
	void SetPrefixName(const std::string& prefixName) { _prefixName = prefixName; }
	
	bool IsSpectralFittingEnabled() {
		return _spectralFittingMode != NoSpectralFitting;
	}
	void SetFitSpectralPol(size_t nTerms) {
		_spectralFittingMode = PolynomialSpectralFitting;
		_spectralFittingTerms = nTerms;
	}
	void SetFitSpectralLogPol(size_t nTerms) {
		_spectralFittingMode = LogPolynomialSpectralFitting;
		_spectralFittingTerms = nTerms;
	}
	void SetDeconvolutionChannels(size_t deconvolutionChannelCount) {
		_requestedDeconvolutionChannelCount = deconvolutionChannelCount;
	}
	
	void InitializeDeconvolutionAlgorithm(const ImagingTable& groupTable, PolarizationEnum psfPolarization, ImageBufferAllocator* imageAllocator, size_t imgWidth, size_t imgHeight, double pixelScaleX, double pixelScaleY, size_t outputChannels, double beamSize, size_t threadCount);
	
	void InitializeImages(class CachedImageSet& residuals, CachedImageSet& models, CachedImageSet& psfs)
	{
		_residualImages = &residuals;
		_modelImages = &models;
		_psfImages = &psfs;
	}
	
	void FreeDeconvolutionAlgorithms();
	
	DeconvolutionAlgorithm& GetAlgorithm()
	{
		return *_cleanAlgorithm;
	}
	
	bool MultiScale() const { return _multiscale; }
	bool FastMultiScale() const { return _fastMultiscale; }
	bool UseMoreSane() const { return _useMoreSane; }
	bool UseIUWT() const { return _useIUWT; }
	bool IsInitialized() const { return _cleanAlgorithm != 0; }
private:
	void performSimpleClean(DynamicSet& residual, DynamicSet& model, const ao::uvector<const double*>& psfs, bool& reachedMajorThreshold, size_t majorIterationNr);
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
	
	double _threshold, _gain, _mGain;
	size_t _nIter;
	bool _allowNegative, _stopOnNegative;
	bool _multiscale, _fastMultiscale;
	double _multiscaleThresholdBias, _multiscaleScaleBias;
	double _cleanBorderRatio;
	std::string _fitsMask, _casaMask;
	bool _useMoreSane, _useIUWT;
	std::string _moreSaneLocation, _moreSaneArgs;
	std::vector<std::string> _moreSaneSigmaLevels;
	std::string _prefixName;
	
	enum SpectralFittingMode _spectralFittingMode;
	size_t _spectralFittingTerms;
	
	/**
	 * The number of channels used during deconvolution. This can be used to
	 * image with more channels than deconvolution. Before deconvolution,
	 * channels are averaged, and after deconvolution they are interpolated.
	 * It is 0 when all channels should be used.
	 */
	size_t _requestedDeconvolutionChannelCount;
	
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
