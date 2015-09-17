#ifndef CLEAN_ALGORITHM_H
#define CLEAN_ALGORITHM_H

#include <string>
#include <cmath>

#include "../polarizationenum.h"
#include "../uvector.h"

namespace ao {
	template<typename T> class lane;
}

class DeconvolutionAlgorithm
{
public:
	enum SpectralFittingMode {
		NoSpectralFitting,
		PolynomialSpectralFitting,
		LogPolynomialSpectralFitting
	};
	
	virtual ~DeconvolutionAlgorithm() { }
	
	void SetMaxNIter(size_t nIter) { _maxIter = nIter; }
	
	void SetThreshold(double threshold) { _threshold = threshold; }
	
	void SetGain(double gain) { _gain = gain; }
	
	void SetMGain(double mGain) { _mGain = mGain; }
	
	void SetAllowNegativeComponents(bool allowNegativeComponents) { _allowNegativeComponents = allowNegativeComponents; }
	
	void SetStopOnNegativeComponents(bool stopOnNegative) { _stopOnNegativeComponent = stopOnNegative; }
	
	void SetCleanBorderRatio(double borderRatio) { _cleanBorderRatio = borderRatio; }
	
	void SetThreadCount(size_t threadCount) { _threadCount = threadCount; }
	
	size_t MaxNIter() const { return _maxIter; }
	double Threshold() const { return _threshold; }
	double Gain() const { return _gain; }
	double MGain() const { return _mGain; }
	double CleanBorderRatio() const { return _cleanBorderRatio; }
	bool AllowNegativeComponents() const { return _allowNegativeComponents; }
	bool StopOnNegativeComponents() const { return _allowNegativeComponents; }
	
	void SetCleanMask(const bool* cleanMask) { _cleanMask = cleanMask; }
	
	size_t IterationNumber() const { return _iterationNumber; }
	
	void SetIterationNumber(size_t iterationNumber) { _iterationNumber = iterationNumber; }
	
	static void ResizeImage(double* dest, size_t newWidth, size_t newHeight, const double* source, size_t width, size_t height);
	
	static void GetModelFromImage(class Model &model, const double* image, size_t width, size_t height, double phaseCentreRA, double phaseCentreDec, double pixelSizeX, double pixelSizeY, double phaseCentreDL, double phaseCentreDM, double spectralIndex, double refFreq, 
																PolarizationEnum polarization = Polarization::StokesI);

	static void RemoveNaNsInPSF(double* psf, size_t width, size_t height);
	
	//static void CalculateFastCleanPSFSize(size_t& psfWidth, size_t& psfHeight, size_t imageWidth, size_t imageHeight);
	
	void CopyConfigFrom(const DeconvolutionAlgorithm& source)
	{
		_threshold = source._threshold;
		_gain = source._gain;
		_mGain = source._mGain;
		_cleanBorderRatio = source._cleanBorderRatio;
		_maxIter = source._maxIter;
		// skip _iterationNumber
		_allowNegativeComponents = source._allowNegativeComponents;
		_stopOnNegativeComponent = source._stopOnNegativeComponent;
		_cleanMask = source._cleanMask;
		_centralFrequencies = source._centralFrequencies;
		_spectralFittingMode = source._spectralFittingMode;
		_spectralFittingTerms = source._spectralFittingTerms;
	}
	
	void SetMultiscaleThresholdBias(double bias)
	{
		_multiscaleThresholdBias = bias;
	}
	void SetMultiscaleScaleBias(double bias)
	{
		_multiscaleScaleBias = bias;
	}
	void SetSpectralFittingMode(SpectralFittingMode mode, size_t nTerms)
	{
		_spectralFittingMode = mode;
		_spectralFittingTerms = nTerms;
	}
	
	void InitializeFrequencies(const ao::uvector<double>& frequencies)
	{
		_centralFrequencies = frequencies;
	}
protected:
	DeconvolutionAlgorithm();
	
	void PerformSpectralFit(double* values);
	
	double _threshold, _gain, _mGain, _cleanBorderRatio;
	double _multiscaleThresholdBias, _multiscaleScaleBias;
	size_t _maxIter, _iterationNumber, _threadCount;
	bool _allowNegativeComponents, _stopOnNegativeComponent;
	const bool* _cleanMask;
	
	ao::uvector<double> _centralFrequencies;
	
	enum SpectralFittingMode _spectralFittingMode;
	size_t _spectralFittingTerms;
};

template<typename ImageSetType>
class TypedDeconvolutionAlgorithm : public DeconvolutionAlgorithm
{
public:
	typedef ImageSetType ImageSet;
	
	virtual ~TypedDeconvolutionAlgorithm() { }
	
	virtual void ExecuteMajorIteration(ImageSetType& dataImage, ImageSetType& modelImage, std::vector<double*> psfImages, size_t width, size_t height, bool& reachedStopGain) = 0;
	
private:
};

class UntypedDeconvolutionAlgorithm : public DeconvolutionAlgorithm
{
public:
	virtual ~UntypedDeconvolutionAlgorithm() { }
	
	virtual void ExecuteMajorIteration(class DynamicSet& dataImage, class DynamicSet& modelImage, const ao::uvector<const double*>& psfImages, size_t width, size_t height, bool& reachedMajorThreshold) = 0;
	
private:
};

#endif
