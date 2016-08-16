#ifndef IUWT_DECONVOLUTION_ALGORITHM_H
#define IUWT_DECONVOLUTION_ALGORITHM_H

#include "../uvector.h"

#include "iuwtdecomposition.h"
#include "imageanalysis.h"

#include <vector>

class IUWTDeconvolutionAlgorithm
{
public:
	IUWTDeconvolutionAlgorithm(size_t width, size_t height, double gain, double mGain, double cleanBorder, bool allowNegativeComponents, const bool* mask, double absoluteThreshold, double thresholdSigmaLevel=4.0, double tolerance=0.75, bool useSNRTest=true);
	
	void PerformMajorIteration(size_t& iterCounter, size_t nIter, class DynamicSet& modelSet, class DynamicSet& dirtySet, const ao::uvector<const double*>& psfs, bool& reachedMajorThreshold);
	
	void Subtract(double* dest, const ao::uvector<double>& rhs);
	void Subtract(ao::uvector<double>& dest, const ao::uvector<double>& rhs)
	{
		Subtract(dest.data(), rhs);
	}
	
private:
	struct ValComponent
	{
		ValComponent() {}
		ValComponent(size_t _x, size_t _y, int _scale, double _val = 0.0) : x(_x), y(_y), scale(_scale), val(_val) { }
		
		std::string ToString() const {
			std::ostringstream str;
			str << x << ',' << y << ", scale " << scale;
			return str.str();
		}
		
		size_t x, y;
		int scale;
		double val;
	};
	
	struct ScaleResponse
	{
		double rms, peakResponse, peakResponseToNextScale, convolvedPeakResponse;
		double bMaj, bMin, bPA;
		size_t convolvedArea;
	};

	double getMaxAbsWithoutMask(const ao::uvector<double>& data, size_t& x, size_t& y, size_t width);
	double getMaxAbsWithMask(const ao::uvector<double>& data, size_t& x, size_t& y, size_t width);
	double getMaxAbs(const ao::uvector<double>& data, size_t& x, size_t& y, size_t width)
	{
		if(_mask == 0)
			return getMaxAbsWithoutMask(data, x, y, width);
		else
			return getMaxAbsWithMask(data, x, y, width);
	}
	
	void measureRMSPerScale(const double* image, const double* convolvedImage, double* scratch, size_t endScale, std::vector<ScaleResponse>& psfResponse);
	
	double mad(const double* dest);
	
	double dotProduct(const ao::uvector<double>& lhs, const ao::uvector<double>& rhs);
	
	void factorAdd(double* dest, const double* rhs, double factor, size_t width, size_t height);
	
	void factorAdd(ao::uvector<double>& dest, const ao::uvector<double>& rhs, double factor);
	
	void boundingBox(size_t& x1, size_t& y1, size_t& x2, size_t& y2, const ao::uvector<double>& image, size_t width, size_t height);
	
	void adjustBox(size_t& x1, size_t& y1, size_t& x2, size_t& y2, size_t width, size_t height, int endScale);
	
	void trim(ao::uvector<double>& dest, const double* source, size_t oldWidth, size_t oldHeight, size_t x1, size_t y1, size_t x2, size_t y2);
	
	void trim(ao::uvector<double>& dest, const ao::uvector<double>& source, size_t oldWidth, size_t oldHeight, size_t x1, size_t y1, size_t x2, size_t y2)
	{
		trim(dest, source.data(), oldWidth, oldHeight, x1, y1, x2, y2);
	}
	
	void trimPsf(ao::uvector<double>& dest, const double* source, size_t oldWidth, size_t oldHeight, size_t newWidth, size_t newHeight)
	{
		trim(dest, source, oldWidth, oldHeight, (oldWidth-newWidth)/2, (oldHeight-newHeight)/2, (oldWidth+newWidth)/2, (oldHeight+newHeight)/2);
	}
	
	void untrim(ao::uvector<double>& image, size_t width, size_t height, size_t x1, size_t y1, size_t x2, size_t y2);
	
	double sum(const ao::uvector<double>& img) const;

	double snr(const IUWTDecomposition& noisyImg, const IUWTDecomposition& model) const;
	
	double rmsDiff(const ao::uvector<double>& a, const ao::uvector<double>& b);
	
	double rms(const ao::uvector<double>& image);
	
	bool runConjugateGradient(IUWTDecomposition& iuwt, const IUWTMask& mask, ao::uvector<double>& maskedDirty, ao::uvector<double>& structureModel, ao::uvector<double>& scratch, const ao::uvector<double>& psfKernel, size_t width, size_t height);
	
	bool fillAndDeconvolveStructure(IUWTDecomposition& iuwt, ao::uvector<double>& dirty, class DynamicSet& structureModelFull, ao::uvector<double>& scratch, const ao::uvector<double>& psf, const ao::uvector<double>& psfKernel, size_t curEndScale, size_t curMinScale, size_t width, size_t height, const ao::uvector<double>& thresholds, const ImageAnalysis::Component& maxComp, bool allowTrimming, const bool* priorMask);
	
	bool findAndDeconvolveStructure(IUWTDecomposition& iuwt, ao::uvector<double>& dirty, const ao::uvector<double>& psf, const ao::uvector<double>& psfKernel, ao::uvector<double>& scratch, class DynamicSet& structureModelFull, size_t curEndScale, size_t curMinScale, double gain, std::vector<ValComponent>& maxComponents);
	
	void performSubImageFitAll(IUWTDecomposition& iuwt, const IUWTMask& mask, const ao::uvector<double>& structureModel, ao::uvector<double>& scratchA, ao::uvector<double>& scratchB, const ImageAnalysis::Component& maxComp, DynamicSet& fittedModel, const double* psf, const ao::uvector<double>& dirty);
	
	void performSubImageFitSingle(IUWTDecomposition& iuwt, const IUWTMask& mask, const ao::uvector<double>& structureModel, ao::uvector<double>& scratchB, const ImageAnalysis::Component& maxComp, const double* psf, ao::uvector<double>& subDirty, double* fittedSubModel, ao::uvector<double>& correctionFactors);
	
	double performSubImageComponentFitBoxed(IUWTDecomposition& iuwt, const IUWTMask& mask, const std::vector<ImageAnalysis::Component2D>& area, ao::uvector<double>& scratch, ao::uvector<double>& maskedDirty, const double* psf, const ao::uvector<double>& psfKernel, size_t x1, size_t y1, size_t x2, size_t y2);
	
	double performSubImageComponentFit(IUWTDecomposition& iuwt, const IUWTMask& mask, const std::vector<ImageAnalysis::Component2D>& area, ao::uvector<double>& scratch, ao::uvector<double>& maskedDirty, const ao::uvector<double>& psfKernel, size_t xOffset, size_t yOffset);
	
	double centralPeak(const ao::uvector<double>& data)
	{
		return data[_width/2 + (_height/2)*_width];
	}
	void constrainedPSFConvolve(double* image, const double* psf, size_t width, size_t height);
	
	bool extractPointSources(const IUWTDecomposition& iuwt, const IUWTMask& mask, const double* dirty, double* model);
		
	size_t _width, _height;
	size_t _curBoxXStart, _curBoxXEnd;
	size_t _curBoxYStart, _curBoxYEnd;
	double _gain, _mGain, _cleanBorder;
	const bool* _mask;
	double _absoluteThreshold, _thresholdSigmaLevel, _tolerance;
	double _psfMaj, _psfMin, _psfPA, _psfVolume;
	ao::uvector<double> _rmses;
	FitsWriter _writer;
	std::vector<ScaleResponse> _psfResponse;
	bool _allowNegativeComponents, _useSNRTest;
	class DynamicSet* _modelSet;
	class DynamicSet* _dirtySet;
	ao::uvector<const double*> _psfs;
	class ThreadPool *_threadPool;
};

#endif
