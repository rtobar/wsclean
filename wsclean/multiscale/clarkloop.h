#ifndef CLARK_LOOP_H
#define CLARK_LOOP_H

#include <cstring>
#include <vector>

#include "../deconvolution/dynamicset.h"

/**
 * In multi-scale, a subminor Clark-optimized loop looks like this:
 * 
 * IterateAndMakeModel():
 * - Make a set S with positions of all the components larger than 'threshold', which are also in the mask
 * - Find the largest component in S
 * Loop {
 * - Measure the largest component per frequency (from S)
 * - Store the model component in S
 * - Subtract this component multiplied with the double convolved PSF and gain from all components in S (per individual image)
 * - Find the new largest component in S
 * }
 * 
 * CorrectResidualDirty():
 * For each individual image {
 * - Put the model components from S onto a full image (using GetFullIndividualModel())
 * - Convolve the model with the SingleConvolvedPSF
 * - Subtract the convolved model from the residual
 * }
 *
 * Finalization:
 * - Put the model components from S onto a full image (using GetFullIndividualModel())
 * - Convolve the model image with the scale kernel
 * - Add the model components to the full model
 */

class ClarkModel
{
public:
	ClarkModel(size_t width, size_t height) :
	_width(width), _height(height)
	{ }
	
	void AddPosition(size_t x, size_t y)
	{ _positions.push_back(std::make_pair(x, y)); }
	
	size_t size() const { return _positions.size(); }
	
	void MakeSets(const DynamicSet& templateSet);
	
	DynamicSet& Residual() { return *_residual; }
	const DynamicSet& Residual() const { return *_residual; }
	
	DynamicSet& Model() { return *_model; }
	const DynamicSet& Model() const { return *_model; }
	
	size_t X(size_t index) const { return _positions[index].first; }
	size_t Y(size_t index) const { return _positions[index].second; }
	size_t FullIndex(size_t index) const { return X(index) + Y(index) * _width; }
	template<bool AllowNegatives>
	size_t GetMaxComponent(double* scratch, double& maxValue) const;
	size_t GetMaxComponent(double* scratch, double& maxValue, bool allowNegatives) const
	{
		if(allowNegatives)
			return GetMaxComponent<true>(scratch, maxValue);
		else
			return GetMaxComponent<false>(scratch, maxValue);
	}
private:
	std::vector<std::pair<size_t,size_t>> _positions;
	std::unique_ptr<DynamicSet> _residual, _model;
	size_t _width, _height;
};

class ClarkLoop
{
public:
	ClarkLoop(size_t width, size_t height) :
	_width(width), _height(height),
	_threshold(0.0), _gain(0.0),
	_currentIteration(0), _maxIterations(0),
	_allowNegativeComponents(true),
	_mask(0), _fitter(0),
	_clarkModel(width, height),
	_fluxCleaned(0.0)
	{ }
	
	void SetThreshold(double threshold)
	{ _threshold = threshold; }
	
	void SetIterationInfo(size_t currentIteration, size_t maxIterations)
	{ _currentIteration = currentIteration; _maxIterations = maxIterations; }
	
	void SetGain(double gain)
	{ _gain = gain; }
	
	void SetAllowNegativeComponents(bool allowNegativeComponents)
	{ _allowNegativeComponents = allowNegativeComponents; }
	
	void SetSpectralFitter(const SpectralFitter* fitter) { _fitter = fitter; }

	size_t CurrentIteration() const { return _currentIteration; }
	
	double FluxCleaned() const { return _fluxCleaned; }
	
	void SetMask(const bool* mask)
	{ _mask = mask; }
	
	void Run(DynamicSet& convolvedResidual, const std::vector<const double*>& doubleConvolvedPsfs);
	
	void CorrectResidualDirty(double* scratchA, double* scratchB, size_t imageIndex, double* residual, const double* singleConvolvedPsf) const;
	
	void GetFullIndividualModel(size_t imageIndex, double* individualModelImg) const;
	
	void UpdateAutoMask(bool* mask) const;
	
private:
	void findPeakPositions(DynamicSet& convolvedResidual);
	
	size_t _width, _height;
	double _threshold, _gain;
	size_t _currentIteration, _maxIterations;
	bool _allowNegativeComponents;
	const bool* _mask;
	const SpectralFitter* _fitter;
	ClarkModel _clarkModel;
	double _fluxCleaned;
};

#endif
