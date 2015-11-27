#ifndef SPECTRAL_FITTER_H
#define SPECTRAL_FITTER_H

#include "../uvector.h"

enum SpectralFittingMode {
	NoSpectralFitting,
	PolynomialSpectralFitting,
	LogPolynomialSpectralFitting
};

class SpectralFitter
{
public:
	SpectralFitter(SpectralFittingMode mode, size_t nTerms) :
		_mode(mode), _nTerms(nTerms)
	{ }
	
	void SetMode(SpectralFittingMode mode, size_t nTerms)
	{
		_mode = mode;
		_nTerms = nTerms;
	}
	
	void FitAndEvaluate(double* values) const;
	
	void Fit(ao::uvector<double>& terms, const double* values) const;
	
	void Evaluate(double* values, const ao::uvector<double>& terms) const;
	
	double Evaluate(const ao::uvector<double>& terms, double frequency) const;
	
	void SetFrequencies(const double* frequencies, size_t n)
	{
		_frequencies.assign(frequencies, frequencies+n);
	}
	
	size_t NTerms() const { return _nTerms; }
	
private:
	double referenceFrequency() const {
		return _frequencies[_frequencies.size()/2];
	}
	enum SpectralFittingMode _mode;
	size_t _nTerms;
	ao::uvector<double> _frequencies;
};

#endif
