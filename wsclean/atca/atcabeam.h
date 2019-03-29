#ifndef ATCA_BEAM_H
#define ATCA_BEAM_H

#include "../uvector.h"

#include <string>

struct VoltagePattern
{
	// Array of size nsamples x frequencies, where the sample index is least significant (fastest changing)
	ao::uvector<double> values;
	ao::uvector<double> frequencies;
	double inverseIncrementRadius;
	double maximumRadiusArcMin;
	
	size_t NSamples() const { return values.size() / frequencies.size(); }
	
	const double* FreqIndexValues(size_t freqIndex) const { return &values[freqIndex * NSamples()]; }
	
	ao::uvector<double> InterpolateValues(double freq) const;
};

class ATCABeam
{
public:
	enum Band {
		ATCA_X,
		ATCA_C,
		ATCA_16,
		ATCA_K,
		ATCA_Q,
		ATCA_W,
		ATCA_L1,
		ATCA_L2,
		ATCA_L3,
		ATCA_S,
		ATCA_C_RI,
		ATCA_Unknown
	};
	
	static enum Band GetBand(double freqGHz);
	static enum Band GetBand(const std::string& str);
	
	static VoltagePattern CalculateVoltagePattern(enum Band band);
	
	static void Calculate(class PrimaryBeamImageSet& beamImages,
		size_t width, size_t height,
		double pixelScaleX, double pixelScaleY, 
		double phaseCentreRA, double phaseCentreDec,
		double phaseCentreDL, double phaseCentreDM,
		double frequencyHz,
		const VoltagePattern& voltagePattern
	);
	
private:
	static void calculateInversePolyNarrowband(VoltagePattern& vPattern, const ao::uvector<double>& coefficients);
	static void calculateInversePolyWideband(VoltagePattern& vPattern, const ao::uvector<double>& coefficients);
	static void calculatePoly(VoltagePattern& vPattern, const ao::uvector<double>& coefficients, double maximumRadiusArcMin);
};

#endif
