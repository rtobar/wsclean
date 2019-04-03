#include "atcabeam.h"

#include <cmath>

#include "../wsclean/logger.h"
#include "../wsclean/primarybeamimageset.h"

enum ATCABeam::Band ATCABeam::GetBand(double freqGHz)
{
	if (freqGHz >= 7.0 && freqGHz < 12.5) {
		return ATCA_X;
	} else if (freqGHz > 4.0 && freqGHz < 7.0) {
		return ATCA_C;  
	} else if (freqGHz > 1.0 && freqGHz < 3.5) {
		return ATCA_16;
	} else if (freqGHz > 15.0 && freqGHz < 26.0) {
		return ATCA_K;
	} else if (freqGHz > 30.0 && freqGHz < 50.0) {
		return ATCA_Q;
	} else if (freqGHz > 80 && freqGHz < 120.0) {
		return ATCA_W;
	} else {
		return ATCA_L1; // ATCA_Unknown ?
	}
}
 
enum ATCABeam::Band ATCABeam::GetBand(const std::string& str)
{
	if (str == "ATCA_L1") {
		return ATCA_L1;
	} else if (str == "ATCA_L2") {
		return ATCA_L2;
	} else if (str == "ATCA_L3") {
		return ATCA_L3;
	} else if (str == "ATCA_16") {
		return ATCA_16;
	} else if (str == "ATCA_S") {
		return ATCA_S;
	} else if (str == "ATCA_C") {
		return ATCA_C;
	} else if (str == "ATCA_C_RI") {
		return ATCA_C_RI;
	} else if (str == "ATCA") {
		return ATCA_Unknown;
	} else if (str == "ATCA_X") {
		return ATCA_X;
	} else if (str == "ATCA_K") {
		return ATCA_K;
	} else if (str == "ATCA_Q") {
		return ATCA_Q;
	} else if (str == "ATCA_W") {
		return ATCA_W;
	}
	return ATCA_Unknown;
}

VoltagePattern ATCABeam::CalculateVoltagePattern(enum Band band)
{
	VoltagePattern vPattern;
	switch(band)
	{
		case ATCA_16:
		{
			Logger::Debug << "Using voltage pattern for ATCA_16 band\n";
			// coef x nfreq
			ao::uvector<double> coefficients({
				1.0, 1.06274e-03, 1.32342e-06, -8.72013e-10, 1.08020e-12,
				1.0, 9.80817e-04, 1.17898e-06, -7.83160e-10, 8.66199e-13,
				1.0, 9.53553e-04, 9.33233e-07, -4.26759e-10, 5.63667e-13,
				1.0, 9.78268e-04, 6.63231e-07, 4.18235e-11, 2.62297e-13,
				1.0, 1.02424e-03, 6.12726e-07, 2.25733e-10, 2.04834e-13,
				1.0, 1.05818e-03, 5.37473e-07, 4.22386e-10, 1.17530e-13,
				1.0, 1.10650e-03, 5.11574e-07, 5.89732e-10, 8.13628e-14
			});
			vPattern.frequencies.resize(7);
			for (size_t i=0; i!=7; ++i) {
				vPattern.frequencies[i] = (1332+i*256)*1.e6;
			}
			vPattern.maximumRadiusArcMin = 53.0;
			calculateInversePolyWideband(vPattern, coefficients);
		} break;
		default:
			throw std::runtime_error("ATCA PB for given spectral band not implemented");
	}
	return vPattern;
}

// from void PBMath1DIPoly::fillPBArray()
void ATCABeam::calculateInversePolyNarrowband(VoltagePattern& vPattern, const ao::uvector<double>& coefficients)
{
	size_t nSamples=10000;
	vPattern.values.resize(nSamples);

	double inverseIncrementRadius = double(nSamples-1) / vPattern.maximumRadiusArcMin;
	for(size_t i=0;i<nSamples;i++) {
		double taper = 0.0;
		double x2 = double(i) / inverseIncrementRadius;
		x2 = x2*x2;
		double y = 1.0;
		for (size_t j=0; j!=coefficients.size(); j++) {
			taper += y * coefficients[j];
			y *= x2;
		}
		if (taper != 0.0) {
			taper = 1.0 / std::sqrt(taper);
		}
		vPattern.values[i] = taper;
  }
};

// from void PBMath1DIPoly::fillPBArray()
void ATCABeam::calculateInversePolyWideband(VoltagePattern& vPattern, const ao::uvector<double>& coefficients)
{
	size_t nSamples=10000;
	size_t nfreq = vPattern.frequencies.size();
	size_t ncoef = coefficients.size() / nfreq;
	vPattern.values.resize(nSamples * nfreq);
	vPattern.inverseIncrementRadius = double(nSamples-1) / vPattern.maximumRadiusArcMin;
	double* output = vPattern.values.data();
	for(size_t n=0; n != nfreq; n++)
	{
		const double* freqcoefficients = &coefficients[n * ncoef];
		for(size_t i=0;i<nSamples;i++)
		{
			double taper = 0.0;
			double x2 = double(i) / vPattern.inverseIncrementRadius;
			x2 = x2*x2;
			double y = 1.0;
			
			for (size_t j=0; j<ncoef; j++) {
				taper += y * freqcoefficients[j];
				y *= x2;
			}
			if (taper != 0.0) {
				taper = 1.0 / std::sqrt(taper);
			}
			*output = taper;
			++output;
		}
	}
};

// From PBMath1DPoly::fillPBArray()
void ATCABeam::calculatePoly(VoltagePattern& vPattern, const ao::uvector<double>& coefficients, double maximumRadiusArcMin)
{
	unsigned nSamples=10000;
	vPattern.values.resize(nSamples);

	double inverseIncrementRadius = double(nSamples-1) / maximumRadiusArcMin;

	for(size_t i=0; i<nSamples; i++)
	{
		double taper = 0.0;
		double x2 = double(i) / inverseIncrementRadius;
		x2 = x2 * x2;
		double y = 1.0;
		for (size_t j=0; j<coefficients.size(); j++)
		{
			taper += y * coefficients[j];
			y *= x2;
		}
		vPattern.values[i] = std::sqrt(taper);
	}
}

ao::uvector<double> VoltagePattern::InterpolateValues(double freq) const
{
	ao::uvector<double> result;
	size_t ifit = 0;
	size_t nFreq = frequencies.size();
	for (ifit=0; ifit!=nFreq; ifit++) {
		if (freq <= frequencies[ifit]) break;
	}
	if (ifit==0) {
		Logger::Debug << "Using voltage pattern coefficients for " << frequencies[0]*1e-6 << " MHz instead of requested " << freq*1e-6 << '\n';
		result.assign(values.begin(), values.begin() + NSamples());
	} else if (ifit==nFreq) {
		Logger::Debug << "Using voltage pattern coefficients for " << frequencies[nFreq-1]*1e-6 << " MHz instead of requested " << freq*1e-6 << '\n';
		result.assign(values.begin() + (nFreq-1) * NSamples(), values.end());
	} else {
		Logger::Debug << "Interpolating voltage pattern coefficients from " << frequencies[ifit-1]*1e-6 << " and " << frequencies[ifit]*1e-6 << " MHz to " << freq*1e-6 << " MHz.\n";
		size_t n = NSamples();
		double l = (freq - frequencies[ifit-1]) / (frequencies[ifit] - frequencies[ifit-1]);
		const double *vpA = FreqIndexValues(ifit-1);
		const double *vpB = FreqIndexValues(ifit);
		result.resize(n);
		for(size_t i=0; i!=n; ++i)
		{
			result[i] = vpA[i]*(1.0-l) + vpB[i]*l;
		}
	}
	return result;
}

void ATCABeam::Calculate(PrimaryBeamImageSet& beamImages,
	size_t width, size_t height,
	double pixelScaleX, double pixelScaleY, 
	double, double,
	double phaseCentreDL, double phaseCentreDM,
	double frequencyHz,
	const VoltagePattern& voltagePattern)
{
	double factor = (180.0 / M_PI) * 60.0 * frequencyHz * 1.0e-9 ;  // arcminutes * GHz
	double rmax2 = voltagePattern.maximumRadiusArcMin / factor;
	rmax2 = rmax2 * rmax2;
	
	ao::uvector<double> interpolatedVals;
	const double *vp;
	if (voltagePattern.frequencies.size() > 1)
	{
		interpolatedVals = voltagePattern.InterpolateValues(frequencyHz);
		vp = interpolatedVals.data();
	}
	else
		vp = voltagePattern.FreqIndexValues(0);
	
	double inverseIncrementRadius=voltagePattern.inverseIncrementRadius;
	
	double x0 = double(width/2) * pixelScaleX - phaseCentreDL;
	double y0 = double(height/2) * pixelScaleY - phaseCentreDM;
	size_t imgIndex = 0;
	Logger::Debug << "Interpolating 1D voltage pattern to output image...\n";
	for(size_t iy=0; iy!=height; ++iy) {
		double ry2 =  pixelScaleY*double(iy) - y0;
		ry2 = ry2*ry2;
		for(size_t ix=0; ix!=width; ++ix) {
			double out;
			double rx = pixelScaleX*double(ix) - x0;
			double r2 = rx*rx + ry2;
			if (r2 > rmax2) {
				out = 0.0;
			}
			else {
				double r = std::sqrt(r2) * factor;
				int indx = int(r*inverseIncrementRadius);
				out = vp[indx];
			}
			
			beamImages[0][imgIndex] = out;
			beamImages[1][imgIndex] = 0.0;
			beamImages[2][imgIndex] = 0.0;
			beamImages[3][imgIndex] = 0.0;
			beamImages[4][imgIndex] = 0.0;
			beamImages[5][imgIndex] = 0.0;
			beamImages[6][imgIndex] = out;
			beamImages[7][imgIndex] = 0.0;
			++imgIndex;
		}
	}
};
