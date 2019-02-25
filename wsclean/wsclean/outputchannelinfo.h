#ifndef OUTPUT_CHANNEL_INFO_H
#define OUTPUT_CHANNEL_INFO_H

#include <cstring>

struct OutputChannelInfo {
	OutputChannelInfo() :
		weight(0.0),
		normalizationFactor(1.0),
		wGridSize(0), visibilityCount(0),
		effectiveVisibilityCount(0.0),
		visibilityWeightSum(0.0),
		beamMaj(0.0), beamMin(0.0), beamPA(0.0),
		theoreticBeamSize(0.0),
		psfNormalizationFactor(1.0)
	{ }
	double weight, normalizationFactor;
	std::size_t wGridSize, visibilityCount;
	double effectiveVisibilityCount, visibilityWeightSum;
	double beamMaj, beamMin, beamPA;
	double theoreticBeamSize, psfNormalizationFactor;
};

#endif
