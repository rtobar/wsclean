#ifndef OUTPUT_CHANNEL_INFO_H
#define OUTPUT_CHANNEL_INFO_H

struct OutputChannelInfo {
	OutputChannelInfo() :
		weight(0.0),
		beamMaj(0.0), beamMin(0.0), beamPA(0.0),
		psfNormalizationFactor(1.0)
	{ }
	double weight, normalizationFactor;
	size_t wGridSize, visibilityCount;
	double effectiveVisibilityCount;
	double beamMaj, beamMin, beamPA;
	double theoreticBeamSize, psfNormalizationFactor;
};

#endif
