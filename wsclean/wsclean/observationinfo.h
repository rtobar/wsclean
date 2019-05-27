#ifndef OBSERVATION_INFO_H
#define OBSERVATION_INFO_H

struct ObservationInfo
{
	double phaseCentreRA, phaseCentreDec;
	double startTime;
	bool hasDenormalPhaseCentre;
	double phaseCentreDL, phaseCentreDM;
	std::string telescopeName;
	std::string observer;
	std::string fieldName;
};

#endif

