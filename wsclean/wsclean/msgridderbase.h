#ifndef MS_GRIDDER_BASE_H
#define MS_GRIDDER_BASE_H

#include "inversionalgorithm.h"
#include "../multibanddata.h"

class MSGridderBase : public InversionAlgorithm
{
public:
	MSGridderBase();
	
	virtual double HighestFrequencyChannel() const { return _freqHigh; }
	virtual double LowestFrequencyChannel() const { return _freqLow; }
	virtual double BandStart() const { return _bandStart; }
	virtual double BandEnd() const { return _bandEnd; }
	virtual double StartTime() const { return _startTime; }
	virtual double PhaseCentreRA() const { return _phaseCentreRA; }
	virtual double PhaseCentreDec() const { return _phaseCentreDec; }
	virtual double PhaseCentreDL() const { return _phaseCentreDL; }
	virtual double PhaseCentreDM() const { return _phaseCentreDM; }
	virtual bool HasDenormalPhaseCentre() const { return _denormalPhaseCentre; }
protected:
	void resetMetaData()
	{
		_hasFrequencies = false;
	}
	
	void initializePhaseCentre(casacore::MeasurementSet& ms, size_t fieldId);
	
	void updateMetaDataForMS(const MultiBandData& selectedBand, double startTime)
	{
		if(_hasFrequencies)
		{
			_freqLow = std::min(_freqLow, selectedBand.LowestFrequency());
			_freqHigh = std::max(_freqHigh, selectedBand.HighestFrequency());
			_bandStart = std::min(_bandStart, selectedBand.BandStart());
			_bandEnd = std::max(_bandEnd, selectedBand.BandEnd());
			_startTime = std::min(_startTime, startTime);
		} else {
			_freqLow = selectedBand.LowestFrequency();
			_freqHigh = selectedBand.HighestFrequency();
			_bandStart = selectedBand.BandStart();
			_bandEnd = selectedBand.BandEnd();
			_startTime = startTime;
			_hasFrequencies = true;
		}
	}

private:
	bool _hasFrequencies;
	double _freqHigh, _freqLow;
	double _bandStart, _bandEnd;
	double _startTime;
	
	double _phaseCentreRA, _phaseCentreDec, _phaseCentreDL, _phaseCentreDM;
	bool _denormalPhaseCentre;
};

#endif
