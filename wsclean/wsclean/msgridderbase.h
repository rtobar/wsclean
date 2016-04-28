#ifndef MS_GRIDDER_BASE_H
#define MS_GRIDDER_BASE_H

#include "inversionalgorithm.h"
#include "../multibanddata.h"

class MSGridderBase : public MeasurementSetGridder
{
public:
	MSGridderBase();
	
	virtual double HighestFrequencyChannel() const override final { return _freqHigh; }
	virtual double LowestFrequencyChannel() const override final { return _freqLow; }
	virtual double BandStart() const override final { return _bandStart; }
	virtual double BandEnd() const override final { return _bandEnd; }
	virtual double StartTime() const override final { return _startTime; }
	virtual double PhaseCentreRA() const override final { return _phaseCentreRA; }
	virtual double PhaseCentreDec() const override final { return _phaseCentreDec; }
	virtual double PhaseCentreDL() const override final { return _phaseCentreDL; }
	virtual double PhaseCentreDM() const override final { return _phaseCentreDM; }
	virtual bool HasDenormalPhaseCentre() const override final { return _denormalPhaseCentre; }
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
