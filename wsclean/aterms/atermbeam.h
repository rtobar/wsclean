#ifndef ATERM_BEAM_H
#define ATERM_BEAM_H

#include "atermbase.h"

class ATermBeam : public ATermBase
{
public:
	ATermBeam() :
	_updateInterval(0),
	_lastATermUpdate(0)
	{ }
	
	virtual bool Calculate(std::complex<float>* buffer, double time, double frequency) final override
	{
		if(time - _lastATermUpdate > _updateInterval)
		{
			_lastATermUpdate = time;
			return calculateBeam(buffer, time + _updateInterval*0.5, frequency);
		}
		else {
			return false;
		}
	}
	
	void SetUpdateInterval(double updateInterval)
	{
		_updateInterval = updateInterval;
		_lastATermUpdate = -_updateInterval - 1;
	}
	
protected:
	virtual bool calculateBeam(std::complex<float>* buffer, double time, double frequency) = 0;
	
private:
	double _updateInterval, _lastATermUpdate;
};

#endif

