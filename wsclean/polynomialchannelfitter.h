#ifndef POLYNOMIAL_CHANNEL_FITTER_H
#define POLYNOMIAL_CHANNEL_FITTER_H

#include "uvector.h"

class PolynomialChannelFitter
{
public:
	void Clear() { _channels.clear(); _dataPoints.clear(); }
	
	void AddChannel(double startFrequency, double endFrequency)
	{
		_channels.push_back(std::make_pair(startFrequency, endFrequency));
	}
	
	void AddDataPoint(size_t channel, double y)
	{
		_dataPoints.push_back(std::make_pair(channel, y));
	}
	
	void Fit(ao::uvector<double>& terms, size_t nTerms);
	
	static double Evaluate(double x, const ao::uvector<double>& terms)
	{
		double val = terms[0];
		double f = 1.0;
		for(size_t i=1; i!=terms.size(); ++i)
		{
			f *= x;
			val += f * terms[i];
		}
		return val;
	}

private:
	/**
	 * Start and end frequencies of the channels
	 */
	ao::uvector<std::pair<double, double>> _channels;
	ao::uvector<std::pair<size_t,double>> _dataPoints;
};

#endif

