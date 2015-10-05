#ifndef POLYNOMIAL_FITTER_H
#define POLYNOMIAL_FITTER_H

#include "uvector.h"

class PolynomialFitter
{
public:
	void AddDataPoint(double x, double y)
	{
		_dataPoints.push_back(std::make_pair(x, y));
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
	ao::uvector<std::pair<double,double>> _dataPoints;
};

#endif
