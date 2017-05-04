#include "polynomialchannelfitter.h"

#include "gsl/gsl_multifit.h"

void PolynomialChannelFitter::Fit(ao::uvector<double>& terms, size_t nTerms)
{
	const size_t nPoints = _dataPoints.size();
	
	gsl_multifit_linear_workspace* work = gsl_multifit_linear_alloc(nPoints, nTerms);
	
	gsl_matrix * x = gsl_matrix_alloc(nPoints, nTerms);  
	gsl_vector * y = gsl_vector_alloc(nPoints);  
	//gsl_vector * w = gsl_vector_alloc(nTerms);
	gsl_vector *c = gsl_vector_calloc(nTerms);
	gsl_matrix *cov = gsl_matrix_calloc(nTerms, nTerms);   
	double chisq;
	
	for(size_t i=0; i!=nPoints; ++i)
	{
		size_t chIndex = _dataPoints[i].first;
		const double nuBegin = _channels[chIndex].first, nuEnd = _channels[chIndex].second;
		double nuBeginTerm = nuBegin, nuEndTerm = nuEnd;
		gsl_matrix_set(x, i, 0, 1.0);
		for(size_t j=1; j!=nTerms; ++j)
		{
			nuBeginTerm *= nuBegin;
			nuEndTerm *= nuEnd;
			// ( mu_i^t - nu_i^t ) / (t [ mu_i - nu_i ] ) with t=j+1
			double val = (nuEndTerm - nuBeginTerm) / (double(j+1)*(nuEnd - nuBegin));
			gsl_matrix_set(x, i, j, val);
		}
		gsl_vector_set(y, i, _dataPoints[i].second);
	}

	int res = gsl_multifit_linear(x, y, c, cov, &chisq, work);
	if(res != 0)
		throw std::runtime_error("Linear fit failed");
	terms.resize(nTerms);
	for(size_t j=0; j!=nTerms; ++j)
	{
		terms[j] = gsl_vector_get(c, j);
	}
}
