#ifndef ATERM_BASE_H
#define ATERM_BASE_H

#include <complex>

class ATermBase
{
public:
	virtual void Calculate(std::complex<float>* buffer, double time, double frequency) = 0;
};

#endif
