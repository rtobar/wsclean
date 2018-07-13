#ifndef ATERM_BASE_H
#define ATERM_BASE_H

#include <complex>

class ATermBase
{
public:
	/**
	 * Calculate the a-terms for the given time and frequency, for all stations.
	 * @param buffer A buffer of size width x height x 4 x nstation, in that order.
	 * @param time The time corresponding to the currently gridded visibilities to
	 * which the aterm will be applied.
	 * @param frequency Frequency of currently gridded visibilities.
	 * @returns @c True when new aterms are calculated. If these aterms are the same as
	 * for the previous call to Calculate(), @c false can be returned and the output
	 * buffer does not need to be updated. The gridder will then make sure to use the
	 * previous aterms, and not reserve extra memory for it etc.
	 */
	virtual bool Calculate(std::complex<float>* buffer, double time, double frequency) = 0;
};

#endif
