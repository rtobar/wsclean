#ifndef FFTW_MULTI_THREAD_ENABLER_H
#define FFTW_MULTI_THREAD_ENABLER_H

#include <cstring>

/**
 * Used to initialize and enable fftw's multithreading. While
 * an instance of this class exists, multi-threading is enabled. When the
 * class is destructed, multi-threading is again disabled.
 * 
 * To make the FFTs in for example the @ref FFTConvolver or @ref FFTResampler
 * multi-threaded, it is enough to construct an instance of this class and keep
 * it until done.
 * For example:
 * 
 * @code
 * void doMultiThreadedFFTActions
 * {
 *   FFTWMultiThreadEnabler fftwMT;
 * 
 *   FFTResampler resampler(..);
 *   ...
 * }
 * @endcode
 */
class FFTWMultiThreadEnabler
{
public:
	/**
	 * Constructor that sets FFTW to use multiple threads.
	 * This will set FFTW to use as many threads as there are cores in the 
	 * system.
	 */
	explicit FFTWMultiThreadEnabler(bool reportNrThreads = true);
	
	/**
	 * Constructor that sets FFTW to use multiple threads.
	 * This will set FFTW to use the given number of threads.
	 */
	explicit FFTWMultiThreadEnabler(size_t nThreads, bool reportNrThreads = true);
	
	/**
	 * Destructor that resets the FFTWs threads.
	 */
	~FFTWMultiThreadEnabler();
};

#endif
