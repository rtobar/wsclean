#include <complex>

#ifndef INTERFACE_H
#define INTERFACE_H

class HighLevelGridderInterface
{
public:
    virtual ~HighLevelGridderInterface() { }

    virtual void set_frequencies(const double* frequencyList, size_t channelCount) = 0;

    virtual void set_stations(const size_t nStations) = 0;

		/**
		 * Set the antialiasing kernel. The kernel should be a kernelSize x kernelSize
		 * array of doubles. The kernel size also sets the size of the aterm kernel.
		 * This method should be called before @ref start_aterm(). The ideal kernel
		 * is a spheroidal kernel.
		 * @param kernelSize Size of the kernel (for IDG: the subgrid size)
		 * @param kernel Array of kernelSize x kernelSize with the kernel.
		 */
    virtual void set_kernel(size_t kernelSize, const double* kernel) = 0;

    virtual void start_w_layer(double layerWInLambda) = 0;

    virtual void finish_w_layer() = 0;

		/**
		 * With K the kernel size, P the number of polarizations(4), S the number of
		 * stations, aterm is an array with dimensions [S][P][K][K].
		 */
    virtual void start_aterm(const std::complex<double>* aterm) = 0;

    virtual void finish_aterm() = 0;

    virtual void set_grid(std::complex<double>* grid) = 0;

    virtual void grid_visibility(
        const std::complex<float>* visibility, // size CH x PL
        const double* uvwInMeters,
        size_t antenna1,
        size_t antenna2,
        size_t timeIndex
    ) = 0;
		
		/**
		 * Performs the FFT on the grid. To be called once all visibilities
		 * have been gridded.
		 */
		virtual void transform_grid_after_gridding() = 0;
		
		/**
		 */
		virtual void transform_grid_before_sampling() = 0;

    virtual void queue_visibility_sampling(
       const double* uvwInMeters,
        size_t antenna1,
        size_t antenna2,
        size_t timeIndex,
        size_t rowId,
        bool& isBufferFull
    ) = 0;

		/**
		 * Prepare the sampling buffer for reading out the sampled visibilities.
		 * This method must be called when @ref queue_visibility_sampling() has
		 * indicated the buffer is full, but can also be called before the buffer
		 * is full. 
		 * Internally, for the IDG this method transposes the data to have the order
		 * as required.
		 */
    virtual void finish_sampled_visibilities() = 0;

    /**
    * Get the visibilities queued for sampling. This fills an
    * CH x PL array. PL is always four, CH is
		* as set in @ref set_frequencies().
 		* This buffer is valid only directly after @ref finish_sampled_visibilities()
		* has been called. A call to @ref queue_visibility_sampling() invalidates the buffer.
		* @returns Pointer to the array with visibilities.
    */
		virtual void get_sampled_visibilities(size_t index, std::complex<float>* data, size_t& rowID) const = 0;
		
		/**
		 * Get the number of rows in the buffer: TI x BL.
		 * Only to be called after @ref finish_sampled_visibilities().
		 * (in principal the client can know this number, but it is
		 * convenient to have it here as well).
		 */
		virtual size_t get_sampling_buffer_size() const = 0;
}; 

#endif
