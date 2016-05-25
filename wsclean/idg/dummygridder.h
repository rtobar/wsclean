#ifndef DUMMY_GRIDDER_H
#define DUMMY_GRIDDER_H

#include <complex>
#include <vector>

#include "interface.h"

class DummyGridder : public HighLevelGridderInterface
{
public:
	virtual ~DummyGridder() { }

	virtual void set_frequencies(const double* frequencyList, size_t channelCount)  
	{
		_channelCount = channelCount;
	}

	virtual void set_stations(const size_t nStations)  
	{ }

	virtual void set_kernel(size_t kernelSize, const double* kernel)  
	{ }

	virtual void start_w_layer(double layerWInLambda)  
	{ }

	virtual void finish_w_layer()  
	{ }

	virtual void start_aterm(const std::complex<double>* aterm)  
	{ }

	virtual void finish_aterm()  
	{ }

	virtual void set_grid(std::complex<double>* grid)  
	{ }

	virtual void grid_visibility(
			const std::complex<float>* visibility, // size CH x PL
			const double* uvwInMeters,
			size_t antenna1,
			size_t antenna2,
			size_t timeIndex
	)  
	{ }
	
	virtual void transform_grid_after_gridding()  
	{ }
	
	virtual void transform_grid_before_sampling()  
	{ }

	virtual void queue_visibility_sampling(
			const double* uvwInMeters,
			size_t antenna1,
			size_t antenna2,
			size_t timeIndex,
			size_t rowId,
			bool& isBufferFull
	)
	{
		rowIds.push_back(rowId);
	}

	virtual void finish_sampled_visibilities()  
	{ }

	virtual void get_sampled_visibilities(size_t index, std::complex<float>* data, size_t& rowID) const  
	{
		for(size_t i=0; i!=_channelCount*4; ++i)
			data[i] = 1.0;
		rowID = rowIds[index];
	}
	
	virtual size_t get_sampling_buffer_size() const  
	{
		return rowIds.size();
	}
private:
	std::vector<size_t> rowIds;
	size_t _channelCount;
};

#endif
