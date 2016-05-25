#ifndef DUMMY_GRIDDER_H
#define DUMMY_GRIDDER_H

#include <complex>
#include <vector>

#include "interface.h"

class DummyGridder : public HighLevelGridderInterface
{
public:
	virtual ~DummyGridder() { }

	virtual void set_frequencies(const double* frequencyList, size_t channelCount) override final
	{
		_channelCount = channelCount;
	}

	virtual void set_stations(const size_t nStations) override final
	{ }

	virtual void set_kernel(size_t kernelSize, const double* kernel) override final
	{ }

	virtual void start_w_layer(double layerWInLambda) override final
	{ }

	virtual void finish_w_layer() override final
	{ }

	virtual void start_aterm(const std::complex<double>* aterm) override final
	{ }

	virtual void finish_aterm() override final
	{ }

	virtual void set_grid(std::complex<double>* grid) override final
	{ }

	virtual void grid_visibility(
			const std::complex<float>* visibility, // size CH x PL
			const double* uvwInMeters,
			size_t antenna1,
			size_t antenna2,
			size_t timeIndex
	) override final
	{ }
	
	virtual void transform_grid_after_gridding() override final
	{ }
	
	virtual void transform_grid_before_sampling() override final
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

	virtual void finish_sampled_visibilities() override final
	{ }

	virtual void get_sampled_visibilities(size_t index, std::complex<float>* data, size_t& rowID) const override final
	{
		for(size_t i=0; i!=_channelCount*4; ++i)
			data[i] = 1.0;
		rowID = rowIds[index];
	}
	
	virtual size_t get_sampling_buffer_size() const override final
	{
		return rowIds.size();
	}
private:
	std::vector<size_t> rowIds;
	size_t _channelCount;
};

#endif
