#ifndef NOISE_MS_ROW_PROVIDER_H
#define NOISE_MS_ROW_PROVIDER_H

#include "directmsrowprovider.h"

#include <random>

class NoiseMSRowProvider : public DirectMSRowProvider
{
public:
	NoiseMSRowProvider(double noiseStdDevJy, const string& msPath, const MSSelection& selection, const std::map<size_t,size_t>& selectedDataDescIds, const std::string &dataColumnName, bool requireModel) :
	DirectMSRowProvider(msPath, selection, selectedDataDescIds, dataColumnName, requireModel),
	_rng(std::random_device{}()),
	_distribution(0.0, noiseStdDevJy)
	{ }
	
	virtual void ReadData(DataArray& data, FlagArray& flags, WeightArray& weights, double& u, double& v, double& w, uint32_t& dataDescId)
	{
		DirectMSRowProvider::ReadData(data, flags, weights, u, v, w, dataDescId);
		for(DataArray::contiter iter = data.cbegin(); iter != data.cend(); ++iter)
		{
			if(std::isfinite(iter->real()) && std::isfinite(iter->imag()))
			{
				iter->real(_distribution(_rng));
				iter->imag(_distribution(_rng));
			}
			else {
				iter->real(std::numeric_limits<float>::quiet_NaN());
				iter->imag(std::numeric_limits<float>::quiet_NaN());
			}
		}
	}
	
private:
	std::mt19937 _rng;
	std::normal_distribution<float> _distribution;
};

#endif
