#include "directmsrowprovider.h"

#include "msprovider.h"

void DirectMSRowProvider::ReadData(MSRowProvider::DataArray& data, MSRowProvider::FlagArray& flags, WeightArray& weights, double& u, double& v, double& w, uint32_t& dataDescId)
{
	u = _currentUVWArray(0);
	v = _currentUVWArray(1);
	w = _currentUVWArray(2);
	dataDescId = _currentDataDescId;
	_dataColumn.get(_currentRow, data);
	_flagColumn.get(_currentRow, flags);
	
	getCurrentWeights(weights);
}

void DirectMSRowProvider::ReadModel(MSRowProvider::DataArray& model)
{
	_modelColumn->get(_currentRow, model);
}
