#include "directmsrowprovider.h"

#include "msprovider.h"

void DirectMSRowProvider::ReadData(MSRowProvider::DataArray& data, MSRowProvider::FlagArray& flags, WeightArray& weights, double& u, double& v, double& w, uint32_t& dataDescId, uint32_t& antenna1, uint32_t& antenna2)
{
	u = _currentUVWArray(0);
	v = _currentUVWArray(1);
	w = _currentUVWArray(2);
	dataDescId = _currentDataDescId;
	_dataColumn.get(_currentRow, data);
	_flagColumn.get(_currentRow, flags);
	antenna1 = _antenna1Column(_currentRow);
	antenna2 = _antenna2Column(_currentRow);
	
	getCurrentWeights(weights);
}

void DirectMSRowProvider::ReadModel(MSRowProvider::DataArray& model)
{
	_modelColumn->get(_currentRow, model);
}
