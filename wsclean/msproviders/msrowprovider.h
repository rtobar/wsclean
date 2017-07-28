#ifndef MS_ROW_PROVIDER_H
#define MS_ROW_PROVIDER_H

#include "../msselection.h"

#include <casacore/casa/Arrays/Array.h>
#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>

#include <cstring>
#include <map>
#include <memory>
#include <vector>

/**
 * An MSRowProvider provides the selected rows of a data set.
 */
class MSRowProvider
{
public:
	MSRowProvider(const string& msPath, const MSSelection& selection, const std::map<size_t,size_t>& selectedDataDescIds, const std::string &dataColumnName, bool requireModel);
	virtual ~MSRowProvider() { }
	
	typedef casacore::Array<std::complex<float>> DataArray;
	typedef casacore::Array<float> WeightArray;
	typedef casacore::Array<bool> FlagArray;
	
	virtual bool AtEnd() const { return _currentRow == _endRow; }
	
	virtual void NextRow();
	
	virtual void ReadData(DataArray& data, FlagArray& flags, WeightArray& weights, double& u, double& v, double& w, uint32_t& dataDescId, uint32_t& antenna1, uint32_t& antenna2, double& time) = 0;
	
	virtual void ReadModel(DataArray& model) = 0;
	
	virtual void OutputStatistics() const { }

	casacore::MeasurementSet& MS() { return _ms; }
	casacore::IPosition DataShape() const { return _dataColumn.shape(0); }
	double StartTime() const { return _timeEpochColumn(_startRow).getValue().get(); }
	double EndTime() const { return _timeEpochColumn(_endRow-1).getValue().get(); }
	size_t StartTimestep() const { return _startTimestep; }
	size_t EndTimestep() const { return _endTimestep; }
	
	size_t CurrentProgress() const { return _currentRow-_startRow; }
	size_t TotalProgress() const { return _endRow-_startRow; }
	
protected:
	casacore::MeasurementSet _ms;
	casacore::ROScalarColumn<int> _antenna1Column;
	casacore::ROScalarColumn<int> _antenna2Column;
	casacore::ROScalarColumn<int> _fieldIdColumn;
	casacore::ROScalarColumn<double> _timeColumn;
	casacore::MEpoch::ROScalarColumn _timeEpochColumn;
	casacore::ROArrayColumn<double> _uvwColumn;
	std::unique_ptr<casacore::ROArrayColumn<float>> _weightSpectrumColumn;
	std::unique_ptr<casacore::ROArrayColumn<float>> _weightScalarColumn;
	casacore::ROArrayColumn<casacore::Complex> _dataColumn;
	casacore::ROArrayColumn<bool> _flagColumn;
	casacore::ROScalarColumn<int> _dataDescIdColumn;
	std::unique_ptr<casacore::ROArrayColumn<casacore::Complex>> _modelColumn;

	void getCurrentWeights(WeightArray& weights);
	const std::map<size_t,size_t>& selectedDataDescIds() const { return _selectedDataDescIds; }
	bool requireModel() const { return _requireModel; }
	
	bool isCurrentRowSelected(int fieldId, int a1, int a2) const;
	
	size_t _currentRow;
	size_t _currentTimestep;
	double _currentTime;
	casacore::Vector<double> _currentUVWArray;
	size_t _currentDataDescId;
	
private:
	std::map<size_t,size_t> _selectedDataDescIds;
	MSSelection _selection;
	bool _msHasWeights;
	bool _requireModel;
	casacore::Array<float> _scratchWeightScalarArray;
	
	size_t _startRow, _endRow;
	size_t _startTimestep, _endTimestep;	
};

#endif
