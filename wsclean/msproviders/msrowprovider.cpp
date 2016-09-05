#include "directmsrowprovider.h"
#include "msprovider.h"

#include "../wsclean/logger.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

MSRowProvider::MSRowProvider(const string& msPath, const MSSelection& selection, const std::map<size_t,size_t>& selectedDataDescIds, const std::string &dataColumnName, bool requireModel) :
	_ms(casacore::MeasurementSet(msPath)),
	_antenna1Column(casacore::ROScalarColumn<int> (_ms, casacore::MS::columnName(casacore::MSMainEnums::ANTENNA1))),
	_antenna2Column(casacore::ROScalarColumn<int>(_ms, casacore::MS::columnName(casacore::MSMainEnums::ANTENNA2))),
	_fieldIdColumn(casacore::ROScalarColumn<int>(_ms, casacore::MS::columnName(casacore::MSMainEnums::FIELD_ID))),
	_timeColumn(casacore::ROScalarColumn<double>(_ms, casacore::MS::columnName(casacore::MSMainEnums::TIME))),
	_timeEpochColumn(casacore::MEpoch::ROScalarColumn(_ms, casacore::MS::columnName(casacore::MSMainEnums::TIME))),
	_uvwColumn(casacore::ROArrayColumn<double>(_ms, casacore::MS::columnName(casacore::MSMainEnums::UVW))),
	_dataColumn(casacore::ROArrayColumn<casacore::Complex>(_ms, dataColumnName)),
	_flagColumn(casacore::ROArrayColumn<bool>(_ms, casacore::MS::columnName(casacore::MSMainEnums::FLAG))),
	_dataDescIdColumn(casacore::ROScalarColumn<int>(_ms, casacore::MS::columnName(casacore::MSMainEnums::DATA_DESC_ID))),
	_selectedDataDescIds(selectedDataDescIds),
	_selection(selection),
	_requireModel(requireModel)
{
	if(requireModel)
		_modelColumn.reset(new casacore::ROArrayColumn<casacore::Complex>(_ms, casacore::MS::columnName(casacore::MSMainEnums::MODEL_DATA)));
	
	const casacore::IPosition shape(DataShape());
	_msHasWeights = MSProvider::openWeightSpectrumColumn(_ms, _weightSpectrumColumn, shape);
	if(!_msHasWeights)
	{
		casacore::IPosition scalarShape(1, shape[0]);
		_scratchWeightScalarArray = casacore::Array<float>(scalarShape);
		_weightScalarColumn.reset(new casacore::ROArrayColumn<float>(_ms, casacore::MS::columnName(casacore::MSMainEnums::WEIGHT)));
	}

	MSProvider::getRowRange(_ms, selection, _startRow, _endRow);
	
	// Determine last timestep
	_startTimestep = selection.HasInterval() ? selection.IntervalStart() : 0;
	size_t timestep = _startTimestep;
	double time = _timeColumn(_startRow);
	for(size_t row=_startRow; row!=_endRow; ++row)
	{
		if(time != _timeColumn(row))
		{
			++timestep;
			time = _timeColumn(row);
		}
	}
	_endTimestep = timestep;
	
	_currentRow = _startRow;
	_currentTimestep = _startTimestep;
	_currentTime = _timeColumn(_startRow);
	_currentUVWArray = _uvwColumn(_startRow);
	_currentDataDescId = _dataDescIdColumn(_startRow);
	
	// If this row is not selected, it is necessary to continue to the first
	// selected row.
	const int
		a1 = _antenna1Column(_currentRow),
		a2 = _antenna2Column(_currentRow),
		fieldId = _fieldIdColumn(_currentRow);
	if(!isCurrentRowSelected(fieldId, a1, a2))
		NextRow();
}

void MSRowProvider::NextRow()
{
	bool isRowSelected;
	do {
		++_currentRow;
		if(_currentRow == _endRow)
		{
			break;
		}
		else {
			const int
				a1 = _antenna1Column(_currentRow),
				a2 = _antenna2Column(_currentRow),
				fieldId = _fieldIdColumn(_currentRow);
			_currentDataDescId = _dataDescIdColumn(_currentRow);
				
			if(_currentTime != _timeColumn(_currentRow))
			{
				++_currentTimestep;
				_currentTime = _timeColumn(_currentRow);
			}
			_currentUVWArray = _uvwColumn(_currentRow);
			isRowSelected = isCurrentRowSelected(fieldId, a1, a2);
		}
	} while(!isRowSelected);
}

bool MSRowProvider::isCurrentRowSelected(int fieldId, int a1, int a2) const
{
	std::map<size_t,size_t>::const_iterator dataDescIdIter = _selectedDataDescIds.find(_currentDataDescId);
	bool isDataDescIdSelected = dataDescIdIter!=_selectedDataDescIds.end();
	return _selection.IsSelected(fieldId, _currentTimestep, a1, a2, _currentUVWArray) && isDataDescIdSelected;
}

void MSRowProvider::getCurrentWeights(WeightArray& weights)
{
	if(_msHasWeights)
		_weightSpectrumColumn->get(_currentRow, weights);
	else {
		_weightScalarColumn->get(_currentRow, _scratchWeightScalarArray);
		MSProvider::expandScalarWeights(_scratchWeightScalarArray, weights);
	}
}
