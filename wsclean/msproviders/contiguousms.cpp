#include "contiguousms.h"
#include "../wsclean/logger.h"
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>

ContiguousMS::ContiguousMS(const string& msPath, const std::string& dataColumnName, const MSSelection& selection, PolarizationEnum polOut, size_t dataDescId, bool includeModel) :
	_timestep(0),
	_time(0.0),
	_dataDescId(dataDescId),
	_isModelColumnPrepared(false),
	_selection(selection),
	_polOut(polOut),
	_msPath(msPath),
	_ms(msPath),
	_antenna1Column(_ms, casacore::MS::columnName(casacore::MSMainEnums::ANTENNA1)),
	_antenna2Column(_ms, casacore::MS::columnName(casacore::MSMainEnums::ANTENNA2)),
	_fieldIdColumn(_ms, casacore::MS::columnName(casacore::MSMainEnums::FIELD_ID)),
	_dataDescIdColumn(_ms, casacore::MS::columnName(casacore::MSMainEnums::DATA_DESC_ID)),
	_timeColumn(_ms, casacore::MS::columnName(casacore::MSMainEnums::TIME)),
	_uvwColumn(_ms, casacore::MS::columnName(casacore::MSMainEnums::UVW)),
	_dataColumnName(dataColumnName),
	_dataColumn(_ms, dataColumnName),
	_flagColumn(_ms, casacore::MS::columnName(casacore::MSMainEnums::FLAG))
{
	Logger::Info << "Opening " << msPath << ", spw " << _dataDescId << " with contiguous MS reader.\n";
	
	_inputPolarizations = GetMSPolarizations(_ms);
 
	const casacore::IPosition shape(_dataColumn.shape(0));
	_dataArray = casacore::Array<std::complex<float>>(shape);
	_weightSpectrumArray = casacore::Array<float>(shape);
	_flagArray = casacore::Array<bool>(shape);
	_bandData = MultiBandData(_ms.spectralWindow(), _ms.dataDescription());
	
	_msHasWeightSpectrum = openWeightSpectrumColumn(_ms, _weightSpectrumColumn, shape);
	if(!_msHasWeightSpectrum)
	{
		casacore::IPosition scalarShape(1, shape[0]);
		_weightScalarArray = casacore::Array<float>(scalarShape);
		_weightScalarColumn.reset(new casacore::ROArrayColumn<float>(_ms, casacore::MS::columnName(casacore::MSMainEnums::WEIGHT)));
	}
	
	getRowRangeAndIDMap(_ms, selection, _startRow, _endRow, std::set<size_t>{dataDescId}, _idToMSRow);
	Reset();
}

void ContiguousMS::Reset()
{
	_row = _startRow - 1;
	_rowId = size_t(-1);
	_time = 0.0;
	if(_selection.HasInterval())
		_timestep = _selection.IntervalStart()-1;
	else
		_timestep = -1;
	NextRow();
}

bool ContiguousMS::CurrentRowAvailable()
{
	if(_row >= _endRow)
		return false;
		
	int fieldId = _fieldIdColumn(_row);
	int a1 = _antenna1Column(_row);
	int a2 = _antenna2Column(_row);
	int dataDescId = _dataDescIdColumn(_row);
	casacore::Vector<double> uvw = _uvwColumn(_row);
	
	while(!_selection.IsSelected(fieldId, _timestep, a1, a2, uvw) || dataDescId != _dataDescId) {
		++_row;
		if(_row >= _endRow)
			return false;
		
		fieldId = _fieldIdColumn(_row);
		a1 = _antenna1Column(_row);
		a2 = _antenna2Column(_row);
		uvw = _uvwColumn(_row);
		dataDescId = _dataDescIdColumn(_row);
		if(_time != _timeColumn(_row))
		{
			++_timestep;
			_time = _timeColumn(_row);
		}
		
		_isMetaRead = false;
		_isDataRead = false;
		_isWeightRead = false;
		_isModelRead = false;
	}
	
	return true;
}

void ContiguousMS::NextRow()
{
	_isMetaRead = false;
	_isDataRead = false;
	_isWeightRead = false;
	_isModelRead = false;
	
	++_rowId;
	int fieldId, a1, a2, dataDescId;
	casacore::Vector<double> uvw;
	do {
		++_row;
		if(_row >= _endRow)
			return;
		
		fieldId = _fieldIdColumn(_row);
		a1 = _antenna1Column(_row);
		a2 = _antenna2Column(_row);
		uvw = _uvwColumn(_row);
		dataDescId = _dataDescIdColumn(_row);
		if(_time != _timeColumn(_row))
		{
			++_timestep;
			_time = _timeColumn(_row);
		}
	} while(!_selection.IsSelected(fieldId, _timestep, a1, a2, uvw) || (dataDescId != _dataDescId) );
}

double ContiguousMS::StartTime()
{
	return casacore::MEpoch::ROScalarColumn(_ms, casacore::MS::columnName(casacore::MS::TIME))(_startRow).getValue().get();
}

void ContiguousMS::ReadMeta(double& u, double& v, double& w, size_t& dataDescId)
{
	readMeta();
	
	casacore::Vector<double> uvwArray = _uvwColumn(_row);
	u = uvwArray(0);
	v = uvwArray(1);
	w = uvwArray(2);
	dataDescId = _dataDescId;
}

void ContiguousMS::ReadMeta(MetaData& metaData)
{
	readMeta();
	
	casacore::Vector<double> uvwArray = _uvwColumn(_row);
	metaData.uInM = uvwArray(0);
	metaData.vInM = uvwArray(1);
	metaData.wInM = uvwArray(2);
	metaData.dataDescId = _dataDescId;
	metaData.antenna1 = _antenna1Column(_row);
	metaData.antenna2 = _antenna2Column(_row);
	metaData.time = _timeColumn(_row);
}

void ContiguousMS::ReadData(std::complex<float>* buffer)
{
	readMeta();
	readData();
	readWeights();
	size_t startChannel, endChannel;
	if(_selection.HasChannelRange())
	{
		startChannel = _selection.ChannelRangeStart();
		endChannel = _selection.ChannelRangeEnd();
	}
	else {
		startChannel = 0;
		endChannel = _bandData[_dataDescId].ChannelCount();
	}
	copyWeightedData(buffer,  startChannel, endChannel, _inputPolarizations, _dataArray, _weightSpectrumArray, _flagArray, _polOut);
}

void ContiguousMS::prepareModelColumn()
{
	initializeModelColumn(_ms);
	
	_modelColumn.reset(new casacore::ArrayColumn<casacore::Complex>(_ms, casacore::MS::columnName(casacore::MSMainEnums::MODEL_DATA)));
	const casacore::IPosition shape(_modelColumn->shape(0));
	_modelArray = casacore::Array<std::complex<float>>(shape);
	_isModelColumnPrepared = true;
}

void ContiguousMS::ReadModel(std::complex<float>* buffer)
{
	if(!_isModelColumnPrepared)
		prepareModelColumn();
	
	readMeta();
	readModel();
	readWeights();
	size_t startChannel, endChannel;
	if(_selection.HasChannelRange())
	{
		startChannel = _selection.ChannelRangeStart();
		endChannel = _selection.ChannelRangeEnd();
	}
	else {
		startChannel = 0;
		endChannel = _bandData[_dataDescId].ChannelCount();
	}
	copyWeightedData(buffer,  startChannel, endChannel, _inputPolarizations, _modelArray, _weightSpectrumArray, _flagArray, _polOut);
}

void ContiguousMS::WriteModel(size_t rowId, std::complex<float>* buffer)
{
	if(!_isModelColumnPrepared)
		prepareModelColumn();
	
	size_t msRowId = _idToMSRow[rowId];
	size_t dataDescId = _dataDescIdColumn(msRowId);
	size_t startChannel, endChannel;
	if(_selection.HasChannelRange())
	{
		startChannel = _selection.ChannelRangeStart();
		endChannel = _selection.ChannelRangeEnd();
	}
	else {
		startChannel = 0;
		endChannel = _bandData[dataDescId].ChannelCount();
	}
	
	_modelColumn->get(msRowId, _modelArray);
	reverseCopyData(_modelArray, startChannel, endChannel, _inputPolarizations, buffer, _polOut);
	_modelColumn->put(msRowId, _modelArray);
}

void ContiguousMS::ReadWeights(std::complex<float>* buffer)
{
	readMeta();
	readData();
	readWeights();
	size_t startChannel, endChannel;
	if(_selection.HasChannelRange())
	{
		startChannel = _selection.ChannelRangeStart();
		endChannel = _selection.ChannelRangeEnd();
	}
	else {
		startChannel = 0;
		endChannel = _bandData[_dataDescId].ChannelCount();
	}
	copyWeights(buffer,  startChannel, endChannel, _inputPolarizations, _dataArray, _weightSpectrumArray, _flagArray, _polOut);
}

void ContiguousMS::ReadWeights(float* buffer)
{
	readMeta();
	readData();
	readWeights();
	size_t startChannel, endChannel;
	if(_selection.HasChannelRange())
	{
		startChannel = _selection.ChannelRangeStart();
		endChannel = _selection.ChannelRangeEnd();
	}
	else {
		startChannel = 0;
		endChannel = _bandData[_dataDescId].ChannelCount();
	}
	copyWeights(buffer,  startChannel, endChannel, _inputPolarizations, _dataArray, _weightSpectrumArray, _flagArray, _polOut);
}

void ContiguousMS::MakeIdToMSRowMapping(vector<size_t>& idToMSRow)
{
	idToMSRow = _idToMSRow;
}
