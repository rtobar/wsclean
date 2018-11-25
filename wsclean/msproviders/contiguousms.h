#ifndef CONTIGUOUSMS_H
#define CONTIGUOUSMS_H

#include "msprovider.h"

#include "../msselection.h"
#include "../multibanddata.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>

#include <memory>

class ContiguousMS : public MSProvider
{
public:
	ContiguousMS(const string& msPath, const std::string& dataColumnName, const MSSelection& selection, PolarizationEnum polOut, size_t dataDescIndex);
	
	ContiguousMS(const ContiguousMS&) = delete;
	
	ContiguousMS& operator=(const ContiguousMS&) = delete;
	
	casacore::MeasurementSet& MS() final override { return _ms; }
	
	const std::string& DataColumnName() final override { return _dataColumnName; }
	
	size_t RowId() const final override { return _rowId; }
	
	bool CurrentRowAvailable() final override;
	
	void NextRow() final override;
	
	void Reset() final override;
	
	void ReadMeta(double& u, double& v, double& w, size_t& dataDescId) final override;
	
	void ReadMeta(MetaData& metaData) final override;
	
	void ReadData(std::complex<float>* buffer) final override;
	
	void ReadModel(std::complex<float>* buffer) final override;
	
	void WriteModel(size_t rowId, std::complex<float>* buffer) final override;
	
	void ReadWeights(std::complex<float>* buffer) final override;
	
	void ReadWeights(float* buffer) final override;
	
	void WriteImagingWeights(size_t rowId, const float* buffer) final override;
	
	void ReopenRW() final override 
	{
		_ms.reopenRW();
	}
	
	double StartTime() final override;
	
	void MakeIdToMSRowMapping(std::vector<size_t>& idToMSRow) final override;
	
	PolarizationEnum Polarization() final override { return _polOut; }
private:
	size_t _row, _rowId;
	size_t _timestep;
	double _time;
	int _dataDescId;
	bool _isMetaRead, _isDataRead, _isModelRead, _isWeightRead;
	bool _isModelColumnPrepared;
	size_t _startRow, _endRow;
	std::vector<size_t> _idToMSRow;
	std::vector<PolarizationEnum> _inputPolarizations;
	MSSelection _selection;
	PolarizationEnum _polOut;
	std::string _msPath;
	casacore::MeasurementSet _ms;
	MultiBandData _bandData;
	bool _msHasWeightSpectrum;

	casacore::ScalarColumn<int> _antenna1Column, _antenna2Column, _fieldIdColumn, _dataDescIdColumn;
	casacore::ScalarColumn<double> _timeColumn;
	casacore::ArrayColumn<double> _uvwColumn;
	std::unique_ptr<casacore::ArrayColumn<float>> _weightSpectrumColumn;
	std::unique_ptr<casacore::ArrayColumn<float>> _weightScalarColumn;
	std::string _dataColumnName;
	casacore::ROArrayColumn<casacore::Complex> _dataColumn;
	casacore::ROArrayColumn<bool> _flagColumn;
	std::unique_ptr<casacore::ArrayColumn<casacore::Complex>> _modelColumn;
	std::unique_ptr<casacore::ArrayColumn<float>> _imagingWeightsColumn;
	
	casacore::Array<std::complex<float>> _dataArray, _modelArray;
	casacore::Array<float> _weightSpectrumArray, _weightScalarArray, _imagingWeightSpectrumArray;
	casacore::Array<bool> _flagArray;
	
	void prepareModelColumn();
	void readMeta()
	{
		if(!_isMetaRead)
		{
			_dataDescId = _dataDescIdColumn(_row);
			_isMetaRead = true;
		}
	}
	void readData()
	{
		if(!_isDataRead)
		{
			_dataColumn.get(_row, _dataArray);
			_isDataRead = true;
		}
	}
	void readWeights()
	{
		if(!_isWeightRead)
		{
			_flagColumn.get(_row, _flagArray);
			if(_msHasWeightSpectrum)
				_weightSpectrumColumn->get(_row, _weightSpectrumArray);
			else {
				_weightScalarColumn->get(_row, _weightScalarArray);
				expandScalarWeights(_weightScalarArray, _weightSpectrumArray);
			}
			_isWeightRead = true;
		}
	}
	void readModel()
	{
		if(!_isModelRead)
		{
			_modelColumn->get(_row, _modelArray);
			_isModelRead = true;
		}
	}
};

#endif
