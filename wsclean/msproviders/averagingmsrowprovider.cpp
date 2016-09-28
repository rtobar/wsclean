#include "averagingmsrowprovider.h"

#include "../multibanddata.h"
#include "../wsclean/logger.h"

#include <casacore/tables/Tables/ArrayColumn.h>

AveragingMSRowProvider::AveragingMSRowProvider(double nWavelengthsAveraging, const string& msPath, const MSSelection& selection, const std::map<size_t, size_t>& selectedDataDescIds, const string& dataColumnName, bool requireModel) :
	MSRowProvider(msPath, selection, selectedDataDescIds, dataColumnName, requireModel)
{
	casacore::MSAntenna antennaTable(_ms.antenna());
	_nAntennae = antennaTable.nrow();
	
	casacore::ROArrayColumn<double> positionColumn(antennaTable, casacore::MSAntenna::columnName(casacore::MSAntennaEnums::POSITION));
	std::vector<Pos> positions(_nAntennae);

	casa::Array<double> posArr(casacore::IPosition(1, 3));
	for(size_t i=0; i!=_nAntennae; ++i)
	{
		positionColumn.get(i, posArr);
		positions[i] = Pos(posArr.data()[0], posArr.data()[1], posArr.data()[2]);
	}
	
	// dataDescId x ant x ant
	_nElements = selectedDataDescIds.size() * _nAntennae * _nAntennae;
	_averagingFactors.assign(_nElements, 0.0);
	_buffers.resize(_nElements);
	MultiBandData bands(_ms.spectralWindow(), _ms.dataDescription());
	
	double dt = (EndTime() - StartTime()) / (EndTimestep() - StartTimestep());
	Logger::Debug << "Assuming integration time of " << dt * (24.0*60.0*60.0) << " seconds.\n";
	
	size_t element = 0;
	size_t averagingSum = 0, minAvgFactor = std::numeric_limits<size_t>::max(), maxAvgFactor = 0;
	for(size_t a1=0; a1!=_nAntennae; ++a1)
	{
		Pos pos1 = positions[a1];
		for(size_t a2=0; a2!=_nAntennae; ++a2)
		{
			Pos pos2 = positions[a2];
			double dx = std::get<0>(pos1) - std::get<0>(pos2);
			double dy = std::get<1>(pos1) - std::get<1>(pos2);
			double dz = std::get<2>(pos1) - std::get<2>(pos2);
			double dist = sqrt(dx*dx + dy*dy + dz*dz);
			for(std::map<size_t, size_t>::const_iterator spwIter=selectedDataDescIds.begin();
					spwIter!=selectedDataDescIds.end(); ++spwIter)
			{
				BandData band = bands[spwIter->first];
				double lambda = band.SmallestWavelength();
				double nWavelengthsPerIntegration = 2.0 * M_PI * dist / lambda * dt;
				_averagingFactors[element] = std::max<size_t>(size_t(floor(nWavelengthsAveraging / nWavelengthsPerIntegration)), 1);
				averagingSum += _averagingFactors[element];
				if(a1 != a2)
				{
					minAvgFactor = std::min<size_t>(minAvgFactor, _averagingFactors[element]);
					maxAvgFactor = std::max<size_t>(maxAvgFactor, _averagingFactors[element]);
				}
				//Logger::Debug << a1 << '\t' << a2 << '\t' << _averagingFactors[element] << '\n';
				++element;
			}
		}
	}
	Logger::Info << "Averaging factor for longest baseline: " << minAvgFactor << " x . For the shortest: " << maxAvgFactor << " x \n";
	
	_spwIndexToDataDescId.resize(selectedDataDescIds.size());
	for(std::map<size_t, size_t>::const_iterator spwIter=selectedDataDescIds.begin();
		spwIter!=selectedDataDescIds.end(); ++spwIter)
	{
		_spwIndexToDataDescId[spwIter->second] = spwIter->first;
	}
	
	_averageFactorSum = 0.0;
	_rowCount = 0;
	_averagedRowCount = 0;

	_currentData = DataArray(DataShape());
	_currentModel = DataArray(DataShape());
	_currentFlags = FlagArray(DataShape());
	_currentWeights = WeightArray(DataShape());
	_averagedDataDescId = _currentDataDescId;
	_flushPosition = 0;
	
	if(!MSRowProvider::AtEnd())
	{
		bool timestepAvailable = processCurrentTimestep();
		if(!timestepAvailable)
			NextRow();
	}
}

bool AveragingMSRowProvider::processCurrentTimestep()
{
	size_t a1 = _antenna1Column(_currentRow);
	size_t a2 = _antenna2Column(_currentRow);
	_averagedDataDescId = _currentDataDescId;
	_averagedAntenna1Index = a1;
	_averagedAntenna2Index = a2;
	
	size_t spwCount = selectedDataDescIds().size();
	size_t elementIndex = spwCount*(a2 + a1*_nAntennae) + selectedDataDescIds().find(_averagedDataDescId)->second;
	size_t avgFactor = _averagingFactors[elementIndex];
	_averageFactorSum += avgFactor;
	++_rowCount;
	
	//if(a1==1 && a2==2)
	//	Logger::Debug << a1 <<'\t' << a2 << '\t' << avgFactor << '\n';
	
	_dataColumn.get(_currentRow, _currentData);
	_flagColumn.get(_currentRow, _currentFlags);
	getCurrentWeights(_currentWeights);
	if(requireModel())
		_modelColumn->get(_currentRow, _currentModel);
	
	if(avgFactor == 1)
		return true;
	else
	{
		size_t bufferSize = DataShape()[0] * DataShape()[1];
		AveragingBuffer& buffer = _buffers[elementIndex];
		if(!buffer.IsInitialized())
			buffer.Initialize(bufferSize, requireModel());
		
		if(requireModel())
			buffer.AddDataAndModel(bufferSize, _currentData.data(), _currentModel.data(), _currentFlags.data(), _currentWeights.data(), _currentUVWArray.data());
		else
			buffer.AddData(bufferSize, _currentData.data(), _currentFlags.data(), _currentWeights.data(), _currentUVWArray.data());
		
		bool foundFullBuffer = (buffer.AveragedDataCount() == avgFactor);
		if(foundFullBuffer)
		{
			if(requireModel())
				buffer.Get(bufferSize, _currentData.data(), _currentModel.data(), _currentFlags.data(), _currentWeights.data(), _currentUVWArray.data());
			else
				buffer.Get(bufferSize, _currentData.data(), _currentFlags.data(), _currentWeights.data(), _currentUVWArray.data());
			buffer.Reset(bufferSize);
		}
		return foundFullBuffer;
	}
}

void AveragingMSRowProvider::NextRow()
{
	++_averagedRowCount;
	if(!MSRowProvider::AtEnd())
	{
		bool foundFullBuffer = false;
		while(!foundFullBuffer)
		{
			MSRowProvider::NextRow();
			if(MSRowProvider::AtEnd())
				break;
			
			foundFullBuffer = processCurrentTimestep();
		}
		// TODO I think this should now say "if(foundFullBuffer) return; "
	}
	
	if(MSRowProvider::AtEnd())
	{
		// There might be residual data in the buffers which have to be read out
		AveragingBuffer* buffer = 0;
		do {
			buffer = &_buffers[_flushPosition];
			++_flushPosition;
		} while(buffer->AveragedDataCount() == 0 && _flushPosition < _nElements);
		
		if(buffer!= 0 && buffer->AveragedDataCount() != 0)
		{
			size_t bufferSize = DataShape()[0] * DataShape()[1];
			if(requireModel())
				buffer->Get(bufferSize, _currentData.data(), _currentModel.data(), _currentFlags.data(), _currentWeights.data(), _currentUVWArray.data());
			else
				buffer->Get(bufferSize, _currentData.data(), _currentFlags.data(), _currentWeights.data(), _currentUVWArray.data());
			
			size_t elementIndex = _flushPosition-1;
			size_t spwCount = selectedDataDescIds().size();
			size_t spwIndex = elementIndex % spwCount;
			size_t subIndex = elementIndex / spwCount;
			_averagedAntenna1Index = subIndex/_nAntennae;
			_averagedAntenna2Index = subIndex%_nAntennae;
			_averagedDataDescId = _spwIndexToDataDescId[spwIndex];
		}
	}
}

void AveragingMSRowProvider::ReadData(MSRowProvider::DataArray& data, MSRowProvider::FlagArray& flags, MSRowProvider::WeightArray& weights, double& u, double& v, double& w, uint32_t& dataDescId, uint32_t& antenna1, uint32_t& antenna2)
{
	size_t bufferSize = DataShape()[0] * DataShape()[1];
	memcpy(data.data(), _currentData.data(), bufferSize*sizeof(std::complex<float>));
	memcpy(flags.data(), _currentFlags.data(), bufferSize*sizeof(bool));
	memcpy(weights.data(), _currentWeights.data(), bufferSize*sizeof(float));
	u = _currentUVWArray.data()[0];
	v = _currentUVWArray.data()[1];
	w = _currentUVWArray.data()[2];
	dataDescId = _averagedDataDescId;
	antenna1 = _averagedAntenna1Index;
	antenna2 = _averagedAntenna2Index;
}

void AveragingMSRowProvider::ReadModel(MSRowProvider::DataArray& model)
{
	size_t bufferSize = DataShape()[0] * DataShape()[1];
	memcpy(model.data(), _currentModel.data(), bufferSize*sizeof(std::complex<float>));
}

void AveragingMSRowProvider::OutputStatistics() const
{
	//Logger::Info << "Average selected integration-time averaging factor after row selection: " << double(_averageFactorSum)/_rowCount << '\n';
	Logger::Info << "Baseline averaging reduced the number of rows to " << 0.1*round(_averagedRowCount*1000.0 / _rowCount) << "%.\n";
}
