#ifndef AVERAGING_MS_ROW_PROVIDER_H
#define AVERAGING_MS_ROW_PROVIDER_H

#include "msrowprovider.h"

#include "../uvector.h"

class AveragingMSRowProvider : public MSRowProvider
{
public:
	AveragingMSRowProvider(double nWavelengthsAveraging, const string& msPath, const MSSelection& selection, const std::map<size_t,size_t>& selectedDataDescIds, const std::string& dataColumnName, bool requireModel);
	
	virtual bool AtEnd() const {
		return _flushPosition >= _nElements;
	}
	
	virtual void NextRow();
	
	virtual void ReadData(DataArray& data, FlagArray& flags, WeightArray& weights, double& u, double& v, double& w, uint32_t& dataDescId);
	
	virtual void ReadModel(DataArray& model);
	
	virtual void OutputStatistics() const;
	
private:
	class AveragingBuffer
	{
	public:
		AveragingBuffer() :
			_data(0),
			_modelData(0),
			_weights(0),
			_averagedDataCount(0),
			_uvwWeight(0)
		{ }
		~AveragingBuffer()
		{
			delete[] _data;
			delete[] _modelData;
			delete[] _weights;
		}
		bool IsInitialized() const { return _data != 0; }
		void Initialize(size_t bufferSize, bool includeModel)
		{
			_data = new std::complex<float>[bufferSize];
			_modelData = includeModel ? new std::complex<float>[bufferSize] : 0;
			_weights = new float[bufferSize];
			Reset(bufferSize);
		}
		void AddData(size_t n, const std::complex<float>* data, const bool* flags, const float* weights, const double* uvw)
		{
			double weightSum = 0.0;
			for(size_t i=0; i!=n; ++i)
			{
				if(!flags[i] && std::isfinite(data[i].real()) && std::isfinite(data[i].imag()))
				{
					_data[i] += data[i] * weights[i];
					_weights[i] += weights[i];
					weightSum += weights[i];
				}
			}
			weightSum /= double(n);
			_uvw[0] += uvw[0] * weightSum;
			_uvw[1] += uvw[1] * weightSum;
			_uvw[2] += uvw[2] * weightSum;
			_uvwWeight += weightSum;
			++_averagedDataCount;
		}
		void AddDataAndModel(size_t n, const std::complex<float>* data, const std::complex<float>* modelData, const bool* flags, const float* weights, const double* uvw)
		{
			double weightSum = 0.0;
			for(size_t i=0; i!=n; ++i)
			{
				if(!flags[i] && std::isfinite(data[i].real()) && std::isfinite(data[i].imag()) &&
					std::isfinite(modelData[i].real()) && std::isfinite(modelData[i].imag()))
				{
					_data[i] += data[i] * weights[i];
					_modelData[i] += modelData[i] * weights[i];
					_weights[i] += weights[i];
					weightSum += weights[i];
				}
			}
			weightSum /= double(n);
			_uvw[0] += uvw[0] * weightSum;
			_uvw[1] += uvw[1] * weightSum;
			_uvw[2] += uvw[2] * weightSum;
			_uvwWeight += weightSum;
			++_averagedDataCount;
		}
		size_t AveragedDataCount() const { return _averagedDataCount; }
		
		void Get(size_t n, std::complex<float>* data, bool* flags, float* weights, double* uvw)
		{
			for(size_t i=0; i!=n; ++i)
			{
				data[i] = _data[i] / _weights[i];
				flags[i] = (weights[i]==0.0);
				weights[i] = _weights[i];
			}
			for(size_t i=0; i!=3; ++i)
				uvw[i] = _uvw[i] / _uvwWeight;
		}
		
		void Get(size_t n, std::complex<float>* data, std::complex<float>* modelData, bool* flags, float* weights, double* uvw)
		{
			for(size_t i=0; i!=n; ++i)
			{
				data[i] = _data[i] / _weights[i];
				modelData[i] = _modelData[i] / weights[i];
				flags[i] = (weights[i]==0.0);
				weights[i] = _weights[i];
			}
			for(size_t i=0; i!=3; ++i)
				uvw[i] = _uvw[i] / _uvwWeight;
		}
		
		void Reset(size_t n)
		{
			for(size_t i=0; i!=n; ++i)
			{
				_data[i] = 0.0;
				_weights[i] = 0.0;
				if(_modelData != 0) _modelData[i] = 0.0;
			}
			_uvw[0] = 0.0; _uvw[1] = 0.0; _uvw[2] = 0.0;
			_averagedDataCount = 0;
			_uvwWeight = 0.0;
		}
	private:
		std::complex<float>* _data;
		std::complex<float>* _modelData;
		double _uvw[3];
		float* _weights;
		size_t _averagedDataCount;
		double _uvwWeight;
	};
	
	bool processCurrentTimestep();
	
	typedef std::tuple<double, double, double> Pos;
	
	ao::uvector<size_t> _averagingFactors;
	std::vector<AveragingBuffer> _buffers;
	ao::uvector<size_t> _spwIndexToDataDescId;
	
	size_t _nAntennae;
	
	DataArray _currentData, _currentModel;
	FlagArray _currentFlags;
	WeightArray _currentWeights;
	size_t _averagedDataDescId;
	size_t _nElements;
	size_t _flushPosition;
	
	// Some statistics
	size_t _averagedRowCount;
	size_t _averageFactorSum;
	size_t _rowCount;
};

#endif
