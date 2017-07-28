#ifndef PARTITIONED_MS
#define PARTITIONED_MS

#include <fstream>
#include <string>
#include <map>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>

#include "../polarization.h"
#include "../uvector.h"
#include "../msselection.h"

#include "msprovider.h"

class PartitionedMS : public MSProvider
{
public:
	class Handle;
	
	struct ChannelRange
	{
		int dataDescId;
		size_t start, end;
		bool operator<(const ChannelRange& rhs) const
		{
			if(dataDescId < rhs.dataDescId) return true;
			if(dataDescId > rhs.dataDescId) return false;
			if(start < rhs.start) return true;
			if(start > rhs.start) return false;
			return end < rhs.end;
		}
	};
	
	PartitionedMS(const Handle& handle, size_t partIndex, PolarizationEnum polarization, size_t bandIndex);
	
	virtual ~PartitionedMS();
	
	PartitionedMS(const PartitionedMS&) = delete;
	PartitionedMS& operator=(const PartitionedMS&) = delete;
	
	virtual casacore::MeasurementSet &MS() final override { openMS(); return *_ms; }
	
	virtual size_t RowId() const final override { return _currentRow; }
	
	virtual bool CurrentRowAvailable() final override;
	
	virtual void NextRow() final override;
	
	virtual void Reset() final override;
	
	virtual void ReadMeta(double& u, double& v, double& w, size_t& dataDescId) final override;
	
	virtual void ReadMeta(MetaData& metaData) final override;
	
	virtual void ReadData(std::complex<float>* buffer) final override;
	
	virtual void ReadModel(std::complex<float>* buffer) final override;
	
	virtual void WriteModel(size_t rowId, std::complex<float>* buffer) final override;
	
	virtual void ReadWeights(float* buffer) final override;
	
	virtual void ReadWeights(std::complex<float>* buffer) final override;
	
	virtual void ReopenRW() final override{ }
	
	virtual double StartTime() final override { return _metaHeader.startTime; }
	
	virtual void MakeIdToMSRowMapping(std::vector<size_t>& idToMSRow) final override;
	
	virtual PolarizationEnum Polarization() final override { return _polarization; }
	
	static Handle Partition(const string& msPath, const std::vector<ChannelRange>& channels, class MSSelection& selection, const string& dataColumnName, bool includeModel, bool initialModelRequired, const class WSCleanSettings& settings);
	
	class Handle {
	public:
		friend class PartitionedMS;
		
		Handle(const Handle& handle) : _data(handle._data)
		{
			++(_data->_referenceCount);
		}
		~Handle() { decrease(); }
		Handle operator=(const Handle& handle)
		{
			if(handle._data != _data)
			{
				decrease();
				_data = handle._data;
				++(_data->_referenceCount);
			}
			return *this;
		}
	private:
		struct HandleData
		{
			HandleData(const std::string& msPath, const string& dataColumnName, const std::string& temporaryDirectory, const std::vector<ChannelRange>& channels, bool initialModelRequired, bool modelUpdateRequired, const std::set<PolarizationEnum>& polarizations, const MSSelection& selection) :
			_msPath(msPath), _dataColumnName(dataColumnName), _temporaryDirectory(temporaryDirectory), _channels(channels), _initialModelRequired(initialModelRequired), _modelUpdateRequired(modelUpdateRequired),
			_polarizations(polarizations), _selection(selection), _referenceCount(1) { }
			
			std::string _msPath, _dataColumnName, _temporaryDirectory;
			std::vector<ChannelRange> _channels;
			bool _initialModelRequired, _modelUpdateRequired;
			std::set<PolarizationEnum> _polarizations;
			MSSelection _selection;
			size_t _referenceCount;
		} *_data;
		
		void decrease();
		Handle(const std::string& msPath, const string& dataColumnName, const std::string& temporaryDirectory, const std::vector<ChannelRange>& channels, bool initialModelRequired, bool modelUpdateRequired, const std::set<PolarizationEnum>& polarizations, const MSSelection& selection) :
			_data(new HandleData(msPath, dataColumnName, temporaryDirectory, channels, initialModelRequired, modelUpdateRequired, polarizations, selection))
		{
		}
	};
private:
	static void unpartition(const Handle& handle);
	
	static void getDataDescIdMap(std::map<size_t,size_t>& dataDescIds, const vector<PartitionedMS::ChannelRange>& channels);
	
	void openMS();
	
	Handle _handle;
	std::string _msPath;
	std::unique_ptr<casacore::MeasurementSet> _ms;
	std::ifstream _metaFile, _weightFile, _dataFile;
	char *_modelFileMap;
	size_t _currentRow;
	bool _readPtrIsOk, _metaPtrIsOk, _weightPtrIsOk;
	ao::uvector<float> _weightBuffer;
	ao::uvector<std::complex<float>> _modelBuffer;
	int _fd;
	PolarizationEnum _polarization;
	size_t _polarizationCountInFile;
	
	struct MetaHeader
	{
		uint64_t selectedRowCount;
		uint32_t filenameLength;
		double startTime;
	} _metaHeader;
	struct MetaRecord
	{
		double u, v, w, time;
		uint16_t antenna1, antenna2, dataDescId;
		static constexpr size_t BINARY_SIZE = 8*4 + 2*3;
		void read(std::istream& str)
		{
			str.read(reinterpret_cast<char*>(&u), sizeof(double));
			str.read(reinterpret_cast<char*>(&v), sizeof(double));
			str.read(reinterpret_cast<char*>(&w), sizeof(double));
			str.read(reinterpret_cast<char*>(&time), sizeof(double));
			str.read(reinterpret_cast<char*>(&antenna1), sizeof(uint16_t));
			str.read(reinterpret_cast<char*>(&antenna2), sizeof(uint16_t));
			str.read(reinterpret_cast<char*>(&dataDescId), sizeof(uint16_t));
		}
		void write(std::ostream& str) const
		{
			str.write(reinterpret_cast<const char*>(&u), sizeof(double));
			str.write(reinterpret_cast<const char*>(&v), sizeof(double));
			str.write(reinterpret_cast<const char*>(&w), sizeof(double));
			str.write(reinterpret_cast<const char*>(&time), sizeof(double));
			str.write(reinterpret_cast<const char*>(&antenna1), sizeof(uint16_t));
			str.write(reinterpret_cast<const char*>(&antenna2), sizeof(uint16_t));
			str.write(reinterpret_cast<const char*>(&dataDescId), sizeof(uint16_t));
		}
	};
	struct PartHeader
	{
		uint64_t channelCount;
		uint64_t channelStart;
		uint32_t dataDescId;
		bool hasModel, hasWeights;
	} _partHeader;
	
	static std::string getPartPrefix(const std::string& msPath, size_t partIndex, PolarizationEnum pol, size_t dataDescId, const std::string& tempDir);
	static std::string getMetaFilename(const std::string& msPath, const std::string& tempDir, size_t dataDescId);
};

#endif
