#include "partitionedms.h"

#include "../multibanddata.h"
#include "../progressbar.h"

#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <errno.h>
#include <fcntl.h>
#include <string.h>

#include <cstdio>
#include <fstream>
#include <sstream>
#include <map>
#include <memory>
#include <vector>

#include <boost/filesystem/path.hpp>

#include <casacore/measures/Measures/MEpoch.h>

#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>

// #define REDUNDANT_VALIDATION 1

PartitionedMS::PartitionedMS(const Handle& handle, size_t partIndex, PolarizationEnum polarization, size_t dataDescId) :
	_metaFile(handle._data->_metaFile),
	_modelFileMap(0),
	_currentRow(0),
	_readPtrIsOk(true),
	_metaPtrIsOk(true),
	_weightPtrIsOk(true)
{
	_metaFile.read(reinterpret_cast<char*>(&_metaHeader), sizeof(MetaHeader));
	std::vector<char> msPath(_metaHeader.filenameLength+1, char(0));
	_metaFile.read(msPath.data(), _metaHeader.filenameLength);
	std::cout << "Opening reordered part " << partIndex << " spw " << dataDescId << " for " << msPath.data() << '\n';
	_ms = casacore::MeasurementSet(msPath.data());
	
	std::string partPrefix = getPartPrefix(msPath.data(), partIndex, polarization, dataDescId, handle._data->_temporaryDirectory);
	
	_dataFile.open(partPrefix+".tmp", std::ios::in);
	if(!_dataFile.good())
		throw std::runtime_error("Error opening temporary data file");
	_dataFile.read(reinterpret_cast<char*>(&_partHeader), sizeof(PartHeader));
	if(!_dataFile.good())
		throw std::runtime_error("Error reading header from file");
	
	if(_partHeader.hasModel)
	{
		_fd = open((partPrefix+"-m.tmp").c_str(), O_RDWR);
		if(_fd == -1)
			throw std::runtime_error("Error opening temporary data file");
		size_t length = _partHeader.channelCount * _metaHeader.selectedRowCount * sizeof(std::complex<float>);
		if(length == 0)
			_modelFileMap = 0;
		else {
			_modelFileMap = reinterpret_cast<char*>( mmap(NULL, length, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_NORESERVE, _fd, 0) );
			if(_modelFileMap == MAP_FAILED)
			{
				int errsv = errno;
				char buffer[1024];
				const char* msg = strerror_r(errsv, buffer, 1024); 
				
				_modelFileMap = 0;
				throw std::runtime_error(std::string("Error creating memory map to temporary model file: mmap() returned MAP_FAILED with error message: ") + msg);
			}
		}
	}
	
	_weightFile.open(partPrefix+"-w.tmp", std::ios::in);
	if(!_weightFile.good())
		throw std::runtime_error("Error opening temporary data file");
	_weightBuffer.resize(_partHeader.channelCount);
	_modelBuffer.resize(_partHeader.channelCount);
}

PartitionedMS::~PartitionedMS()
{
	if(_modelFileMap!=0)
	{
		size_t length = _partHeader.channelCount * _metaHeader.selectedRowCount * sizeof(std::complex<float>);
		if(length != 0)
			munmap(_modelFileMap, length);
	}
	if(_partHeader.hasModel)
		close(_fd);
}

void PartitionedMS::Reset()
{
	_currentRow = 0;
	_metaFile.seekg(sizeof(MetaHeader) + _metaHeader.filenameLength, std::ios::beg);
	_dataFile.seekg(sizeof(PartHeader), std::ios::beg);
	_weightFile.seekg(0, std::ios::beg);
	_readPtrIsOk = true;
	_metaPtrIsOk = true;
	_weightPtrIsOk = true;
}

bool PartitionedMS::CurrentRowAvailable()
{
	return _currentRow < _metaHeader.selectedRowCount;
}

void PartitionedMS::NextRow()
{
	++_currentRow;
	if(_currentRow < _metaHeader.selectedRowCount)
	{
		if(_readPtrIsOk)
			_dataFile.seekg(_partHeader.channelCount * sizeof(std::complex<float>), std::ios::cur);
		else
			_readPtrIsOk = true;
		
		if(_metaPtrIsOk)
			_metaFile.seekg(sizeof(MetaRecord), std::ios::cur);
		else
			_metaPtrIsOk = true;
		
		if(_weightPtrIsOk && _partHeader.hasWeights)
			_weightFile.seekg(_partHeader.channelCount * sizeof(float), std::ios::cur);
		_weightPtrIsOk = true;
	}
}

void PartitionedMS::ReadMeta(double& u, double& v, double& w, size_t& dataDescId)
{
	if(!_metaPtrIsOk)
		_metaFile.seekg(-sizeof(MetaRecord), std::ios::cur);
	_metaPtrIsOk = false;
	
	MetaRecord record;
	_metaFile.read(reinterpret_cast<char*>(&record), sizeof(MetaRecord));
	u = record.u;
	v = record.v;
	w = record.w;
	dataDescId = record.dataDescId;
}

void PartitionedMS::ReadData(std::complex<float>* buffer)
{
	if(!_readPtrIsOk)
	{
		_dataFile.seekg(-_partHeader.channelCount * sizeof(std::complex<float>), std::ios::cur);
	}
#ifdef REDUNDANT_VALIDATION
	size_t pos = size_t(_dataFile.tellg()) - sizeof(PartHeader);
	size_t fact = _partHeader.hasModel ? 2 : 1;
	if(pos != fact * _currentRow * _partHeader.channelCount * sizeof(std::complex<float>))
	{
		std::ostringstream s;
		s << "Not on right pos: " << pos << " instead of " << fact * _currentRow * _partHeader.channelCount * sizeof(std::complex<float>) <<
			" (row " << (pos / (fact * _partHeader.channelCount * sizeof(std::complex<float>))) << " instead of " << _currentRow << ")";
		throw std::runtime_error(s.str());
	}
#endif
	_dataFile.read(reinterpret_cast<char*>(buffer), _partHeader.channelCount * sizeof(std::complex<float>));
	_readPtrIsOk = false;
}

void PartitionedMS::ReadModel(std::complex<float>* buffer)
{
#ifdef REDUNDANT_VALIDATION
	if(!_partHeader.hasModel)
		throw std::runtime_error("Partitioned MS initialized without model");
#endif
	size_t rowLength = _partHeader.channelCount * sizeof(std::complex<float>);
	memcpy(reinterpret_cast<char*>(buffer), _modelFileMap + rowLength*_currentRow, rowLength);
}

void PartitionedMS::WriteModel(size_t rowId, std::complex<float>* buffer)
{
#ifdef REDUNDANT_VALIDATION
	if(!_partHeader.hasModel)
		throw std::runtime_error("Partitioned MS initialized without model");
#endif
	_weightFile.seekg(_partHeader.channelCount * sizeof(float) * rowId, std::ios::beg);
	_weightFile.read(reinterpret_cast<char*>(_weightBuffer.data()), _partHeader.channelCount * sizeof(float));
	for(size_t i=0; i!=_partHeader.channelCount; ++i)
		buffer[i] *= _weightBuffer[i];
	
	size_t rowLength = _partHeader.channelCount * sizeof(std::complex<float>);
	std::complex<float>* modelWritePtr = reinterpret_cast<std::complex<float>*>(_modelFileMap + rowLength*rowId);
	
	// In case the value was not sampled in this pass, it will be set to infinite and should not overwrite the current
	// value in the set.
	for(size_t i=0; i!=_partHeader.channelCount; ++i)
	{
		if(std::isfinite(buffer[i].real()))
			modelWritePtr[i] = buffer[i];
	}
}

void PartitionedMS::ReadWeights(std::complex<float>* buffer)
{
	if(!_weightPtrIsOk)
		_weightFile.seekg(-_partHeader.channelCount * sizeof(float), std::ios::cur);
	float* displacedBuffer = reinterpret_cast<float*>(buffer)+_partHeader.channelCount;
	_weightFile.read(reinterpret_cast<char*>(displacedBuffer), _partHeader.channelCount * sizeof(float));
	_weightPtrIsOk = false;
	copyRealToComplex(buffer, displacedBuffer, _partHeader.channelCount);
}

void PartitionedMS::ReadWeights(float* buffer)
{
	if(!_weightPtrIsOk)
		_weightFile.seekg(-_partHeader.channelCount * sizeof(float), std::ios::cur);
	_weightFile.read(reinterpret_cast<char*>(buffer), _partHeader.channelCount * sizeof(float));
	_weightPtrIsOk = false;
}

std::string PartitionedMS::getPartPrefix(const std::string& msPathStr, size_t partIndex, PolarizationEnum pol, size_t dataDescId, const std::string& tempDir)
{
	boost::filesystem::path
		msPath(msPathStr),
		prefixPath;
	if(tempDir.empty())
		prefixPath = msPath;
	else
		prefixPath = boost::filesystem::path(tempDir) / msPath.filename();
		
	std::string prefix(prefixPath.string());
	while(!prefix.empty() && *prefix.rbegin() == '/')
		prefix.resize(prefix.size()-1);
	
	std::ostringstream partPrefix;
	partPrefix << prefix << "-part";
	if(partIndex < 1000) partPrefix << '0';
	if(partIndex < 100) partPrefix << '0';
	if(partIndex < 10) partPrefix << '0';
	partPrefix << partIndex;
	partPrefix << "-";
	partPrefix << Polarization::TypeToShortString(pol);
	partPrefix << "-b" << dataDescId;
	return partPrefix.str();
}

string PartitionedMS::getMetaFilename(const string& msPathStr, const std::string& tempDir)
{
	boost::filesystem::path
		msPath(msPathStr),
		prefixPath;
	if(tempDir.empty())
		prefixPath = msPath;
	else
		prefixPath = boost::filesystem::path(tempDir) / msPath.filename();
	std::string prefix(prefixPath.string());
	while(!prefix.empty() && *prefix.rbegin() == '/')
		prefix.resize(prefix.size()-1);
	return prefix + "-parted-meta.tmp";
}

// should be private but is not allowed on older compilers
struct PartitionFiles
{
	std::ofstream
		*data,
		*weight,
		*model;
};


/*
 * When partitioned:
 * One global file stores:
 * - Metadata:
 *   * Number of selected rows
 *   * Filename length + string
 *   * [ UVW, dataDescId ]
 * The binary parts store the following information:
 * - Number of channels
 * - Start channel in MS
 * - Total weight in part
 * - Data    (single polarization, as requested)
 * - Weights (single, only needed when imaging PSF)
 * - Model, optionally
 */
PartitionedMS::Handle PartitionedMS::Partition(const string& msPath, const std::vector<ChannelRange>& channels, MSSelection& selection, const string& dataColumnName, bool includeWeights, bool includeModel, bool initialModelRequired, bool modelUpdateRequired, const std::set<PolarizationEnum>& polsOut, const std::string& temporaryDirectory)
{
	size_t channelParts = channels.size();
	casacore::MeasurementSet ms(msPath);

	// Ordered as files[pol x channelpart]
	std::vector<PartitionFiles>
		files(channelParts*polsOut.size());
	
	size_t fileIndex = 0;
	for(size_t part=0; part!=channelParts; ++part)
	{
		for(std::set<PolarizationEnum>::const_iterator p=polsOut.begin(); p!=polsOut.end(); ++p)
		{
			PartitionFiles& f = files[fileIndex];
			std::string partPrefix = getPartPrefix(msPath, part, *p, channels[part].dataDescId, temporaryDirectory);
			f.data = new std::ofstream(partPrefix + ".tmp");
			if(includeWeights)
				f.weight = new std::ofstream(partPrefix + "-w.tmp");
			if(initialModelRequired)
				f.model = new std::ofstream(partPrefix + "-m.tmp");
			f.data->seekp(sizeof(PartHeader), std::ios::beg);
			
			++fileIndex;
		}
	}
	std::vector<PolarizationEnum> msPolarizations = GetMSPolarizations(ms);
	
	MultiBandData band(ms.spectralWindow(), ms.dataDescription());
	casacore::ROScalarColumn<int> antenna1Column(ms, casacore::MS::columnName(casacore::MSMainEnums::ANTENNA1));
	casacore::ROScalarColumn<int> antenna2Column(ms, casacore::MS::columnName(casacore::MSMainEnums::ANTENNA2));
	casacore::ROScalarColumn<int> fieldIdColumn(ms, casacore::MS::columnName(casacore::MSMainEnums::FIELD_ID));
	casacore::ROScalarColumn<double> timeColumn(ms, casacore::MS::columnName(casacore::MSMainEnums::TIME));
	casacore::MEpoch::ROScalarColumn timeEpochColumn(ms, casacore::MS::columnName(casacore::MSMainEnums::TIME));
	casacore::ROArrayColumn<double> uvwColumn(ms, casacore::MS::columnName(casacore::MSMainEnums::UVW));
	std::unique_ptr<casacore::ROArrayColumn<float>> weightColumn;
	casacore::ROArrayColumn<casacore::Complex> dataColumn(ms, dataColumnName);
	casacore::ROArrayColumn<bool> flagColumn(ms, casacore::MS::columnName(casacore::MSMainEnums::FLAG));
	casacore::ROScalarColumn<int> dataDescIdColumn(ms, ms.columnName(casacore::MSMainEnums::DATA_DESC_ID));

	std::unique_ptr<casacore::ROArrayColumn<casacore::Complex>> modelColumn;
	if(initialModelRequired)
	{
		modelColumn.reset(new casacore::ROArrayColumn<casacore::Complex>(ms, ms.columnName(casacore::MSMainEnums::MODEL_DATA)));
	}
	
	const casacore::IPosition shape(dataColumn.shape(0));
	const size_t polarizationCount = shape[0];
	size_t channelCount = shape[1];
	
	bool isWeightDefined;
	if(ms.isColumn(casacore::MSMainEnums::WEIGHT_SPECTRUM))
	{
		weightColumn.reset(new casacore::ROArrayColumn<float>(ms, casacore::MS::columnName(casacore::MSMainEnums::WEIGHT_SPECTRUM)));
		isWeightDefined = weightColumn->isDefined(0);
	} else {
		isWeightDefined = false;
	}
	bool msHasWeights = false;
	casacore::Array<float> weightArray(shape);
	if(isWeightDefined)
	{
		casacore::IPosition modelShape = weightColumn->shape(0);
		msHasWeights = (modelShape == shape);
	}
	if(!msHasWeights)
	{
		weightArray.set(1);
		std::cout << "WARNING: This measurement set has no or an invalid WEIGHT_SPECTRUM column; all visibilities are assumed to have equal weight.\n";
	}
	
	size_t startRow, endRow;
	getRowRange(ms, selection, startRow, endRow);
	
	// Count selected rows
	uint64_t selectedRowCount = 0;
	size_t timestep = selection.HasInterval() ? selection.IntervalStart() : 0;
	double time = timeColumn(startRow);
	for(size_t row=startRow; row!=endRow; ++row)
	{
		const int
			a1 = antenna1Column(row), a2 = antenna2Column(row),
			fieldId = fieldIdColumn(row);
		casacore::Vector<double> uvw = uvwColumn(row);
		if(time != timeColumn(row))
		{
			++timestep;
			time = timeColumn(row);
		}
		if(selection.IsSelected(fieldId, timestep, a1, a2, uvw))
			++selectedRowCount;
	}
	std::cout << "Reordering " << msPath << " (" << selectedRowCount << " selected rows) into " << channelParts << " x " << polsOut.size() << " parts.\n";

	// Write header of meta file
	std::string metaFilename = getMetaFilename(msPath, temporaryDirectory);
	std::ofstream metaFile(metaFilename);
	MetaHeader metaHeader;
	memset(&metaHeader, 0, sizeof(MetaHeader));
	metaHeader.selectedRowCount = selectedRowCount;
	metaHeader.filenameLength = msPath.size();
	metaHeader.startTime = timeEpochColumn(startRow).getValue().get();
	metaFile.write(reinterpret_cast<char*>(&metaHeader), sizeof(metaHeader));
	metaFile.write(msPath.c_str(), msPath.size());
	
	// Write actual data
	timestep = selection.HasInterval() ? selection.IntervalStart() : 0;
	time = timeColumn(startRow);
	
	std::vector<std::complex<float>> dataBuffer(polarizationCount * channelCount);
	std::vector<float> weightBuffer(polarizationCount * channelCount);
	
	casacore::Array<std::complex<float>> dataArray(shape), modelArray(shape);
	casacore::Array<bool> flagArray(shape);
	ProgressBar progress1("Reordering");
	for(size_t row=startRow; row!=endRow; ++row)
	{
		progress1.SetProgress(row-startRow, endRow-startRow);
		const int
			a1 = antenna1Column(row), a2 = antenna2Column(row),
			fieldId = fieldIdColumn(row);
			
		if(time != timeColumn(row))
		{
			++timestep;
			time = timeColumn(row);
		}
		casacore::Vector<double> uvwArray = uvwColumn(row);
		if(selection.IsSelected(fieldId, timestep, a1, a2, uvwArray))
		{
			size_t dataDescId = dataDescIdColumn(row);
			MetaRecord meta;
			memset(&meta, 0, sizeof(MetaRecord));
			meta.u = uvwArray(0);
			meta.v = uvwArray(1);
			meta.w = uvwArray(2);
			meta.dataDescId = dataDescId;
			metaFile.write(reinterpret_cast<char*>(&meta), sizeof(MetaRecord));
			if(!metaFile.good())
				throw std::runtime_error("Error writing to temporary file");
				
			dataColumn.get(row, dataArray);
			if(initialModelRequired)
				modelColumn->get(row, modelArray);
			if(msHasWeights)
				weightColumn->get(row, weightArray);
			flagColumn.get(row, flagArray);
			
			fileIndex = 0;
			for(size_t part=0; part!=channelParts; ++part)
			{
				size_t
					dataDescId = channels[part].dataDescId, // TODO use it
					partStartCh = channels[part].start,
					partEndCh = channels[part].end;
				
				for(std::set<PolarizationEnum>::const_iterator p=polsOut.begin(); p!=polsOut.end(); ++p)
				{
					PartitionFiles& f = files[fileIndex];
					copyWeightedData(dataBuffer.data(), partStartCh, partEndCh, msPolarizations, dataArray, weightArray, flagArray, *p);
					f.data->write(reinterpret_cast<char*>(dataBuffer.data()), (partEndCh - partStartCh) * sizeof(std::complex<float>));
					if(!f.data->good())
						throw std::runtime_error("Error writing to temporary data file");
					
					if(initialModelRequired)
					{
						copyWeightedData(dataBuffer.data(), partStartCh, partEndCh, msPolarizations, modelArray, weightArray, flagArray, *p);
						f.model->write(reinterpret_cast<char*>(dataBuffer.data()), (partEndCh - partStartCh) * sizeof(std::complex<float>));
						if(!f.model->good())
							throw std::runtime_error("Error writing to temporary data file");
					}
					
					if(includeWeights)
					{
						copyWeights(weightBuffer.data(), partStartCh, partEndCh, msPolarizations, dataArray, weightArray, flagArray, *p);
						f.weight->write(reinterpret_cast<char*>(weightBuffer.data()), (partEndCh - partStartCh) * sizeof(float));
						if(!f.weight->good())
							throw std::runtime_error("Error writing to temporary weights file");
					}
					++fileIndex;
				}
			}
		}
	}
	progress1.SetProgress(ms.nrow(), ms.nrow());
	
	// Write header to parts and write empty model files (if requested)
	PartHeader header;
	memset(&header, 0, sizeof(PartHeader));
	header.hasModel = includeModel;
	header.hasWeights = includeWeights;
	fileIndex = 0;
	dataBuffer.assign(channelCount, 0.0);
	std::unique_ptr<ProgressBar> progress2;
	if(includeModel && !initialModelRequired)
		progress2.reset(new ProgressBar("Initializing model visibilities"));
	for(size_t part=0; part!=channelParts; ++part)
	{
		header.channelStart = channels[part].start,
		header.channelCount = channels[part].end - header.channelStart;
		header.dataDescId = channels[part].dataDescId;
		for(std::set<PolarizationEnum>::const_iterator p=polsOut.begin(); p!=polsOut.end(); ++p)
		{
			PartitionFiles& f = files[fileIndex];
			f.data->seekp(0, std::ios::beg);
			f.data->write(reinterpret_cast<char*>(&header), sizeof(PartHeader));
			if(!f.data->good())
				throw std::runtime_error("Error writing to temporary data file");
			
			delete f.data;
			if(includeWeights)
				delete f.weight;
			if(initialModelRequired)
				delete f.model;
			++fileIndex;
			
			// If model is requested, fill model file with zeros
			if(includeModel && !initialModelRequired)
			{
				std::string partPrefix = getPartPrefix(msPath, part, *p, header.dataDescId, temporaryDirectory);
				std::ofstream modelFile(partPrefix + "-m.tmp");
				for(size_t i=0; i!=selectedRowCount; ++i)
				{
					modelFile.write(reinterpret_cast<char*>(dataBuffer.data()), header.channelCount * sizeof(std::complex<float>));
					progress2->SetProgress(part*selectedRowCount + i, channelParts*selectedRowCount);
				}
			}
		}
	}
	progress2.reset();
	
	return Handle(metaFilename, msPath, dataColumnName, temporaryDirectory, channels, modelUpdateRequired, polsOut, selection);
}

void PartitionedMS::unpartition(const PartitionedMS::Handle& handle)
{
	const std::set<PolarizationEnum> pols = handle._data->_polarizations;
	std::ifstream metaFile(handle._data->_metaFile);
	MetaHeader metaHeader;
	metaFile.read(reinterpret_cast<char*>(&metaHeader), sizeof(MetaHeader));
	std::vector<char> msPath(metaHeader.filenameLength+1, char(0));
	metaFile.read(msPath.data(), metaHeader.filenameLength);
	
	ChannelRange firstRange = handle._data->_channels[0];
	std::ifstream firstDataFile(getPartPrefix(msPath.data(), 0, *pols.begin(), firstRange.dataDescId, handle._data->_temporaryDirectory)+".tmp", std::ios::in);
	if(!firstDataFile.good())
		throw std::runtime_error("Error opening temporary data file");
	PartHeader firstPartHeader;
	firstDataFile.read(reinterpret_cast<char*>(&firstPartHeader), sizeof(PartHeader));
	if(!firstDataFile.good())
		throw std::runtime_error("Error reading from temporary data file");
	
	if(firstPartHeader.hasModel)
	{
		const size_t channelParts = handle._data->_channels.size();
		std::vector<std::ifstream*> modelFiles(channelParts*pols.size()), weightFiles(channelParts*pols.size());
		size_t fileIndex = 0;
		for(size_t part=0; part!=channelParts; ++part)
		{
			size_t dataDescId = handle._data->_channels[part].dataDescId;
			for(std::set<PolarizationEnum>::const_iterator p=pols.begin(); p!=pols.end(); ++p)
			{
				std::string partPrefix = getPartPrefix(msPath.data(), part, *p, dataDescId, handle._data->_temporaryDirectory);
				modelFiles[fileIndex] = new std::ifstream(partPrefix + "-m.tmp");
				if(firstPartHeader.hasWeights)
					weightFiles[fileIndex] = new std::ifstream(partPrefix + "-w.tmp");
				++fileIndex;
			}
		}
		
		casacore::MeasurementSet ms(msPath.data(), casacore::Table::Update);
		const std::vector<PolarizationEnum> msPolarizations = GetMSPolarizations(ms);
		initializeModelColumn(ms);
		casacore::ROScalarColumn<int> antenna1Column(ms, casacore::MS::columnName(casacore::MSMainEnums::ANTENNA1));
		casacore::ROScalarColumn<int> antenna2Column(ms, casacore::MS::columnName(casacore::MSMainEnums::ANTENNA2));
		casacore::ROScalarColumn<int> fieldIdColumn(ms, casacore::MS::columnName(casacore::MSMainEnums::FIELD_ID));
		casacore::ROScalarColumn<double> timeColumn(ms, casacore::MS::columnName(casacore::MSMainEnums::TIME));
		casacore::ROArrayColumn<casacore::Complex> dataColumn(ms, handle._data->_dataColumnName);
		casacore::ArrayColumn<casacore::Complex> modelColumn(ms, casacore::MS::columnName(casacore::MSMainEnums::MODEL_DATA));
		casacore::ROArrayColumn<double> uvwColumn(ms, casacore::MS::columnName(casacore::MSMainEnums::UVW));
		
		const casacore::IPosition shape(dataColumn.shape(0));
		size_t channelCount = shape[1];
		
		std::vector<std::complex<float>> modelDataBuffer(channelCount);
		std::vector<float> weightBuffer(channelCount);
		casacore::Array<std::complex<float>> modelDataArray(shape);
	
		ProgressBar progress(std::string("Writing changed model back to ") + msPath.data());
		size_t startRow, endRow;
		getRowRange(ms, handle._data->_selection, startRow, endRow);
		size_t timestep = handle._data->_selection.HasInterval() ? handle._data->_selection.IntervalStart() : 0;
		double time = timeColumn(startRow);
		for(size_t row=startRow; row!=endRow; ++row)
		{
			progress.SetProgress(row - startRow, startRow - endRow);
			const int
				a1 = antenna1Column(row), a2 = antenna2Column(row),
				fieldId = fieldIdColumn(row);
			casacore::Vector<double> uvw = uvwColumn(row);
				
			if(time != timeColumn(row))
			{
				++timestep;
				time = timeColumn(row);
			}
			if(handle._data->_selection.IsSelected(fieldId, timestep, a1, a2, uvw))
			{
				modelColumn.get(row, modelDataArray);
				size_t fileIndex = 0;
				for(size_t part=0; part!=channelParts; ++part)
				{
					size_t
						dataDescId = handle._data->_channels[part].dataDescId, //TODO take this into account
						partStartCh = handle._data->_channels[part].start,
						partEndCh = handle._data->_channels[part].end;
					
					for(std::set<PolarizationEnum>::const_iterator p=pols.begin(); p!=pols.end(); ++p)
					{
						modelFiles[fileIndex]->read(reinterpret_cast<char*>(modelDataBuffer.data()), (partEndCh - partStartCh) * sizeof(std::complex<float>));
						if(!modelFiles[fileIndex]->good())
							throw std::runtime_error("Error reading from temporary model data file");
						if(firstPartHeader.hasWeights)
						{
							weightFiles[fileIndex]->read(reinterpret_cast<char*>(weightBuffer.data()), (partEndCh - partStartCh) * sizeof(float));
							if(!weightFiles[fileIndex]->good())
								throw std::runtime_error("Error reading from temporary weight data file");
							for(size_t i=0; i!=partEndCh - partStartCh; ++i)
								modelDataBuffer[i] /= weightBuffer[i];
						}
						reverseCopyData(modelDataArray, partStartCh, partEndCh, msPolarizations, modelDataBuffer.data(), *p);
						
						++fileIndex;
					}
				}
				modelColumn.put(row, modelDataArray);
			}
		}
		progress.SetProgress(ms.nrow(),ms.nrow());
		
		fileIndex = 0;
		for(size_t part=0; part!=channelParts; ++part)
		{
			for(std::set<PolarizationEnum>::const_iterator p=pols.begin(); p!=pols.end(); ++p)
			{
				delete modelFiles[fileIndex];
				if(firstPartHeader.hasWeights)
					delete weightFiles[fileIndex];
				++fileIndex;
			}
		}
	}
}

void PartitionedMS::Handle::decrease()
{
	--(_data->_referenceCount);
	if(_data->_referenceCount == 0)
	{
		if(_data->_modelUpdateRequired)
			PartitionedMS::unpartition(*this);
		
		std::cout << "Cleaning up temporary files...\n";
		
		//MetaHeader metaHeader;
		//std::ifstream metaFile(_metaFile);
		//metaFile.read(reinterpret_cast<char*>(&metaHeader), sizeof(MetaHeader));
		
		for(size_t part=0; part!=_data->_channels.size(); ++part)
		{
			for(std::set<PolarizationEnum>::const_iterator p=_data->_polarizations.begin(); p!=_data->_polarizations.end(); ++p)
			{
				std::string prefix = getPartPrefix(_data->_msPath, part, *p, _data->_channels[part].dataDescId, _data->_temporaryDirectory);
				std::remove((prefix + ".tmp").c_str());
				std::remove((prefix + "-w.tmp").c_str());
				std::remove((prefix + "-m.tmp").c_str());
			}
		}
		std::remove(_data->_metaFile.c_str());
		delete _data;
	}
}

void PartitionedMS::MakeMSRowToRowIdMapping(std::vector<size_t>& msToId, const MSSelection& selection)
{
	const size_t nRow = _ms.nrow();
	casacore::ROArrayColumn<double> uvwColumn(_ms, casacore::MS::columnName(casacore::MSMainEnums::UVW));
	casacore::ROScalarColumn<int> antenna1Column(_ms, casacore::MS::columnName(casacore::MSMainEnums::ANTENNA1));
	casacore::ROScalarColumn<int> antenna2Column(_ms, casacore::MS::columnName(casacore::MSMainEnums::ANTENNA2));
	casacore::ROScalarColumn<int> fieldIdColumn(_ms, casacore::MS::columnName(casacore::MSMainEnums::FIELD_ID));
	casacore::ROScalarColumn<double> timeColumn(_ms, casacore::MS::columnName(casacore::MSMainEnums::TIME));
	size_t startRow, endRow;
	getRowRange(_ms, selection, startRow, endRow);
	
	msToId.assign(startRow, 0);
	size_t currentRowId = 0;
	size_t timestep = selection.HasInterval() ? selection.IntervalStart() : 0;
	double time = timeColumn(startRow);
	for(size_t row=startRow; row!=endRow; ++row)
	{
		msToId.push_back(currentRowId);
		const int
			a1 = antenna1Column(row), a2 = antenna2Column(row),
			fieldId = fieldIdColumn(row);
		casacore::Vector<double> uvw = uvwColumn(row);
		if(time != timeColumn(row))
		{
			++timestep;
			time = timeColumn(row);
		}
		if(selection.IsSelected(fieldId, timestep, a1, a2, uvw))
			++currentRowId;
	}
	for(size_t i=0; i!=nRow-endRow; ++i)
		msToId.push_back(0);
}
