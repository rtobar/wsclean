#include "idgmsgridder.h"

#include <cmath>

//#include "interface.h"
//#include "dummygridder.h"
#include <idg-api.h>

#include "../msproviders/msprovider.h"

#include <boost/thread/thread.hpp>

#include "../wsclean/logger.h"

IdgMsGridder::IdgMsGridder() :
	_inversionLane(1024),
	_predictionCalcLane(1024),
	_predictionWriteLane(1024),
	_outputProvider(nullptr)
{ 
	Logger::Info << "number of MS: " << MeasurementSetCount() << "\n";;
}

IdgMsGridder::~IdgMsGridder()
{ }

void IdgMsGridder::Invert()
{
	const size_t width = TrimWidth(), height = TrimHeight();
	
	// Stokes I is always the first requested pol. So, only when Stokes I
	// is requested, do the actual inversion. Since all pols are produced at once by IDG,
	// when Stokes Q/U/V is requested, the result of earlier gridding is returned.
	// For the first inversion, WSClean will ask for the PSF of Stokes I. The next run
	// is the dirty of Stokes I, which should thus overwrite all images.
	if(Polarization() == Polarization::StokesI)
	{
		std::vector<MSData> msDataVector;
		initializeMSDataVector(msDataVector, 4);
		
		resetVisibilityCounters();

		_image.assign(4 * width * height, 0.0);
		
		for(size_t i=0; i!=MeasurementSetCount(); ++i)
		{
			// Adds the gridding result to _image member
			gridMeasurementSet(msDataVector[i]);
		}
		
		std::cout << "total weight: " << totalWeight() << std::endl;
		
		
		// Normalize by total weight
		
		for(size_t ii=0; ii != 4 * width * height; ++ii)
		{
			_image[ii] = _image[ii]/totalWeight();
		}
		
		// result is now in _image member
		// Can be accessed by subsequent calls to ImageRealResult()
		
	}
	else if(_image.empty()) {
		throw std::runtime_error("IdgMsGridder::Invert() was called out of sequence");
	}
}

void IdgMsGridder::gridMeasurementSet(MSGridderBase::MSData& msData)
{
	const size_t width = TrimWidth();
	const float max_w = msData.maxW;
	_selectedBands = msData.SelectedBand();

	// TODO for now we map the ms antennas directly to the gridder's antenna,
	// including non-selected antennas. Later this can be made more efficient.
	casacore::MeasurementSet& ms = msData.msProvider->MS();
	size_t nStations = ms.antenna().nrow();

	std::vector<std::vector<double>> bands;
	for(size_t i=0; i!=_selectedBands.BandCount(); ++i)
	{
		bands.push_back(std::vector<double>(_selectedBands[i].begin(), _selectedBands[i].end()));
	}

	_bufferset = std::unique_ptr<idg::api::BufferSet>(idg::api::BufferSet::create(
			idg::api::Type::CPU_OPTIMIZED, 
//			idg::api::Type::HYBRID_CUDA_CPU_OPTIMIZED,
			128, bands, nStations, 
			width, _actualPixelSizeX, max_w, idg::api::BufferSetType::gridding));
	
	casacore::ScalarColumn<int> antenna1Col(ms, casacore::MeasurementSet::columnName(casacore::MSMainEnums::ANTENNA1));
	casacore::ScalarColumn<int> antenna2Col(ms, casacore::MeasurementSet::columnName(casacore::MSMainEnums::ANTENNA2));
	casacore::ScalarColumn<double> timeCol(ms, casacore::MeasurementSet::columnName(casacore::MSMainEnums::TIME));
	
	_inversionLane.clear();
	
	ao::uvector<float> weightBuffer(_selectedBands.MaxChannels()*4);
	ao::uvector<std::complex<float>> modelBuffer(_selectedBands.MaxChannels()*4);
	ao::uvector<bool> isSelected(_selectedBands.MaxChannels()*4, true);
	std::vector<size_t> idToMSRow;
	msData.msProvider->MakeIdToMSRowMapping(idToMSRow);
	// The gridder doesn't need to know the absolute time index; this value indexes relatively to where we
	// start in the measurement set, and only increases when the time changes.
	int timeIndex = -1;
	double currentTime = -1.0;
	for(msData.msProvider->Reset() ; msData.msProvider->CurrentRowAvailable() ; msData.msProvider->NextRow())
	{
		size_t rowIndex = idToMSRow[msData.msProvider->RowId()];
		if(currentTime != timeCol(rowIndex))
		{
			currentTime = timeCol(rowIndex);
			timeIndex++;
		}
		size_t dataDescId;
		double uInMeters, vInMeters, wInMeters;
		msData.msProvider->ReadMeta(uInMeters, vInMeters, wInMeters, dataDescId);
		const BandData& curBand(_selectedBands[dataDescId]);
		IDGInversionRow rowData;
		
		rowData.data = new std::complex<float>[curBand.ChannelCount()*4];
		std::complex<float> *data4 = new std::complex<float>[curBand.ChannelCount()*4];

		rowData.uvw[0] = uInMeters;
		rowData.uvw[1] = vInMeters;
		rowData.uvw[2] = wInMeters;
		
		rowData.antenna1 = antenna1Col(rowIndex);
		rowData.antenna2 = antenna2Col(rowIndex);
		rowData.timeIndex = timeIndex;
		rowData.dataDescId = dataDescId;
		
		readAndWeightVisibilities<4>(*msData.msProvider, rowData, curBand, weightBuffer.data(), modelBuffer.data(), isSelected.data());

		rowData.uvw[1] = -vInMeters;  // DEBUG vdtol, flip axis
		rowData.uvw[2] = -wInMeters;  //

		_bufferset->get_gridder(rowData.dataDescId)->grid_visibilities(timeIndex, antenna1Col(rowIndex), antenna2Col(rowIndex), rowData.uvw, rowData.data);
		delete[] data4;
		delete[] rowData.data;
	}
	_inversionLane.write_end();
	
	// TODO needs to add, not replace, because gridMeasurementSet is called in a loop over measurement sets
	_bufferset->finished();
	_bufferset->get_image(_image.data());
	//_bufferset.reset();
}

void IdgMsGridder::Predict(double* image)
{
	const size_t width = TrimWidth(), height = TrimHeight();

	if (Polarization() == Polarization::StokesI)
	{
		_image.assign(4 * width * height, 0.0);
	}

	size_t polIndex = Polarization::StokesToIndex(Polarization());
	for(size_t i=0; i != width * height; ++i)
		_image[i + polIndex*width*height] = image[i];

	// Stokes V is always the last requested pol. So, only when Stokes V
	// is requested, do the actual prediction. Since all pols are predicted at once by IDG,
	// when Stokes I/Q/U prediction is requested, one waits until V is also predicted.
	// Stokes V is always requested last
	
	// !! DOES NOT WORK because a new IdgMsGridder object is created for each polarization 
	// (when wsclean is called with -predict option, not sure about predict step while cleaning)

	// So now the  actual predict can happen
	if (Polarization() == Polarization::StokesI)
	{
		// Do actual predict
		std::vector<MSData> msDataVector;
		initializeMSDataVector(msDataVector, 4);

		for(size_t i=0; i!=MeasurementSetCount(); ++i)
		{
			predictMeasurementSet(msDataVector[i]);
		}
	}
}

void IdgMsGridder::predictMeasurementSet(MSGridderBase::MSData& msData)
{
	const size_t width = TrimWidth();
	const float max_w = msData.maxW;


	msData.msProvider->ReopenRW();

	_selectedBands = msData.SelectedBand();
	_outputProvider = msData.msProvider;
	
	casacore::MeasurementSet& ms = msData.msProvider->MS();
	size_t nr_stations = ms.antenna().nrow();

	std::vector<std::vector<double>> bands;
	for(size_t i=0; i!=_selectedBands.BandCount(); ++i)
	{
		bands.push_back(std::vector<double>(_selectedBands[i].begin(), _selectedBands[i].end()));
	}

	_bufferset = std::unique_ptr<idg::api::BufferSet>(idg::api::BufferSet::create(
			idg::api::Type::CPU_OPTIMIZED,
//			idg::api::Type::HYBRID_CUDA_CPU_OPTIMIZED, 
			128, bands, nr_stations, 
			width, _actualPixelSizeX, max_w, idg::api::BufferSetType::degridding));
	_bufferset->set_image(_image.data());

	casacore::ScalarColumn<int> antenna1Col(ms, casacore::MeasurementSet::columnName(casacore::MSMainEnums::ANTENNA1));
	casacore::ScalarColumn<int> antenna2Col(ms, casacore::MeasurementSet::columnName(casacore::MSMainEnums::ANTENNA2));
	casacore::ScalarColumn<double> timeCol(ms, casacore::MeasurementSet::columnName(casacore::MSMainEnums::TIME));
	
	_predictionCalcLane.clear();
	_predictionWriteLane.clear();
	boost::mutex mutex;
	boost::thread predictCalcThread(&IdgMsGridder::predictCalcThreadFunction, this);
	boost::thread predictWriteThread(&IdgMsGridder::predictWriteThreadFunction, this, &mutex);

	ao::uvector<std::complex<float>> buffer(_selectedBands.MaxChannels()*4);
	std::vector<size_t> idToMSRow;
	msData.msProvider->MakeIdToMSRowMapping(idToMSRow);
	int timeIndex = -1;
	double currentTime = -1.0;
	// The mutex needs to be locked during the while-guard and other IO access
	boost::mutex::scoped_lock lock(mutex);
	for(msData.msProvider->Reset() ; msData.msProvider->CurrentRowAvailable() ; msData.msProvider->NextRow())
	{
		size_t provRowId = msData.msProvider->RowId();
		size_t msRowIndex = idToMSRow[provRowId];
		if(currentTime != timeCol(msRowIndex))
		{
			currentTime = timeCol(msRowIndex);
			timeIndex++;
		}
		size_t dataDescId;
		double uInMeters, vInMeters, wInMeters;
		msData.msProvider->ReadMeta(uInMeters, vInMeters, wInMeters, dataDescId);
		lock.unlock();
		
		IDGPredictionRow row;
		row.uvw[0] = uInMeters;
		row.uvw[1] = -vInMeters;
		row.uvw[2] = -wInMeters;
		row.antenna1 = antenna1Col(msRowIndex);
		row.antenna2 = antenna2Col(msRowIndex);
		row.timeIndex = timeIndex;
		row.dataDescId = dataDescId;
		row.rowId = provRowId;
		
		// (For debug) comment this out to don't let IDG hang
		_predictionCalcLane.write(row);
		
		lock.lock(); // lock for guard
	}
	lock.unlock();
	
	_predictionCalcLane.write_end();
	_predictionWriteLane.write_end();
	
	predictWriteThread.join();
	predictCalcThread.join();
	_bufferset.reset();
}

void IdgMsGridder::predictCalcThreadFunction()
{
	IDGPredictionRow row;
	
	// read the first sample, exit if none is there
	if (!_predictionCalcLane.read(row)) return;
	
	//casacore::MeasurementSet ms = _outputProvider->MS();
	//casacore::ArrayColumn<casacore::Complex> modelColumn(ms, casacore::MS::columnName(casacore::MSMainEnums::MODEL_DATA));
	//const casacore::IPosition shape(modelColumn.shape(0));
	vector<size_t> idToMSRow;
	_outputProvider->MakeIdToMSRowMapping(idToMSRow);

	while(true)
	{
		bool isBufferFull = _bufferset->get_degridder(row.dataDescId)->request_visibilities(row.rowId, row.timeIndex, row.antenna1, row.antenna2, row.uvw);

		// compute if the buffer is full, or if there are no more samples
		// if the buffer is full, no new sample will be read
		if (isBufferFull || !_predictionCalcLane.read(row))
		{
			auto available_row_ids = _bufferset->get_degridder(row.dataDescId)->compute();
			std::cout << "computed " << available_row_ids.size() << " rows" << std::endl;
			for(auto i : available_row_ids)
			{
				//size_t rowid(i.first);
				//std::complex<float>* dataptr(i.second);
				//casacore::Array<std::complex<float>> modeldataArray(shape, dataptr, casacore::SHARE);
				//modelColumn.put(idToMSRow[rowid], modeldataArray);
				
				_outputProvider->WriteModel(i.first, i.second);
			}
			_bufferset->get_degridder(row.dataDescId)->finished_reading();
			
			// If the buffer is not full
			// we were computing because there were no more samples, return.
			if (!isBufferFull) return;
		}
	}
}

void IdgMsGridder::predictWriteThreadFunction(boost::mutex* mutex)
{
	IDGRowForWriting row;
	while(_predictionWriteLane.read(row))
	{
		boost::mutex::scoped_lock lock(*mutex);
		// TODO we should not write visibilities that were outside the w-range
// 		_outputProvider->WriteModel(row.rowId, row.data);
		delete[] row.data;
	}
}

void IdgMsGridder::Predict(double* real, double* imaginary)
{
	throw std::runtime_error("IDG gridder cannot make complex images");
}

double* IdgMsGridder::ImageRealResult()
{
	const size_t width = TrimWidth(), height = TrimHeight();
	size_t polIndex = Polarization::StokesToIndex(Polarization());
	return _image.data() + height*width*polIndex;
}

double* IdgMsGridder::ImageImaginaryResult()
{
	throw std::runtime_error("IDG gridder cannot make complex images");
}

void IdgMsGridder::GetGriddingCorrectionImage(double* image) const
{
	const size_t width = TrimWidth(), height = TrimHeight();
	for(size_t i=0; i!=width*height; ++i)
		image[i] = 1.0;
}

bool IdgMsGridder::HasGriddingCorrectionImage() const
{
	return false; // For now (TODO)
}

