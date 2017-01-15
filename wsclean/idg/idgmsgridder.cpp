#include "idgmsgridder.h"

//#include "interface.h"
//#include "dummygridder.h"
#include <idg.h>

#include "../msproviders/msprovider.h"

#include <boost/thread/thread.hpp>

IdgMsGridder::IdgMsGridder() :
	_kernelSize(32),
	_inversionLane(1024),
	_predictionCalcLane(1024),
	_predictionWriteLane(1024),
	_outputProvider(nullptr)
{ }

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
	if(Polarization() == ::Polarization::StokesI)
	{
		std::vector<MSData> msDataVector;
		initializeMSDataVector(msDataVector, 4);
		
		resetVisibilityCounters();
		
		// TODO use untrimmed dimensions and perform trimming
		_grid.assign(width * height * 4, 0.0);
		
		for(size_t i=0; i!=MeasurementSetCount(); ++i)
		{
			gridMeasurementSet(msDataVector[i]);
		}
		
		idg::GridderPlan* plan = new idg::GridderPlan();
		plan->set_grid(4, height, width, _grid.data());
		plan->set_spheroidal(_kernelSize, _kernelSize, _kernel.data());
		plan->set_cell_size(_actualPixelSizeX, _actualPixelSizeY);
		int kernel_size = (_kernelSize / 4) + 1;
		plan->set_w_kernel(kernel_size);
		plan->internal_set_subgrid_size(_kernelSize);
		plan->transform_grid();
		plan->bake();
		delete plan;
	}
	else if(_grid.empty()) {
		throw std::runtime_error("IdgMsGridder::Invert() was called out of sequence");
	}
	
	// Copy the requested polarization result into the output image
	_image.resize(width * height);
	size_t polIndex = Polarization::StokesToIndex(Polarization());
	for(size_t i=0; i != width * height; ++i)
		_image[i] = _grid[i + polIndex*width*height].real();
}

void IdgMsGridder::constructGridders(const MultiBandData& selectedBands, size_t nStations, bool constructDegridders)
{
	// TODO This should be the prolate spheroidal in image domain
	_kernel.assign(_kernelSize * _kernelSize, 1.0);  // TODO this is NN gridding
		
	// A gridder per band is needed
	if(constructDegridders)
		_degridderPlans.resize(selectedBands.BandCount());
	else
		_gridderPlans.resize(selectedBands.BandCount());
	ao::uvector<std::complex<double>> aterm(nStations * 4 * _kernelSize * _kernelSize, 0.0); // TODO fill with sensible data
	//aterm[0] = 1.0; // TODO
	aterm[(_kernelSize / 2) * (_kernelSize+1) * 4] = 1.0; // TODO
	aterm[(_kernelSize / 2) * (_kernelSize+1) * 4+1] = 0.0; // TODO
	aterm[(_kernelSize / 2) * (_kernelSize+1) * 4+2] = 0.0; // TODO
	aterm[(_kernelSize / 2) * (_kernelSize+1) * 4+3] = 1.0; // TODO
	const size_t width = TrimWidth(), height = TrimHeight();
	for(size_t i=0; i!=selectedBands.BandCount(); ++i)
	{
		ao::uvector<double> frequencyList(selectedBands[i].begin(), selectedBands[i].end());
		idg::Scheme* plan;
		if(constructDegridders)
		{
			_degridderPlans[i] = new idg::DegridderPlan(idg::Type::CPU_REFERENCE, 1024);
			plan = _degridderPlans[i];
		}
		else {
			_gridderPlans[i] = new idg::GridderPlan(idg::Type::CPU_REFERENCE, 1024*10);
			plan = _gridderPlans[i];
		}
		plan->set_frequencies(selectedBands[i].ChannelCount(), frequencyList.data());
		plan->set_stations(nStations);
		plan->set_grid(4, height, width, _grid.data());
		plan->set_spheroidal(_kernelSize, _kernelSize, _kernel.data());
		plan->set_cell_size(_actualPixelSizeX, _actualPixelSizeY);
//		plan->set_image_size(_actualPixelSizeX*width);
		int kernel_size = (_kernelSize / 4) + 1;
		plan->set_w_kernel(kernel_size);
		plan->internal_set_subgrid_size(_kernelSize);
		plan->bake();
		
		plan->start_w_layer(0.0); //TODO
		
		// This doesn't seem to have effect yet:
		//plan->start_aterm(nStations, _kernelSize, _kernelSize, 4, aterm.data());
	}
}

void IdgMsGridder::gridMeasurementSet(MSGridderBase::MSData& msData)
{
	MultiBandData selectedBands = msData.SelectedBand();

	// TODO for now we map the ms antennas directly to the gridder's antenna,
	// including non-selected antennas. Later this can be made more efficient.
	casa::MeasurementSet& ms = msData.msProvider->MS();
	size_t nStations = ms.antenna().nrow();

	constructGridders(selectedBands, nStations, false);
	
	casacore::ScalarColumn<int> antenna1Col(ms, casacore::MeasurementSet::columnName(casacore::MSMainEnums::ANTENNA1));
	casacore::ScalarColumn<int> antenna2Col(ms, casacore::MeasurementSet::columnName(casacore::MSMainEnums::ANTENNA2));
	casacore::ScalarColumn<double> timeCol(ms, casacore::MeasurementSet::columnName(casacore::MSMainEnums::TIME));
	
	_inversionLane.clear();
	boost::thread gridThread(&IdgMsGridder::gridThreadFunction, this);
	
	ao::uvector<float> weightBuffer(selectedBands.MaxChannels()*4);
	ao::uvector<std::complex<float>> modelBuffer(selectedBands.MaxChannels()*4);
	ao::uvector<bool> isSelected(selectedBands.MaxChannels()*4, true); // TODO make w-pass-aware
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
		const BandData& curBand(selectedBands[dataDescId]);
		IDGInversionRow rowData;
		rowData.data = new std::complex<float>[curBand.ChannelCount()*4];
		rowData.uvw[0] = uInMeters;
		rowData.uvw[1] = vInMeters;
		rowData.uvw[2] = wInMeters;
		rowData.antenna1 = antenna1Col(rowIndex);
		rowData.antenna2 = antenna2Col(rowIndex);
		rowData.timeIndex = timeIndex;
		rowData.dataDescId = dataDescId;
		
		readAndWeightVisibilities<4>(*msData.msProvider, rowData, curBand, weightBuffer.data(), modelBuffer.data(), isSelected.data());
		
		_inversionLane.write(rowData);
	}
	_inversionLane.write_end();
	gridThread.join();
	
	for(idg::GridderPlan* iface : _gridderPlans)
	{
		iface->finish_aterm();
		iface->finish_w_layer();
		iface->finished();
		
		delete iface;
	}
}

void IdgMsGridder::gridThreadFunction()
{
	IDGInversionRow row;
	while(_inversionLane.read(row))
	{
		idg::GridderPlan* plan = _gridderPlans[row.dataDescId];
		plan->grid_visibilities(row.timeIndex, row.antenna1, row.antenna2, row.uvw, row.data);
		delete[] row.data;
	}
}

void IdgMsGridder::Predict(double* real)
{
	const size_t width = TrimWidth(), height = TrimHeight();
	
	// Fill the grid with the image information for this polarization
	_grid.resize(width * height * 4);
	size_t polIndex = Polarization::StokesToIndex(Polarization());
	for(size_t i=0; i != width * height; ++i)
		_grid[i + polIndex*width*height].real(_image[i]);

	// Stokes V is always the last requested pol. So, only when Stokes V
	// is requested, do the actual prediction. Since all pols are predicted at once by IDG,
	// when Stokes I/Q/U prediction is requested, one waits until V is also predicted.
	if(Polarization() == ::Polarization::StokesI)
	{
		_grid.resize(width * height * 4);
		
		// use a temporary gridder to transform the grid
		std::unique_ptr<idg::GridderPlan> iface( new idg::GridderPlan() ); // TODO;
		iface->set_grid(4, height, width, _grid.data());
		iface->set_spheroidal(_kernelSize, _kernel.data());
		iface->transform_grid();
		
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
	_selectedBands = msData.SelectedBand();
	_outputProvider = msData.msProvider;
	
	casa::MeasurementSet& ms = msData.msProvider->MS();
	size_t nStations = ms.antenna().nrow();

	constructGridders(_selectedBands, nStations, true);
	
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
		row.uvw[1] = vInMeters;
		row.uvw[2] = wInMeters;
		row.antenna1 = antenna1Col(msRowIndex);
		row.antenna2 = antenna2Col(msRowIndex);
		row.timeIndex = timeIndex;
		row.dataDescId = dataDescId;
		row.rowId = provRowId;
		
		_predictionCalcLane.write(row);
		
		lock.lock(); // lock for guard
	}
	lock.unlock();
	
	_predictionCalcLane.write_end();
	_predictionWriteLane.write_end();
	
	predictWriteThread.join();
	predictCalcThread.join();
}

void IdgMsGridder::predictCalcThreadFunction()
{
	IDGPredictionRow row;
	while(_predictionCalcLane.read(row))
	{
		idg::DegridderPlan* plan = _degridderPlans[row.dataDescId];
		// TODO was the buffer before or after the request full?
		bool isBufferFull = plan->request_visibilities(row.rowId, row.timeIndex, row.antenna1, row.antenna2, row.uvw);
		if(isBufferFull)
		{
			// Get the band belonging to the gridder that is full
			const BandData& curBand(_selectedBands[row.dataDescId]);
			
			auto available_row_ids = plan->compute();
			for(size_t i : available_row_ids)
			{
				IDGRowForWriting writeRow;
				writeRow.data = new std::complex<float>[curBand.ChannelCount()*4];
				plan->read_visibilities(i, writeRow.data); 
				writeRow.rowId = i;
				_predictionWriteLane.write(writeRow);
			}
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
		_outputProvider->WriteModel(row.rowId, row.data);
	}
}

void IdgMsGridder::Predict(double* real, double* imaginary)
{
	throw std::runtime_error("IDG gridder cannot make complex images");
}

double* IdgMsGridder::ImageRealResult()
{
	return _image.data();
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
