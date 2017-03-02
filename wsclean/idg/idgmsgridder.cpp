#include "idgmsgridder.h"

//#include "interface.h"
//#include "dummygridder.h"
#include <idg.h>
#include <idg-utility.h>

#include "../msproviders/msprovider.h"

#include <boost/thread/thread.hpp>

#include "../wsclean/logger.h"

IdgMsGridder::IdgMsGridder() :
	_subgridSize(24),
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
    
    _taper_subgrid.resize(_subgridSize);
    _taper_grid.resize(width);
    idg::init_optimal_taper_1D(_subgridSize, width, 7.0, 1.25, _taper_subgrid.data(), _taper_grid.data());
    

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
		std::cout << "total weight: " << totalWeight() << std::endl;
        
		// TODO use untrimmed dimensions and perform trimming
		_grid.assign(width * height * 4, 0.0);
		
		for(size_t i=0; i!=MeasurementSetCount(); ++i)
		{
			gridMeasurementSet(msDataVector[i]);
		}
		
		idg::GridderPlan* plan = new idg::GridderPlan(idg::Type::CPU_OPTIMIZED);
		plan->set_grid(4, height, width, _grid.data());
// 		plan->set_cell_size(_actualPixelSizeX, _actualPixelSizeY);
// 		plan->set_w_kernel(7);
// 		plan->internal_set_subgrid_size(_subgridSize);
                
		plan->transform_grid();
		plan->bake();
        
        std::cout << "total weight: " << totalWeight() << std::endl;
        
        float inv_spheroidal[width];
        for(size_t i=0; i < width ; ++i)
        {
            float y = _taper_grid[i];
            inv_spheroidal[i] = 0.0;
            if (y > 0.001) inv_spheroidal[i] = 1.0/y;
        }
        
        for(size_t ii=0; ii != width * height; ++ii)
        {
            size_t i = ii % width;
            size_t j = ii / width;
            _grid[ii] *= 0.5*inv_spheroidal[i]*inv_spheroidal[j]/totalWeight()*width*height;
        }

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
	const size_t width = TrimWidth(), height = TrimHeight();

    ao::uvector<float> taper;
    taper.resize(_subgridSize * _subgridSize);
    for(int i = 0; i < _subgridSize; i++)
        for(int j = 0; j < _subgridSize; j++)
            taper[i*_subgridSize + j] = _taper_subgrid[i] * _taper_subgrid[j];
		
	// A gridder per band is needed
    if(constructDegridders)
		_degridderPlans.resize(selectedBands.BandCount());
	else
		_gridderPlans.resize(selectedBands.BandCount());
	ao::uvector<std::complex<double>> aterm(nStations * 4 * _subgridSize * _subgridSize, 0.0); // TODO fill with sensible data
	//aterm[0] = 1.0; // TODO
	aterm[(_subgridSize / 2) * (_subgridSize+1) * 4] = 1.0; // TODO
	aterm[(_subgridSize / 2) * (_subgridSize+1) * 4+1] = 0.0; // TODO
	aterm[(_subgridSize / 2) * (_subgridSize+1) * 4+2] = 0.0; // TODO
	aterm[(_subgridSize / 2) * (_subgridSize+1) * 4+3] = 1.0; // TODO
	for(size_t i=0; i!=selectedBands.BandCount(); ++i)
	{
		ao::uvector<double> frequencyList(selectedBands[i].begin(), selectedBands[i].end());
		idg::Scheme* plan;
		if(constructDegridders)
		{
//			_degridderPlans[i] = new idg::DegridderPlan(idg::Type::CUDA_GENERIC, 128);
 			_degridderPlans[i] = new idg::DegridderPlan(idg::Type::CPU_OPTIMIZED, 256);
			plan = _degridderPlans[i];
		}
		else {
//			_gridderPlans[i] = new idg::GridderPlan(idg::Type::CUDA_GENERIC, 128);
 			_gridderPlans[i] = new idg::GridderPlan(idg::Type::CPU_OPTIMIZED, 256);
			plan = _gridderPlans[i];
		}
		plan->set_frequencies(selectedBands[i].ChannelCount(), frequencyList.data());
		plan->set_stations(nStations);
		plan->set_grid(4, height, width, _grid.data());
		plan->set_spheroidal(_subgridSize, _subgridSize, taper.data());
		plan->set_cell_size(_actualPixelSizeX, _actualPixelSizeY);
//		plan->set_image_size(_actualPixelSizeX*width);
		int w_kernel_size = (3*_subgridSize / 4) + 1;
		plan->set_w_kernel(w_kernel_size);
		plan->internal_set_subgrid_size(_subgridSize);
		plan->bake();
		
		plan->start_w_layer(0.0); //TODO
                Logger::Info << "start_w_layer\n";
		
		// This doesn't seem to have effect yet:
		//plan->start_aterm(nStations, _subgridSize, _subgridSize, 4, aterm.data());
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
		rowData.uvw[1] = -vInMeters;  // DEBUG vdtol, flip axis
		rowData.uvw[2] = -wInMeters;
		rowData.antenna1 = antenna1Col(rowIndex);
		rowData.antenna2 = antenna2Col(rowIndex);
		rowData.timeIndex = timeIndex;
		rowData.dataDescId = dataDescId;
		
		readAndWeightVisibilities<4>(*msData.msProvider, rowData, curBand, weightBuffer.data(), modelBuffer.data(), isSelected.data());
        idg::GridderPlan* plan = _gridderPlans[rowData.dataDescId];
				// For debug, comment out next line to don't let IDG hang
        plan->grid_visibilities(timeIndex, antenna1Col(rowIndex), antenna2Col(rowIndex), rowData.uvw, rowData.data);
        delete[] rowData.data;
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

void IdgMsGridder::Predict(double* image)
{
	const size_t width = TrimWidth(), height = TrimHeight();
    _taper_subgrid.resize(_subgridSize);
    _taper_grid.resize(width);
    idg::init_optimal_taper_1D(_subgridSize, width, 7.0, 1.25, _taper_subgrid.data(), _taper_grid.data());
	
        
	// Stokes V is always the last requested pol. So, only when Stokes V
	// is requested, do the actual prediction. Since all pols are predicted at once by IDG,
	// when Stokes I/Q/U prediction is requested, one waits until V is also predicted.
	if(Polarization() == Polarization::StokesI)
	{
      _grid.resize(width * height * 4);

      // Fill the grid with the image information for this polarization
      _grid.resize(width * height * 4);
        
        float inv_spheroidal[width];
        for(size_t i=0; i < width ; ++i)
        {
            float y = _taper_grid[i];
            inv_spheroidal[i] = 0.0;
            if (y > 0.001) inv_spheroidal[i] = 1.0/y;
        }
        
        for(size_t ii=0; ii != width * height; ++ii)
        {
            size_t i = ii % width;
            size_t j = ii / width;
            _grid[ii].real(image[ii]*inv_spheroidal[i]*inv_spheroidal[j]);
            _grid[ii + 3*width*height].real(image[ii]*inv_spheroidal[i]*inv_spheroidal[j]);
        }
        
        
		// use a temporary gridder to transform the grid

		std::unique_ptr<idg::DegridderPlan> iface( new idg::DegridderPlan(idg::Type::CPU_OPTIMIZED) );

		iface->set_grid(4, height, width, _grid.data());
                
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
        msData.msProvider->ReopenRW();

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
}

void IdgMsGridder::predictCalcThreadFunction()
{
	IDGPredictionRow row;
        // read the first sample, exit if none is there
        if (!_predictionCalcLane.read(row)) return;
        
        
        casa::MeasurementSet ms = _outputProvider->MS();
        casacore::ArrayColumn<casacore::Complex> modelColumn(ms, casacore::MS::columnName(casacore::MSMainEnums::MODEL_DATA));
        const casacore::IPosition shape(modelColumn.shape(0));
        vector<size_t> idToMSRow;
        _outputProvider->MakeIdToMSRowMapping(idToMSRow);
        
	while(true)
	{
		idg::DegridderPlan* plan = _degridderPlans[row.dataDescId];
		bool isBufferFull = plan->request_visibilities(row.rowId, row.timeIndex, row.antenna1, row.antenna2, row.uvw);
                
                // compute if the buffer is full, or if there are no more samples
                // if the buffer is full, no new sample will be read
		if (isBufferFull || !_predictionCalcLane.read(row))
		{
			auto available_row_ids = plan->compute();
//                         std::cout << "computed " << available_row_ids.size() << " rows" << std::endl;
			for(auto i : available_row_ids)
			{
                          size_t rowid(i.first);
                          std::complex<float>* dataptr(i.second);
                          casacore::Array<std::complex<float>> modeldataArray(shape, dataptr, casacore::SHARE);

                          modelColumn.put(idToMSRow[rowid], modeldataArray);
                          //                             _outputProvider->WriteModel(i.first, i.second);
//                             for(int ii = -4;;)
//                             {
//                                 std::cout << i.second[ii];
//                                 if (++ii == 24) break;
//                                 std::cout << ", ";
//                             }
//                             std::cout << std::endl;
			}
			plan->finished_reading();
                        
                        // The buffer was not full
                        // So, we were computing because there were no more samples
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
