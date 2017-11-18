// kate: space-indent off; tab-width 2; indent-width 2; replace-tabs off; eol unix;

#include "idgmsgridder.h"

#include <cmath>
#include <fstream>
#include <thread>

//#include "interface.h"
//#include "dummygridder.h"
#include <idg-api.h>

#include "../msproviders/msprovider.h"

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include "../wsclean/logger.h"
#include "../wsclean/wscleansettings.h"

#include "../aterms/lofarbeamterm.h"

IdgMsGridder::IdgMsGridder(const WSCleanSettings& settings) :
	_outputProvider(nullptr),
	_settings(settings),
	_proxyType(idg::api::Type::CPU_OPTIMIZED),
	_buffersize(256)
{
	readConfiguration();
	setIdgType();
	_bufferset = std::unique_ptr<idg::api::BufferSet>(
		idg::api::BufferSet::create(_proxyType));
}

IdgMsGridder::~IdgMsGridder()
{ }

void IdgMsGridder::Invert()
{
	const size_t untrimmedWidth = ImageWidth();
	const size_t width = TrimWidth(), height = TrimHeight();

	assert(width == height);
	assert(untrimmedWidth == ImageHeight());

	_options["padded_size"] = untrimmedWidth;

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

		double max_w = 0;
		for(size_t i=0; i!=MeasurementSetCount(); ++i)
		{
			max_w = std::max(max_w, msDataVector[i].maxW);
		}

		_bufferset->init(width, _actualPixelSizeX, max_w, _options);
		for(size_t i=0; i!=MeasurementSetCount(); ++i)
		{
			// Adds the gridding result to _image member
			gridMeasurementSet(msDataVector[i]);
		}
		
		std::cout << "total weight: " << totalWeight() << std::endl;
		_image.assign(4 * width * height, 0.0);
		_bufferset->get_image(_image.data());
		
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
	_selectedBands = msData.SelectedBand();

	// TODO for now we map the ms antennas directly to the gridder's antenna,
	// including non-selected antennas. Later this can be made more efficient.
	casacore::MeasurementSet& ms = msData.msProvider->MS();
	size_t nr_stations = ms.antenna().nrow();
	
	std::vector<std::vector<double>> bands;
	for(size_t i=0; i!=_selectedBands.BandCount(); ++i)
	{
		bands.push_back(std::vector<double>(_selectedBands[i].begin(), _selectedBands[i].end()));
	}
	const float max_baseline = msData.maxBaselineInM;

	_bufferset->init_buffers(_buffersize, bands, nr_stations, max_baseline, _options, idg::api::BufferSetType::gridding);
	
	std::unique_ptr<LofarBeamTerm> aTermMaker;
	ao::uvector<std::complex<float>> aTermBuffer;
	double aTermUpdateInterval = 900.0; // 15min
	if(_settings.gridWithBeam)
	{
		size_t subgridsize = _bufferset->get_subgridsize();
		double dl = _actualPixelSizeX * TrimWidth() / subgridsize, dm = _actualPixelSizeY * TrimHeight() / subgridsize;
		double pdl = PhaseCentreDL(), pdm = PhaseCentreDM();
		aTermMaker.reset(new LofarBeamTerm(ms, subgridsize, subgridsize, dl, dm, pdl, pdm, _settings.useDifferentialLofarBeam));
		aTermBuffer.resize(subgridsize*subgridsize*4*nr_stations);
	}

	ao::uvector<float> weightBuffer(_selectedBands.MaxChannels()*4);
	ao::uvector<std::complex<float>> modelBuffer(_selectedBands.MaxChannels()*4);
	ao::uvector<bool> isSelected(_selectedBands.MaxChannels()*4, true);
	ao::uvector<std::complex<float>> dataBuffer(_selectedBands.MaxChannels()*4);
	// The gridder doesn't need to know the absolute time index; this value indexes relatively to where we
	// start in the measurement set, and only increases when the time changes.
	int timeIndex = -1;
	double currentTime = -1.0, lastATermUpdate = -aTermUpdateInterval-1;
	for(msData.msProvider->Reset() ; msData.msProvider->CurrentRowAvailable() ; msData.msProvider->NextRow())
	{
		MSProvider::MetaData metaData;
		msData.msProvider->ReadMeta(metaData);
		if(currentTime != metaData.time)
		{
			currentTime = metaData.time;
			timeIndex++;
		}
		const BandData& curBand(_selectedBands[metaData.dataDescId]);
		IDGInversionRow rowData;
		
		rowData.data = dataBuffer.data();
		rowData.uvw[0] = metaData.uInM;
		rowData.uvw[1] = metaData.vInM;
		rowData.uvw[2] = metaData.wInM;
		
		rowData.antenna1 = metaData.antenna1;
		rowData.antenna2 = metaData.antenna2;
		rowData.timeIndex = timeIndex;
		rowData.dataDescId = metaData.dataDescId;
		readAndWeightVisibilities<4>(*msData.msProvider, rowData, curBand, weightBuffer.data(), modelBuffer.data(), isSelected.data());

		rowData.uvw[1] = -metaData.vInM;  // DEBUG vdtol, flip axis
		rowData.uvw[2] = -metaData.wInM;  //

		if(_settings.gridWithBeam && currentTime - lastATermUpdate > aTermUpdateInterval)
		{
			Logger::Debug << "Calculating a-terms for timestep " << timeIndex << "\n";
			aTermMaker->Calculate(aTermBuffer.data(), currentTime + aTermUpdateInterval*0.5, _selectedBands.CentreFrequency());
			_bufferset->get_gridder(rowData.dataDescId)->set_aterm(timeIndex, aTermBuffer.data());
			lastATermUpdate = currentTime;
		}

		_bufferset->get_gridder(rowData.dataDescId)->grid_visibilities(timeIndex, metaData.antenna1, metaData.antenna2, rowData.uvw, rowData.data);
	}
	
	_bufferset->finished();
}

void IdgMsGridder::Predict(double* image)
{
	const size_t untrimmedWidth = ImageWidth();
	const size_t width = TrimWidth(), height = TrimHeight();

	assert(width == height);
	assert(untrimmedWidth == ImageHeight());

	_options["padded_size"] = untrimmedWidth;

	if (Polarization() == Polarization::StokesI)
	{
		_image.assign(4 * width * height, 0.0);
	}

	size_t polIndex = Polarization::StokesToIndex(Polarization());
	for(size_t i=0; i != width * height; ++i)
		_image[i + polIndex*width*height] = image[i];

	// Stokes V is the last requested pol, unless only Stokes I is imaged. Only when the last
	// polarization is requested, do the actual prediction.
	if (Polarization() == Polarization::StokesV || (Polarization() == Polarization::StokesI && _settings.polarizations.size()==1))
	{
		// Do actual predict
		std::vector<MSData> msDataVector;
		initializeMSDataVector(msDataVector, 4);

		double max_w = 0;
		for(size_t i=0; i!=MeasurementSetCount(); ++i)
		{
			max_w = std::max(max_w, msDataVector[i].maxW);
		}

		_bufferset->init(width, _actualPixelSizeX, max_w, _options);
		_bufferset->set_image(_image.data());

		for(size_t i=0; i!=MeasurementSetCount(); ++i)
		{
			predictMeasurementSet(msDataVector[i]);
		}
	}
}

void IdgMsGridder::setIdgType()
{
	switch(_settings.idgMode)
	{
		default:
			return;
		case WSCleanSettings::IDG_CPU:
			_proxyType = idg::api::Type::CPU_OPTIMIZED;
			return;
		case WSCleanSettings::IDG_GPU:
			_proxyType = idg::api::Type::CUDA_GENERIC;
			return;
		case WSCleanSettings::IDG_HYBRID:
			_proxyType = idg::api::Type::HYBRID_CUDA_CPU_OPTIMIZED;
			return;
	}
}


void IdgMsGridder::predictMeasurementSet(MSGridderBase::MSData& msData)
{
	const float max_baseline = msData.maxBaselineInM;

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

	_bufferset->init_buffers(_buffersize, bands, nr_stations, max_baseline, _options, idg::api::BufferSetType::degridding);

	std::unique_ptr<LofarBeamTerm> aTermMaker;
	ao::uvector<std::complex<float>> aTermBuffer;
	double aTermUpdateInterval = 900.0; // 15min
	if(_settings.gridWithBeam)
	{
		size_t subgridsize = _bufferset->get_subgridsize();
		double dl = _actualPixelSizeX * TrimWidth() / subgridsize, dm = _actualPixelSizeY * TrimHeight() / subgridsize;
		double pdl = PhaseCentreDL(), pdm = PhaseCentreDM();
		aTermMaker.reset(new LofarBeamTerm(ms, subgridsize, subgridsize, dl, dm, pdl, pdm, _settings.useDifferentialLofarBeam));
		aTermBuffer.resize(subgridsize*subgridsize*4*nr_stations);
	}
	
	ao::uvector<std::complex<float>> buffer(_selectedBands.MaxChannels()*4);
	int timeIndex = -1;
	double currentTime = -1.0, lastATermUpdate = -aTermUpdateInterval-1;
	for(msData.msProvider->Reset() ; msData.msProvider->CurrentRowAvailable() ; msData.msProvider->NextRow())
	{
		MSProvider::MetaData metaData;
		msData.msProvider->ReadMeta(metaData);
		size_t provRowId = msData.msProvider->RowId();
		if(currentTime != metaData.time)
		{
			currentTime = metaData.time;
			timeIndex++;
		}
		
		IDGPredictionRow row;
		row.uvw[0] = metaData.uInM;
		row.uvw[1] = -metaData.vInM;
		row.uvw[2] = -metaData.wInM;
		row.antenna1 = metaData.antenna1;
		row.antenna2 = metaData.antenna2;
		row.timeIndex = timeIndex;
		row.dataDescId = metaData.dataDescId;
		row.rowId = provRowId;
		predictRow(row);

		if(_settings.gridWithBeam && currentTime - lastATermUpdate > aTermUpdateInterval)
		{
			Logger::Debug << "Calculating a-terms for timestep " << timeIndex << "\n";
			aTermMaker->Calculate(aTermBuffer.data(), currentTime + aTermUpdateInterval*0.5, _selectedBands.CentreFrequency());
			_bufferset->get_degridder(row.dataDescId)->set_aterm(timeIndex, aTermBuffer.data());
			lastATermUpdate = currentTime;
		}
  }
	
	for(size_t d=0; d!=_selectedBands.DataDescCount(); ++d)
		computePredictionBuffer(d);
}

void IdgMsGridder::predictRow(IDGPredictionRow& row)
{
	bool isBufferFull = _bufferset->get_degridder(row.dataDescId)->request_visibilities(row.rowId, row.timeIndex, row.antenna1, row.antenna2, row.uvw);

	if(isBufferFull)
	{
		computePredictionBuffer(row.dataDescId);
	}
}

void IdgMsGridder::computePredictionBuffer(size_t dataDescId)
{
	auto available_row_ids = _bufferset->get_degridder(dataDescId)->compute();
	Logger::Debug << "Computed " << available_row_ids.size() << " rows.\n";
	for(auto i : available_row_ids)
	{
		_outputProvider->WriteModel(i.first, i.second);
	}
	_bufferset->get_degridder(dataDescId)->finished_reading();
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

void IdgMsGridder::readConfiguration()
{
	namespace po = boost::program_options; 
	po::options_description desc("Options"); 
	desc.add_options() 
	("proxy", "idg proxy")
	("max_nr_w_layers", po::value<int>(), "")
	("buffersize", po::value<int>(), ""); 

	po::variables_map vm;
	std::cout << "trying to open config file" << std::endl;
    std::ifstream ifs("idg.conf");
	if (ifs.fail())
	{
		std::cout << "could not open config file" << std::endl;
	}
	else
	{
		std::cout << "reading config file" << std::endl;
		try 
		{ 
			po::store(po::parse_config_file(ifs, desc), vm);
		}
		catch(po::error& e) 
		{ 
		} 
	}
	if (vm.count("proxy")) 
	{
		std::string proxy(vm["proxy"].as<string>());
		boost::to_lower(proxy);
		std::cout << "proxy = " << proxy << std::endl;
		if (proxy == "cpu-optimized") _proxyType = idg::api::Type::CPU_OPTIMIZED;
		if (proxy == "cpu-reference") _proxyType = idg::api::Type::CPU_REFERENCE;
		if (proxy == "cuda-generic") _proxyType = idg::api::Type::CUDA_GENERIC;
		if (proxy == "hybrid-cuda-cpu-optimized") _proxyType = idg::api::Type::HYBRID_CUDA_CPU_OPTIMIZED;
	}
	if (vm.count("buffersize")) 
	{
		_buffersize = vm["buffersize"].as<int>();
	}
	if (vm.count("max_nr_w_layers")) 
	{
		_options["max_nr_w_layers"] = vm["max_nr_w_layers"].as<int>();
	}
}
