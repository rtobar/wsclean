// kate: space-indent off; tab-width 2; indent-width 2; replace-tabs off; eol unix;

#include "idgmsgridder.h"

#include <cmath>
#include <thread>

#include <idg-api.h>

#include "../msproviders/msprovider.h"

#include "../wsclean/imagefilename.h"
#include "../wsclean/imagingtable.h"
#include "../wsclean/logger.h"
#include "../wsclean/wscleansettings.h"

#include "../fitsreader.h"
#include "../fitswriter.h"

#include "../aterms/atermconfig.h"
#include "../aterms/lofarbeamterm.h"

#include "idgconfiguration.h"

IdgMsGridder::IdgMsGridder(const WSCleanSettings& settings) :
	_averageBeam(nullptr),
	_outputProvider(nullptr),
	_settings(settings),
	_proxyType(idg::api::Type::CPU_OPTIMIZED),
	_buffersize(256)
{
	IdgConfiguration::Read(_proxyType, _buffersize, _options);
	setIdgType();
	_bufferset = std::unique_ptr<idg::api::BufferSet>(
		idg::api::BufferSet::create(_proxyType));
	if(_settings.gridWithBeam) _options["a_term_kernel_size"] = float(5.0);
}

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
		if (!_metaDataCache->average_beam) _metaDataCache->average_beam.reset(new AverageBeam());
		_averageBeam = static_cast<AverageBeam*>(_metaDataCache->average_beam.get());

		std::vector<MSData> msDataVector;
		initializeMSDataVector(msDataVector);
		

		double max_w = 0;
		for(size_t i=0; i!=MeasurementSetCount(); ++i)
		{
			max_w = std::max(max_w, msDataVector[i].maxWWithFlags);
		}

		_bufferset->init(width, _actualPixelSizeX, max_w+1.0, _options);

		if (DoImagePSF())
		{
			// Computing the PSF
			// For the PSF the aterm is computed but not applied
			// The aterm is computed so that the average beam can be computed
			_bufferset->set_apply_aterm(false);
			_bufferset->unset_matrix_inverse_beam();
			_bufferset->init_compute_avg_beam(idg::api::compute_flags::compute_and_grid);
			resetVisibilityCounters();
			for(size_t i=0; i!=MeasurementSetCount(); ++i)
			{
				// Adds the gridding result to _image member
				gridMeasurementSet(msDataVector[i]);
			}
			_bufferset->finalize_compute_avg_beam();
			_averageBeam->SetScalarBeam(_bufferset->get_scalar_beam());
			_averageBeam->SetMatrixInverseBeam(_bufferset->get_matrix_inverse_beam());
			_image.assign(4 * width * height, 0.0);
			_bufferset->get_image(_image.data());

			// Normalize by total weight
			Logger::Debug << "Total weight: " << totalWeight() << '\n';

			double center_pixel_value = _image[height/2 * width + width/2]; // TODO check memory layout, is this correct? for now it does not matter, because width == height

			if (center_pixel_value)
			{
				for(size_t ii=0; ii != 4 * width * height; ++ii)
				{
					_image[ii] /= center_pixel_value;
				}
			}

		}
		else {
			// Compute a dirty/residual image
			// with application of the a term
			_bufferset->set_apply_aterm(true);

			// Because compensation for the average beam happens at subgrid level
			// it needs to be known in advance.
			// If it is not in the cache it needs to be computed first
			if (!_averageBeam->Empty())
			{
				// Set avg beam from cache
				Logger::Debug << "Using average beam from cache.\n";
				_bufferset->set_scalar_beam(_averageBeam->ScalarBeam());
				_bufferset->set_matrix_inverse_beam(_averageBeam->MatrixInverseBeam());
			}
			else {
				// Compute avg beam
				Logger::Debug << "Computing average beam.\n";
				_bufferset->init_compute_avg_beam(idg::api::compute_flags::compute_only);
				for(size_t i=0; i!=MeasurementSetCount(); ++i)
				{
					gridMeasurementSet(msDataVector[i]);
				}
				_bufferset->finalize_compute_avg_beam();
				Logger::Debug << "Finished computing average beam.\n";
				_averageBeam->SetScalarBeam(_bufferset->get_scalar_beam());
				_averageBeam->SetMatrixInverseBeam(_bufferset->get_matrix_inverse_beam());
			}

			resetVisibilityCounters();
			for(size_t i=0; i!=MeasurementSetCount(); ++i)
			{
				// Adds the gridding result to _image member
				gridMeasurementSet(msDataVector[i]);
			}
			_image.assign(4 * width * height, 0.0);
			_bufferset->get_image(_image.data());
		}

		// result is now in _image member
		// Can be accessed by subsequent calls to ImageRealResult()
		
	}
	else if(_image.empty()) {
		throw std::runtime_error("IdgMsGridder::Invert() was called out of sequence");
	}
}

std::unique_ptr<class ATermBase> IdgMsGridder::getATermMaker(MSGridderBase::MSData& msData)
{
	casacore::MeasurementSet& ms = msData.msProvider->MS();
	size_t nr_stations = ms.antenna().nrow();
	std::unique_ptr<ATermBase> aTermMaker;
	ao::uvector<std::complex<float>> aTermBuffer;
	if(!_settings.atermConfigFilename.empty() || _settings.gridWithBeam)
	{
		// IDG uses a flipped coordinate system which is moved by half a pixel:
		double dl = -_bufferset->get_subgrid_pixelsize();
		double dm = -_bufferset->get_subgrid_pixelsize();
		double pdl = PhaseCentreDL() - 0.5*dl, pdm = PhaseCentreDM() + 0.5*dm;
		size_t subgridsize = _bufferset->get_subgridsize();
		if(!_settings.atermConfigFilename.empty())
		{
			std::unique_ptr<ATermConfig> config(new ATermConfig(ms, nr_stations, subgridsize, subgridsize, dl, dm, pdl, pdm));
			config->Read(_settings.atermConfigFilename);
			return std::move(config);
		}
		else
			return std::unique_ptr<LofarBeamTerm>(new LofarBeamTerm(ms, subgridsize, subgridsize, dl, dm, pdl, pdm, _settings.beamAtermUpdateTime, _settings.useDifferentialLofarBeam));
	}
	else {
		return std::unique_ptr<ATermBase>();
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
	
	std::unique_ptr<ATermBase> aTermMaker = getATermMaker(msData);
	double aTermUpdateInterval = _settings.beamAtermUpdateTime; // in seconds
	ao::uvector<std::complex<float>> aTermBuffer;
	size_t subgridsize = _bufferset->get_subgridsize();
	if(aTermMaker)
		aTermBuffer.resize(subgridsize*subgridsize*4*nr_stations);

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
			
			if(aTermMaker && (!_settings.gridWithBeam || currentTime - lastATermUpdate > aTermUpdateInterval))
			{
				if(aTermMaker->Calculate(aTermBuffer.data(), currentTime + aTermUpdateInterval*0.5, _selectedBands.CentreFrequency()))
				{
					_bufferset->get_gridder(metaData.dataDescId)->set_aterm(timeIndex, aTermBuffer.data());
					Logger::Debug << "Calculated a-terms for timestep " << timeIndex << "\n";
				}
				lastATermUpdate = currentTime;
			}
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

		_bufferset->get_gridder(rowData.dataDescId)->grid_visibilities(timeIndex, metaData.antenna1, metaData.antenna2, rowData.uvw, rowData.data, weightBuffer.data());
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
		if (!_metaDataCache->average_beam)
		{
			Logger::Info << "no average_beam in cache, creating an empty one.\n";
			_metaDataCache->average_beam.reset(new AverageBeam());
		}
		_averageBeam = static_cast<AverageBeam*>(_metaDataCache->average_beam.get());
	}

	size_t polIndex = Polarization::StokesToIndex(Polarization());
	for(size_t i=0; i != width * height; ++i)
		_image[i + polIndex*width*height] = image[i];

	// Stokes V is the last requested pol, unless only Stokes I is imaged. Only when the last
	// polarization is given, do the actual prediction.
	if (Polarization() == Polarization::StokesV || (Polarization() == Polarization::StokesI && _settings.polarizations.size()==1))
	{
		// Do actual predict

		bool do_scale = false;
		if (!_averageBeam->Empty())
		{
			// Set avg beam from cache
			Logger::Debug << "Average beam is already in cache.\n";
			_bufferset->set_scalar_beam(_averageBeam->ScalarBeam());
			do_scale = true;
		}

		std::vector<MSData> msDataVector;
		initializeMSDataVector(msDataVector);

		double max_w = 0;
		for(size_t i=0; i!=MeasurementSetCount(); ++i)
		{
			max_w = std::max(max_w, msDataVector[i].maxWWithFlags);
		}

		_bufferset->init(width, _actualPixelSizeX, max_w+1.0, _options);
		_bufferset->set_image(_image.data(), do_scale);

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

	double aTermUpdateInterval = _settings.beamAtermUpdateTime; // in seconds
	std::unique_ptr<ATermBase> aTermMaker = getATermMaker(msData);
	ao::uvector<std::complex<float>> aTermBuffer;
	size_t subgridsize = _bufferset->get_subgridsize();
	if(aTermMaker)
		aTermBuffer.resize(subgridsize*subgridsize*4*nr_stations);
	
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
			
			if(aTermMaker && (!_settings.gridWithBeam || currentTime - lastATermUpdate > aTermUpdateInterval))
			{
				if(aTermMaker->Calculate(aTermBuffer.data(), currentTime + aTermUpdateInterval*0.5, _selectedBands.CentreFrequency()))
				{
					_bufferset->get_degridder(metaData.dataDescId)->set_aterm(timeIndex, aTermBuffer.data());
					Logger::Debug << "Calculated new a-terms for timestep " << timeIndex << "\n";
				}
				lastATermUpdate = currentTime;
			}
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
  }
	
	for(size_t d=0; d!=_selectedBands.DataDescCount(); ++d)
		computePredictionBuffer(d);
}

void IdgMsGridder::predictRow(IDGPredictionRow& row)
{
	while (_bufferset->get_degridder(row.dataDescId)->request_visibilities(row.rowId, row.timeIndex, row.antenna1, row.antenna2, row.uvw))
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

void IdgMsGridder::Predict(double* /*real*/, double* /*imaginary*/)
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
	return false;
}

void IdgMsGridder::SaveBeamImage(const ImagingTableEntry& entry, ImageFilename& filename) const
{
	if (!_averageBeam || _averageBeam->Empty())
	{
		throw std::runtime_error("IDG gridder can not save beam image. Beam has not been computed yet.");
	}
	FitsWriter writer;
	writer.SetImageDimensions(_settings.trimmedImageWidth, _settings.trimmedImageHeight, PhaseCentreRA(), PhaseCentreDec(), _settings.pixelScaleX, _settings.pixelScaleY);
	writer.SetPhaseCentreShift(PhaseCentreDL(), PhaseCentreDM());
	ImageFilename polName(filename);
	polName.SetPolarization(Polarization::StokesI);
	writer.SetPolarization(Polarization::StokesI);
	writer.SetFrequency(entry.CentralFrequency(), entry.bandEndFrequency - entry.bandStartFrequency);
	writer.Write(polName.GetBeamPrefix(_settings) + ".fits", _averageBeam->ScalarBeam()->data());
}

void IdgMsGridder::SavePBCorrectedImages(FitsWriter& writer, ImageFilename& filename, const std::string& filenameKind, ImageBufferAllocator& allocator) const
{
	ImageBufferAllocator::Ptr image;
	for(size_t polIndex = 0; polIndex != 4; ++polIndex)
	{
		PolarizationEnum pol = Polarization::IndexToStokes(polIndex);
		ImageFilename name(filename);
		name.SetPolarization(pol);
		FitsReader reader(name.GetPrefix(_settings) + "-" + filenameKind + ".fits");
		if(!image)
			allocator.Allocate(reader.ImageWidth() * reader.ImageHeight(), image);
		reader.Read(image.data());
		
		const float* beam = _averageBeam->ScalarBeam()->data();
		for(size_t i=0; i!=reader.ImageWidth() * reader.ImageHeight(); ++i)
		{
			image[i] /= beam[i];
		}
		
		writer.SetPolarization(pol);
		writer.Write(name.GetPrefix(_settings) + "-" + filenameKind + "-pb.fits", image.data());
	}
}
