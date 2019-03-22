#include "griddingtaskmanager.h"

#include "wscleansettings.h"
#include "wsmsgridder.h"
#include "directmsgridder.h"

#include "../idg/idgmsgridder.h"

GriddingTaskManager::GriddingTaskManager(const class WSCleanSettings& settings, ImageBufferAllocator& allocator) :
	_taskList(settings.parallelGridding),
	_settings(settings),
	_allocator(allocator),
	_gridder(createGridder())
{
	prepareGridder(*_gridder);
}

GriddingTaskManager::~GriddingTaskManager()
{
	if(!_threadList.empty())
		Finish();
}

std::unique_ptr<MSGridderBase> GriddingTaskManager::createGridder() const
{
	if(_settings.useIDG)
		return std::unique_ptr<MSGridderBase>(new IdgMsGridder(_settings, _allocator));
	else if(_settings.directFT)
	{
		switch(_settings.directFTPrecision) {
		case DirectFTPrecision::Half:
			throw std::runtime_error("Half precision is not implemented");
			//return std::unique_ptr<MSGridderBase>(new DirectMSGridder<half_float::half>(&_imageAllocator, _settings.threadCount));
			break;
		case DirectFTPrecision::Float:
			return std::unique_ptr<MSGridderBase>(new DirectMSGridder<float>(&_allocator, _settings.threadCount));
			break;
		default:
		case DirectFTPrecision::Double:
			return std::unique_ptr<MSGridderBase>(new DirectMSGridder<double>(&_allocator, _settings.threadCount));
			break;
		case DirectFTPrecision::LongDouble:
			return std::unique_ptr<MSGridderBase>(new DirectMSGridder<long double>(&_allocator, _settings.threadCount));
			break;
		}
	}
	else
		return std::unique_ptr<MSGridderBase>(new WSMSGridder(&_allocator, _settings.threadCount, _settings.memFraction, _settings.absMemLimit));
}

void GriddingTaskManager::Run(GriddingTask& task, std::function<void (GriddingResult &)> finishCallback)
{
	if(_settings.useMPI)
	{
		// TODO
	}
	else if(_settings.parallelGridding == 1)
	{
		GriddingResult result = runDirect(task, *_gridder);
		finishCallback(result);
	}
	else {
		// Start an extra thread if not maxed out already
		if(_threadList.size() < _settings.parallelGridding)
			_threadList.emplace_back(&GriddingTaskManager::processQueue, this);
		else
			_taskList.wait_for_empty(); // if all threads are busy, block until one available (in order not to stack too many tasks)
		
		std::lock_guard<std::mutex> lock(_mutex);
		while(!_readyList.empty())
		{
			// Call callbacks for any finished tasks
			_readyList.back().second(_readyList.back().first);
			_readyList.pop_back();
		}
			
		_taskList.write(std::pair<GriddingTask, std::function<void (GriddingResult &)>>(std::move(task), finishCallback));
	}
}

void GriddingTaskManager::Finish()
{
	if(_settings.useMPI)
	{
		// TODO
	}
	else if(_settings.parallelGridding != 1)
	{
		_taskList.write_end();
		for(std::thread& t : _threadList)
			t.join();
		_threadList.clear();
		_taskList.clear();
		while(!_readyList.empty())
		{
			// Call callbacks for any finished tasks
			_readyList.back().second(_readyList.back().first);
			_readyList.pop_back();
		}
	}
}

GriddingResult GriddingTaskManager::runDirect(GriddingTask& task, MSGridderBase& gridder)
{
	gridder.ClearMeasurementSetList();
	for(auto& p : task.msList)
		gridder.AddMeasurementSet(p.first.get(), p.second);
	gridder.SetPolarization(task.polarization);
	gridder.SetIsComplex(task.polarization == Polarization::XY || task.polarization == Polarization::YX);
	gridder.SetVerbose(task.verbose);
	gridder.SetMetaDataCache(task.cache);
	gridder.SetPrecalculatedWeightInfo(task.precalculatedWeightInfo);
	if(task.operation == GriddingTask::Invert)
	{
		gridder.SetDoImagePSF(task.imagePSF);
		gridder.SetDoSubtractModel(task.subtractModel);
		gridder.SetStoreImagingWeights(task.storeImagingWeights);
		gridder.Invert();
	}
	else {
		gridder.SetAddToModel(task.addToModel);
		if(task.polarization == Polarization::XY || task.polarization == Polarization::YX)
			gridder.Predict(std::move(task.modelImageReal), std::move(task.modelImageImaginary));
		else
			gridder.Predict(std::move(task.modelImageReal));
	}
	
	GriddingResult result;
	result.imageRealResult = gridder.ImageRealResult();
	if(gridder.IsComplex())
		result.imageImaginaryResult = gridder.ImageImaginaryResult();
	result.phaseCentreRA = gridder.PhaseCentreRA();
	result.phaseCentreDec = gridder.PhaseCentreDec();
	result.beamSize = gridder.BeamSize();
	result.imageWeight = gridder.ImageWeight();
	result.normalizationFactor = gridder.NormalizationFactor();
	result.actualWGridSize = gridder.ActualWGridSize();
	result.griddedVisibilityCount = gridder.GriddedVisibilityCount();
	result.effectiveGriddedVisibilityCount = gridder.EffectiveGriddedVisibilityCount();
	result.visibilityWeightSum = gridder.VisibilityWeightSum();
	result.actualInversionWidth = gridder.ActualInversionWidth();
	result.actualInversionHeight = gridder.ActualInversionHeight();
	result.hasGriddingCorrectionImage = gridder.HasGriddingCorrectionImage();
	return result;
}

void GriddingTaskManager::processQueue()
{
	std::unique_ptr<MSGridderBase> gridder(createGridder());
	prepareGridder(*gridder);
	gridder->SetCasacoreMutex(_casacoreMutex);
	
	std::pair<GriddingTask, std::function<void (GriddingResult &)>> taskPair;
	while(_taskList.read(taskPair))
	{
		GriddingResult result = runDirect(taskPair.first, *gridder);
		
		std::lock_guard<std::mutex> lock(_mutex);
		_readyList.emplace_back(std::move(result), taskPair.second);
	}
}

void GriddingTaskManager::prepareGridder(MSGridderBase& gridder)
{
	gridder.SetGridMode(_settings.gridMode);
	gridder.SetImageWidth(_settings.paddedImageWidth);
	gridder.SetImageHeight(_settings.paddedImageHeight);
	gridder.SetTrimSize(_settings.trimmedImageWidth, _settings.trimmedImageHeight);
	gridder.SetNWSize(_settings.widthForNWCalculation, _settings.heightForNWCalculation);
	gridder.SetNWFactor(_settings.nWLayersFactor);
	gridder.SetPixelSizeX(_settings.pixelScaleX);
	gridder.SetPixelSizeY(_settings.pixelScaleY);
	if(_settings.nWLayers != 0)
		gridder.SetWGridSize(_settings.nWLayers);
	else
		gridder.SetNoWGridSize();
	gridder.SetAntialiasingKernelSize(_settings.antialiasingKernelSize);
	gridder.SetOverSamplingFactor(_settings.overSamplingFactor);
	gridder.SetDataColumnName(_settings.dataColumnName);
	gridder.SetWeighting(_settings.weightMode);
	gridder.SetWLimit(_settings.wLimit/100.0);
	gridder.SetSmallInversion(_settings.smallInversion);
	gridder.SetVisibilityWeightingMode(_settings.visibilityWeightingMode);
}

