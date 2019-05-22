#ifndef GRIDDING_OPERATION_H
#define GRIDDING_OPERATION_H

#include <cstring>
#include <condition_variable>
#include <functional>
#include <thread>
#include <vector>

#include "msgridderbase.h"
#include "measurementsetgridder.h"
#include "imagebufferallocator.h"

#include "../lane.h"
#include "../imageweights.h"
#include "../polarization.h"
#include "../msselection.h"

#include "../msproviders/msprovider.h"

struct GriddingResult
{
	ImageBufferAllocator::Ptr imageRealResult;
	ImageBufferAllocator::Ptr imageImaginaryResult;
	double phaseCentreRA, phaseCentreDec;
	double beamSize;
	double imageWeight;
	double normalizationFactor;
	size_t actualWGridSize;
	size_t griddedVisibilityCount;
	double effectiveGriddedVisibilityCount;
	double visibilityWeightSum;
	size_t actualInversionWidth, actualInversionHeight;
	bool hasGriddingCorrectionImage;
};

class GriddingTask
{
public:
	enum { Invert, Predict } operation;
	bool imagePSF;
	bool subtractModel;
	PolarizationEnum polarization;
	bool verbose;
	MSGridderBase::MetaDataCache* cache;
	bool storeImagingWeights;
	
	std::shared_ptr<ImageWeights> precalculatedWeightInfo;
	std::vector<std::pair<std::unique_ptr<MSProvider>, MSSelection>> msList;
	
	// For prediction
	bool addToModel;
	ImageBufferAllocator::Ptr modelImageReal;
	ImageBufferAllocator::Ptr modelImageImaginary;
};

class GriddingTaskManager
{
public:
	GriddingTaskManager(const class WSCleanSettings& settings, class ImageBufferAllocator& allocator);
	~GriddingTaskManager();
	
	void Run(GriddingTask& task, std::function<void(GriddingResult&)> finishCallback);
	
	void Finish();
	
	MSGridderBase* Gridder() { return _gridder.get(); }
	
private:
	GriddingResult runDirect(GriddingTask& task, MSGridderBase& gridder);
	
	std::unique_ptr<MSGridderBase> createGridder() const;
	void prepareGridder (MSGridderBase& gridder);
	void processQueue();
	
	std::mutex _mutex, _casacoreMutex;
	std::vector<std::thread> _threadList;
	ao::lane<std::pair<GriddingTask, std::function<void(GriddingResult&)>>> _taskList;
	std::vector<std::pair<GriddingResult, std::function<void(GriddingResult&)>>> _readyList;
	
	const class WSCleanSettings& _settings;
	ImageBufferAllocator& _allocator;
	std::unique_ptr<MSGridderBase> _gridder;
};

#endif
