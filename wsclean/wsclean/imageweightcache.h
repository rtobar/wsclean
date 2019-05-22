#ifndef IMAGE_WEIGHT_CACHE_H
#define IMAGE_WEIGHT_CACHE_H

#include "measurementsetgridder.h"
#include "logger.h"

#include "../imageweights.h"
#include "../weightmode.h"

#include <limits>
#include <mutex>

class ImageWeightCache
{
public:
	ImageWeightCache(const WeightMode& weightMode, size_t imageWidth, size_t imageHeight, double pixelScaleX, double pixelScaleY, double minUVInLambda, double maxUVInLambda, double rankFilterLevel, size_t rankFilterSize, bool weightsAsTaper) :
		_weightMode(weightMode),
		_imageWidth(imageWidth),
		_imageHeight(imageHeight),
		_pixelScaleX(pixelScaleX),
		_pixelScaleY(pixelScaleY),
		_minUVInLambda(minUVInLambda),
		_maxUVInLambda(maxUVInLambda),
		_rankFilterLevel(rankFilterLevel),
		_rankFilterSize(rankFilterSize),
		_gaussianTaperBeamSize(0),
		_tukeyTaperInLambda(0), _tukeyInnerTaperInLambda(0),
		_edgeTaperInLambda(0),
		_edgeTukeyTaperInLambda(0),
		_weightsAsTaper(weightsAsTaper),
		_currentWeightChannel(std::numeric_limits<size_t>::max()),
		_currentWeightInterval(std::numeric_limits<size_t>::max())
	{
	}
	
	void SetTaperInfo(double gaussianTaperBeamSize, double tukeyTaperInLambda, double tukeyInnerTaperInLambda, double edgeTaperInLambda, double edgeTukeyTaperInLambda)
	{
		_gaussianTaperBeamSize = gaussianTaperBeamSize;
		_tukeyTaperInLambda = tukeyTaperInLambda;
		_tukeyInnerTaperInLambda = tukeyInnerTaperInLambda;
		_edgeTaperInLambda = edgeTaperInLambda;
		_edgeTukeyTaperInLambda = edgeTukeyTaperInLambda;
	}
	
	std::shared_ptr<ImageWeights> Get(const std::vector<std::pair<std::unique_ptr<MSProvider>, MSSelection>>& msList, size_t outChannelIndex, size_t outIntervalIndex)
	{
		std::unique_lock<std::mutex> lock(_mutex);
		if(outChannelIndex != _currentWeightChannel || outIntervalIndex != _currentWeightInterval)
		{
			lock.unlock();
			std::unique_ptr<ImageWeights> weights = recalculateWeights(msList);
			lock.lock();
			_currentWeightChannel = outChannelIndex;
			_currentWeightInterval = outIntervalIndex;
			_cachedWeights = std::move(weights);
		}
		return _cachedWeights;
	}
	
	std::shared_ptr<ImageWeights> Get(casacore::MeasurementSet& ms, MSSelection& selection)
	{
		std::unique_ptr<ImageWeights> weights = MakeEmptyWeights();
		weights->Grid(ms, selection);
		initializeWeightTapers(*weights);
		return std::move(weights);
	}
	
	std::unique_ptr<ImageWeights> MakeEmptyWeights() const
	{
		return std::unique_ptr<ImageWeights>(new ImageWeights(_weightMode, _imageWidth, _imageHeight, _pixelScaleX, _pixelScaleY, _weightsAsTaper, _weightMode.SuperWeight()));
	};
	
	std::shared_ptr<ImageWeights> GetMFWeights() const
	{
		return _cachedWeights;
	}
	
	void SetMFWeights(std::unique_ptr<ImageWeights> weights)
	{
		initializeWeightTapers(*weights);
		std::unique_lock<std::mutex> lock(_mutex);
		_cachedWeights = std::move(weights);
		_currentWeightChannel = std::numeric_limits<size_t>::max();
		_currentWeightInterval = std::numeric_limits<size_t>::max();
	}
	
private:
	std::unique_ptr<ImageWeights> recalculateWeights(const std::vector<std::pair<std::unique_ptr<MSProvider>, MSSelection>>& msList)
	{
		Logger::Info << "Precalculating weights for " << _weightMode.ToString() << " weighting... ";
		Logger::Info.Flush();
		std::unique_ptr<ImageWeights> weights = MakeEmptyWeights();
		for(size_t i=0; i!=msList.size(); ++i)
		{
			weights->Grid(*msList[i].first, msList[i].second);
			if(msList.size() > 1)
				(Logger::Info << i << ' ').Flush();
		}
		weights->FinishGridding();
		initializeWeightTapers(*weights);
		Logger::Info << "DONE\n";
		return weights;
	}
	
	void initializeWeightTapers(ImageWeights& weights)
	{
		if(_rankFilterLevel >= 1.0)
			weights.RankFilter(_rankFilterLevel, _rankFilterSize);
		
		if(_gaussianTaperBeamSize != 0.0)
			weights.SetGaussianTaper(_gaussianTaperBeamSize);
		
		if(_tukeyInnerTaperInLambda != 0.0)
			weights.SetTukeyInnerTaper(_tukeyInnerTaperInLambda, _minUVInLambda);
		else if(_minUVInLambda!=0.0)
			weights.SetMinUVRange(_minUVInLambda);
		
		if(_tukeyTaperInLambda != 0.0)
			weights.SetTukeyTaper(_tukeyTaperInLambda, _maxUVInLambda);
		else if(_maxUVInLambda!=0.0)
			weights.SetMaxUVRange(_maxUVInLambda);
		
		if(_edgeTukeyTaperInLambda != 0.0)
			weights.SetEdgeTukeyTaper(_edgeTukeyTaperInLambda, _edgeTaperInLambda);
		else if(_edgeTaperInLambda != 0.0)
			weights.SetEdgeTaper(_edgeTaperInLambda);
	}

	std::shared_ptr<ImageWeights> _cachedWeights;
	const WeightMode _weightMode;
	size_t _imageWidth, _imageHeight;
	double _pixelScaleX, _pixelScaleY;
	double _minUVInLambda, _maxUVInLambda;
	double _rankFilterLevel;
	size_t _rankFilterSize;
	double _gaussianTaperBeamSize;
	double _tukeyTaperInLambda, _tukeyInnerTaperInLambda;
	double _edgeTaperInLambda;
	double _edgeTukeyTaperInLambda;
	bool _weightsAsTaper;
	std::mutex _mutex;
	
	size_t _currentWeightChannel, _currentWeightInterval;
};

#endif
