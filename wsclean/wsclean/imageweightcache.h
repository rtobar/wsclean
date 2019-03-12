#ifndef IMAGE_WEIGHT_CACHE_H
#define IMAGE_WEIGHT_CACHE_H

#include "measurementsetgridder.h"
#include "logger.h"

#include "../imageweights.h"
#include "../weightmode.h"

#include <limits>

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
	
	void Update(const std::vector<std::pair<std::unique_ptr<MSProvider>, MSSelection>>& msList, size_t outChannelIndex, size_t outIntervalIndex)
	{
		if(outChannelIndex != _currentWeightChannel || outIntervalIndex != _currentWeightInterval)
		{
			_currentWeightChannel = outChannelIndex;
			_currentWeightInterval = outIntervalIndex;
			
			recalculateWeights(msList);
		}
	}
	
	void ResetWeights()
	{
		_imageWeights.reset(new ImageWeights(_weightMode, _imageWidth, _imageHeight, _pixelScaleX, _pixelScaleY, _weightsAsTaper, _weightMode.SuperWeight()));
	};
	
	ImageWeights& Weights()
	{
		return *_imageWeights;
	}
	
	const ImageWeights& Weights() const
	{
		return *_imageWeights;
	}
	
	void InitializeWeightTapers()
	{
		if(_rankFilterLevel >= 1.0)
			_imageWeights->RankFilter(_rankFilterLevel, _rankFilterSize);
		
		if(_gaussianTaperBeamSize != 0.0)
			_imageWeights->SetGaussianTaper(_gaussianTaperBeamSize);
		
		if(_tukeyInnerTaperInLambda != 0.0)
			_imageWeights->SetTukeyInnerTaper(_tukeyInnerTaperInLambda, _minUVInLambda);
		else if(_minUVInLambda!=0.0)
			_imageWeights->SetMinUVRange(_minUVInLambda);
		
		if(_tukeyTaperInLambda != 0.0)
			_imageWeights->SetTukeyTaper(_tukeyTaperInLambda, _maxUVInLambda);
		else if(_maxUVInLambda!=0.0)
			_imageWeights->SetMaxUVRange(_maxUVInLambda);
		
		if(_edgeTukeyTaperInLambda != 0.0)
			_imageWeights->SetEdgeTukeyTaper(_edgeTukeyTaperInLambda, _edgeTaperInLambda);
		else if(_edgeTaperInLambda != 0.0)
			_imageWeights->SetEdgeTaper(_edgeTaperInLambda);
	}

private:
	void recalculateWeights(const std::vector<std::pair<std::unique_ptr<MSProvider>, MSSelection>>& msList)
	{
		Logger::Info << "Precalculating weights for " << _weightMode.ToString() << " weighting... ";
		Logger::Info.Flush();
		ResetWeights();
		for(size_t i=0; i!=msList.size(); ++i)
		{
			_imageWeights->Grid(*msList[i].first, msList[i].second);
			if(msList.size() > 1)
				(Logger::Info << i << ' ').Flush();
		}
		_imageWeights->FinishGridding();
		InitializeWeightTapers();
		Logger::Info << "DONE\n";
	}
	
	std::unique_ptr<ImageWeights> _imageWeights;
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
	
	size_t _currentWeightChannel, _currentWeightInterval;
};

#endif
