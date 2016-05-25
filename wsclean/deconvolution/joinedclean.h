#ifndef JOINED_CLEAN_H
#define JOINED_CLEAN_H

#include "deconvolutionalgorithm.h"
#include "imageset.h"
#include "simpleclean.h"

namespace ao {
	template<typename T> class lane;
}

template<typename ImageSetType>
class JoinedClean : public TypedDeconvolutionAlgorithm<ImageSetType>
{
public:
	virtual void ExecuteMajorIteration(ImageSetType& dataImage, ImageSetType& modelImage, const ao::uvector<const double*>& psfImages, size_t width, size_t height, bool& reachedStopGain)  ;
	
private:
	size_t _width, _height;
	ao::uvector<double> _curPeakValues;
	
	struct CleanTask
	{
		size_t cleanCompX, cleanCompY;
	};
	struct CleanResult
	{
		CleanResult() : nextPeakX(0), nextPeakY(0), peakLevelUnnormalized(0.0)
		{ }
		size_t nextPeakX, nextPeakY;
		double peakLevelUnnormalized;
	};
	struct CleanThreadData
	{
		size_t startY, endY;
		ImageSetType* dataImage;
		ao::uvector<const double*> psfImages;
	};

	void findPeak(const ImageSetType& image, size_t& x, size_t& y) const
	{
		if(this->_cleanMask == 0)
			findPeak(image, x, y, 0, _height);
		else
			findPeak(image, x, y, 0, _height, this->_cleanMask);
	}
	void findPeak(const ImageSetType& image, size_t& x, size_t& y, size_t startY, size_t stopY) const;
	void findPeak(const ImageSetType& image, size_t& x, size_t& y, size_t startY, size_t stopY, const bool* mask) const;
	
	std::string peakDescription(const ImageSetType& image, size_t& x, size_t& y);
	
	void cleanThreadFunc(ao::lane<CleanTask>* taskLane, ao::lane<CleanResult>* resultLane, CleanThreadData cleanData);
	
	void subtractImage(double *image, const double *psf, size_t x, size_t y, double factor, size_t startY, size_t endY) const
	{
		SimpleClean::PartialSubtractImage(image, _width, _height, psf, _width, _height, x, y, factor, startY, endY);
	}
};

#endif
