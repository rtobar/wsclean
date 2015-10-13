#ifndef SIMPLE_CLEAN_H
#define SIMPLE_CLEAN_H

#include <string>
#include <cmath>
#include <limits>

#include "deconvolutionalgorithm.h"
#include "imageset.h"

//#define FORCE_NON_AVX 1

namespace ao {
	template<typename T> class lane;
}

class SimpleClean : public TypedDeconvolutionAlgorithm<deconvolution::SingleImageSet>
{
	public:
		static double FindPeakSimple(const double *image, size_t width, size_t height, size_t &x, size_t &y, bool allowNegativeComponents, size_t startY, size_t endY, double borderRatio);
		
		static double FindPeak(const double *image, size_t width, size_t height, size_t &x, size_t &y, bool allowNegativeComponents, const bool* cleanMask);

#if defined __AVX__ && !defined FORCE_NON_AVX
		template<bool AllowNegativeComponent>
		static double FindPeakAVX(const double *image, size_t width, size_t height, size_t &x, size_t &y, size_t startY, size_t endY, double borderRatio);
		
		static double FindPeakAVX(const double *image, size_t width, size_t height, size_t &x, size_t &y, bool allowNegativeComponents, size_t startY, size_t endY, double borderRatio)
		{
			if(allowNegativeComponents)
				return FindPeakAVX<true>(image, width, height, x, y, startY, endY, borderRatio);
			else
				return FindPeakAVX<false>(image, width, height, x, y, startY, endY, borderRatio);
		}
		static double FindPeak(const double *image, size_t width, size_t height, size_t &x, size_t &y, bool allowNegativeComponents, size_t startY, size_t endY, double borderRatio)
		{
			return FindPeakAVX(image, width, height, x, y, allowNegativeComponents, startY, endY, borderRatio);
		}
#else
		static double FindPeak(const double *image, size_t width, size_t height, size_t &x, size_t &y, bool allowNegativeComponents, size_t startY, size_t endY, double borderRatio)
		{
			return FindPeakSimple(image, width, height, x, y, allowNegativeComponents, startY, endY, borderRatio);
		}
#endif

		static double FindPeak(const double *image, size_t width, size_t height, size_t &x, size_t &y, bool allowNegativeComponents, size_t startY, size_t endY, const bool* cleanMask, double borderRatio);
		
		static void SubtractImage(double *image, const double *psf, size_t width, size_t height, size_t x, size_t y, double factor);
		
		static void PartialSubtractImage(double *image, const double *psf, size_t width, size_t height, size_t x, size_t y, double factor, size_t startY, size_t endY);
		
		static void PartialSubtractImage(double *image, size_t imgWidth, size_t imgHeight, const double *psf, size_t psfWidth, size_t psfHeight, size_t x, size_t y, double factor, size_t startY, size_t endY);
		
#ifdef __AVX__
		static void PartialSubtractImageAVX(double *image, size_t imgWidth, size_t imgHeight, const double *psf, size_t psfWidth, size_t psfHeight, size_t x, size_t y, double factor, size_t startY, size_t endY);
#endif
		
		/**
		 * Single threaded implementation -- just for reference.
		 */
		void ExecuteMajorIterationST(double *dataImage, double *modelImage, const double *psfImage, size_t width, size_t height);
		
    virtual void ExecuteMajorIteration(ImageSet& dataImage, ImageSet& modelImage, const ao::uvector<const double*>& psfImages, size_t width, size_t height, bool& reachedStopGain)
		{
			ExecuteMajorIteration(dataImage.GetImage(0), modelImage.GetImage(0), psfImages[0], width, height, reachedStopGain);
		}
		
		void ExecuteMajorIteration(double* dataImage, double* modelImage, const double* psfImage, size_t width, size_t height, bool& reachedStopGain);
	private:
		struct CleanTask
		{
			size_t cleanCompX, cleanCompY;
			double peakLevel;
		};
		struct CleanResult
		{
			CleanResult() : nextPeakX(0), nextPeakY(0), peakLevel(0.0)
			{ }
			size_t nextPeakX, nextPeakY;
			double peakLevel;
		};
		struct CleanThreadData
		{
			size_t startY, endY;
			double *dataImage;
			size_t imgWidth, imgHeight;
			const double *psfImage;
			size_t psfWidth, psfHeight;
		};
		void cleanThreadFunc(ao::lane<CleanTask>* taskLane, ao::lane<CleanResult>* resultLane, CleanThreadData cleanData);
};

#endif
