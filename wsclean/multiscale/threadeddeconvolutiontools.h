#ifndef THREADED_DECONVOLUTION_TOOLS_H
#define THREADED_DECONVOLUTION_TOOLS_H

#include <vector>

#include "../lane.h"
#include "../uvector.h"

#include <boost/thread/thread.hpp>

class ThreadedDeconvolutionTools
{
public:
	ThreadedDeconvolutionTools(size_t threadCount);
	~ThreadedDeconvolutionTools();
	
	struct PeakData
	{
		double value, rms;
		size_t x, y;
	};
	
	void SubtractImage(double *image, const double *psf, size_t width, size_t height, size_t x, size_t y, double factor);
	
	// This one is for many transforms of the same scale
	void MultiScaleTransform(class MultiScaleTransforms* msTransforms, const ao::uvector<double*>& images, double* scratch, double scale);
	
	// This one is for transform of different scales
	void MultiScaleTransform(class MultiScaleTransforms* msTransforms, class ImageBufferAllocator* allocator, const ao::uvector<double*>& images, ao::uvector<double> scales);
	
	void FindMultiScalePeak(class MultiScaleTransforms* msTransforms, class ImageBufferAllocator* allocator, const double* image, const ao::uvector<double>& scales, std::vector<PeakData>& results, bool allowNegativeComponents, const bool* mask, double borderRatio, bool calculateRMS);
	
	static double RMS(const double* image, size_t n)
	{
		double result = 0.0;
		for(size_t i=0; i!=n; ++i)
			result += image[i] * image[i];
		return sqrt(result/double(n));
	}
	
private:
	struct ThreadResult {
	};
	struct FindMultiScalePeakResult : public ThreadResult {
		double value, rms;
		size_t x, y;
	};
	
	struct ThreadTask {
		virtual ThreadResult* operator()() = 0;
		virtual ~ThreadTask() { }
	};
	struct SubtractionTask : public ThreadTask {
		virtual ThreadResult* operator()();
		
		double *image;
		const double *psf;
		size_t width, height, x, y;
		double factor;
		size_t startY, endY;
	};
	struct FinishMultiScaleTransformTask : public ThreadTask {
		virtual ThreadResult* operator()();
		
		class MultiScaleTransforms* msTransforms;
		double* image;
		double* kernel;
	};
	struct MultiScaleTransformTask : public ThreadTask {
		virtual ThreadResult* operator()();
		
		class MultiScaleTransforms* msTransforms;
		double* image;
		double* scratch;
		double scale;
	};
	struct FindMultiScalePeakTask : public ThreadTask {
		virtual ThreadResult* operator()();
		
		class MultiScaleTransforms* msTransforms;
		double* image;
		double* scratch;
		double scale;
		bool allowNegativeComponents;
		const bool* mask;
		double borderRatio;
		bool calculateRMS;
	};
	
	std::vector<ao::lane<ThreadTask*>*> _taskLanes;
	std::vector<ao::lane<ThreadResult*>*> _resultLanes;
	size_t _threadCount;
	boost::thread_group _threadGroup;
	
	void threadFunc(ao::lane<ThreadTask*>* taskLane, ao::lane<ThreadResult*>* resultLane)
	{
		ThreadTask* task;
		while(taskLane->read(task))
		{
			ThreadResult* result = (*task)();
			resultLane->write(result);
			delete task;
		}
	}
};

#endif
