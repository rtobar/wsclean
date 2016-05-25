#include "threadeddeconvolutiontools.h"
#include "multiscaletransforms.h"
#include "../deconvolution/simpleclean.h"

ThreadedDeconvolutionTools::ThreadedDeconvolutionTools(size_t threadCount) :
	_taskLanes(threadCount),
	_resultLanes(threadCount),
	_threadCount(threadCount)
{
	for(size_t i=0; i!=_threadCount; ++i)
	{
		_taskLanes[i] = new ao::lane<ThreadTask*>(1);
		_resultLanes[i] = new ao::lane<ThreadResult*>(1);
		_threadGroup.add_thread(new boost::thread(&ThreadedDeconvolutionTools::threadFunc, this, _taskLanes[i], _resultLanes[i]));
	}
}

ThreadedDeconvolutionTools::~ThreadedDeconvolutionTools()
{
	for(size_t i=0; i!=_threadCount; ++i)
	{
		_taskLanes[i]->write_end();
	}
	
	_threadGroup.join_all();
	
	for(size_t i=0; i!=_threadCount; ++i)
	{
		delete _taskLanes[i];
		delete _resultLanes[i];
	}
}

void ThreadedDeconvolutionTools::SubtractImage(double* image, const double* psf, size_t width, size_t height, size_t x, size_t y, double factor)
{
	for(size_t thr=0; thr!=_threadCount; ++thr)
	{
		SubtractionTask* task = new SubtractionTask();
		task->image = image;
		task->psf = psf;
		task->width = width;
		task->height = height;
		task->x = x;
		task->y = y;
		task->factor = factor;
		task->startY = height*thr/_threadCount;
		task->endY = height*(thr+1)/_threadCount;
		if(thr == _threadCount-1)
		{
			ThreadResult* result = (*task)();
			delete result;
		}
		else {
			_taskLanes[thr]->write(task);
		}
	}
	for(size_t thr=0; thr!=_threadCount-1; ++thr)
	{
		ThreadResult* result = 0;
		_resultLanes[thr]->read(result);
		delete result;
	}
}

ThreadedDeconvolutionTools::ThreadResult* ThreadedDeconvolutionTools::SubtractionTask::operator()()
{
	SimpleClean::PartialSubtractImage(image, psf, width, height, x, y, factor, startY, endY);
	return 0;
}

void ThreadedDeconvolutionTools::MultiScaleTransform(MultiScaleTransforms* msTransforms, const ao::uvector<double*>& images, double* scratch, double scale)
{
	size_t imageIndex = 0;
	size_t nextThread = 0;
	msTransforms->PrepareTransform(scratch, scale);
	while(imageIndex < images.size())
	{
		FinishMultiScaleTransformTask* task = new FinishMultiScaleTransformTask();
		task->msTransforms = msTransforms;
		task->image = images[imageIndex];
		task->kernel = scratch;
		_taskLanes[nextThread]->write(task);
		
		++nextThread;
		if(nextThread == _threadCount)
		{
			for(size_t thr=0; thr!=nextThread; ++thr)
			{
				ThreadResult* result = 0;
				_resultLanes[thr]->read(result);
				delete result;
			}
			nextThread = 0;
		}
		++imageIndex;
	}
	for(size_t thr=0; thr!=nextThread; ++thr)
	{
		ThreadResult* result = 0;
		_resultLanes[thr]->read(result);
		delete result;
	}
}

ThreadedDeconvolutionTools::ThreadResult* ThreadedDeconvolutionTools::FinishMultiScaleTransformTask::operator()()
{
	msTransforms->FinishTransform(image, kernel);
	return 0;
}

void ThreadedDeconvolutionTools::MultiScaleTransform(MultiScaleTransforms* msTransforms, ImageBufferAllocator* allocator, const ao::uvector<double*>& images, ao::uvector<double> scales)
{
	size_t imageIndex = 0;
	size_t nextThread = 0;
	
	size_t scratchCount = std::min(images.size(), _threadCount);
	std::unique_ptr<ImageBufferAllocator::Ptr[]> scratchImages(
		new ImageBufferAllocator::Ptr[scratchCount]);
	for(size_t i=0; i!=scratchCount; ++i)
		allocator->Allocate(msTransforms->Width() * msTransforms->Height(), scratchImages[i]);
	
	while(imageIndex < images.size())
	{
		MultiScaleTransformTask* task = new MultiScaleTransformTask();
		task->msTransforms = msTransforms;
		task->image = images[imageIndex];
		task->scratch = scratchImages[nextThread].data();
		task->scale = scales[imageIndex];
		_taskLanes[nextThread]->write(task);
		
		++nextThread;
		if(nextThread == _threadCount)
		{
			for(size_t thr=0; thr!=nextThread; ++thr)
			{
				ThreadResult* result = 0;
				_resultLanes[thr]->read(result);
				delete result;
			}
			nextThread = 0;
		}
		++imageIndex;
	}
	for(size_t thr=0; thr!=nextThread; ++thr)
	{
		ThreadResult* result = 0;
		_resultLanes[thr]->read(result);
		delete result;
	}
}

ThreadedDeconvolutionTools::ThreadResult* ThreadedDeconvolutionTools::MultiScaleTransformTask::operator()()
{
	msTransforms->Transform(image, scratch, scale);
	return 0;
}

void ThreadedDeconvolutionTools::FindMultiScalePeak(MultiScaleTransforms* msTransforms, ImageBufferAllocator* allocator, const double* image, const ao::uvector<double>& scales, std::vector<ThreadedDeconvolutionTools::PeakData>& results, bool allowNegativeComponents, const bool* mask, double borderRatio, bool calculateRMS)
{
	size_t imageIndex = 0;
	size_t nextThread = 0;
	size_t resultIndex = 0;
	
	results.resize(scales.size());
	const size_t dataSize = msTransforms->Width() * msTransforms->Height();
	
	size_t size = std::min(scales.size(), _threadCount);
	std::unique_ptr<ImageBufferAllocator::Ptr[]> imageData(
		new ImageBufferAllocator::Ptr[size]);
	std::unique_ptr<ImageBufferAllocator::Ptr[]> scratchData(
		new ImageBufferAllocator::Ptr[size]);
	for(size_t i=0; i!=size; ++i)
	{
		allocator->Allocate(dataSize, imageData[i]);
		allocator->Allocate(dataSize, scratchData[i]);
	}
	
	while(imageIndex < scales.size())
	{
		FindMultiScalePeakTask* task = new FindMultiScalePeakTask();
		task->msTransforms = msTransforms;
		memcpy(imageData[nextThread].data(), image, dataSize*sizeof(double));
		task->image = imageData[nextThread].data();
		task->scratch = scratchData[nextThread].data();
		task->scale = scales[imageIndex];
		task->allowNegativeComponents = allowNegativeComponents;
		task->mask = mask;
		task->borderRatio = borderRatio;
		task->calculateRMS = calculateRMS;
		_taskLanes[nextThread]->write(task);
		
		++nextThread;
		if(nextThread == _threadCount)
		{
			for(size_t thr=0; thr!=nextThread; ++thr)
			{
				ThreadResult* result = 0;
				_resultLanes[thr]->read(result);
				results[resultIndex].value = static_cast<FindMultiScalePeakResult*>(result)->value;
				results[resultIndex].x = static_cast<FindMultiScalePeakResult*>(result)->x;
				results[resultIndex].y = static_cast<FindMultiScalePeakResult*>(result)->y;
				results[resultIndex].rms = static_cast<FindMultiScalePeakResult*>(result)->rms;
				delete result;
				++resultIndex;
			}
			nextThread = 0;
		}
		++imageIndex;
	}
	for(size_t thr=0; thr!=nextThread; ++thr)
	{
		ThreadResult* result = 0;
		_resultLanes[thr]->read(result);
		results[resultIndex].value = static_cast<FindMultiScalePeakResult*>(result)->value;
		results[resultIndex].x = static_cast<FindMultiScalePeakResult*>(result)->x;
		results[resultIndex].y = static_cast<FindMultiScalePeakResult*>(result)->y;
		results[resultIndex].rms = static_cast<FindMultiScalePeakResult*>(result)->rms;
		delete result;
		++resultIndex;
	}
}

ThreadedDeconvolutionTools::ThreadResult* ThreadedDeconvolutionTools::FindMultiScalePeakTask::operator()()
{
	msTransforms->Transform(image, scratch, scale);
	size_t width = msTransforms->Width(), height = msTransforms->Height();
	FindMultiScalePeakResult* result = new FindMultiScalePeakResult();
	if(calculateRMS)
		result->rms = RMS(image, width*height);
	else
		result->rms=-1.0;
	if(mask == 0)
		result->value = SimpleClean::FindPeak(image, width, height, result->x, result->y, allowNegativeComponents, 0, height, borderRatio);
	else
		result->value = SimpleClean::FindPeak(image, width, height, result->x, result->y, allowNegativeComponents, 0, height, mask, borderRatio);
	return result;
}
