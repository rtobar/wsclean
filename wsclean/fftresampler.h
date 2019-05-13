#ifndef FFT_RESAMPLE_H
#define FFT_RESAMPLE_H

#include "lane.h"

#include <vector>

#include <fftw3.h>
#include <boost/thread/thread.hpp>

#include "uvector.h"
#include "windowfunction.h"

class FFTResampler
{
private:
	struct Task
	{
		double *input, *output;
	};
	
public:
	FFTResampler(size_t inWidth, size_t inHeight, size_t outWidth, size_t outHeight, size_t cpuCount, bool verbose=true);
	
	~FFTResampler();
	
	void AddTask(double* input, double* output)
	{
		Task task;
		task.input = input;
		task.output = output;
		_tasks.write(task);
	}
	
	void Start()
	{
		for(size_t i=0; i!=_tasks.capacity(); ++i)
		{
			_threads.add_thread(new boost::thread(&FFTResampler::runThread, this));
		}
	}
	
	void Finish()
	{
		_tasks.write_end();
		_threads.join_all();
		_tasks.clear();
	}
	
	void Resample(double* input, double* output)
	{
		Task task;
		task.input = input;
		task.output = output;
		runSingle(task, false);
	}
	
	void SingleFT(const double* input, double* realOutput, double* imaginaryOutput);
	
	/**
	 * Only to be used with SingleFT (it makes resampling thread unsafe!)
	 */
	void SetTukeyWindow(double insetSize, bool correctWindow)
	{
		_windowFunction = WindowFunction::Tukey;
		_tukeyInsetSize = insetSize;
		_correctWindow = correctWindow;
		_windowRowIn.clear();
		_windowColIn.clear();
		_windowOut.clear();
	}
	
	void SetWindowFunction(WindowFunction::Type window, bool correctWindow)
	{
		_windowFunction = window;
		_correctWindow = correctWindow;
		_windowRowIn.clear();
		_windowColIn.clear();
		_windowOut.clear();
	}
	
private:
	void runThread();
	void runSingle(const Task& task, bool skipWindow) const;
	void applyWindow(double* data) const;
	void unapplyWindow(double* data) const;
	void makeWindow(ao::uvector<double>& data, size_t width) const;
	void makeTukeyWindow(ao::uvector<double>& data, size_t width) const;
	
	size_t _inputWidth, _inputHeight;
	size_t _outputWidth, _outputHeight;
	size_t _fftWidth, _fftHeight;
	WindowFunction::Type _windowFunction;
	double _tukeyInsetSize;
	mutable ao::uvector<double> _windowRowIn;
	mutable ao::uvector<double> _windowColIn;
	mutable ao::uvector<double> _windowOut;
	bool _correctWindow;
	
	fftw_plan _inToFPlan, _fToOutPlan;
	
	ao::lane<Task> _tasks;
	boost::thread_group _threads;
	bool _verbose;
};

#endif
