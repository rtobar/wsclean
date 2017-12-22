#ifndef PARALLEL_FOR_H
#define PARALLEL_FOR_H

#include <cstring>
#include <thread>
#include <mutex>
#include <vector>

#include <sched.h>

namespace ao {

	template<typename Iter>
	class parallel_for
	{
	public:
		template<typename Function>
 		parallel_for(Iter start, Iter end, Function function, bool singleThreaded=false) :
			_current(start), _end(end)
		{
			unsigned nthreads = cpus();
			if(nthreads > 0)
				--nthreads;
			std::vector<std::thread> threads;
			if(!singleThreaded)
			{
				threads.reserve(nthreads);
				for(unsigned t=0; t!=nthreads; ++t)
					threads.push_back(std::thread(&parallel_for::run<Function>, this, function));
			}
			run<Function>(function);
			for(std::thread& thr : threads)
				thr.join();
		}
		
	private:
		template<typename Function>
		void run(Function function)
		{
			Iter iter;
			while(next(iter)) {
				function(iter);
			}
		}
		
		bool next(Iter& iter)
		{
			std::unique_lock<std::mutex> lock(_mutex);
			if(_current == _end)
				return false;
			else {
				iter = _current;
				++_current;
				return true;
			}
		}
		
		static unsigned cpus()
		{
	#ifdef __APPLE__
			return sysconf(_SC_NPROCESSORS_ONLN);
	#else
			cpu_set_t cs;
			CPU_ZERO(&cs);
			sched_getaffinity(0, sizeof cs , &cs);

			int count = 0;
			for (int i = 0; i < CPU_SETSIZE; i++)
			{
				if (CPU_ISSET(i, &cs))
					++count;
			}

			return count;
	#endif
		}
		
		Iter _current, _end;
		std::mutex _mutex;
	};
}

#endif
