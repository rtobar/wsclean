#ifndef MSIOSYSTEM_H
#define MSIOSYSTEM_H

#include <casacore/casa/OS/HostInfo.h>

#include <stdio.h>
#include <unistd.h>
#include <sched.h>

class System
{
	public:
		static long TotalMemory()
		{
			return casacore::HostInfo::memoryTotal()*1024;
		}
		
		static unsigned ProcessorCount()
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
				count++;
			}

			return count;
#endif
		}
};

#endif //MSIOSYSTEM_H
