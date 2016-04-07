#ifndef SYSTEM_H
#define SYSTEM_H

#include <stdio.h>
#include <unistd.h>
#include <sched.h>

class System {
public:
	static unsigned ProcessorCount()
	{
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
	}
};

#endif
