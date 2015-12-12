#include "progressbar.h"

#include "wsclean/logger.h"

#include <iostream>

ProgressBar::ProgressBar(const std::string& taskDescription) :
	_taskDescription(taskDescription),
	_displayedDots(0)
{
	Logger::Info << taskDescription << ":";
	if(taskDescription.size() < 40)
		Logger::Info << " 0%";
	else
		Logger::Info << "\n 0%";
	Logger::Info.Flush();
}

ProgressBar::~ProgressBar()
{
	SetProgress(1,1);
}

void ProgressBar::SetProgress(size_t taskIndex, size_t taskCount)
{
	unsigned progress = (taskIndex * 100 / taskCount);
	unsigned dots = progress / 2;
	
	if(dots > _displayedDots)
	{
		while(dots != _displayedDots)
		{
			++_displayedDots;
			if(_displayedDots % 5 == 0)
				Logger::Info << ((_displayedDots/5)*10) << '%';
			else
				Logger::Info << '.';
		}
		if(progress == 100)
			Logger::Info << '\n';
		Logger::Info.Flush();
	}
}
