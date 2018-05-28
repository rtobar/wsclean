#include "progressbar.h"

#include "wsclean/logger.h"

#include <iostream>

ProgressBar::ProgressBar(const std::string& taskDescription) :
	_taskDescription(taskDescription),
	_displayedDots(-1)
{
}

ProgressBar::~ProgressBar()
{
	SetProgress(1,1);
}

ProgressBar& ProgressBar::operator=(ProgressBar&& rhs)
{
	SetProgress(1, 1);
	_displayedDots = rhs._displayedDots;
	_taskDescription = std::move(rhs._taskDescription);
	rhs._displayedDots = 50;
	return *this;
}

void ProgressBar::SetProgress(size_t taskIndex, size_t taskCount)
{
	if(_displayedDots==-1) {
		Logger::Info << _taskDescription << ":";
		if(_taskDescription.size() < 40)
			Logger::Info << " 0%";
		else
			Logger::Info << "\n 0%";
		_displayedDots=0;
		Logger::Info.Flush();
	}
	int progress = (taskIndex * 100 / taskCount);
	int dots = progress / 2;
	
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
