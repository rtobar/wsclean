#include "wscleansettings.h"

#include <sstream>

void WSCleanSettings::Validate()
{
	// antialiasingKernelSize should be odd
	if(antialiasingKernelSize%2 == 0)
	{
		std::stringstream s;
		s << "Bad anti-aliasing kernel size given of " << antialiasingKernelSize << ". The kernel size has to be odd.";
		throw std::runtime_error(s.str());
	}
}
