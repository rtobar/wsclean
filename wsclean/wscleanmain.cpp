#include "wsclean/commandline.h"
#include "wsclean/logger.h"
#include "wsclean/wsclean.h"

#include <exception>
#include <iostream>

int main(int argc, char *argv[])
{
	try {
		WSClean wsclean;
		if(CommandLine::Parse(wsclean, argc, argv))
			CommandLine::Run(wsclean);
		return 0;
	} catch(std::exception& e)
	{
		Logger::Error
			<< "+ + + + + + + + + + + + + + + + + + +\n"
			<< "+ An exception occured:\n"
			<< "+ >>> " << e.what() << "\n"
			<< "+ + + + + + + + + + + + + + + + + + +\n";
		return -1;
	}
}
