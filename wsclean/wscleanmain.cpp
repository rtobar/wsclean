#include "wsclean/commandline.h"
#include "wsclean/logger.h"

#include <exception>
#include <iostream>

int main(int argc, char *argv[])
{
	try {
		int returnValue = CommandLine::Run(argc, argv);
		return returnValue;
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
