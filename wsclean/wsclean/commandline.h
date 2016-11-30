#ifndef WSCLEAN_COMMAND_LINE_H
#define WSCLEAN_COMMAND_LINE_H

#include <cstring>

class CommandLine
{
public:
	static int Run(int argc, char *argv[]);
	
private:
	static void printHeader();
	static void printHelp();
	static size_t parse_size_t(const char* param, const char* name);
};

#endif

