#ifndef WSCLEAN_COMMAND_LINE_H
#define WSCLEAN_COMMAND_LINE_H

#include <string>
#include <cstring>

#include "wscleansettings.h"

class CommandLine
{
public:
	static bool Parse(class WSClean& wsclean, int argc, char *argv[]);
	static void Run(class WSClean& wsclean);
	
private:
	static void deprecated(const std::string& param, const std::string& replacement);
	static void printHeader();
	static void printHelp();
	static size_t parse_size_t(const char* param, const char* name);
	static double parse_double(const char* param, double lowerLimit, const char* name, bool inclusive=true);
	static double parse_double(const char* param, const char* name);
};

#endif

