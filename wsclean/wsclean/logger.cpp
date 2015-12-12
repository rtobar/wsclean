#include "logger.h"

enum Logger::LoggerLevel Logger::_coutLevel = Logger::InfoLevel;

Logger::LogWriter<Logger::DebugLevel> Logger::Debug;

Logger::LogWriter<Logger::InfoLevel> Logger::Info;

Logger::LogWriter<Logger::WarningLevel> Logger::Warn;

Logger::LogWriter<Logger::ErrorLevel> Logger::Error;

Logger::LogWriter<Logger::FatalLevel> Logger::Fatal;

Logger::LogWriter<Logger::NoLevel, true> Logger::Progress;

void Logger::SetVerbosity(VerbosityLevel verbosityLevel)
{
	switch(verbosityLevel)
	{
		case QuietVerbosity:
			_coutLevel = NoLevel;
			break;
		case NormalVerbosity:
			_coutLevel = InfoLevel;
		break;
		case VerboseVerbosity:
			_coutLevel = DebugLevel;
			break;
	}
}
