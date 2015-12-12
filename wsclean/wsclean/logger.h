#ifndef AOFLAGGER
#define AOFLAGGER

#include <sstream>
#include <iostream>

#include <boost/thread/mutex.hpp>

class Logger
{
	public:
		enum LoggerLevel { NoLevel=5, FatalLevel=4, ErrorLevel=3, WarningLevel=2, InfoLevel=1, DebugLevel=0 };
		
		enum VerbosityLevel { QuietVerbosity, NormalVerbosity, VerboseVerbosity };

		template<enum LoggerLevel Level, bool ToStdErr=false>
		class LogWriter
		{
			public:
				LogWriter &operator<<(const std::string &str)
				{
					boost::mutex::scoped_lock lock(_mutex);
					ToStdOut(str);
					return *this;
				}
				LogWriter &operator<<(const char *str)
				{
					(*this) << std::string(str);
					return *this;
				}
				LogWriter &operator<<(const char c)
				{
					boost::mutex::scoped_lock lock(_mutex);
					ToStdOut(c);
					return *this;
				}
				template<typename S>
				LogWriter &operator<<(const S &str)
				{
					boost::mutex::scoped_lock lock(_mutex);
					ToStdOut(str);
					return *this;
				}
				void Flush()
				{
					boost::mutex::scoped_lock lock(_mutex);
					std::cout.flush();
				}
			private:
				std::stringstream _buffer;
				boost::mutex _mutex;

				template<typename S>
				void ToStdOut(const S &str)
				{
					if((int) _coutLevel <= (int) Level)
					{
						if(ToStdErr)
							std::cerr << str;
						else
							std::cout << str;
					}
				}
		};

		static void SetVerbosity(VerbosityLevel verbosityLevel);
		
		static bool IsVerbose() { return _coutLevel==DebugLevel; }

		static class LogWriter<DebugLevel> Debug;
		static class LogWriter<InfoLevel> Info;
		static class LogWriter<WarningLevel> Warn;
		static class LogWriter<ErrorLevel> Error;
		static class LogWriter<FatalLevel> Fatal;
		static class LogWriter<NoLevel, true> Progress;
	private:
		Logger()
		{
		}

		static enum LoggerLevel _coutLevel;
};

#endif
