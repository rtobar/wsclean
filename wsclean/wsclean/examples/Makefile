wsgridderexample:	wspredictionexample.cpp ../wstackinggridder.cpp ../logger.cpp ../../fftwmultithreadenabler.cpp
	g++ -Wall -o wspredictionexample -std=c++11 -DAVOID_CASACORE wspredictionexample.cpp ../wstackinggridder.cpp ../logger.cpp ../../fftwmultithreadenabler.cpp -lfftw3 -lfftw3_threads -lboost_date_time -lboost_thread -lboost_system

channelfitexample:	channelfitexample.cpp
	g++ -Wall -o channelfitexample -std=c++11 channelfitexample.cpp ../../polynomialchannelfitter.cpp ../../polynomialfitter.cpp -lgsl -lgslcblas
