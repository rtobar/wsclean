cmake_minimum_required(VERSION 2.8.8)

project(wsclean)

#if(NOT CMAKE_BUILD_TYPE)
#  set(CMAKE_BUILD_TYPE "Release")
#endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -O3 -Wall -DNDEBUG")

if(PORTABLE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native") 
endif(PORTABLE)

# Casacore has a separate CMake file in this directory
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake)

SET(CASACORE_MAKE_REQUIRED_EXTERNALS_OPTIONAL TRUE)
find_package(Casacore REQUIRED COMPONENTS casa ms tables measures fits)

find_library(FFTW3_LIB fftw3 REQUIRED)
find_library(FFTW3_THREADS_LIB fftw3_threads REQUIRED)
#Prevent accidentally finding old BoostConfig.cmake file from casapy
set(Boost_NO_BOOST_CMAKE ON)
find_package(Boost COMPONENTS date_time filesystem thread system REQUIRED)
find_library(PTHREAD_LIB pthread REQUIRED)
find_library(FITSIO_LIB cfitsio REQUIRED)
find_library(DL_LIB dl REQUIRED)
find_path(FITSIO_INCLUDE_DIR NAMES fitsio.h)
find_path(FFTW3_INCLUDE_DIR NAMES fftw3.h)
find_library(GSL_LIB NAMES gsl)
find_path(GSL_INCLUDE_DIR NAMES gsl/gsl_version.h)
find_library(GSL_CBLAS_LIB NAMES gslcblas)

set(CASA_LIBS ${CASA_MS_LIB} ${CASA_TABLES_LIB} ${CASA_MEASURES_LIB} ${CASA_CASA_LIB} ${CASA_FITS_LIB} ${DL_LIB})

find_path(CASA_INCLUDE_DIR NAMES casacore/ms/MeasurementSets/MeasurementSet.h PATHS /usr/local/include/)

find_path(LOFAR_STATION_RESPONSE_DIR NAMES StationResponse/Station.h)
find_library(LOFAR_STATION_RESPONSE_LIB stationresponse)
if(LOFAR_STATION_RESPONSE_DIR AND LOFAR_STATION_RESPONSE_LIB)
  include_directories(${LOFAR_STATION_RESPONSE_DIR})
  set(LBEAM_LIBS ${LOFAR_STATION_RESPONSE_LIB})
  set(LBEAM_FILES lofar/lbeamevaluator.cpp lofar/lmspredicter.cpp)
  # Some versions of the LOFAR code require "HAVE_BOOST" to be set explicitly
  # (This is supposed to be fixed in newer LOFAR repositories, AO 2015-11-27)
  add_definitions(-DHAVE_LOFAR_BEAM -DHAVE_BOOST)
  message("LOFAR beam library found.")
else()
  set(LBEAM_LIBS "")
  set(LBEAM_FILES lofar/lmspredicter.cpp)
  message("LOFAR beam library not found: LOFAR beam will not be available in WSClean. This is not an error, it means that LOFAR beam correction and the 'lofar-dft' mode will be unavailable.")
endif(LOFAR_STATION_RESPONSE_DIR AND LOFAR_STATION_RESPONSE_LIB)

include_directories(${CASACORE_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIR})
include_directories(${FITSIO_INCLUDE_DIR})
include_directories(${FFTW3_INCLUDE_DIR})
include_directories(${GSL_INCLUDE_DIR})

include(CheckCXXSourceCompiles)

check_cxx_source_compiles(
"#include \"${CMAKE_CURRENT_SOURCE_DIR}/aocommon/uvector_11.h\"
int main(int argc, char* argv[])
{
return 0;
}" UVECTOR_COMPILES)
check_cxx_source_compiles(
"#include \"${CMAKE_CURRENT_SOURCE_DIR}/aocommon/lane_11.h\"
int main(int argc, char* argv[])
{
return 0;
}" LANE11_COMPILES)

if(UVECTOR_COMPILES)
  add_definitions(-DHAVE_UVECTOR)
endif(UVECTOR_COMPILES)
if(LANE11_COMPILES)
  add_definitions(-DHAVE_LANE11)
endif(LANE11_COMPILES)

# GSL is required for WSClean, so always available
add_definitions(-DHAVE_GSL)

# The following stuff will set the "rpath" correctly, so that
# LD_LIBRARY_PATH doesn't have to be set.

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

add_library(wsclean-object OBJECT
  casamaskreader.cpp dftpredictionalgorithm.cpp fftconvolver.cpp fftresampler.cpp fftwmultithreadenabler.cpp fitsiochecker.cpp fitsreader.cpp fitswriter.cpp image.cpp imageweights.cpp modelrenderer.cpp multibanddata.cpp nlplfitter.cpp polynomialfitter.cpp progressbar.cpp stopwatch.cpp
  deconvolution/deconvolution.cpp deconvolution/deconvolutionalgorithm.cpp deconvolution/genericclean.cpp deconvolution/imageset.cpp deconvolution/moresane.cpp deconvolution/simpleclean.cpp deconvolution/spectralfitter.cpp
  interface/wscleaninterface.cpp
  iuwt/imageanalysis.cpp iuwt/iuwtdecomposition.cpp iuwt/iuwtdeconvolutionalgorithm.cpp iuwt/iuwtmask.cpp
  lofar/lbeamimagemaker.cpp
  model/model.cpp
  msproviders/averagingmsrowprovider.cpp msproviders/contiguousms.cpp msproviders/directmsrowprovider.cpp msproviders/msprovider.cpp msproviders/msrowprovider.cpp msproviders/partitionedms.cpp
  multiscale/clarkloop.cpp multiscale/multiscalealgorithm.cpp multiscale/multiscaletransforms.cpp multiscale/threadeddeconvolutiontools.cpp
  wsclean/commandline.cpp wsclean/imagingtable.cpp wsclean/logger.cpp wsclean/msgridderbase.cpp wsclean/wscfitswriter.cpp wsclean/wsclean.cpp
  wsclean/wscleansettings.cpp wsclean/wsmsgridder.cpp wsclean/wstackinggridder.cpp
	${LBEAM_FILES})

set_property(TARGET wsclean-object PROPERTY POSITION_INDEPENDENT_CODE 1)

set(WSCLEANFILES $<TARGET_OBJECTS:wsclean-object>)

add_library(wsclean-lib STATIC ${WSCLEANFILES})
set_target_properties(wsclean-lib PROPERTIES OUTPUT_NAME wsclean)

add_library(wsclean-shared SHARED ${WSCLEANFILES})
target_link_libraries(wsclean-shared)

add_executable(wsclean wscleanmain.cpp)
target_link_libraries(wsclean wsclean-lib ${CASACORE_LIBRARIES} ${FFTW3_LIB} ${FFTW3_THREADS_LIB} ${Boost_DATE_TIME_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${FITSIO_LIB} ${GSL_LIB} ${GSL_CBLAS_LIB} ${PTHREAD_LIB} ${LBEAM_LIBS})

add_executable(interfaceexample EXCLUDE_FROM_ALL interface/interfaceexample.c)
target_link_libraries(interfaceexample wsclean-lib ${CASACORE_LIBRARIES} ${FFTW3_LIB} ${FFTW3_THREADS_LIB} ${Boost_DATE_TIME_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${FITSIO_LIB} ${GSL_LIB} ${GSL_CBLAS_LIB} ${LBEAM_LIBS})

set_target_properties(wsclean-object PROPERTIES COMPILE_FLAGS "-std=c++0x")
set_target_properties(wsclean PROPERTIES COMPILE_FLAGS "-std=c++0x")
set_target_properties(wsclean-lib PROPERTIES COMPILE_FLAGS "-std=c++0x")
set_target_properties(wsclean-shared PROPERTIES COMPILE_FLAGS "-std=c++0x")

install(TARGETS wsclean DESTINATION bin)
install(TARGETS wsclean-lib DESTINATION lib)
install(FILES interface/wscleaninterface.h DESTINATION include)

# add target to generate API documentation with Doxygen
find_package(Doxygen)

if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM)

endif(DOXYGEN_FOUND)

find_package(Boost COMPONENTS date_time thread system filesystem unit_test_framework)
if(Boost_FOUND)
  add_executable(runtest EXCLUDE_FROM_ALL
		tests/test.cpp tests/testclean.cpp tests/testradeccoord.cpp
		tests/testpolynomialfitter.cpp tests/testfitsdateobstime.cpp
		tests/testbaselinedependentaveraging.cpp
		tests/testimageset.cpp
		tests/testimage.cpp
		${WSCLEANFILES})
  target_link_libraries(runtest ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ${CASACORE_LIBRARIES} ${FFTW3_LIB} ${FFTW3_THREADS_LIB} ${Boost_DATE_TIME_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${CFITSIO_LIBRARY} ${GSL_LIB} ${GSL_CBLAS_LIB} ${LBEAM_LIBS} ${IDG_OPT_LIBS})
  add_test(runtest runtest)
  add_custom_target(check COMMAND runtest -l unit_scope DEPENDS runtest)
else()
  message("Boost testing framework not found (not required for wsclean: only required for running make check).")
endif()

message(STATUS "Flags passed to C++ compiler: " ${CMAKE_CXX_FLAGS})
