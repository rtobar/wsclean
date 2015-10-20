cmake_minimum_required(VERSION 2.6)

project(wsclean)

if(PORTABLE)
  set(CMAKE_CXX_FLAGS "-O3 -Wall -DNDEBUG -ggdb")
else()
  set(CMAKE_CXX_FLAGS "-O3 -Wall -DNDEBUG -march=native -msse4.1 -ggdb") 
endif(PORTABLE)
  
find_library(CASA_MS_LIB casa_ms REQUIRED)
find_library(CASA_TABLES_LIB casa_tables REQUIRED)
find_library(CASA_CASA_LIB casa_casa REQUIRED)
find_library(CASA_MEASURES_LIB casa_measures REQUIRED)
find_library(CASA_FITS_LIB casa_fits REQUIRED)
find_library(FFTW3_LIB fftw3 REQUIRED)
find_library(FFTW3_THREADS_LIB fftw3_threads REQUIRED)
find_package(Boost COMPONENTS filesystem thread system REQUIRED)
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
  add_definitions(-DHAVE_LOFAR_BEAM)
  message("LOFAR beam library found.")
else()
  set(LBEAM_LIBS "")
  set(LBEAM_FILES lofar/lmspredicter.cpp)
  message("LOFAR beam library not found: LOFAR beam will not be available in WSClean. This is not an error, it means the 'lofar-dft' mode will be unavailable, which is not used often. Hence, this message can probably be ignored.")
endif(LOFAR_STATION_RESPONSE_DIR AND LOFAR_STATION_RESPONSE_LIB)

include_directories(${CASA_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIR})
include_directories(${FITSIO_INCLUDE_DIR})
include_directories(${FFTW3_INCLUDE_DIR})
include_directories(${GSL_INCLUDE_DIR})

set(CMAKE_REQUIRED_FLAGS "-std=c++0x")
include(CheckCXXSourceCompiles)

check_cxx_source_compiles(
"#include \"${CMAKE_CURRENT_SOURCE_DIR}/aocommon/uvector.h\"
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
  casamaskreader.cpp dftpredictionalgorithm.cpp fftconvolver.cpp fftresampler.cpp fftwmultithreadenabler.cpp fitsiochecker.cpp fitsreader.cpp fitswriter.cpp imageweights.cpp nlplfitter.cpp modelrenderer.cpp polynomialfitter.cpp progressbar.cpp stopwatch.cpp
  deconvolution/deconvolution.cpp deconvolution/deconvolutionalgorithm.cpp deconvolution/dynamicset.cpp deconvolution/fastmultiscaleclean.cpp deconvolution/joinedclean.cpp deconvolution/moresane.cpp deconvolution/simpleclean.cpp deconvolution/spectralfitter.cpp
  interface/wscleaninterface.cpp
  iuwt/imageanalysis.cpp iuwt/iuwtdecomposition.cpp iuwt/iuwtdeconvolutionalgorithm.cpp iuwt/iuwtmask.cpp
  model/model.cpp
  msproviders/contiguousms.cpp msproviders/msprovider.cpp msproviders/partitionedms.cpp
  multiscale/multiscalealgorithm.cpp multiscale/multiscaletransforms.cpp multiscale/threadeddeconvolutiontools.cpp
  wsclean/imagingtable.cpp wsclean/wsclean.cpp wsclean/wsmsgridder.cpp wsclean/wstackinggridder.cpp
	${LBEAM_FILES})

set_property(TARGET wsclean-object PROPERTY POSITION_INDEPENDENT_CODE 1)

set(WSCLEANFILES $<TARGET_OBJECTS:wsclean-object>)

add_library(wsclean-lib STATIC ${WSCLEANFILES})
set_target_properties(wsclean-lib PROPERTIES OUTPUT_NAME wsclean)

add_library(wsclean-shared SHARED ${WSCLEANFILES})
target_link_libraries(wsclean-shared)

add_executable(wsclean wscleanmain.cpp)
target_link_libraries(wsclean wsclean-lib ${CASA_LIBS} ${FFTW3_LIB} ${FFTW3_THREADS_LIB} ${Boost_FILESYSTEM_LIBRARY} ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${FITSIO_LIB} ${GSL_LIB} ${GSL_CBLAS_LIB} ${PTHREAD_LIB} ${LBEAM_LIBS})

add_executable(interfaceexample EXCLUDE_FROM_ALL interface/interfaceexample.c)
target_link_libraries(interfaceexample wsclean-lib ${CASA_LIBS} ${FFTW3_LIB} ${FFTW3_THREADS_LIB} ${Boost_FILESYSTEM_LIBRARY} ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${FITSIO_LIB} ${GSL_LIB} ${GSL_CBLAS_LIB} ${LBEAM_LIBS})

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
