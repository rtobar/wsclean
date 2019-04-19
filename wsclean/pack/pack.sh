NCPU=16
if [[ "$1" == "" ]] ; then
    echo Syntax: ./pack.sh \<version\>
else
    VERSION="$1"
    WORKDIR=`pwd`
    rm -rf /tmp/wsclean /tmp/wsclean-${VERSION}
    mkdir /tmp/wsclean
    mkdir /tmp/wsclean/aocommon
    mkdir /tmp/wsclean/atca
    mkdir /tmp/wsclean/aterms
    mkdir /tmp/wsclean/CMake
    mkdir /tmp/wsclean/deconvolution
    mkdir /tmp/wsclean/idg
    mkdir /tmp/wsclean/interface
    mkdir /tmp/wsclean/iuwt
    mkdir /tmp/wsclean/lofar
    mkdir /tmp/wsclean/model
    mkdir /tmp/wsclean/msproviders
    mkdir /tmp/wsclean/multiscale
    mkdir /tmp/wsclean/mwa
    mkdir /tmp/wsclean/tests
    mkdir /tmp/wsclean/units
    mkdir /tmp/wsclean/wsclean
    mkdir /tmp/wsclean/wsclean/examples
    cd ..
    cp -v CMakeLists.txt CMakeVersionInfo.txt wscversion.h.in Doxyfile.in application.* areaset.* banddata.* buffered_lane.* casamaskreader.* dftpredictionalgorithm.* fftconvolver.* fftresampler.* fftwmanager.* fitsiochecker.* fitsreader.* fitswriter.* gaussianfitter.* image.* imageweights.* lane.* multibanddata.* nlplfitter.* numberlist.* matrix2x2.* modelrenderer.* msselection.* ndppp.* parsetreader.* polarization.* polynomialchannelfitter.* polynomialfitter.* progressbar.* rmsimage.* stopwatch.* system.* threadpool.* uvector.* weightmode.* wscleanmain.cpp /tmp/wsclean/
    cp -v aocommon/barrier.h aocommon/cloned_ptr.h aocommon/lane.h aocommon/lane_11.h aocommon/parallelfor.h aocommon/uvector_11.h /tmp/wsclean/aocommon/
    cp -v atca/*.{h,cpp} /tmp/wsclean/atca/
    cp -v aterms/*.{h,cpp} /tmp/wsclean/aterms/
    cp -v CMake/*.cmake /tmp/wsclean/CMake/
    cp -v deconvolution/*.{h,cpp} /tmp/wsclean/deconvolution/
    cp -v idg/*.{h,cpp} /tmp/wsclean/idg/
    cp -v interface/*.{c,h,cpp,py} /tmp/wsclean/interface/
    cp -v iuwt/*.{h,cpp} /tmp/wsclean/iuwt/
    cp -v lofar/*.{h,cpp} /tmp/wsclean/lofar/
    cp -v model/*.{h,cpp} /tmp/wsclean/model/
    cp -v msproviders/*.{h,cpp} /tmp/wsclean/msproviders
    cp -v multiscale/*.{h,cpp} /tmp/wsclean/multiscale
    cp -v mwa/*.{h,cpp} /tmp/wsclean/mwa
    cp -v tests/*.cpp /tmp/wsclean/tests
    cp -v units/*.h /tmp/wsclean/units
    cp -v wsclean/*.{h,cpp} /tmp/wsclean/wsclean
    cp -v wsclean/examples/{Makefile,*.cpp} /tmp/wsclean/wsclean/examples/
    cd /tmp
    mv wsclean wsclean-${VERSION}
    tar -cjvf ${WORKDIR}/wsclean-${VERSION}.tar.bz2 wsclean-${VERSION}/
    rm -rf /tmp/wsclean-${VERSION}
    tar -xjvf ${WORKDIR}/wsclean-${VERSION}.tar.bz2
    cd /tmp/wsclean-${VERSION}
    
    echo Building WITHOUT LOFAR lib
    mkdir build
    cd build
    cmake ../
    make -j ${NCPU} && make check -j ${NCPU} && ./wsclean --version|tee ../runs.txt
    cd ..
    
    echo Building WITH LOFAR lib
    rm build -rf
    mkdir build
    cd build
    cmake ../ -DCMAKE_PREFIX_PATH="/home/anoko/Software/LOFARBeam-install"
    make -j ${NCPU} && make check -j ${NCPU} && ./wsclean --version|tee -a ../runs.txt
    cd ..
    
    echo Building WITH LOFAR lib AND IDG
    rm build -rf
    mkdir build
    cd build
    cmake ../ -DCMAKE_PREFIX_PATH="/home/anoko/Software/LOFARBeam-install;/home/anoko/Software/idg-install"
    make -j ${NCPU} && make check -j ${NCPU} && ./wsclean --version|tee -a ../runs.txt
    cd ..
    
    cat runs.txt
    rm -f runs.txt
    cat CMakeVersionInfo.txt
fi
