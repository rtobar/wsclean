if [[ "$1" == "" ]] ; then
		echo Syntax: ./pack.sh \<version\>
else
		VERSION="$1"
		WORKDIR=`pwd`
		rm -rf /tmp/wsclean /tmp/wsclean-${VERSION}
		mkdir /tmp/wsclean
		mkdir /tmp/wsclean/aocommon
		mkdir /tmp/wsclean/CMake
		mkdir /tmp/wsclean/deconvolution
		mkdir /tmp/wsclean/idg
		mkdir /tmp/wsclean/interface
		mkdir /tmp/wsclean/iuwt
		mkdir /tmp/wsclean/lofar
		mkdir /tmp/wsclean/model
		mkdir /tmp/wsclean/msproviders
		mkdir /tmp/wsclean/multiscale
		mkdir /tmp/wsclean/wsclean
		mkdir /tmp/wsclean/wsclean/examples
		cd ..
		cp -v CMakeLists.txt Doxyfile.in angle.h application.* areaset.* banddata.* buffered_lane.* casamaskreader.* dftpredictionalgorithm.* fftconvolver.* fftresampler.* fftwmultithreadenabler.* fitsiochecker.* fitsreader.* fitswriter.* gaussianfitter.h imagecoordinates.* imageoperations.* imageweights.* lane.* multibanddata.* nlplfitter.* numberlist.* matrix2x2.* modelrenderer.* msselection.* ndppp.* polarizationenum.* polynomialfitter.* progressbar.* radeccoord.* stopwatch.* system.* threadpool.* uvector.* weightmode.* wscleanmain.cpp wscversion.h /tmp/wsclean/
		cp -v aocommon/lane.h aocommon/lane_03.h aocommon/lane_11.h aocommon/uvector.h aocommon/uvector_03.h /tmp/wsclean/aocommon/
		cp -v CMake/*.cmake /tmp/wsclean/CMake/
		cp -v deconvolution/*.{h,cpp} /tmp/wsclean/deconvolution/
		cp -v idg/*.{h,cpp} /tmp/wsclean/idg/
		cp -v interface/*.{c,h,cpp,py} /tmp/wsclean/interface/
		cp -v iuwt/*.{h,cpp} /tmp/wsclean/iuwt/
		cp -v lofar/*.{h,cpp} /tmp/wsclean/lofar/
		cp -v model/*.{h,cpp} /tmp/wsclean/model/
		cp -v msproviders/*.{h,cpp} /tmp/wsclean/msproviders
		cp -v multiscale/*.{h,cpp} /tmp/wsclean/multiscale
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
		make -j 4
		cat ../wscversion.h
		cd ..

		echo Building WITH LOFAR lib
		rm build -rf
		mkdir build
		cd build
		cmake ../ -DCMAKE_PREFIX_PATH="/home/anoko/Software/LOFAR-install"
		make -j 4
		cd ..
fi
