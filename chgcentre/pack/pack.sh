WORKDIR=`pwd`
rm -rf /tmp/chgcentre
mkdir /tmp/chgcentre
cd ..
cp -vr CMakeLists.txt CMake/ *.{h,cpp} /tmp/chgcentre/
cd /tmp
tar -cjvf ${WORKDIR}/chgcentre.tar.bz2 chgcentre
rm -rf /tmp/chgcentre
tar -xjvf ${WORKDIR}/chgcentre.tar.bz2
cd /tmp/chgcentre
mkdir build
cd build
cmake ../
make -j 12
