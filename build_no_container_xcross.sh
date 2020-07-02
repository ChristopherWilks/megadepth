#!/usr/bin/env bash

set -e

export APPLE_VER=12
OSXCROSS_ROOT=/opt/osxcross/target/bin
export PATH=$OSXCROSS_ROOT:$PATH
export OSX_ARCH=x86_64
export OSX_CC=o64-gcc
export OSX_CXX=o64-g++
export CC=${OSX_CC}
export CXX=${OSX_CXX}
export AR=${OSXCROSS_ROOT}/$OSX_ARCH-apple-darwin${APPLE_VER}-ar
export RANLIB=${OSXCROSS_ROOT}/$OSX_ARCH-apple-darwin${APPLE_VER}-ranlib

#build dynamic by default
build_type=$1
bc=`perl -e '$bt="'$build_type'"; if($bt=~/static/i) { print "megadepth_static"; } elsif($bt=~/statlib/i) { print "megadepth_statlib"; } else { print "megadepth_dynamic"; }'`

if [[ ! -s libdeflate ]] ; then
    ./get_libdeflate.sh
fi

if [[ ! -e htslib/libhts.a ]] ; then
    export CPPFLAGS="-I../libdeflate"
    export LDFLAGS="-L../libdeflate -ldeflate"
    ./get_htslib.sh $OSX_ARCH-apple-darwin${APPLE_VER}-gcc
    export CPPFLAGS=
    export LDFLAGS=
fi

if [[ ! -s libBigWig ]] ; then
    ./get_libBigWig.sh
fi

set -x

#compile a no-curl static version of libBigWig
if [[ $bc == 'megadepth_static' ]]; then
    pushd libBigWig
    make clean
    make -f Makefile.nocurl lib-static
    popd
fi

#compile a position indenpendent code static version of libBigWig
#but *with* curl calls present, to be resolved later in the
#dynamic linking of megadepth
if [[ $bc == 'megadepth_statlib' ]]; then
    pushd libBigWig
    make clean
    make -f Makefile.fpic lib-static
    popd
fi

#compile a the original, dynamic version of libBigWig
if [[ $bc == 'megadepth_dynamic' ]]; then
    pushd libBigWig
    make clean
    make -f Makefile.orig lib-shared
    popd
fi

export LD_LIBRARY_PATH=./htslib:./libBigWig:$LD_LIBRARY_PATH

DR=build-release-temp
mkdir -p ${DR}
pushd ${DR}
#CC=${OSX_CC} CXX=${OSX_CXX} cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_AR=/opt/osxcross/target/bin/$OSX_ARCH-apple-darwin${APPLE_VER}-ar -DCMAKE_RANLIB=/opt/osxcross/target/bin/$OSX_ARCH-apple-darwin${APPLE_VER}-ranlib -D CMAKE_C_COMPILER=${OSX_CC} -D CMAKE_CXX_COMPILER=${OSX_CXX} ..
#CC=${OSX_CC} CXX=${OSX_CXX} AR=/opt/osxcross/target/bin/$OSX_ARCH-apple-darwin${APPLE_VER}-ar RANLIB=/opt/osxcross/target/bin/$OSX_ARCH-apple-darwin${APPLE_VER}-ranlib make ${bc}
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_AR=${AR} -DCMAKE_RANLIB=${RANLIB} -D CMAKE_C_COMPILER=${OSX_CC} -D CMAKE_CXX_COMPILER=${OSX_CXX} ..
make ${bc}
popd
cp ${DR}/${bc} ./
ln -fs ./$bc megadepth
./megadepth --version
rm -rf ${DR}

DR=build-debug-temp
mkdir -p ${DR}
pushd ${DR}
cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_AR=${AR} -DCMAKE_RANLIB=${RANLIB} -D CMAKE_C_COMPILER=${OSX_CC} -D CMAKE_CXX_COMPILER=${OSX_CXX} ..
make ${bc}
popd
cp ${DR}/${bc} ./${bc}_debug
ln -fs ./${bc}_debug megadepth_debug
./megadepth_debug --version
rm -rf ${DR}
