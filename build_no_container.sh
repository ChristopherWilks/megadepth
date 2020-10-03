#!/usr/bin/env bash

set -ex

#build dynamic by default
build_type=$1
bc=`perl -e '$bt="'$build_type'"; if($bt=~/static/i) { print "megadepth_static"; } elsif($bt=~/statlib/i) { print "megadepth_statlib"; } else { print "megadepth_dynamic"; }'`

ln -fs CMakeLists.txt.ci CMakeLists.txt

#clear symlink main lib dirs
rm -rf zlib htslib libBigWig libdeflate build-release-temp

if [[ ! -s zlib_ci ]] ; then
    ./get_zlib.sh
    mv zlib zlib_ci
fi
ln -fs zlib_ci zlib

if [[ ! -s libdeflate_ci ]] ; then
    ./get_libdeflate.sh
    mv libdeflate libdeflate_ci
fi
ln -fs libdeflate_ci libdeflate

if [[ ! -s htslib_ci ]] ; then
    export CPPFLAGS="-I../libdeflate"
    export LDFLAGS="-L../libdeflate -ldeflate"
    ./get_htslib.sh linux
    export CPPFLAGS=
    export LDFLAGS=
    mv htslib htslib_ci
fi
ln -fs htslib_ci htslib

if [[ ! -s libBigWig_ci ]] ; then
    ./get_libBigWig.sh
    mv libBigWig libBigWig_ci
fi
ln -fs libBigWig_ci libBigWig

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
cmake -DCMAKE_BUILD_TYPE=Release ..
make ${bc}
popd
cp ${DR}/${bc} ./${bc}_ci
ln -fs ./${bc}_ci megadepth
./megadepth --version
rm -rf ${DR}

DR=build-debug-temp
mkdir -p ${DR}
pushd ${DR}
cmake -DCMAKE_BUILD_TYPE=Debug ..
make ${bc}
popd
cp ${DR}/${bc} ./${bc}_debug
ln -fs ./${bc}_debug megadepth_debug
./megadepth_debug --version
rm -rf ${DR}
