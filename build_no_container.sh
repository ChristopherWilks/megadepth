#!/usr/bin/env bash

set -e

#build dynamic by default
build_type=$1
bc=`perl -e '$bt="'$build_type'"; if($bt=~/static/i) { print "megadepth_static"; } elsif($bt=~/statlib/i) { print "megadepth_statlib"; } else { print "megadepth_dynamic"; }'`

if [[ ! -s zlib ]] ; then
    ./get_zlib.sh
fi

if [[ ! -s libdeflate ]] ; then
    ./get_libdeflate.sh
fi

if [[ ! -s htslib ]] ; then
    export CPPFLAGS="-I../libdeflate"
    export LDFLAGS="-L../libdeflate -ldeflate"
    ./get_htslib.sh
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
cmake -DCMAKE_BUILD_TYPE=Release ..
make ${bc}
popd
cp ${DR}/${bc} ./
ln -fs ./$bc megadepth
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
