#!/usr/bin/env bash

set -e

export PATH=/opt/rh/devtoolset-8/root/usr/bin:$PATH
export CFLAGS="-g -O2 -fvisibility=hidden -I/hbb_shlib/include -DCURL_STATICLIB -fPIC"
export LDFLAGS="-L/hbb_shlib/lib -static-libstdc++"
export STATICLIB_CFLAGS="-g -O2 -fvisibility=hidden -I/hbb_shlib/include -fPIC"

#build dynamic by default
build_type=$1
bc=`perl -e '$bt="'$build_type'"; if($bt=~/static/i) { print "megadepth_static"; } elsif($bt=~/statlib/i) { print "megadepth_statlib2"; } else { print "megadepth_dynamic"; }'`

#if [[ ! -s zlib ]] ; then
#    ./get_zlib.sh
#fi

if [[ ! -s libdeflate ]] ; then
    ./get_libdeflate.sh
fi

if [[ ! -s htslib ]] ; then
    export CPPFLAGS="$CPPFLAGS -I../libdeflate"
    export LDFLAGS="$LDFLAGS -L../libdeflate -ldeflate"
    ./get_htslib.sh
    #export CPPFLAGS=
    #export LDFLAGS=
fi

if [[ ! -s libBigWig ]] ; then
    ./get_libBigWig.sh
    pushd libBigWig
    #export CFLAGS="-I../libcurl_windows/include -I../zlib_windows -g -Wall -O3 -Wsign-compare -DCURL_STATICLIB"
    #export CFLAGS="$CFLAGS -DCURL_STATICLIB"
    make clean
    make -f Makefile.fpic lib-static
    popd
fi

set -x

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
