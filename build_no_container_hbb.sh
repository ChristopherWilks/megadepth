#!/usr/bin/env bash

set -ex
yum install cmake -yy

working_dir=$(dirname $0)
pushd $working_dir

#clear symlink main lib dirs
rm -rf zlib htslib libBigWig libdeflate build-release-temp

export PATH=/opt/rh/devtoolset-8/root/usr/bin:$PATH
export CFLAGS="-g -O2 -fvisibility=hidden -DCURL_STATICLIB -fPIC"
export CPPFLAGS="-g -O2 -fvisibility=hidden -I/hbb_shlib/include"
export LDFLAGS="-L/hbb_shlib/lib -static-libstdc++"
export SHLIB_LDFLAGS="-static-libstdc++"
export STATICLIB_CFLAGS="-g -O2 -fvisibility=hidden -fPIC"
export STATICLIB_CPPFLAGS="-g -O2 -fvisibility=hidden -I/hbb_shlib/include"

build_type=$1
bc=`perl -e '$bt="'$build_type'"; if($bt=~/static/i) { print "megadepth_static"; } else { print "megadepth"; }'`

if [[ "$bc" == "megadepth_static" ]]; then
    ln -fs CMakeLists.txt.static CMakeLists.txt
else
    ln -fs CMakeLists.txt.statlib CMakeLists.txt
fi

#dont need our own zlib, since it's already statically compiled in HBB

if [[ ! -s libdeflate_hbb ]] ; then
    ./get_libdeflate.sh
    mv libdeflate libdeflate_hbb
fi
ln -fs libdeflate_hbb libdeflate

if [[ ! -s htslib_hbb ]] ; then
    export CPPFLAGS="$CPPFLAGS -I../libdeflate"
    #for staticlly linking libcurl
    export LDFLAGS="-static-libstdc++ -L/hbb_shlib/lib -L../libdeflate"
    export LIBS="-lm -ldl -lssl -lcrypto -lz -ldeflate -lrt -pthread"
    ./get_htslib.sh linux hbb
    export CPPFLAGS="-g -O2 -fvisibility=hidden -I/hbb_shlib/include"
    export LDFLAGS="-L/hbb_shlib/lib -static-libstdc++"
    export LIBS=
    #reset env vars
    mv htslib htslib_hbb
fi
ln -fs htslib_hbb htslib

if [[ ! -s libBigWig_hbb ]] ; then
    ./get_libBigWig.sh
    pushd libBigWig
    export CFLAGS="-O2 -I/hbb_shlib/include -DCURL_STATICLIB -fPIC"
    make clean
    make -f Makefile.fpic lib-static
    popd
    mv libBigWig libBigWig_hbb
fi
ln -fs libBigWig_hbb libBigWig
export CFLAGS="-g -O2 -fvisibility=hidden -I/hbb_shlib/include -DCURL_STATICLIB -fPIC"
export LDFLAGS="-L/hbb_shlib/lib -static-libstdc++"

set -x

export LD_LIBRARY_PATH=./htslib:./libBigWig:$LD_LIBRARY_PATH

DR=build-release-temp
mkdir -p ${DR}
pushd ${DR}
cmake -DCMAKE_BUILD_TYPE=Release ..
make ${bc}
popd
cp ${DR}/${bc} ./
#ln -fs ./$bc megadepth
./megadepth --version
rm -rf ${DR}
mv megadepth megadepth.full
strip -s megadepth.full -o megadepth

DR=build-debug-temp
mkdir -p ${DR}
pushd ${DR}
cmake -DCMAKE_BUILD_TYPE=Debug ..
make ${bc}
popd
cp ${DR}/${bc} ./${bc}_debug
#ln -fs ./${bc}_debug megadepth_debug
./megadepth_debug --version
rm -rf ${DR}
