#!/usr/bin/env bash

set -e

if [[ ! -d /code ]] ; then
    echo "Must mount megadepth directory at /code inside container"
fi 
cd /code

rm -rf htslib libBigWig
./get_htslib.sh
./get_libBigWig.sh

set -x

DR=build-release-temp
mkdir -p ${DR}
pushd ${DR}
cmake -DCMAKE_BUILD_TYPE=Release ..
make
ldd -v megadepth
./megadepth --help
popd
cp ${DR}/megadepth ./megadepth
rm -rf ${DR}

DR=build-debug-temp
mkdir -p ${DR}
pushd ${DR}
cmake -DCMAKE_BUILD_TYPE=Debug ..
make
ldd -v megadepth
./megadepth --help
popd
cp ${DR}/megadepth ./megadepth-debug
rm -rf ${DR}

zip megadepth.zip megadepth megadepth-debug
#rm -f megadepth megadepth-debug
