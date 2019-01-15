#!/usr/bin/env bash

set -e

if [[ ! -d /code ]] ; then
    echo "Must mount bamcount directory at /code inside container"
fi 
cd /code

./get_htslib.sh
./get_libBigWig.sh

set -x

export LD_LIBRARY_PATH=/code/libBigWig:/code/htslib

DR=build-release-temp
mkdir -p ${DR}
pushd ${DR}
cmake -DCMAKE_BUILD_TYPE=Release ..
make
ldd -v bamcount
./bamcount --help
popd
cp ${DR}/bamcount ./bamcount
rm -rf ${DR}

DR=build-debug-temp
mkdir -p ${DR}
pushd ${DR}
cmake -DCMAKE_BUILD_TYPE=Debug ..
make
ldd -v bamcount
./bamcount --help
popd
cp ${DR}/bamcount ./bamcount-debug
rm -rf ${DR}

zip bamcount.zip bamcount bamcount-debug
rm -f bamcount bamcount-debug
