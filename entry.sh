#!/usr/bin/env bash

set -ex

if [[ ! -d /code ]] ; then
    echo "Must mount bamcount directory at /code inside container"
fi 
cd /code

DR=build-release-temp
mkdir -p ${DR}
pushd ${DR}
cmake -DCMAKE_BUILD_TYPE=Release ..
make
ldd -v bamcount
./bamcount
popd
rm -rf ${DR}

DR=build-debug-temp
mkdir -p ${DR}
pushd ${DR}
cmake -DCMAKE_BUILD_TYPE=Debug ..
make
ldd -v bamcount
./bamcount
popd
rm -rf ${DR}

cp build-release/bamcount ./bamcount
cp build-debug/bamcount ./bamcount-debug
zip bamcount.zip bamcount bamcount-debug
rm -f bamcount bamcount-debug
