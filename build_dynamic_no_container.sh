#!/usr/bin/env bash
#only builds the dynamically linked version, w/o using a Docker container

set -e

set -x

DR=build-release-temp
mkdir -p ${DR}
pushd ${DR}
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON ..
make bamcount-dyn
ldd -v bamcount-dyn
#./bamcount --help
popd
cp ${DR}/bamcount-dyn ./bamcount
#rm -rf ${DR}

DR=build-debug-temp
mkdir -p ${DR}
pushd ${DR}
cmake -DCMAKE_BUILD_TYPE=Debug ..
make bamcount-dyn
ldd -v bamcount-dyn
#./bamcount --help
popd
cp ${DR}/bamcount-dyn ./bamcount-debug
#rm -rf ${DR}
