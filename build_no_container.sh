#!/usr/bin/env bash

set -e

#build dynamic by default
build_type=$1
bc=`perl -e '$bt="'$build_type'"; if($bt=~/static/i) { print "bamcount"; } else { print "bamcount-dyn"; }'`

if [[ ! -s htslib ]] ; then
    ./get_htslib.sh
fi

if [[ ! -s libBigWig ]] ; then
    ./get_libBigWig.sh
fi

set -x

export LD_LIBRARY_PATH=./htslib:./libBigWig:$LD_LIBRARY_PATH

DR=build-release-temp
mkdir -p ${DR}
pushd ${DR}
cmake -DCMAKE_BUILD_TYPE=Release ..
make ${bc}
popd
cp ${DR}/${bc} ./bamcount
./bamcount --version
rm -rf ${DR}

DR=build-debug-temp
mkdir -p ${DR}
pushd ${DR}
cmake -DCMAKE_BUILD_TYPE=Debug ..
make ${bc}
popd
cp ${DR}/${bc} ./bamcount-debug
./bamcount-debug --version
rm -rf ${DR}
