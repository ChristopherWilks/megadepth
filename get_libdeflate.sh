#!/usr/bin/env bash

set -ex

compiler=$1
platform=$2

target_dir=libdeflate
if [[ -n $platform ]]; then
    target_dir="libdeflate_"${platform}
fi

VER=1.6
TARGZ=${VER}.tar.gz
FN=libdeflate-${TARGZ}
DIR=libdeflate-${VER}
curl -L https://github.com/ebiggers/libdeflate/archive/v${TARGZ} > $FN
tar -zxvf $FN
rm -f ${FN}
mv $DIR $target_dir
pushd $target_dir
target="libdeflate.a"
if [[ -z $compiler ]]; then
    #from https://github.com/samtools/htslib/issues/688
    make CFLAGS="$CFLAGS -fPIC -O3" $target
else
    if [[ "$platform" == "windows" ]]; then
        make CC=${compiler}-gcc CFLAGS='-O3' libdeflatestatic.lib
        ln -fs libdeflatestatic.lib libdeflate.a
    else
        make CC=${compiler}-gcc CFLAGS='-fPIC -O3' $target
    fi
fi
popd
