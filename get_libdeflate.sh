#!/usr/bin/env bash

set -ex

VER=1.6
AR=v${VER}.tar.gz
FN=libdeflate-${AR}
DIR=libdeflate-${VER}
curl -L https://github.com/ebiggers/libdeflate/archive/${AR} > $FN
tar -zxvf $FN
rm -f ${FN}
pushd $DIR
#from https://github.com/samtools/htslib/issues/688
make CFLAGS='-fPIC -O3' libdeflate.a
popd
mv $DIR libdeflate
