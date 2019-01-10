#!/usr/bin/env bash

set -ex

VER=0.4.2
AR=${VER}.tar.gz
FN=libBigWig-${AR}
DIR=libBigWig-${VER}
curl -L https://github.com/dpryan79/libBigWig/archive/${AR} > $FN
tar -zxvf $FN
rm -f ${FN}
pushd $DIR
make
popd
mv $DIR libBigWig
