#!/usr/bin/env bash

set -ex

VER=1.2.11
AR=zlib-${VER}.tar.gz
wget https://www.zlib.net/${AR}
gzip -dc ${AR} | tar xvf - 
rm -f ${AR}
pushd zlib-${VER}
./configure
make
popd
mv zlib-${VER} zlib
