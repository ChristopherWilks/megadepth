#!/usr/bin/env bash

set -ex

VER=1.2.11
ar=zlib-${VER}.tar.gz
wget https://www.zlib.net/${ar}
gzip -dc ${ar} | tar xvf - 
rm -f ${ar}
pushd zlib-${VER}
./configure
make
popd
mv zlib-${VER} zlib
