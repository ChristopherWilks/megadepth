#!/usr/bin/env bash

set -ex

VER=1.2.13
if [[ -z $SUBMODULE ]]; then
    ar=zlib-${VER}.tar.gz
    wget https://www.zlib.net/${ar}
    gzip -dc ${ar} | tar xvf - 
    rm -f ${ar}
    pushd zlib-${VER}
else 
    pushd zlib
fi
./configure
make
popd
if [[ -z $SUBMODULE ]]; then
    mv zlib-${VER} zlib
fi
