#!/usr/bin/env bash

set -ex

VER=1.9
AR=htslib-${VER}.tar.bz2
wget https://github.com/samtools/htslib/releases/download/${VER}/${AR}
tar jxvf ${AR}
rm -f ${AR}
pushd htslib-${VER}
autoheader
autoconf
./configure --disable-bz2 --disable-lzma --disable-libcurl --without-libdeflate
make
popd
mv htslib-${VER} htslib
