#!/usr/bin/env bash

set -ex

#macos:
#   x86_64-apple-darwin12 (CC=o64-gcc)
#OR
#   i386-apple-darwin1 (CC=o32-gcc)
#windows:
#   x86_64-w64-mingw32 (64bit)
#OR
#   i686-w64-mingw32 (32bit)

compiler=$1
#"macos" or "windows" (or nothing for normal linux build)
platform=$2

target_dir=htslib
if [[ -n $platform ]]; then
    target_dir="htslib_"${platform}
fi

#VER=1.9
VER=1.10.2
ar=htslib-${VER}.tar.bz2
if [[ ! -s $target_dir ]] ; then
    curl -OL https://github.com/samtools/htslib/releases/download/${VER}/${ar}
    bzip2 -dc ${ar} | tar xvf - 
    rm -f ${ar}
    mv htslib-${VER} $target_dir
fi
pushd $target_dir

#autoheader
#autoconf
make clean

if [[ -z $compiler ]]; then
    ./configure --disable-bz2 --disable-lzma --with-libdeflate
    make
else
    ./configure --disable-bz2 --disable-lzma --with-libdeflate --host=$compiler
    if [[ "$platform" == "macos" ]]; then
        #inherit CC, AR, and RANLIB from build_no_container_xcross.sh
        echo $CC
        echo $AR
        echo $RANLIB
        #only make static lib for cross-compilation for now
        #export CC=/opt/osxcross/target/bin/${compiler}-gcc
        #export AR=/opt/osxcross/target/bin/${compiler}-ar
        #export RANLIB=/opt/osxcross/target/bin/${compiler}-ranlib
    else # windows
        export CC=${compiler}-gcc
        export AR=${compiler}-ar
        export RANLIB=${compiler}-ranlib
    fi
    make CC=$CC AR=$AR RANLIB=$RANLIB libhts.a
fi
popd
