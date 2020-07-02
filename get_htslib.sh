#!/usr/bin/env bash

set -ex

#e.g. i386-apple-darwin1 (CC=o32-gcc)
#or x86_64-apple-darwin12-gcc (CC=o64-gcc)
#CC must be set as well
xcross=$1

export CPPFLAGS="-I../libdeflate"
export LDFLAGS="-L../libdeflate -ldeflate"

VER=1.9
AR=htslib-${VER}.tar.bz2
if [[ ! -s htslib ]] ; then
    curl -OL https://github.com/samtools/htslib/releases/download/${VER}/${AR}
    bzip2 -dc ${AR} | tar xvf - 
    rm -f ${AR}
    pushd htslib-${VER}
    MOVE=1
else
    pushd htslib
fi
autoheader
autoconf
make clean
if [[ -z $xcross ]]; then
    ./configure --disable-bz2 --disable-lzma --disable-libcurl --with-libdeflate
    make
else
    ./configure --disable-bz2 --disable-lzma --disable-libcurl --with-libdeflate --host=$xcross
    #only make static lib for cross-compilation for now
    export AR=/opt/osxcross/target/bin/x86_64-apple-darwin12-ar
    export RANLIB=/opt/osxcross/target/bin/x86_64-apple-darwin12-ranlib
    cp Makefile Makefile.orig
    cat Makefile.orig | perl -ne 'chomp; $f=$_; $f=~s!^(AR\s*=\s*ar\s*)$!AR='$AR'!; $f=~s!^(RANLIB\s*=\s*ar\s*)$!AR='$RANLIB'!; print "$f\n";' > Makefile
    make libhts.a
fi

popd

if [[ -n "$MOVE" ]]; then
    mv htslib-${VER} htslib
fi
