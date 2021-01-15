#!/usr/bin/env bash

set -ex

platform=$1

target_dir=libBigWig

#in case we're not using git submodules to pull dependencies
if [[ -z $SUBMODULE ]]; then
    if [[ -n $platform ]]; then
        target_dir="libBigWig_"${platform}
    fi
    VER=0.4.4
    TARGZ=${VER}.tar.gz
    FN=libBigWig-${TARGZ}
    DIR=libBigWig-${VER}
    curl -L https://github.com/dpryan79/libBigWig/archive/${TARGZ} > $FN
    tar -zxvf $FN
    rm -f ${FN}
    mv $DIR $target_dir
fi

pushd $target_dir

cp Makefile Makefile.orig
cat Makefile.orig | perl -ne 'chomp; $f=$_; $f=~s/^(CFLAGS.+)$/$1 -DNOCURL/; $f=~s/^(LIBS =.*)(-lcurl)/$1/; $f=~s/^(LDFLAGS.+)=/$1\?=/; print "$f\n";' > Makefile.nocurl
cat Makefile.orig | perl -ne 'chomp; $f=$_; $f=~s/^(CFLAGS \?= )$/$1 -fPIC /; $f=~s/^(LDFLAGS.+)=/$1\?=/; print "$f\n";' > Makefile.fpic
popd
