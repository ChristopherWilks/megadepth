#!/usr/bin/env bash

set -ex

VER=0.4.4
AR=${VER}.tar.gz
FN=libBigWig-${AR}
DIR=libBigWig-${VER}
curl -L https://github.com/dpryan79/libBigWig/archive/${AR} > $FN
tar -zxvf $FN
rm -f ${FN}
pushd $DIR
cp Makefile Makefile.bak
cat Makefile.bak | perl -ne 'chomp; $f=$_; $f=~s/^(CFLAGS.+)$/$1 -DNOCURL/; $f=~s/^(LIBS =.*)(-lcurl)/$1/; print "$f\n";' > Makefile
make
popd
mv $DIR libBigWig
