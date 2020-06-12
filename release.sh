#!/usr/bin/env bash
ver=`cat VERSION`
rm -rf megadepth-${ver} *.zip megadepth_*
mkdir megadepth-${ver}
rsync -av README.md *.sh Dockerfile CMakeLists.txt megadepth.cpp megadepth-${ver}/
for f in static dynamic statlib; do
    /bin/bash -x build_no_container.sh $f
done
ln -fs megadepth_statlib megadepth
ln -fs megadepth_statlib_debug megadepth_debug
zip megadepth-${ver}.zip megadepth megadepth_debug megadepth_dynamic megadepth_dynamic_debug megadepth_static megadepth_static_debug  megadepth_statlib megadepth_statlib_debug
rm -rf megadepth-${ver}/megadepth-${ver}
zip -r ${ver}.zip megadepth-${ver}
tar -zcvf ${ver}.tar.gz megadepth-${ver}
