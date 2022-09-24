#for interactive debugging switch to just bash shell
#older image ids: de2c52df2597, 555a36ef79dd
docker run --ulimit core=-1 -it --rm --volume `pwd`:/build:rw --entrypoint=/build/build_no_container_hbb.sh docker.pkg.github.com/phusion/holy-build-box/hbb-64cmake:latest -i
mv -f megadepth megadepth_hbb
ln -fs megadepth_hbb megadepth
