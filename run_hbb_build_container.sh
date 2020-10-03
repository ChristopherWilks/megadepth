#for interactive debugging switch to just bash shell
#docker run --ulimit core=-1 -it --rm --volume `pwd`:/build:rw --entrypoint=/bin/bash de2c52df2597 -i
docker run --ulimit core=-1 -it --rm --volume `pwd`:/build:rw --entrypoint=/build/build_no_container_hbb.sh de2c52df2597 -i
