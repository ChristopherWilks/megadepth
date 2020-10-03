#docker run --ulimit core=-1 -it --volume `pwd`:/build:rw --entrypoint=/bin/bash a3939c35343f -i
docker run --ulimit core=-1 -it --volume `pwd`:/build:rw --entrypoint=/build/build_no_container_windows.sh a3939c35343f -i
