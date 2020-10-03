#docker run --ulimit core=-1 -it --volume `pwd`:/build:rw --entrypoint=/bin/bash osxcross_gcc -i
docker run --ulimit core=-1 -it --volume `pwd`:/build:rw --entrypoint=/build/build_no_container_macos.sh osxcross_gcc -i
