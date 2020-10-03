#!/usr/bin/env bash
#build "statlib" version (or main linux version) of megadepth
./run_hbb_build_container.sh
#build latest Docker runner image:
/bin/bash -x create_docker_to_run_megadepth.sh
