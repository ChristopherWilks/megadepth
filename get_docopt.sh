#!/usr/bin/env bash

# bwam.cpp uses docopt and CMakeLists.txt expects to find its .h .cpp
# file in the docopt subdir

set -ex

git clone https://github.com/docopt/docopt.cpp.git
mv docopt.cpp docopt
