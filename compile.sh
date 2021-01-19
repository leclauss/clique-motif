#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# update submodules
git submodule update --init --recursive

# algorithms
cd $DIR/algorithms

# clique algorithms
cd clique

# lmc
cd lmc
make all
cd ..

# end of clique algorithms
cd ..

# graph algorithms
cd graph

# scamp
cd scamp
mkdir build
cd build
cmake ..
cmake --build . --config Release
cd ../..

# end of graph algorithms
cd ..

