#!/bin/bash

# Compiles the sub-projects.
# Needs the following packages to be installed: git make build-essential gcc g++ cmake pkg-config

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# update submodules
echo "Updating submodules ..."
git submodule update --init --recursive

# algorithms
cd $DIR/algorithms

# clique algorithms
cd clique

# lmc
echo "Compiling LMC ..."
cd lmc
make all
cd ..

# end of clique algorithms
cd ..

# graph algorithms
cd graph

# scamp
echo "Compiling SCAMP ..."
cd scamp
mkdir build
cd build
cmake ..
cmake --build . --config Release
cd ../..

# end of graph algorithms

# tsgenerator
echo "Compiling tsgenerator-dev ..."
cd $DIR/tsgenerator
cd tsgenerator-dev
mkdir build
cd build
cmake ..
make
echo "Creating tsgenerator-dev package ..."
make package
echo "Installing tsgenerator-dev ..."
sudo apt install ./tsgenerator-dev-3.0.0.deb
sudo ldconfig
echo "Compiling tsgenerator ..."
cd ../../tsgenerator
mkdir build
cd build
cmake ..
make

