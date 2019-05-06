#!/usr/bin/env bash

mkdir -p "$1"
cd "$1"

# load local configuration file not checked into git
THISDIR="${BASH_SOURCE%/*}/"
source "$THISDIR/make-build-dir.sh.local"

cmake -D CMAKE_C_COMPILER=$CC -D CMAKE_CXX_COMPILER=$CXX -D CMAKE_BUILD_TYPE=Release -D PARALLELIZE=ON -D FLOATING_POINT_MULTIPLICITY=OFF ../..
