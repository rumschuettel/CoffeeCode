#!/usr/bin/env bash

mkdir -p "$1"
cd "$1"

# load local configuration file not checked into git
THISDIR="${BASH_SOURCE%/*}/"
source "$THISDIR/make-build-dir.sh.local"

cmake -D CMAKE_BUILD_TYPE=Release -D PARALLELIZE=ON -D FLOATING_POINT_MULTIPLICITY=OFF -D PATH_BOOST=/home/johannes_bausch_gmail_com/boost_1_69_0 -D PATH_NAUTY=/home/johannes_bausch_gmail_com/nauty27rc2 -D SYMMETRIC_SOLVER=ON -D REDUCE_LAMBDA_IF_POSSIBLE=ON ../..
