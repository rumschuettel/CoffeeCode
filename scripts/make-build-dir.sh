#!/usr/bin/env bash

mkdir -p "$1"
cd "$1"

# sub-shells don't automatically have conda activated
# freaking new way of activating conda goes here
source /home/jkrb2/opt/anaconda5/etc/profile.d/conda.sh
conda activate gcc8

cmake -D CMAKE_BUILD_TYPE=Release -D OPTIMIZE_FOR_DEPOLARIZING=ON -D FLOATING_POINT_MULTIPLICITY=OFF ../..
