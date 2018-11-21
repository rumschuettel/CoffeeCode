#!/usr/bin/env bash

HOSTNAME=`hostname`
mkdir -p "build/$HOSTNAME"

for ((K_SYS=1; K_SYS<=20; K_SYS++))
do
    for ((K_ENV=1; K_ENV<=$K_SYS; K_ENV++))
    do
        cmd="g++ -DNDEBUG -std=c++14 -DK_SYS=$K_SYS -DK_ENV=$K_ENV -O3 -march=native CoffeeCode.cpp -o build/$HOSTNAME/CoffeeCode.$K_SYS.$K_ENV.out"
        echo $cmd
        $cmd
    done
done