#!/bin/bash -e

#GSL_RNG_SEED=$1

# e.g. ./run.sh 1337 1 -> don't compile
# ... sh 0 -> compile
cd gauss
GSL_RNG_SEED=$1 ./build.sh $2

# cd gauss_pool_or_not/
# 	./build.sh $1 $2


