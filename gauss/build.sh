#!/bin/bash -e

# gcc -o Gaussian_poolOrNot Gaussian_poolOrNot.c GSLfun.c -lgsl -lgslcblas -lm && ./Gaussian_poolOrNot
# this if statement works
# if $1 < 1: compile
# else: don't compile
if [[ $1 -lt 1 ]]; then
echo "compiled";
gcc -o Gaussian_poolOrNot Gaussian_poolOrNot.c GSLfun.c -lgsl -lgslcblas -lm 
fi

./Gaussian_poolOrNot # "GSL_RNG_SEED" $1 


