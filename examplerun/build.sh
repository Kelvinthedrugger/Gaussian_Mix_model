#!/bin/bash -e
g++ -std=c++17 -Wall -pedantic gsl_test.cpp -o gsl_test -lgsl -lgslcblas && ./gsl_test
