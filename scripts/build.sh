#!/usr/bin/env bash
set -e
shopt -s expand_aliases
alias time='date; time'

scriptdir=$(cd $(dirname $0); pwd -P)
sourcedir=$(cd $scriptdir/..; pwd -P)

cmake -S . -B build # mkdir build && cd build && cmake ..
cmake --build build
(time ./build/src/main) |& tee results.txt
