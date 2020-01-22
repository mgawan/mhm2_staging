#!/usr/bin/env bash

set -e

srcdir=`pwd`/src
echo $srcdir

if [ "$1" == "clean" ]; then
    #cmake --build $srcdir --target clean
    rm -rf build/*
else
    if hash module 2>/dev/null; then
        module load upcxx
    fi
    cd `pwd`/build
    if [ "$1" == "Debug" ]; then
        cmake $srcdir -DCMAKE_BUILD_TYPE=Debug
    elif [ "$1" == "Release" ]; then
        cmake $srcdir -DCMAKE_BUILD_TYPE=Release
    fi
    make -j
fi

