#!/bin/bash

set -e

srcdir=`pwd`/src
echo $srcdir

if [ "$1" == "clean" ]; then
    #cmake --build $srcdir --target clean
    rm -rf build/*
else
    cd `pwd`/build
    if [ "$1" == "Debug" ]; then
        cmake $srcdir -DCMAKE_BUILD_TYPE=Debug
    elif [ "$1" == "Release" ]; then
        cmake $srcdir -DCMAKE_BUILD_TYPE=Release
    fi
    make -j
fi

