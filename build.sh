#!/usr/bin/env bash

set -e

rootdir=`pwd`
srcdir=$rootdir/src
echo $srcdir

if [ "$1" == "clean" ]; then
    #cmake --build $srcdir --target clean
    rm -rf build/*
else
#    if hash module 2>/dev/null; then
#        module load upcxx
#    fi
    cd $rootdir/build
    if [ "$1" == "Debug" ]; then
        cmake $srcdir -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=$rootdir/bin
    elif [ "$1" == "Release" ]; then
        cmake $srcdir -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$rootdir/bin
    fi
    make -j install
    # ensure python wrapper is in build
    #cp -v $srcdir/mhmxx.py $srcdir/../build
fi

