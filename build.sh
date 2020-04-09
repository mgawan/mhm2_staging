#!/usr/bin/env bash

set -e

SECONDS=0

rootdir=`pwd`

INSTALL_PATH=${MHMXX_INSTALL_PATH:=$rootdir/install}

if [ "$1" == "clean" ]; then
    rm -rf build/*
    rm -rf $INSTALL_PATH
else
    cd $rootdir/build
    if [ "$1" == "Debug" ] || [ "$1" == "Release" ]; then
        rm -rf $INSTALL_PATH
        cmake $rootdir -DCMAKE_BUILD_TYPE=$1 -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH 
    fi
    make -j install
fi

echo "Build took $((SECONDS))s"
