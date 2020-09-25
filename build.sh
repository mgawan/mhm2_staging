#!/usr/bin/env bash

if [ -n "$MHM2_BUILD_ENV" ]; then
    source $MHM2_BUILD_ENV
fi

set -e

SECONDS=0

rootdir=`pwd`

INSTALL_PATH=${MHM2_INSTALL_PATH:=$rootdir/install}

rm -rf $INSTALL_PATH/bin/mhm2

if [ "$1" == "clean" ]; then
    rm -rf .build/*
    # if this isn't removed then the the rebuild will not work
    rm -rf $INSTALL_PATH/cmake
    exit 0
else
    mkdir -p $rootdir/.build
    cd $rootdir/.build
    if [ "$1" == "Debug" ] || [ "$1" == "Release" ] || [ "$1" == "RelWithDebInfo" ]; then
        rm -rf *
        rm -rf $INSTALL_PATH/cmake
        cmake $rootdir -DCMAKE_BUILD_TYPE=$1 -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH $MHM2_CMAKE_EXTRAS
    fi
    make -j 12 install
fi

echo "Build took $((SECONDS))s"
