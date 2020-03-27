#!/usr/bin/env bash

set -e

rootdir=`pwd`
srcdir=$rootdir/src

INSTALL_PATH=${MHMXX_INSTALL_PATH:=$rootdir/bin}

if [ "$1" == "clean" ]; then
    rm -rf build/*
else
    cd $rootdir/build
    if [ "$1" == "Debug" ]; then
        cmake $srcdir -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH
    elif [ "$1" == "Release" ]; then
        cmake $srcdir -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH
    fi
    make -j install
fi

