#!/usr/bin/env bash

if [ -n "$MHMXX_BUILD_ENV" ]; then
    source $MHMXX_BUILD_ENV
fi

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
        rm -rf build/*
        rm -rf $INSTALL_PATH
        cmake $rootdir -DCMAKE_BUILD_TYPE=$1 -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH 
    fi
    make -j install
fi

if [ -n "$MHMXX_BUILD_ENV" ]; then
    env_id=`echo $MHMXX_BUILD_ENV|cut -d '.' -f2|cut -d '-' -f2-`
    cd $INSTALL_PATH/bin
    mv mhmxx mhmxx.$env_id
    ln -s mhmxx.$env_id mhmxx
fi

echo "Build took $((SECONDS))s"
