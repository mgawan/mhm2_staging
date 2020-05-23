#!/usr/bin/env bash

if [ -n "$MHMXX_BUILD_ENV" ]; then
    source $MHMXX_BUILD_ENV
fi

set -e

SECONDS=0

rootdir=`pwd`

INSTALL_PATH=${MHMXX_INSTALL_PATH:=$rootdir/install}

rm -rf $INSTALL_PATH/bin/mhmxx

if [ "$1" == "clean" ]; then
    rm -rf .build/*
    # if this isn't removed then the the rebuild will not work
    rm -rf $INSTALL_PATH/cmake
    exit 0
else
    mkdir -p $rootdir/.build
    cd $rootdir/.build
    if [ "$1" == "Debug" ] || [ "$1" == "Release" ]; then
        rm -rf .build/*
        rm -rf $INSTALL_PATH/cmake
        cmake $rootdir -DCMAKE_BUILD_TYPE=$1 -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH
        #cmake $rootdir -DCMAKE_BUILD_TYPE=$1 -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH -GNinja
    fi
    make -j install
    # ninja is about 20% faster on my 80 core server, but is not as widely supported (eg. not on Cori)
    #ninja -j 0
fi

if [ -n "$MHMXX_BUILD_ENV" ]; then
    env_id=`echo $MHMXX_BUILD_ENV|cut -d '.' -f2|cut -d '-' -f2-`
    cd $INSTALL_PATH/bin
    rm -f mhmxx.$env_id
    mv mhmxx mhmxx.$env_id
    ln -s mhmxx.$env_id mhmxx
fi

echo "Build took $((SECONDS))s"
