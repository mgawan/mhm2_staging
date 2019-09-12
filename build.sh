#!/bin/bash

set -e

cd `pwd`/src/build

if [ "$2" == "clean" ] || [ "$1" == "clean" ]; then
    cmake --build . --target clean
fi

if [ "$1" == "Debug" ]; then
    cmake . -DCMAKE_BUILD_TYPE=Debug
elif [ "$1" == "Release" ]; then
    cmake . -DCMAKE_BUILD_TYPE=Release
fi

make -j

