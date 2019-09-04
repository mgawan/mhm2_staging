#!/bin/bash

set -e

cd `pwd`/src/build

build_type=`cmake -L . | grep CMAKE_BUILD_TYPE`

if [ "$1" == "Debug" ] && [ "$build_type" == "CMAKE_BUILD_TYPE:STRING=Release" ]; then
    cmake . -DCMAKE_BUILD_TYPE=Debug
elif [ "$1" == "Release" ] && [ "$build_type" == "CMAKE_BUILD_TYPE:STRING=Debug" ]; then
    cmake . -DCMAKE_BUILD_TYPE=Release
fi

if [ "$2" == "clean" ]; then
    cmake --build . --target clean
fi

cmake --build . -- -j

