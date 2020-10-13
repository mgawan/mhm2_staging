#!/bin/bash

USAGE="USAGE: $0 [MHM2-tag/version/branch ]
   Execute *either* with a valid git tag/version/branch with Dockerfile in 
            or in the top level of mhm2 source to build a development version

   Optionally set the env variable DOCKER_FILE to use a different build than the default
"

MHM2_VER=${MHM2_VER:=$1}
DOCKER_FILE=${DOCKER_FILE:=Dockerfile}
BUILD_ARG="--build-arg=MHM2_VER=${MHM2_VER}"
if [ -z "$MHM2_VER" ]
then
  if [ -f VERSION -a -d .git -a -f src/mhm2.py -a -f upcxx-utils/VERSION ]
  then
    # top level mhm2 directory okay
    if [ ! -f ${DOCKER_FILE} ] 
    then
      DOCKER_FILE=contrib/docker/Dockerfile
    fi
    MHM2_VER=$(git describe --tags)
    BUILD_ARG=""
    echo "Building using the current source ($MHM2_VER)"
  else 
    echo "$USAGE"
    exit 1
  fi
else
  echo "Building using the tagged version $MHM2_VER"
fi

docker build -f ${DOCKER_FILE} ${BUILD_ARG} --tag=mhm2:${MHM2_VER} .


