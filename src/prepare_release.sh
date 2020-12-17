#!/bin/bash

USAGE="$0 mhm2-version-tag"

ver=$1
if [ -z "$ver" ]
then
  echo "$USAGE"
  exit 1
fi

set -e
set -x
pkg=$(mktemp -d)
if [ ! -d "$pkg" ]
then
  echo "Could not mktemp!"
  exit 1
fi

cd $pkg

git clone --recurse-submodules https://bitbucket.org/berkeleylab/mhm2.git mhm2-$ver
( cd mhm2-$ver && git checkout v${ver} && git submodule update && rm -rf .git)

tar -cvzf mhm2-$ver.tar.gz mhm2-$ver

cd -
mv $pkg/mhm2-$ver.tar.gz .
rm -rf $pkg



