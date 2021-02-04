#!/bin/bash

UPCXXVER=${UPCXXVER:=2020.10.0}

getcores()
{
  if lscpu
  then
     :
  fi 2>/dev/null | awk '/^CPU\(s\):/ {print $2}'
  if sysctl -a
  then
     :
  fi 2>/dev/null | awk '/^machdep.cpu.core_count/ {print $2}'
}
BUILD_THREADS=${BUILD_THREADS:=$(getcores)}

conduit=$1
shift
if [ -n "$1" ]
then
  shm=$1
else
  shm=posix
fi

if [ -n "$2" ]
then
  installdir=${2}
else
  installdir=$HOME/install
fi

if [ -n "$3" ]
then
  builddir=${3}
else
  builddir=${TMPDIR:=/dev/shm}
fi
[ -d "$builddir" ] && [ -w "$builddir" ] || builddir=/tmp

if [ -n "$4" ]
then
  codedir=${4}
else
  codedir=$HOME
fi

CC=${CC:=$(which gcc)}
CXX=${CXX:=$(which g++)}

[ -d $builddir ] || builddir=/tmp

USAGE="$0 Conduit(smp|ibv|ofi|aries|mpi|udp) [ SHARED_MEMORY(posix|xpmem|sysv|file|none) [ INSTALL_DIR($installdir) [ BUILD_DIR($builddir) [ CODE_DIR(${codedir}) ] ] ] ]"

if [ -z "$conduit" ] 
then
  echo $USAGE
  echo "You must choose a default networking conduit: smp, ibv, ofi, aries, mpi or udp"
  echo
  exit 0
fi

oops()
{
  echo "uh oh, something bad happened!"
  exit 1
}

echo "Building bupc and upcxx, conduit=$conduit shared_memory=$shm installdir=$installdir builddir=$builddir codedir=$codedir"

trap oops 0

set -e
set -x


cd $codedir
builddir=$builddir/bupc-$USER-hipmer-builds
rm -rf ${builddir}
mkdir -p ${builddir}

if [ ! -x $installdir/bin/upcxx ]
then

  upcxx_opts=${upcxx_opts:=}
  echo "Building UPC++ -- ${upcxx_opts}"
  # build upcxx
  cd $codedir

  UPCXXDIR=upcxx-${UPCXXVER}
  UPCXXTAR=${UPCXXDIR}.tar.gz
  UPCXXURLBASE=https://bitbucket.org/berkeleylab/upcxx/downloads
  if [ "$UPCXXVER" == "latest" ]
  then
    UPCXXURLBASE=https://upcxx.lbl.gov/third-party/hipmer
  fi


  rm -f $UPCXXTAR
  [ -f $UPCXXTAR ] || curl -LO $UPCXXURLBASE/${UPCXXTAR}
  cd $builddir
  [ -d ${UPCXXDIR} ] || tar -xvzf $codedir/${UPCXXTAR}
  [ -d ${UPCXXDIR} ] || mv upcxx-*/ ${UPCXXDIR}
  cd ${UPCXXDIR}
set -x
if [ -z "${with_pmirun_cmd}" ]
then
  CC=$CC CXX=$CXX ./configure --prefix=$installdir ${upcxx_opts} && (make -j ${BUILD_THREADS} || make) && make install
else
  CC=$CC CXX=$CXX ./configure --prefix=$installdir ${upcxx_opts} --with-pmirun-cmd="${with_pmirun_cmd}"  && (make -j ${BUILD_THREADS} || make) && make install
fi
#  ./install $installdir
else
  echo "upcxx is already installed $installdir/bin/upcxx"
fi

set -x
trap "" 0

