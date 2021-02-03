module rm upcxx
module rm upcxx-gpu
module rm PrgEnv-intel
module rm PrgEnv-cray
module rm PrgEnv-gnu

module rm craype
module rm craype-mic-knl
module rm craype-haswell
module rm craype-x86-skylake

module load cmem

module load PrgEnv-gnu
module load craype
#module load craype-x86-skylake

module load cmake/3.18.2
module load git

# use custom build of upcxx
# upcxx_opts="--enable-ibv-multirail --with-ibv-max-hcas=4 --without-cross" ./install_upcxx.sh ibv posix /global/common/software/m2865/upcxx-2020.10.0-ibv-cmem
export PATH=/global/common/software/m2865/upcxx-2020.10.0-ibv-cmem/bin:$PATH

# use ibv and do not probe for shared heap size
export UPCXX_NETWORK=ibv
export GASNET_PHYSMEM_MAX=16G
export GASNET_PHYSMEM_PROBE=0
export GASNET_TMPDIR=$SCRATCH/tmp
mkdir -p $GASNET_TMPDIR

module list

export MHM2_CMAKE_EXTRAS="-DCMAKE_CXX_COMPILER=$(which g++) -DCMAKE_C_COMPILER=$(which gcc) "

