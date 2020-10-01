module load cmake/3.17.0
module load git
module load python/3.7

module load gcc/9.3.0
# This is necessary for proper spawning through upcxx-run and may need to be in the login scripts
export LD_LIBRARY_PATH=/global/software/sl-7.x86_64/modules/langs/gcc/9.3.0/lib64:${LD_LIBRARY_PATH}

export GASNET_ODP_VERBOSE=0
export GASNET_PHYSMEM_MAX='10 GB'

export UPCXX_NETWORK=ibv
export MHM2_CMAKE_EXTRAS="-DCMAKE_CXX_COMPILER=$(which g++) -DCMAKE_C_COMPILER=$(which gcc)"

