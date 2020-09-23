module rm upcxx
module rm upcxx-gpu
module rm PrgEnv-intel
module rm PrgEnv-cray
module rm PrgEnv-gnu

module rm craype
module rm craype-mic-knl
module rm craype-haswell
module rm craype-x86-skylake

module load esslurm

module load PrgEnv-gnu
module load craype
module load craype-x86-skylake

module load cuda
module load cmake
module load git
module load upcxx-gpu

module list

export OMP_NUM_THREADS=1
export MHMXX_CMAKE_EXTRAS="-DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc"
