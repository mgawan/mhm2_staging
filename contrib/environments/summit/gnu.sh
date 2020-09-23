module use /gpfs/alpine/world-shared/csc296/summit/modulefiles


module rm xl
module load gcc/7.4.0
module load cuda/10.1.243
module load git/2.20.1
module load cmake/3.17.3

module load upcxx-cuda/2020.3.2

MHM2_CMAKE_EXTRAS="-DCMAKE_C_COMPILER=$(which mpicc) -DCMAKE_CXX_COMPILER=$(which mpicxx) -DCMAKE_CUDA_COMPILER=$(which nvcc)"
