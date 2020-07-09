module use /gpfs/alpine/world-shared/csc296/summit/modulefiles


module rm xl
module load gcc/9.1.0
module load git/2.20.1
module load cmake/3.17.3

module load upcxx/2020.3.0

# cmake -DCMAKE_C_COMPILER=$(which mpicc) -DCMAKE_CXX_COMPILER=$(which mpicxx) 
