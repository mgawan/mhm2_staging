module rm upcxx
module rm PrgEnv-intel
module rm PrgEnv-cray
module load PrgEnv-gnu

module rm craype-mic-knl
module rm craype-haswell

# if this is loaded then the adept-sw library ends up static which causes the link to fail
#module rm craype

module load esslurm
module load cuda
module load cmake
module load git
module load upcxx-gpu

module list

export OMP_NUM_THREADS=1
export MHMXX_CMAKE_EXTRAS="-DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc"
