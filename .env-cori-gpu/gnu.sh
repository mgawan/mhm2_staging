module rm upcxx
module rm PrgEnv-intel
module rm PrgEnv-cray
module load PrgEnv-gnu

module rm craype-mic-knl
module load craype-haswell

module load cuda
module rm cmake
module load cmake
module load git
module load upcxx-gpu
#module load upcxx

module list

export OMP_NUM_THREADS=1
