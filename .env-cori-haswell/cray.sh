module rm upcxx
module rm PrgEnv-gnu
module rm PrgEnv-intel
module load PrgEnv-cray
module rm craype-mic-knl
module load craype-haswell

module load cmake
module load git
module load upcxx

module list
