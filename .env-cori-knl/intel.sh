module rm upcxx
module rm PrgEnv-gnu
module rm PrgEnv-cray
module load PrgEnv-intel
module rm craype-haswell
module load craype-mic-knl

module load cmake
module load git
module load upcxx

module list
