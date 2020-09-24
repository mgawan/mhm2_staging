module rm upcxx
module rm PrgEnv-intel
module rm PrgEnv-gnu
module load PrgEnv-cray
module rm craype-haswell
module load craype-mic-knl

module load cmake
module load git
module load upcxx

module list
