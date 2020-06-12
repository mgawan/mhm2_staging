module rm upcxx
module rm PrgEnv-intel
module rm PrgEnv-cray
module load PrgEnv-gnu
module rm craype-haswell
module load craype-mic-knl

module load cmake
module load git
module load upcxx/2020.3.0

module list
