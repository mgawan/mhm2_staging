module load craype
module rm upcxx
module rm upcxx-gpu
module rm PrgEnv-intel
module rm PrgEnv-cray
module load PrgEnv-gnu
module rm craype-haswell
module load craype-mic-knl

module load cmake
module load git
module load upcxx

module list
