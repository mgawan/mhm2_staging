module load craype
module rm upcxx
module rm upcxx-gpu
module rm PrgEnv-cray
module rm PrgEnv-intel
module load PrgEnv-gnu
module rm craype-mic-knl
module load craype-haswell

module load cmake
module load git
module load upcxx

module list
