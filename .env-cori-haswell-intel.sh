module rm PrgEnv-cray
module rm PrgEnv-gnu
module load PrgEnv-intel
module rm craype-mic-knl
module load craype-haswell

module load cmake
module load git
module load upcxx/2020.3.0

module list
