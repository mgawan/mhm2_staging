#pragma once

#ifdef __cplusplus
extern "C" {
#endif

extern const char *MHMXX_VERSION;
extern const char *MHMXX_VERSION_DATE;
extern const char *MHMXX_BUILD_DATE;

#define MAX_BUILD_KMER_STR "Maximum kmer len=MAX_BUILD_KMER"

#ifdef USE_CUDA
#define WITH_CUDA "CUDA + CPU SW"
#else
#define WITH_CUDA "CPU SW Only"
#endif



#ifdef __cplusplus
}
#endif


