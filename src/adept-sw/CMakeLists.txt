add_library(ADEPT_SW_LIBRARY driver.cpp kernel.cpp gpu_alns.cpp utils_gpu.cpp)
set_source_files_properties(driver.cpp kernel.cpp gpu_alns.cpp utils_gpu.cpp PROPERTIES LANGUAGE CUDA)
#target_compile_options(ADEPT_SW_LIBRARY PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-arch=compute_70 -code=compute_70>)

