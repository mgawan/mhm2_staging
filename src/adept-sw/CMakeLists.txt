if(NOT ENABLE_CUDA)
    message(FATAL_ERROR "Trying to build ADEPT-SW but CUDA is not enabled")
endif()

find_package(OpenMP REQUIRED)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler=-Wall,${OpenMP_CXX_FLAGS} -gencode arch=compute_70,code=sm_70")
message(STATUS "CMAKE_CUDA_FLAGS=${CMAKE_CUDA_FLAGS}")

option(ADEPT_SW_SHARED "Adept-SW shared library" ON)
option(ADEPT_SW_STATIC "Adept-SW static library" ON)

if (ADEPT_SW_SHARED)
  add_library(ADEPT_SW_LIBRARY_shared SHARED driver.cpp kernel.cpp gpu_alns.cpp utils_gpu.cpp)
  target_link_libraries(ADEPT_SW_LIBRARY_shared INTERFACE OpenMP::OpenMP_CXX) 
  install(TARGETS ADEPT_SW_LIBRARY_shared LIBRARY DESTINATION lib)
  if (NOT ADEPT_SW_STATIC)
    add_library(ADEPT_SW_LIBRARY ALIAS ADEPT_SW_LIBRARY_shared)
  endif()
endif()
if (ADEPT_SW_STATIC)
  add_library(ADEPT_SW_LIBRARY_static STATIC driver.cpp kernel.cpp gpu_alns.cpp utils_gpu.cpp)
  target_link_libraries(ADEPT_SW_LIBRARY_static INTERFACE OpenMP::OpenMP_CXX)
  install(TARGETS ADEPT_SW_LIBRARY_static ARCHIVE DESTINATION lib)
  add_library(ADEPT_SW_LIBRARY ALIAS ADEPT_SW_LIBRARY_static)
endif()


set_source_files_properties(driver.cpp kernel.cpp gpu_alns.cpp utils_gpu.cpp PROPERTIES LANGUAGE CUDA LINKER_LANGUAGE CUDA)

