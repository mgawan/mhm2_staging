if(NOT ENABLE_CUDA)
    message(FATAL_ERROR "Trying to build GPU_LOC_ASSEM but CUDA is not enabled")
endif()

set(CMAKE_CUDA_STANDARD 14 CACHE STRING "")
set(CMAKE_CUDA_STANDARD_REQUIRED ON CACHE BOOL "")
set(CMAKE_CUDA_EXTENSIONS OFF CACHE BOOL "")
enable_language(CUDA)

find_package(CUDA REQUIRED)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)

if (${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.18)
  set(CMAKE_CUDA_ARCHITECTURES 70)
else()
  # CUDA_SELECT_NVCC_ARCH_FLAGS(ARCH_FLAGS_LIST)
  if(ARCH_FLAGS_LIST)
    STRING(REPLACE ";" " " ARCH_FLAGS "${ARCH_FLAGS_LIST}")
  else()
    foreach(_CMAKE_CUDA_COMPUTE_CAPABILITY 70)
      set(ARCH_FLAGS "${ARCH_FLAGS} -gencode arch=compute_${_CMAKE_CUDA_COMPUTE_CAPABILITY},code=sm_${_CMAKE_CUDA_COMPUTE_CAPABILITY}")
    endforeach()
    if(NOT CMAKE_CUDA_COMPUTE_CAPABILITY)
      set(CMAKE_CUDA_COMPUTE_CAPABILITY 70 CACHE STRING "")
    endif()
    message(WARNING "Could not auto-detect the GPU arch flags building for all gpu compute capabilities 3 - 7 and PTX for ${CMAKE_CUDA_COMPUTE_CAPABILITY}")
    set(ARCH_FLAGS "${ARCH_FLAGS} -gencode arch=compute_${CMAKE_CUDA_COMPUTE_CAPABILITY},code=sm_${CMAKE_CUDA_COMPUTE_CAPABILITY}")
    set(ARCH_FLAGS "${ARCH_FLAGS} -gencode arch=compute_${CMAKE_CUDA_COMPUTE_CAPABILITY},code=compute_${CMAKE_CUDA_COMPUTE_CAPABILITY}")
  endif()
endif()

message(STATUS "ARCH_FLAGS=${ARCH_FLAGS}")


set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler=-Wall ${ARCH_FLAGS}")
message(STATUS "CMAKE_CUDA_FLAGS=${CMAKE_CUDA_FLAGS}")

add_library(GPU_LOC_obj OBJECT kernel.cpp driver.cpp)

set_source_files_properties(kernel.cpp driver.cpp PROPERTIES LANGUAGE CUDA LINKER_LANGUAGE CUDA)
option(GPU_LOC_SHARED "GPU_LOC shared library" OFF)
option(GPU_LOC_STATIC "GPU_LOC static library" ON)

if (GPU_LOC_SHARED)
  add_library(GPU_LOC_SHARED SHARED $<TARGET_OBJECTS:GPU_LOC_obj>)
  set_property(TARGET GPU_LOC_SHARED PROPERTY CUDA_RESOLVE_DEVICE_SYMBOLS ON)
  install(TARGETS GPU_LOC_SHARED LIBRARY DESTINATION lib)
endif()
if (GPU_LOC_STATIC)
  add_library(GPU_LOC_STATIC STATIC $<TARGET_OBJECTS:GPU_LOC_obj>)
  set_property(TARGET GPU_LOC_STATIC PROPERTY CUDA_RESOLVE_DEVICE_SYMBOLS ON)
  install(TARGETS GPU_LOC_STATIC ARCHIVE DESTINATION lib)
endif()

