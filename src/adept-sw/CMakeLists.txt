if(NOT ENABLE_CUDA)
    message(FATAL_ERROR "Trying to build ADEPT-SW but CUDA is not enabled")
endif()

set(CMAKE_CUDA_STANDARD 14 CACHE STRING "")
set(CMAKE_CUDA_STANDARD_REQUIRED ON CACHE BOOL "")
set(CMAKE_CUDA_EXTENSIONS OFF CACHE BOOL "")
enable_language(CUDA)

find_package(OpenMP REQUIRED)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler=-Wall,${OpenMP_CXX_FLAGS} -gencode arch=compute_70,code=sm_70")
message(STATUS "CMAKE_CUDA_FLAGS=${CMAKE_CUDA_FLAGS}")

add_library(ADEPT_SW_LIBRARY_obj OBJECT driver.cpp kernel.cpp gpu_alns.cpp utils_gpu.cpp)
if (${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.13)
  target_link_libraries(ADEPT_SW_LIBRARY_obj INTERFACE OpenMP::OpenMP_CXX)
endif()

set_source_files_properties(driver.cpp kernel.cpp gpu_alns.cpp utils_gpu.cpp PROPERTIES LANGUAGE CUDA LINKER_LANGUAGE CUDA)
option(ADEPT_SW_SHARED "Adept-SW shared library" ON)
option(ADEPT_SW_STATIC "Adept-SW static library" ON)

if (ADEPT_SW_SHARED)
  add_library(ADEPT_SW_LIBRARY_shared SHARED $<TARGET_OBJECTS:ADEPT_SW_LIBRARY_obj>)
  set_property(TARGET ADEPT_SW_LIBRARY_shared PROPERTY CUDA_RESOLVE_DEVICE_SYMBOLS ON)
  target_link_libraries(ADEPT_SW_LIBRARY_shared INTERFACE OpenMP::OpenMP_CXX) 
  install(TARGETS ADEPT_SW_LIBRARY_shared LIBRARY DESTINATION lib)
endif()
if (ADEPT_SW_STATIC)
  add_library(ADEPT_SW_LIBRARY_static STATIC $<TARGET_OBJECTS:ADEPT_SW_LIBRARY_obj>)
  target_link_libraries(ADEPT_SW_LIBRARY_static INTERFACE OpenMP::OpenMP_CXX)
  set_property(TARGET ADEPT_SW_LIBRARY_static PROPERTY CUDA_RESOLVE_DEVICE_SYMBOLS ON)
  # necessary to link static lib with broken dependency calculations for cmake_device_link objects
  set(ADEPT_SW_BROKEN_DEVICE_LINK_HACK ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/ADEPT_SW_LIBRARY_static.dir/cmake_device_link.o PARENT_SCOPE)
  install(TARGETS ADEPT_SW_LIBRARY_static ARCHIVE DESTINATION lib)
endif()

