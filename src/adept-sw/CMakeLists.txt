if(NOT ENABLE_CUDA)
    message(FATAL_ERROR "Trying to build ADEPT-SW but CUDA is not enabled")
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler=-Wall -gencode arch=compute_70,code=sm_70")

add_library(ADEPT_SW_LIBRARY SHARED driver.cpp kernel.cpp gpu_alns.cpp utils_gpu.cpp)
target_link_libraries(ADEPT_SW_LIBRARY)
install(TARGETS ADEPT_SW_LIBRARY ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

set_source_files_properties(driver.cpp kernel.cpp gpu_alns.cpp utils_gpu.cpp PROPERTIES LANGUAGE CUDA LINKER_LANGUAGE CUDA)


