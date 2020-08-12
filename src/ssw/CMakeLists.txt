# all c++ must be at least c++14

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


add_library(SSW_LIBRARY ssw.cpp ssw_core.cpp)
set_target_properties(SSW_LIBRARY PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS OFF)
