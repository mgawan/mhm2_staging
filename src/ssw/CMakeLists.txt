# all c++ must be at least c++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_COMPILER upcxx)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O")

add_library(SSW_LIBRARY ssw.cpp ssw_core.cpp)
