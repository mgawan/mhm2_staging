cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(MHM)

FIND_PACKAGE(ZLIB 1.2.3 REQUIRED)
include_directories(BEFORE ${ZLIB_INCLUDE_DIRECTORIES})
GET_FILENAME_COMPONENT(ZLIB_BASEPATH ${ZLIB_INCLUDE_DIRS} PATH)
LINK_DIRECTORIES(${ZLIB_BASEPATH}/lib64 ${ZLIB_BASEPATH}/lib)
MESSAGE("Adding ZLIB for -I${ZLIB_INCLUDE_DIRS} and -L${ZLIB_BASEPATH}/lib64 -L${ZLIB_BASEPATH}/lib)")

# all c++ must be at least c++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_COMPILER upcxx)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O")

add_definitions(-DCONFIG_USE_COLORS)

# include mateidavid's zstr c++ library
include_directories("/src/zstr/src")

add_executable(mhm main.cpp merge_reads.cpp)

target_link_libraries(mhm ${ZLIB_LIBRARIES})
#install(TARGETS mhm DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
