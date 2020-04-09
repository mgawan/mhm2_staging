# The name of our project is "MHMxx". CMakeLists files in this project can
# refer to the root source directory of the project as ${MHMXX_SOURCE_DIR} and
# to the root binary directory of the project as ${MHMXX_BINARY_DIR}
cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(MHMxx)
#----------------------------------------------------------------------------------------#
#   project specification
#----------------------------------------------------------------------------------------#
file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/VERSION MHMXX_BASE_VERSION)

# use PROJECT_SOURCE_DIR, not CMAKE_SOURCE_DIR
set(MHMXX_SOURCE_DIR ${PROJECT_SOURCE_DIR})
set(MHMXX_BINARY_DIR ${PROJECT_BINARY_DIR})

# Build options
if (NOT MAX_BUILD_KMER)
  set(MAX_BUILD_KMER 160 CACHE STRING "The maximum kmer allowed to use (must be multiple of 32)")
endif()
message(STATUS "Building ${PROJECT_NAME} with a maximum kmer size of ${MAX_BUILD_KMER}")
add_definitions(-DMAX_BUILD_KMER=${MAX_BUILD_KMER})

# include modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/Modules)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
            FORCE)
endif()

message("Building ${CMAKE_BUILD_TYPE} version")
if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    add_definitions(-DDEBUG)
endif()

# check compiler and set appropriate linker flags
include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)
include(CheckCXXSourceCompiles)
include(CheckCXXSourceRuns)
include(GetGitVersion)


# all c++ must be at least c++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

if(NOT APPLE)
    add_definitions(-D_XOPEN_SOURCE=700)
endif()

# verify submodules
include(CheckSubmodules)
CHECK_SUBMODULES(zstr upcxx-utils)# todo klign

# set version
include_directories(include)
GET_GIT_VERSION(MHMXX)


# add upcxx-utils
if (DEFINED ENV{UPCXX_UTILS_DIR} AND EXISTS $ENV{UPCXX_UTILS_DIR})
  set(CMAKE_MODULE_PATH "$ENV{UPCXX_UTILS_DIR}/cmake;${CMAKE_MODULE_PATH}")
  message(STATUS "Using externally provided upcxx-utils installation: $ENV{UPCXX_UTILS_DIR}")
  message(STATUS "Prepended $ENV{UPCXX_UTILS_DIR}/cmake to CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
else()
  message(STATUS "No external upcxx-utils package provided")
  message(STATUS "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
endif()
find_package(UPCXX_UTILS 0.3.0 QUIET)
if(NOT UPCXX_UTILS_FOUND)

  message(STATUS "Did not find an installation of upcxx_utils, including submodule in src/upcxx_utils")
  include_directories(upcxx-utils/include)
  add_subdirectory(upcxx-utils)

else()
  message(STATUS "Found exßisting upcxx-utils install: ${UPCXX_UTILS_FOUND} at ${UPCXX_UTILS_DIR}")
  set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${UPCXX_UTILS_DIR}")
  message(STATUS "Appended the UPCXX_UTILS_DIR to the CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
endif()
message(STATUS "added upcxx_utils (${UPCXX_UTILS_VERSION})")

# include mateidavid's zstr c++ library
include_directories("zstr/src")

add_subdirectory(src)
