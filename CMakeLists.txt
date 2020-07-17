 # HipMer v 2.0, Copyright (c) 2020, The Regents of the University of California,
 # through Lawrence Berkeley National Laboratory (subject to receipt of any required
 # approvals from the U.S. Dept. of Energy).  All rights reserved."
 
 # Redistribution and use in source and binary forms, with or without modification,
 # are permitted provided that the following conditions are met:
 
 # (1) Redistributions of source code must retain the above copyright notice, this
 # list of conditions and the following disclaimer.
 
 # (2) Redistributions in binary form must reproduce the above copyright notice,
 # this list of conditions and the following disclaimer in the documentation and/or
 # other materials provided with the distribution.
 
 # (3) Neither the name of the University of California, Lawrence Berkeley National
 # Laboratory, U.S. Dept. of Energy nor the names of its contributors may be used to
 # endorse or promote products derived from this software without specific prior
 # written permission.
 
 # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
 # EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 # OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
 # SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 # INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 # TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 # BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 # CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 # ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
 # DAMAGE.
 
 # You are under no obligation whatsoever to provide any bug fixes, patches, or upgrades
 # to the features, functionality or performance of the source code ("Enhancements") to
 # anyone; however, if you choose to make your Enhancements available either publicly,
 # or directly to Lawrence Berkeley National Laboratory, without imposing a separate
 # written license agreement for such Enhancements, then you hereby grant the following
 # license: a  non-exclusive, royalty-free perpetual license to install, use, modify,
 # prepare derivative works, incorporate into other computer software, distribute, and
 # sublicense such enhancements or derivative works thereof, in binary and source code
 # form.

# The name of our project is "MHMxx". CMakeLists files in this project can
# refer to the root source directory of the project as ${MHMXX_SOURCE_DIR} and
# to the root binary directory of the project as ${MHMXX_BINARY_DIR}
cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
#project(MHMxx LANGUAGE C CXX CUDA)
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

find_package(CpuArch
    REQUIRED
    COMPONENTS         
    OPTIONAL_COMPONENTS mmx avx sse sse2 fma ssse3 avx2 avx512f powerpc64 altivec vsx power8-fusion power8-vector )

LIST(APPEND CpuArch_FEATURES vsx powerpc64 power8-fusion power8-vector) # CpuArch sometimes does not find it

foreach(_ARCH ${CpuArch_FEATURES})
   check_c_compiler_flag(-m${_ARCH} HAVE_C_OPT_${_ARCH})
   if (HAVE_C_OPT_${_ARCH})
     message(STATUS "Compiling C with -m${_ARCH} option")
     string(APPEND CMAKE_C_FLAGS " -m${_ARCH}")
   endif()
   check_cxx_compiler_flag(-m${_ARCH} HAVE_CXX_OPT_${_ARCH})
   if (HAVE_CXX_OPT_${_ARCH})
     message(STATUS "Compiling CXX with -m${_ARCH} option")
     string(APPEND CMAKE_CXX_FLAGS " -m${_ARCH}")
   endif()
endforeach()

include(CMakeDefinitions.txt)

# verify submodules
include(CheckSubmodules)
CHECK_SUBMODULES(zstr upcxx-utils)# todo klign

# set version
include_directories(include)
GET_GIT_VERSION(MHMXX)

# add upcxx-utils
if (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/upcxx-utils/CMakeLists.txt")
  message("Did not find upcxx-utils/CMakeLists.txt")
  if (DEFINED ENV{UPCXX_UTILS_DIR} AND EXISTS $ENV{UPCXX_UTILS_DIR})
    set(CMAKE_MODULE_PATH "$ENV{UPCXX_UTILS_DIR}/cmake;${CMAKE_MODULE_PATH}")
    message(STATUS "Using externally provided upcxx-utils installation: $ENV{UPCXX_UTILS_DIR}")
    message(STATUS "Prepended $ENV{UPCXX_UTILS_DIR}/cmake to CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
  else()
    message(STATUS "No external upcxx-utils package provided")
    message(STATUS "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
  endif()
  find_package(UPCXX_UTILS 0.3.0 REQUIRED)
else()
  message(STATUS "Using upcxx-utils submodule in source tree")
  set(UPCXX_UTILS_AS_SUBMODULE TRUE)
endif()

if(UPCXX_UTILS_AS_SUBMODULE OR NOT UPCXX_UTILS_FOUND)

  message(STATUS "Did not find an installation of upcxx_utils, including submodule in src/upcxx_utils")
  include_directories(upcxx-utils/include)
  add_subdirectory(upcxx-utils)
  set(UPCXX_UTILS_FOUND TRUE)
  set(UPCXX_UTILS_DIR upcxx-utils)
  add_dependencies(MHMXX_VERSION UPCXX_UTILS_VERSION)

else()
  message(STATUS "Found existing upcxx-utils install: ${UPCXX_UTILS_FOUND} at ${UPCXX_UTILS_DIR}")
endif()
message(STATUS "added upcxx_utils (${UPCXX_UTILS_VERSION})")

message(STATUS "finding UPCXX for MHMxx build")
find_package(UPCXX 2019.9.0.1 REQUIRED)
if(NOT UPCXX_FOUND)
  message(FATAL_ERROR "UPCXX is required to build and it should have been found in the upcxx-utils submodule under src")
endif()

include(CheckLanguage)
check_language(CUDA)
if (CMAKE_CUDA_COMPILER)
  enable_language(CUDA)
endif()


# include mateidavid's zstr c++ library
#include_directories("zstr/src")

add_subdirectory(src)

