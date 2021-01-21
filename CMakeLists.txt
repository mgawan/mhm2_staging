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

# The name of our project is "MHM2". CMakeLists files in this project can
# refer to the root source directory of the project as ${MHM2_SOURCE_DIR} and
# to the root binary directory of the project as ${MHM2_BINARY_DIR}
cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(MHM2)

message(STATUS "Building ${CMAKE_PROJECT_NAME} with CMake ${CMAKE_VERSION}")
#----------------------------------------------------------------------------------------#
#   project specification
#----------------------------------------------------------------------------------------#
file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/VERSION MHM2_BASE_VERSION)

# use PROJECT_SOURCE_DIR, not CMAKE_SOURCE_DIR
set(MHM2_SOURCE_DIR ${PROJECT_SOURCE_DIR})
set(MHM2_BINARY_DIR ${PROJECT_BINARY_DIR})

# include modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/Modules ${PROJECT_SOURCE_DIR}/upcxx-utils/cmake/Modules)

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

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif()

include(CheckLanguage)
check_language(CUDA)
if(CMAKE_CUDA_COMPILER)
    option(ENABLE_CUDA "Enable CUDA" ON)
    if (ENABLE_CUDA)
      SET(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
      SET(CMAKE_CUDA_RESOLVE_DEVICE_SYMBOLS OFF) # only build device link objects for GPU targets
      enable_language(CUDA)
    endif()
else()
    message(STATUS "No CUDA environment detected")
    set(ENABLE_CUDA OFF CACHE BOOL "Enable CUDA" FORCE)
endif()
if(ENABLE_CUDA)
  message(STATUS "Building for GPU with CUDA")
else()
  message(STATUS "Building for CPU only")
endif()

option(MHM2_VECTORS "Test for optimal arch-specific compiler flags" ON)
if (MHM2_VECTORS)
  find_package(CpuArch
    REQUIRED
    COMPONENTS
    OPTIONAL_COMPONENTS mmx avx sse sse2 fma ssse3 avx2 avx512f powerpc64 altivec vsx power8-fusion power8-vector )

  if (NOT APPLE)
    LIST(APPEND CpuArch_FEATURES vsx powerpc64 power8-fusion power8-vector) # CpuArch sometimes does not find it
  endif()

  option(MHM2_DEBUG_VECTORS "Turn on vectorization compile flags for Debug build" OFF)
  foreach(_ARCH ${CpuArch_FEATURES})
    check_c_compiler_flag(-m${_ARCH} HAVE_C_OPT_${_ARCH})
    if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug" AND NOT MHM2_DEBUG_VECTORS)
      message("NOT enabling -m${_ARCH} compile option in this Debug build")
    else()
      if (HAVE_C_OPT_${_ARCH})
        message(STATUS "Compiling C with -m${_ARCH} option")
        string(APPEND CMAKE_C_FLAGS " -m${_ARCH}")
      endif()
      check_cxx_compiler_flag(-m${_ARCH} HAVE_CXX_OPT_${_ARCH})
      if (HAVE_CXX_OPT_${_ARCH})
        message(STATUS "Compiling CXX with -m${_ARCH} option")
        string(APPEND CMAKE_CXX_FLAGS " -m${_ARCH}")
      endif()
    endif()
  endforeach()
endif()

option(ENABLE_GASNET_STATS "Turn on gasnet stats recording" OFF)
message("Building ${CMAKE_BUILD_TYPE} version")
if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
  add_definitions(-DDEBUG)
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
  #set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_C_FLAGS_DEBUG}")
  set(ENABLE_GASNET_STATS ON)
elseif("${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo")
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
  #set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_C_FLAGS_RELWITHDEBINFO}")
  set(ENABLE_GASNET_STATS ON)
  set(ENV{UPCXX_CODEMODE} "debug")  # build with all upcxx assertions turned on (may still be necessary to build with in older CMake versions)
elseif("${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebug")
  # suggested by https://stackoverflow.com/questions/48754619/what-are-cmake-build-type-debug-release-relwithdebinfo-and-minsizerel
  # RelWithDebug [`-O2 -g -DDEBUG`] - aka RelWithDebInfo but `-DNDEBUG` removed and `-DDEBUG` added.
  # Additionally DBG logs are disabled BUT dbg files might be made for all ranks to capture all higher level log messages
  string(REGEX REPLACE "-DNDEBUG" "" CMAKE_CXX_FLAGS_RELWITHDEBUG "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DDEBUG" )
  string(REGEX REPLACE "-DNDEBUG" "" CMAKE_C_FLAGS_RELWITHDEBUG "${CMAKE_C_FLAGS_RELWITHDEBINFO} -DDEBUG" )
  add_definitions(-DNO_DBG_LOGS)
  set(ENABLE_GASNET_STATS ON)
  set(ENV{UPCXX_CODEMODE} "debug")  # build with all upcxx assertions turned on
  set(CMAKE_CXX_FLAGS_RELWITHDEBUG "${CMAKE_CXX_FLAGS_RELWITHDEBUG}" CACHE STRING "CXX flags in RelWithDebug build mode" FORCE)
  set(CMAKE_C_FLAGS_RELWITHDEBUG "${CMAKE_C_FLAGS_RELWITHDEBUG}" CACHE STRING "C flags in RelWithDebug build mode" FORCE)
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELWITHDEBUG}")
  #set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_C_FLAGS_RELWITHDEBUG}")
elseif("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
  #set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_C_FLAGS_RELEASE}")
endif()

if(ENABLE_GASNET_STATS)
  add_definitions(-DENABLE_GASNET_STATS)
endif()

option(UPCXX_VERBOSE "" ON)
set(MIN_UPCXX_VERSION "2020.3" CACHE STRING "Minimum version of UPCXX required for MHM2")
find_package(UPCXX ${MIN_UPCXX_VERSION})
if(NOT UPCXX_FOUND)
  message(FATAL_ERROR "UPCXX (version ${MIN_UPCXX_VERSION}) is required to build.
  Please verify the installation of upcxx
  and verify that upcxx is in your PATH: which upcxx && upcxx --version
  and has been installed in the correct location: find $(dirname $(which upcxx))/../share -name 'UPCXX*' -type f
    which should discover UPCXXConfig.cmake and UPCXXConfigVersion.cmake")
endif()
find_package(UPCXX ${MIN_UPCXX_VERSION} REQUIRED)

message(STATUS "Building for ${CMAKE_BUILD_TYPE}: CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}")

find_package(Threads REQUIRED)
find_package(ZLIB 1.2.3 REQUIRED)
include_directories(BEFORE ${ZLIB_INCLUDE_DIRS})
get_filename_component(ZLIB_BASEPATH ${ZLIB_INCLUDE_DIRS} PATH)
link_directories(${ZLIB_BASEPATH}/lib64 ${ZLIB_BASEPATH}/lib)
message("Adding ZLIB for -I${ZLIB_INCLUDE_DIRS} and -L${ZLIB_BASEPATH}/lib64 -L${ZLIB_BASEPATH}/lib)")

# Build options
if (NOT MAX_BUILD_KMER)
  set(MAX_BUILD_KMER 160 CACHE STRING "The maximum kmer allowed to use (must be multiple of 32)")
endif()
math(EXPR TEST_MAX_KMER_SIZE "((${MAX_BUILD_KMER} + 31) / 32) * 32")
if (NOT ${MAX_BUILD_KMER} EQUAL ${TEST_MAX_KMER_SIZE})
  message(FATAL_ERROR "MAX_BUILD_KMER must be a multiple of 32.  MAX_BUILD_KMER=${MAX_BUILD_KMER} is not")
endif()
add_definitions(-DMAX_BUILD_KMER=${MAX_BUILD_KMER})

# enumerate the kmer lengths that will be built
set(MHM2_KMER_LENGTHS)
set(tmplen 32)
while( NOT ${MAX_BUILD_KMER} LESS ${tmplen} )
    set(MHM2_KMER_LENGTHS ${MHM2_KMER_LENGTHS} ${tmplen} )
    math(EXPR tmplen "(${tmplen} + 32)")
endwhile()

message(STATUS "Building ${PROJECT_NAME} with a maximum kmer size of ${MAX_BUILD_KMER}: ${MHM2_KMER_LENGTHS}")


if(NOT APPLE)
  add_definitions(-D_XOPEN_SOURCE=700)
endif()


include(CMakeDefinitions.txt)

# verify submodules
include(CheckSubmodules)
CHECK_SUBMODULES(zstr upcxx-utils test/googletest)# todo klign

# set version
include_directories(include src)
GET_GIT_VERSION(MHM2)

# add upcxx-utils
if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/upcxx-utils/CMakeLists.txt")
  message("Did not find upcxx-utils/CMakeLists.txt")
  if (DEFINED ENV{UPCXX_UTILS_DIR} AND EXISTS $ENV{UPCXX_UTILS_DIR})
    set(CMAKE_MODULE_PATH "$ENV{UPCXX_UTILS_DIR}/cmake;${CMAKE_MODULE_PATH}")
    message(STATUS "Using externally provided upcxx-utils installation: $ENV{UPCXX_UTILS_DIR}")
    message(STATUS "Prepended $ENV{UPCXX_UTILS_DIR}/cmake to CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
  else()
    message(STATUS "No external upcxx-utils package provided")
    message(STATUS "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
  endif()
  find_package(UPCXX_UTILS 0.3.4 REQUIRED)
else()
  message(STATUS "Using upcxx-utils submodule in source tree")
  set(UPCXX_UTILS_AS_SUBMODULE TRUE)
endif()

if(UPCXX_UTILS_AS_SUBMODULE OR NOT UPCXX_UTILS_FOUND)
  message(STATUS "Did not find an installation of upcxx_utils, including submodule in src/upcxx_utils")
  option(UPCXX_ENABLE_TESTING OFF)
  include_directories(upcxx-utils/include)
  add_subdirectory(upcxx-utils)
  set(UPCXX_UTILS_FOUND TRUE)
  set(UPCXX_UTILS_DIR upcxx-utils)
  add_dependencies(MHM2_VERSION UPCXX_UTILS_VERSION)
else()
  message(STATUS "Found existing upcxx-utils install: ${UPCXX_UTILS_FOUND} at ${UPCXX_UTILS_DIR}")
endif()
message(STATUS "added upcxx_utils (${UPCXX_UTILS_VERSION})")

# include mateidavid's zstr c++ library
include_directories("zstr/src")

add_subdirectory(src)
install(PROGRAMS upcxx-run-summit DESTINATION ${CMAKE_INSTALL_PREFIX}/bin RENAME mhm2-upcxx-run-summit)
install(PROGRAMS ci/check_asm_quality.py ci/ci_asm_qual_test.sh DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
install(FILES ci/good-arctic-sample0.txt ci/good-arctic-sample0-k63.txt ci/good-arctic-sample0-k31k63.txt ci/good-arcticsynth.txt DESTINATION ${CMAKE_INSTALL_PREFIX}/share)

option(MHM2_ENABLE_TESTING "Enable testing in MHM2" ON)
if (MHM2_ENABLE_TESTING)
  add_subdirectory(test)
endif()
