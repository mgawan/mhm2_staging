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

cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

if (${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.13 AND DEFINED UPCXX_LIBRARIES)
  message(STATUS "MHM2 is using the UPCXX::upcxx library interface: ${UPCXX_LIBRARIES}")
  find_package(UPCXX)
else()
  find_program(UPCXX_EXEC upcxx)
  set(CMAKE_CXX_COMPILER ${UPCXX_EXEC})
  message(STATUS "MHM2 is using upcxx directly as the UPCXX::upcxx library interface is not available in this low version of cmake: ${CMAKE_VERSION}"    )
  unset(UPCXX_LIBRARIES)
endif()

message("Building ${CMAKE_BUILD_TYPE} version")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-reorder -Wno-sign-compare -Wno-unused-function -Wno-catch-value")
find_package(Threads REQUIRED)

find_package(ZLIB 1.2.3 REQUIRED)
include_directories(BEFORE ${ZLIB_INCLUDE_DIRECTORIES})
get_filename_component(ZLIB_BASEPATH ${ZLIB_INCLUDE_DIRS} PATH)
link_directories(${ZLIB_BASEPATH}/lib64 ${ZLIB_BASEPATH}/lib)
message("Adding ZLIB for -I${ZLIB_INCLUDE_DIRS} and -L${ZLIB_BASEPATH}/lib64 -L${ZLIB_BASEPATH}/lib)")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")


if(ENABLE_CUDA)
    add_definitions(-DENABLE_GPUS)
    message(STATUS "Building for GPUs with CUDA")
endif()

# a bit of jiggery-pokery to get a list of all compile time definitions in the code
get_directory_property(dir_defs COMPILE_DEFINITIONS)
set(all_defs "")
message("Definitions:")
foreach(d ${dir_defs})
    message("  " ${d})
    set(all_defs "${all_defs} ${d}")
endforeach()

add_definitions(-DALL_DEFNS=\"${all_defs}\")

# include ssw library
include_directories("ssw")
add_subdirectory(ssw)

set(MHM2_TARGET_OBJECTS)
add_library(hash_funcs OBJECT hash_funcs.c)
list(APPEND MHM2_TARGET_OBJECTS $<TARGET_OBJECTS:hash_funcs>)

# build all the extern template souce code instantiation explosion
set(EXTERN_TEMPLATE_FILES)
foreach(KMER_LENGTH ${MHM2_KMER_LENGTHS})
  foreach(TEMPLATE "kmer_dht" "kcount" "kmer" "klign" "scaffolding" "contigging" "post_assembly")
    set(TEMPLATE_FILE "${TEMPLATE}-extern-template-${KMER_LENGTH}")
    configure_file("${TEMPLATE}-extern-template.in.cpp" "${TEMPLATE_FILE}.cpp" @ONLY)
    list(APPEND EXTERN_TEMPLATE_FILES ${TEMPLATE_FILE})
  endforeach()
endforeach()

#foreach(KMER_LENGTH ${MHM2_KMER_LENGTHS})
#  foreach(TEMPLATE "kcount" "klign" "scaffolding" "contigging" "post_assembly")
#    set(TEMPLATE_FILE "${TEMPLATE}-extern-template-${KMER_LENGTH}")
#    string(TOUPPER ${TEMPLATE} TEMPLATE_UPPER)
#    configure_file("extern_template.in.cpp" "${TEMPLATE_FILE}.cpp" @ONLY)
#    list(APPEND EXTERN_TEMPLATE_FILES ${TEMPLATE_FILE})
#  endforeach()
#endforeach()


foreach(tgt fastq alignments contigs packed_reads utils bloom ctg_graph options
            main_static_objs merge_reads dbjg_traversal klign cgraph build_ctg_graph
            walk_ctg_graph spanner splinter localassm histogrammer aln_depths
            ${EXTERN_TEMPLATE_FILES}
        )
  add_library(${tgt} OBJECT ${tgt}.cpp)
  if (${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.13 AND DEFINED UPCXX_LIBRARIES)
    target_link_libraries(${tgt} ${UPCXX_LIBRARIES} ${UPCXX_UTILS_LIBRARIES})
  endif()
  list(APPEND MHM2_TARGET_OBJECTS $<TARGET_OBJECTS:${tgt}>)
endforeach()

set(MHM2_TARGET_OBJECTS ${MHM2_TARGET_OBJECTS} PARENT_SCOPE)
add_executable(mhm2 main.cpp ${MHM2_TARGET_OBJECTS})

set(MHM2_LINK_LIBRARIES Threads::Threads SSW_LIBRARY ${ZLIB_LIBRARIES} MHM2_VERSION_LIB)

if(ENABLE_CUDA)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR})
    include_directories("adept-sw")
    add_subdirectory(adept-sw)
    set_property(TARGET mhm2 PROPERTY CUDA_RESOLVE_DEVICE_SYMBOLS OFF) # just to be sure. It should be off in parent dir
    set(MHM2_LINK_LIBRARIES ${MHM2_LINK_LIBRARIES} ADEPT_SW_LIBRARY_static)
    target_link_libraries(mhm2
        ${MHM2_LINK_LIBRARIES}
        ${UPCXX_LIBRARIES} ${UPCXX_UTILS_LIBRARIES})

else()
    target_link_libraries(mhm2
        ${MHM2_LINK_LIBRARIES}
        ${UPCXX_LIBRARIES} ${UPCXX_UTILS_LIBRARIES})
endif()

set(MHM2_LINK_LIBRARIES ${MHM2_LINK_LIBRARIES} PARENT_SCOPE)

install(TARGETS mhm2 DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
  SET(MHM2PY_RUNTIME_ENV "GASNET_BACKTRACE='1', ")
elseif("${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo")
  SET(MHM2PY_RUNTIME_ENV "GASNET_BACKTRACE='1', ")
endif()
configure_file(mhm2.py ${CMAKE_CURRENT_BINARY_DIR}/mhm2.py )
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/mhm2.py DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
