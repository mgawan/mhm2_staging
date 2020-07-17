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
    message(STATUS "UPCXX_UTILS is using the UPCXX::upcxx library interface: ${UPCXX_LIBRARIES}")
else()
    find_program(UPCXX_EXEC upcxx)
    set(CMAKE_CXX_COMPILER ${UPCXX_EXEC})
    message(STATUS "MHMxx is using upcxx directly as the UPCXX::upcxx library interface is not available in this low version of cmake: ${CMAKE_VERSION}"        )
    unset(UPCXX_LIBRARIES)
endif()

message("Building ${CMAKE_BUILD_TYPE} version")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-reorder -Wno-sign-compare -Wno-unused-function -Wno-catch-value")

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

include_directories("adept-sw")
if(ENABLE_CUDA)
    add_subdirectory(adept-sw)
endif()

add_executable(mhmxx
    main.cpp merge_reads.cpp kcount.cpp dbjg_traversal.cpp hash_funcs.c klign.cpp cgraph.cpp build_ctg_graph.cpp
    walk_ctg_graph.cpp spanner.cpp splinter.cpp localassm.cpp histogrammer.cpp aln_depths.cpp)

if(ENABLE_CUDA)
    target_link_libraries(mhmxx
        Threads::Threads
        SSW_LIBRARY
        ${ZLIB_LIBRARIES}
        ${UPCXX_LIBRARIES}
        ${UPCXX_UTILS_LIBRARIES}
        MHMXX_VERSION_LIB
        ADEPT_SW_LIBRARY)
else()
    target_link_libraries(mhmxx
        Threads::Threads
        SSW_LIBRARY
        ${ZLIB_LIBRARIES}
        ${UPCXX_LIBRARIES}
        ${UPCXX_UTILS_LIBRARIES}
        MHMXX_VERSION_LIB)
endif()

install(PROGRAMS mhmxx.py DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
install(TARGETS mhmxx DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
install(TARGETS mhmxx EXPORT ADEPT_SW_LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
