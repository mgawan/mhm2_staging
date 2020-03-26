cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(MHMXX)

message("Building ${CMAKE_BUILD_TYPE} version")

#find_package(UPCXX REQUIRED)
find_package(ZLIB 1.2.3 REQUIRED)
include_directories(BEFORE ${ZLIB_INCLUDE_DIRECTORIES})
get_filename_component(ZLIB_BASEPATH ${ZLIB_INCLUDE_DIRS} PATH)
link_directories(${ZLIB_BASEPATH}/lib64 ${ZLIB_BASEPATH}/lib)
message("Adding ZLIB for -I${ZLIB_INCLUDE_DIRS} and -L${ZLIB_BASEPATH}/lib64 -L${ZLIB_BASEPATH}/lib)")

# For some reason these don't work on NERSC's Cori or Linux. So We just explicitly include the -std=c++17 in the CXX flags below
#set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_COMPILER upcxx)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -network=udp")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O")
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    # ignore annoying warning on icc
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd597")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wall -Wno-reorder -Wno-sign-compare -Wno-unused-function -pthread")

include(CMakeDefinitions.txt)

# a bit of jiggery-pokery to get a list of all compile time definitions in the code
get_directory_property(dir_defs COMPILE_DEFINITIONS)
set(all_defs "")
message("Definitions:")
foreach(d ${dir_defs})
    message("  " ${d})
    set(all_defs "${all_defs} ${d}")
endforeach()

add_definitions(-DALL_DEFNS=\"${all_defs}\")

# include mateidavid's zstr c++ library
include_directories("${CMAKE_SOURCE_DIR}/zstr/src")
# include ssw library
include_directories("${CMAKE_SOURCE_DIR}/ssw")
# all the other headers
include_directories("${CMAKE_SOURCE_DIR}/../include")

add_subdirectory(ssw)

add_executable(mhmxx
    main.cpp merge_reads.cpp kcount.cpp dbjg_traversal.cpp hash_funcs.c klign.cpp cgraph.cpp build_ctg_graph.cpp 
    walk_ctg_graph.cpp spanner.cpp splinter.cpp localassm.cpp)

target_link_libraries(mhmxx SSW_LIBRARY ${ZLIB_LIBRARIES})

#install(TARGETS mhmxx DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

# Set up version
find_package(Git)
if (GIT_FOUND)
    # Get the current working branch
    execute_process(
        COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    # Get the latest abbreviated commit hash of the working branch
    execute_process(
        COMMAND ${GIT_EXECUTABLE} log -1 --format=%h
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_COMMIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    execute_process(
        COMMAND bash -c "${GIT_EXECUTABLE} diff --quiet --exit-code || echo -dirty"
        OUTPUT_VARIABLE GIT_DIFF
        OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    execute_process(
        COMMAND ${GIT_EXECUTABLE} tag
        OUTPUT_VARIABLE GIT_TAG
        OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    string(TIMESTAMP _timestamp)
    add_definitions("-DMHMXX_VERSION=\"${GIT_TAG}-${GIT_COMMIT_HASH}${GIT_DIFF}-${GIT_BRANCH} (${CMAKE_BUILD_TYPE}) built on ${_timestamp}\"")
else()
    string(TIMESTAMP _timestamp)
    add_definitions("-DMHMXX_VERSION=\"unknown version (${CMAKE_BUILD_TYPE}) build on ${_timestamp}\"")
endif()    
