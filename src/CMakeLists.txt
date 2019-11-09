cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(MHM)

message("Building ${CMAKE_BUILD_TYPE} version")

find_package(ZLIB 1.2.3 REQUIRED)
include_directories(BEFORE ${ZLIB_INCLUDE_DIRECTORIES})
get_filename_component(ZLIB_BASEPATH ${ZLIB_INCLUDE_DIRS} PATH)
link_directories(${ZLIB_BASEPATH}/lib64 ${ZLIB_BASEPATH}/lib)
message("Adding ZLIB for -I${ZLIB_INCLUDE_DIRS} and -L${ZLIB_BASEPATH}/lib64 -L${ZLIB_BASEPATH}/lib)")

# all c++ must be at least c++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_COMPILER upcxx)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O")

add_definitions(-DMAX_KMER_LEN=128)
# quality encoding is logarithmic, so the probability of an error is:
# cutoff 10: 10% chance
# cutoff 20: 1% chance
# cutoff 30: 0.1% chance
add_definitions(-DQUAL_CUTOFF=20)
add_definitions(-DMAX_KMER_MAPPINGS=10)
add_definitions(-DBLOOM_FP=0.05)
add_definitions(-DMAX_FILE_PATH=255)
add_definitions(-DMAX_RANKS_PER_DIR=1024)
add_definitions(-DMAX_RPCS_IN_FLIGHT=8192) 
# 4MB in flight
add_definitions(-DMIN_INFLIGHT_BYTES=41944304)
add_definitions(-DUSE_COLORS)
add_definitions(-DUSE_BYTELL)
add_definitions(-DMIN_CTG_PRINT_LEN=200)
add_definitions(-DDEPTH_DIFF_THRES=0.1)
add_definitions(-DUSE_PROGBAR)
add_definitions(-DALN_EXTRA=0)
add_definitions(-DMAX_DEPTH_DIFF=15)
add_definitions(-DMIN_DEPTH_DIFF=5)
add_definitions(-DMAX_SEARCH_LEVEL=5)
add_definitions(-DMAX_QUEUE_SIZE=100)
add_definitions(-DALN_SLOP=0)
add_definitions(-DASSM_CLEN_THRES=500)
#FIXME: these should not be defined here!!
add_definitions(-DBREAK_SCAFF_NS=20)
add_definitions(-DMAX_CTG_GRAPH_DEGREE=100)
add_definitions(-DWALK_SUPPORT_THRES=3)
add_definitions(-DMAX_ALNS_PER_READ=20)
add_definitions(-DALIGN_EXPAND_BASES=5)
add_definitions(-DFIVE_PRIME_WIGGLE_ROOM=5)
add_definitions(-DTHREE_PRIME_WIGGLE_ROOM=5)
add_definitions(-DTMP_READ_LEN=150)

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

add_executable(mhm
    main.cpp merge_reads.cpp kcount.cpp dbjg_traversal.cpp hash_funcs.c klign.cpp cgraph.cpp build_ctg_graph.cpp walk_ctg_graph.cpp
    spanner.cpp)

target_link_libraries(mhm SSW_LIBRARY)
target_link_libraries(mhm ${ZLIB_LIBRARIES})


#install(TARGETS mhm DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

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
    add_definitions("-DMHM_VERSION=\"${GIT_TAG}-${GIT_COMMIT_HASH}${GIT_DIFF}-${GIT_BRANCH} (${CMAKE_BUILD_TYPE}) ${_timestamp}\"")
else()
    string(TIMESTAMP _timestamp)
    add_definitions("-DMHM_VERSION=\"unknown version (${CMAKE_BUILD_TYPE}) ${_timestamp}\"")
endif()    
