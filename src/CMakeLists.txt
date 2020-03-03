cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(MHM)

message("Building ${CMAKE_BUILD_TYPE} version")

#find_package(UPCXX REQUIRED)
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

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -network=udp")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O")
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    # ignore annoying warning on icc
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd597")
endif()

add_definitions(-DMAX_KMER_LEN=160)
# quality encoding is logarithmic, so the probability of an error is:
# cutoff 10: 10% chance
# cutoff 20: 1% chance
# cutoff 30: 0.1% chance
add_definitions(-DQUAL_CUTOFF=20)
add_definitions(-DFILTER_BAD_QUAL_IN_READ)
add_definitions(-DMAX_KMER_MAPPINGS=10)
add_definitions(-DBLOOM_FP=0.05)
add_definitions(-DMAX_FILE_PATH=255)
add_definitions(-DMAX_RANKS_PER_DIR=1024)
add_definitions(-DMAX_RPCS_IN_FLIGHT=8192) 
# 4MB in flight
add_definitions(-DMIN_INFLIGHT_BYTES=41944304)
add_definitions(-DUSE_COLORS)
add_definitions(-DUSE_BYTELL)
add_definitions(-DUSE_PROGBAR)
add_definitions(-DMIN_CTG_PRINT_LEN=500)
add_definitions(-DDEPTH_DIFF_THRES=0.1)
add_definitions(-DMAX_DEPTH_DIFF=15)
add_definitions(-DMIN_DEPTH_DIFF=5)
add_definitions(-DMAX_SEARCH_LEVEL=5)
add_definitions(-DMAX_QUEUE_SIZE=100)
add_definitions(-DASSM_CLEN_THRES=500)
add_definitions(-DALN_SEED_SPACE=8)
add_definitions(-DUSE_ALNS_FOR_GAP_CORRECTION)
# This should be a parameter - it is a good way of tuning errors vs ctgy - the higher it is, the higher the errors and ctgy
add_definitions(-DMISMATCH_THRES=1)
# another good candidate for errors v ctgy tradeoff. THis applies to SSW alignments in gap closing
add_definitions(-DSSW_MISMATCH_TOLERANCE=3)
# FIXME: these should not be defined here!!
add_definitions(-DBREAK_SCAFF_NS=10)
add_definitions(-DMAX_CTG_GRAPH_DEGREE=100)
add_definitions(-DWALK_SUPPORT_THRES=3)
add_definitions(-DMAX_ALNS_PER_READ=20)
add_definitions(-DALIGN_EXPAND_BASES=5)
add_definitions(-DALN_WIGGLE=1)
add_definitions(-DMAX_NUM_GET_CTGS=100000)
# FIXME: these should be parameters
add_definitions(-DLASSM_MAX_KMER_LEN=121)
add_definitions(-DLASSM_MIN_KMER_LEN=21)
add_definitions(-DLASSM_SHIFT_SIZE=8)
add_definitions(-DLASSM_MIN_QUAL=10)
add_definitions(-DLASSM_MIN_HI_QUAL=20)
add_definitions(-DLASSM_MIN_VIABLE_DEPTH=0.2)
add_definitions(-DLASSM_MIN_EXPECTED_DEPTH=0.3)
add_definitions(-DLASSM_RATING_THRES=0)
add_definitions(-DLASSM_MAX_COUNT_MERS_READS=3000)


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
    main.cpp merge_reads.cpp kcount.cpp dbjg_traversal.cpp hash_funcs.c klign.cpp cgraph.cpp build_ctg_graph.cpp walk_ctg_graph.cpp
    spanner.cpp splinter.cpp localassm.cpp)

target_link_libraries(mhmxx SSW_LIBRARY)
target_link_libraries(mhmxx ${ZLIB_LIBRARIES})


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
    add_definitions("-DMHM_VERSION=\"${GIT_TAG}-${GIT_COMMIT_HASH}${GIT_DIFF}-${GIT_BRANCH} (${CMAKE_BUILD_TYPE}) ${_timestamp}\"")
else()
    string(TIMESTAMP _timestamp)
    add_definitions("-DMHM_VERSION=\"unknown version (${CMAKE_BUILD_TYPE}) ${_timestamp}\"")
endif()    
