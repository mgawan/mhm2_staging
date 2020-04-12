cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

message(STATUS "finding UPCXX for MHMxx build")
find_package(UPCXX 2019.9.0.1 REQUIRED)
if (${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.13 AND DEFINED UPCXX_LIBRARIES)
  message(STATUS "UPCXX_UTILS is using the UPCXX::upcxx library interface: ${UPCXX_LIBRARIES}")
  find_package(UPCXX REQUIRED)
else()
  find_program(UPCXX_EXEC upcxx)
  set(CMAKE_CXX_COMPILER ${UPCXX_EXEC})
  message(STATUS "UPCXX_UTILS is using upcxx directly as the UPCXX::upcxx library interface is not available in this low version of cmake: ${CMAKE_VERSION}"    )
  unset(UPCXX_LIBRARIES)
endif()

message("Building ${CMAKE_BUILD_TYPE} version")

find_package(Threads REQUIRED)

find_package(ZLIB 1.2.3 REQUIRED)
include_directories(BEFORE ${ZLIB_INCLUDE_DIRECTORIES})
get_filename_component(ZLIB_BASEPATH ${ZLIB_INCLUDE_DIRS} PATH)
link_directories(${ZLIB_BASEPATH}/lib64 ${ZLIB_BASEPATH}/lib)
message("Adding ZLIB for -I${ZLIB_INCLUDE_DIRS} and -L${ZLIB_BASEPATH}/lib64 -L${ZLIB_BASEPATH}/lib)")

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

add_executable(mhmxx
    main.cpp merge_reads.cpp kcount.cpp dbjg_traversal.cpp hash_funcs.c klign.cpp cgraph.cpp build_ctg_graph.cpp 
    walk_ctg_graph.cpp spanner.cpp splinter.cpp localassm.cpp histogrammer.cpp)

target_link_libraries(mhmxx Threads::Threads SSW_LIBRARY ${ZLIB_LIBRARIES} ${UPCXX_LIBRARIES} ${UPCXX_UTILS_LIBRARIES} MHMXX_VERSION_LIB)

install(PROGRAMS mhmxx.py DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
install(TARGETS mhmxx DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
