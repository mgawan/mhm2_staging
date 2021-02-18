# MHM2 testing

#
# The intent here is to include unit tests, like endian correctness for plain old C/C++
# and not depend on upcxx / upcxx-run to execute them, so make test can work on all platforms
# including multi node without calling upcxx::init() and upcxx::finalize()
#
enable_testing()
add_subdirectory(googletest)

if(ENABLE_CUDA)
    add_definitions(-DENABLE_GPUS)
    message(STATUS "Building tests for GPUs with CUDA")
    include_directories(${CMAKE_SOURCE_DIR}/src)
endif()

set(BINARY mhm2_test)

file(GLOB TEST_SOURCES LIST_DIRECTORIES false *.hpp *.cpp)

set(SOURCES ${TEST_SOURCES})

include_directories("${CMAKE_SOURCE_DIR}/src/ssw")

add_executable(${BINARY} ${TEST_SOURCES} ${MHM2_TARGET_OBJECTS})

add_test(NAME ${BINARY} COMMAND ${BINARY})

message("MHM2_LINK_LIBRARIES: ${MHM2_LINK_LIBRARIES}")

target_link_libraries(${BINARY} PUBLIC ${MHM2_LINK_LIBRARIES} ${UPCXX_UTILS_LIBRARIES} ${UPCXX_LIBRARIES}  gtest)
install(TARGETS ${BINARY} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

# enable "make check"
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} -V)
add_dependencies(check ${BINARY})