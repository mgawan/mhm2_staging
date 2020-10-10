# MHM2 testing

#
# The intent here is to include unit tests, like endian correctness for plain old C/C++
# and not depend on upcxx / upcxx-run to execute them, so make test can work on all platforms
# including multi node without calling upcxx::init() and upcxx::finalize()
#

add_subdirectory(googletest)

set(BINARY ${CMAKE_PROJECT_NAME}_test)

file(GLOB_RECURSE TEST_SOURCES LIST_DIRECTORIES false *.hpp *.cpp)

set(SOURCES ${TEST_SOURCES})

add_executable(${BINARY} ${TEST_SOURCES} ${MHM2_TARGET_OBJECTS})

add_test(NAME ${BINARY} COMMAND ${BINARY})

message("MHM2_LINK_LIBRARIES: ${MHM2_LINK_LIBRARIES}")

target_link_libraries(${BINARY} PUBLIC ${MHM2_LINK_LIBRARIES} ${UPCXX_UTILS_LIBRARIES} ${UPCXX_LIBRARIES}  gtest)
