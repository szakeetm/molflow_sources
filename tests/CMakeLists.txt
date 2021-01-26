cmake_minimum_required(VERSION 3.12.2 FATAL_ERROR)

project(molflow_tests)

# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
        ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
        EXCLUDE_FROM_ALL)

# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
if (CMAKE_VERSION VERSION_LESS 2.8.11)
    include_directories("${gtest_SOURCE_DIR}/include")
endif()

link_directories(${LINK_DIR_1} ${LINK_DIR_2})

set(EXTRA_SRC
        ../src/ParameterParser.cpp)
# Now simply link against gtest or gtest_main as needed. Eg
add_executable(testsuite simulation_testsuite.cpp ${EXTRA_SRC})
target_link_libraries(testsuite gtest_main molflowSim)

target_include_directories(testsuite PRIVATE ../include ../src_shared)
target_compile_features(testsuite PRIVATE cxx_std_17)
add_test(NAME example_test COMMAND testsuite)

file(COPY ./simulation.cfg ./test_lr1000_pipe.xml DESTINATION ${CMAKE_EXECUTABLE_OUTPUT_DIRECTORY})