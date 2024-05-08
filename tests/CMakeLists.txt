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

# Get the current working branch
execute_process(
        COMMAND git rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE)

# Get the latest commit hash
execute_process(
        COMMAND git rev-parse HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_COMMIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE
        )



set(EXTRA_SRC
        ${CMAKE_HOME_DIRECTORY}/src/IO/CSVExporter.cpp
        ${CMAKE_HOME_DIRECTORY}/src_shared/File.cpp
        ${CMAKE_HOME_DIRECTORY}/src_shared/FlowMPI.h #contains templates
        ${CMAKE_HOME_DIRECTORY}/src_shared/FlowMPI.cpp
        )

# Now simply link against gtest or gtest_main as needed. Eg
add_executable(testsuite testsuite.cpp ${EXTRA_SRC})

set_target_properties(testsuite PROPERTIES EXECUTABLE_OUTPUT_DIRECTORY ${CMAKE_EXECUTABLE_OUTPUT_DIRECTORY}/tests)
set_target_properties(testsuite PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/tests)
target_link_libraries(testsuite PRIVATE gtest_main molflow_core)
target_compile_definitions(testsuite PRIVATE
        "-DGIT_COMMIT_HASH=\"${GIT_COMMIT_HASH}\"")

target_include_directories(testsuite PRIVATE ${CMAKE_HOME_DIRECTORY}/include ${CMAKE_HOME_DIRECTORY}/src_shared ${CMAKE_HOME_DIRECTORY}/src)

target_link_directories(testsuite PRIVATE ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

target_compile_features(testsuite PRIVATE cxx_std_17)
#add_test(NAME example_test COMMAND testsuite)

file(COPY
        ./TestCases
        ${CMAKE_HOME_DIRECTORY}/copy_to_build/parameter_catalog
        DESTINATION ${CMAKE_EXECUTABLE_OUTPUT_DIRECTORY}/tests)