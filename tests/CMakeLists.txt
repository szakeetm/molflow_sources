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
        OUTPUT_STRIP_TRAILING_WHITESPACE)


link_directories(${LINK_DIR_1} ${LINK_DIR_2})

set(EXTRA_SRC
        ../src/IO/CSVExporter.cpp
        ../src_shared/File.cpp
        ../src_shared/FlowMPI.h #contains templates
        ../src_shared/FlowMPI.cpp
        )

# Now simply link against gtest or gtest_main as needed. Eg
add_executable(testsuite simulation_testsuite.cpp ${EXTRA_SRC})


#[[set(CMAKE_EXECUTABLE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${OUTPUT_BIN_REL}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY    "${CMAKE_BINARY_DIR}/${OUTPUT_BIN_REL}")]]
set_target_properties(testsuite PROPERTIES EXECUTABLE_OUTPUT_DIRECTORY ${CMAKE_EXECUTABLE_OUTPUT_DIRECTORY}/tests)
set_target_properties(testsuite PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/tests)
target_link_libraries(testsuite PRIVATE gtest_main molflowSim)
target_compile_definitions(testsuite PRIVATE
        "-DGIT_COMMIT_HASH=\"${GIT_COMMIT_HASH}\"")

target_include_directories(testsuite PRIVATE ../include ../src_shared ../src)
# set the path to the library folder
IF (WIN32)

    set(LINK_DIR_1 ../lib/win/${MY_BUILD_TYPE})
    set(LINK_DIR_2 ../lib_external/win/${MY_BUILD_TYPE})
ELSEIF(APPLE)
    set(LINK_DIR_1 ../lib_external/mac)
ELSE()
    IF(os_version_suffix MATCHES "\\.el[1-9]")
        set(LINK_DIR_1 ../lib_external/linux_fedora)
    ELSE()
        set(LINK_DIR_1 ../lib_external/linux_debian)
    ENDIF()
ENDIF()
target_link_directories(testsuite PRIVATE ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
target_link_directories(testsuite PRIVATE ${LINK_DIR_1} ${LINK_DIR_2})
if(MSVC)
    # zip extraction
    target_link_libraries(testsuite PRIVATE
            #bzip2.lib
            #lzma.lib
            #zlib.lib
            )
endif(MSVC)

target_compile_features(testsuite PRIVATE cxx_std_17)
#add_test(NAME example_test COMMAND testsuite)

file(COPY ./simulation.cfg
        ./TestCases
        ../copy_to_build/parameter_catalog
        DESTINATION ${CMAKE_EXECUTABLE_OUTPUT_DIRECTORY}/tests)