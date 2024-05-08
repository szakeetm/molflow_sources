cmake_minimum_required(VERSION 3.12.2 FATAL_ERROR)

project(molflow_proj)

message("-- CMAKE_SYSTEM_INFO_FILE: ${CMAKE_SYSTEM_INFO_FILE}")
message("-- CMAKE_SYSTEM_NAME:      ${CMAKE_SYSTEM_NAME}")
message("-- CMAKE_SYSTEM_PROCESSOR: ${CMAKE_SYSTEM_PROCESSOR}")
message("-- CMAKE_SYSTEM:           ${CMAKE_SYSTEM}")

#Use new BVH by default
option(USE_KDTREE "simulation will use kd structure instead of bvh" OFF)
option(USE_CLANG "build application with clang" OFF) # OFF is the default
option(USE_PROFILING "disable optimisation for profiling the application" OFF) # OFF is the default
option(ENABLE_IMGUI_TESTS "build ImGui graphical tests" OFF)
option(NO_INTERFACE "only build CLI binary and dependencies" OFF) # OFF is the default
if(NO_INTERFACE)
    message(STATUS "Will not build GUI related libraries and executables.")
endif()

# Define Release by default.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
    message(STATUS "Build type not specified, will use Release by default.")
endif(NOT CMAKE_BUILD_TYPE)
if(CMAKE_BUILD_TYPE MATCHES Debug OR CMAKE_BUILD_TYPE MATCHES RelWithDebInfo) # case sensitive "Debug" and "RelWithDebInfo"
    ADD_DEFINITIONS(-DDEBUG)
    ADD_DEFINITIONS(-D_DEBUG)
    MESSAGE("Added DEBUG preprocessor macro.")
endif()
if(ENABLE_IMGUI_TESTS)
    ADD_DEFINITIONS(-DENABLE_IMGUI_TESTSS)
    MESSAGE("Building ImGui tests.")
endif()
message(STATUS "Build type used: " ${CMAKE_BUILD_TYPE})
set(CMAKE_VERBOSE_MAKEFILE OFF)

################# Flags ################
# Defines Flags for Windows and Linux. #
########################################

if(USE_KDTREE)
    MESSAGE("Will use K-D Tree (adding USE_KDTREE macro)")
    ADD_DEFINITIONS(-DUSE_KDTREE)
else()
    MESSAGE("Will use new BVH (default)")
endif()

set(CMAKE_CXX_STANDARD 17)

# https://docs.microsoft.com/en-us/cpp/build/reference/compiler-options-listed-alphabetically?view=vs-2019
if(MSVC)
    if(CMAKE_BUILD_TYPE MATCHES Debug)
        string(REGEX REPLACE "/Zi" "/ZI" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}") #Edit and continue (/ZI) instead of default /Zi
    endif()
else() #not MSVC
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        string(APPEND CMAKE_CXX_FLAGS " -stdlib=libc++") #instead of stdlibc++
    endif()
endif()

if(USE_PROFILING)
    MESSAGE("Setting profiling flags.")
    string(APPEND CMAKE_C_FLAGS " -O0 -pg")
    string(APPEND CMAKE_CXX_FLAGS " -O0 -pg")
    string(APPEND CMAKE_EXE_LINKER_FLAGS " -pg")
    string(APPEND CMAKE_SHARED_LINKER_FLAGS " -pg")
endif(USE_PROFILING)

set(CMAKE_EXPORT_COMPILE_COMMANDS OFF) # Don't generate .json files

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake/")
include(CMake/Molflow.cmake)

# Project-wide folders
set(CPP_DIR_SRC ${CMAKE_HOME_DIRECTORY}/src)
set(SIMU_DIR ${CMAKE_HOME_DIRECTORY}/src/Simulation)
set(IO_DIR ${CMAKE_HOME_DIRECTORY}/src/IO)
set(CPP_DIR_SRC_SHARED ${CMAKE_HOME_DIRECTORY}/src_shared)
set(GLAPP_DIR ${CMAKE_HOME_DIRECTORY}/src_shared/GLApp)
set(HELPER_DIR ${CMAKE_HOME_DIRECTORY}/src_shared/Helper)

set(HEADER_DIR_SRC ${CMAKE_HOME_DIRECTORY}/src)
set(HEADER_DIR_SRC_SHARED ${CMAKE_HOME_DIRECTORY}/src_shared)
set(HEADER_DIR_INCLUDE ${CMAKE_HOME_DIRECTORY}/include)

add_subdirectory(src_shared EXCLUDE_FROM_ALL)
add_subdirectory(cmake_modules/molflow_core)
add_subdirectory(cmake_modules/molflow_cli)

if(NOT NO_INTERFACE)
    add_subdirectory(cmake_modules/molflow_gui)
    add_dependencies(molflow src_shared_core src_shared_gui compress)
endif()

add_dependencies(molflowCLI molflow_core)

option(USE_TESTS "Build the tests" OFF)
if(USE_TESTS)
    enable_testing()
    include(GoogleTest)
    add_subdirectory(tests)
    add_dependencies(testsuite molflow_core molflowCLI)
endif()