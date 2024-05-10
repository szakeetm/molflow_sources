cmake_minimum_required(VERSION 3.12.2 FATAL_ERROR)

project(molflow_proj) #does not build anything on its own, includes molflowCLI or molflow projects from the cmake_modules folder

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
    ADD_DEFINITIONS(-DENABLE_IMGUI_TESTS)
    MESSAGE("Building ImGui tests.")
endif()
message(STATUS "Build type used: " ${CMAKE_BUILD_TYPE})
set(CMAKE_VERBOSE_MAKEFILE OFF)


if(USE_KDTREE)
    MESSAGE("Will use K-D Tree (adding USE_KDTREE macro)")
    ADD_DEFINITIONS(-DUSE_KDTREE)
else()
    MESSAGE("Will use new BVH (default)")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake/")

include(cmake_modules/molflow_proj.cmake) # Technical details offloaded to keep this root CMakeLists.txt file easily readable

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