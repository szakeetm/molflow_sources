cmake_minimum_required(VERSION 3.12.2 FATAL_ERROR)

project(molflow_proj)

message("-- CMAKE_SYSTEM_INFO_FILE: ${CMAKE_SYSTEM_INFO_FILE}")
message("-- CMAKE_SYSTEM_NAME:      ${CMAKE_SYSTEM_NAME}")
message("-- CMAKE_SYSTEM_PROCESSOR: ${CMAKE_SYSTEM_PROCESSOR}")
message("-- CMAKE_SYSTEM:           ${CMAKE_SYSTEM}")

string (REGEX MATCH "\\.el[1-9]" os_version_suffix ${CMAKE_SYSTEM})
message("-- os_version_suffix:      ${os_version_suffix}")

#list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

#Use new BVH by default
option(USE_KDTREE "simulation will use kd structure instead of bvh" OFF)

option(USE_CLANG "build application with clang" OFF) # OFF is the default
option(USE_PROFILING "disable optimisation for profiling the application" OFF) # OFF is the default

option(NO_INTERFACE "only build CLI binary and dependencies" OFF) # OFF is the default
if(NO_INTERFACE)
    message(STATUS "Will not build GUI related libraries and executables.")
endif()

# Define Release by default.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
    message(STATUS "Build type not specified, will use Release by default.")
endif(NOT CMAKE_BUILD_TYPE)
if(CMAKE_BUILD_TYPE MATCHES Debug) #case sensitive "Debug"!
    ADD_DEFINITIONS(-DDEBUG)
    ADD_DEFINITIONS(-D_DEBUG)
    MESSAGE("Debug build mode, added DEBUG preprocessor macro.")
endif()
message(STATUS "Build type used: " ${CMAKE_BUILD_TYPE})
set(CMAKE_VERBOSE_MAKEFILE OFF)

################# Flags ################
# Defines Flags for Windows and Linux. #
########################################

if(USE_KDTREE)
    MESSAGE("Simulation with KDTREE...")
    ADD_DEFINITIONS(-DUSE_KDTREE)
else()
    MESSAGE("Simulation with new BVH (default)...")
    #ADD_DEFINITIONS(-DUSE_KDTREE)
endif()

# https://docs.microsoft.com/en-us/cpp/build/reference/compiler-options-listed-alphabetically?view=vs-2019
if(MSVC)
    if(CMAKE_BUILD_TYPE MATCHES Debug)
        string(REGEX REPLACE "/Zi" "/ZI" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}") #Edit and continue (/ZI) instead of default /Zi
    endif()
else() #not MSVC
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        string(APPEND CMAKE_CXX_FLAGS " -stdlib=libc++")
    endif()
endif()

if(USE_PROFILING)
    string(APPEND CMAKE_C_FLAGS " -O0 -pg")
    string(APPEND CMAKE_CXX_FLAGS " -O0 -pg")
    string(APPEND CMAKE_EXE_LINKER_FLAGS " -pg")
    string(APPEND CMAKE_SHARED_LINKER_FLAGS " -pg")
endif(USE_PROFILING)

set(CMAKE_EXPORT_COMPILE_COMMANDS OFF)

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
add_subdirectory(CMake/molflow_core)
add_subdirectory(CMake/molflow_cli)


if(NOT NO_INTERFACE)
    add_subdirectory(CMake/molflow_gui)
    add_dependencies(molflow src_shared_simulator src_shared_interface compress)
endif()

add_dependencies(molflow_core ziplib)
add_dependencies(molflowCLI molflow_core)

option(USE_TESTS "Build the tests" OFF)
if(USE_TESTS)
    enable_testing()
    include(GoogleTest)
    add_subdirectory(tests)
    add_dependencies(testsuite molflow_core molflowCLI)
endif()

# Allow to install molflowCLI into Unix OS folders (e.g. /usr/local/bin)
# molflow GUI not supported right now, as extra folders are needed

#SET(MOLFLOW_INSTALL_DIR ${CMAKE_INSTALL_BINDIR})
SET(MOLFLOW_INSTALL_DIR molflow)

include(GNUInstallDirs)

# if no path is provided, use user home path
# e.g. /home/user/molflow
# cmake --install . --prefix /usr/local/bin/
# will install to /usr/local/bin/molflow
# cmake -DCMAKE_INSTALL_PREFIX:PATH=/Users/pascal/molflow_test/ .
#
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "$ENV{HOME}" CACHE PATH "..." FORCE)
endif()

install(CODE "MESSAGE(\"Installing in home directory ${CMAKE_INSTALL_PREFIX}\")")
install(DIRECTORY copy_to_build/images
        copy_to_build/desorption_yields
        copy_to_build/parameter_catalog
        DESTINATION ${MOLFLOW_INSTALL_DIR}
        #USE_SOURCE_PERMISSIONS
        PATTERN "images/*"
        PATTERN "desorption_yields/*"
        PATTERN "parameter_catalog/*"
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ GROUP_WRITE
        )

if(NOT NO_INTERFACE)
    install(TARGETS molflow
            CONFIGURATIONS Release
            RUNTIME
            DESTINATION ${MOLFLOW_INSTALL_DIR}
    )
endif()
install(TARGETS molflowCLI
        CONFIGURATIONS Release
        RUNTIME
        DESTINATION ${MOLFLOW_INSTALL_DIR}
        )
