cmake_minimum_required(VERSION 3.12.2 FATAL_ERROR)

project(molflow_proj)

message("-- CMAKE_SYSTEM_INFO_FILE: ${CMAKE_SYSTEM_INFO_FILE}")
message("-- CMAKE_SYSTEM_NAME:      ${CMAKE_SYSTEM_NAME}")
message("-- CMAKE_SYSTEM_PROCESSOR: ${CMAKE_SYSTEM_PROCESSOR}")
message("-- CMAKE_SYSTEM:           ${CMAKE_SYSTEM}")

string (REGEX MATCH "\\.el[1-9]" os_version_suffix ${CMAKE_SYSTEM})
message("-- os_version_suffix:      ${os_version_suffix}")

#list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
option(USE_CLANG "build application with clang" OFF) # OFF is the default
option(USE_PROFILING "disable optimisation for profiling the application" OFF) # OFF is the default
option(WITH_GPU "activate GPU build" ON) # OFF is the default

if(WIN32 AND MSVC)
    message("Creating GPU build")
    set(WITH_GPU ON)
else()
    message("GPU build not available with this compiler and/or OS")
    #set(WITH_GPU OFF)
endif(WIN32 AND MSVC)

# Define Release by default.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
    message(STATUS "Build type not specified: Use Release by default.")
endif(NOT CMAKE_BUILD_TYPE)
if(CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    MESSAGE("Adding Debug flag...")
    ADD_DEFINITIONS(-DDEBUG)
    ADD_DEFINITIONS(-D_DEBUG)
endif()
message(STATUS "Build type used: " ${CMAKE_BUILD_TYPE})
set(CMAKE_VERBOSE_MAKEFILE OFF)

################# Flags ################
# Defines Flags for Windows and Linux. #
########################################

# https://docs.microsoft.com/en-us/cpp/build/reference/compiler-options-listed-alphabetically?view=vs-2019
if(MSVC)
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++17")

    #set(CMAKE_GENERATOR_PLATFORM x64)
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /mwindows")

    #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mwindows")
    # main vs WinMain
    # https://stackoverflow.com/questions/2752792/whats-the-equivalent-of-gccs-mwindows-option-in-cmake
    #set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /subsystem:windows /ENTRY:mainCRTStartup")
    #set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /W3 /EHsc")
    #set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /W3 /EHsc")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /W4 /MDd /Od /EHsc")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO  "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} /W4 /GL /O2 /Ob2 /EHsc")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /W4 /GL /O2 /EHsc")
endif(MSVC)
if(NOT MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
	#Flag to set debug symbols
	set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
	
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
    endif()
endif(NOT MSVC)

if(USE_PROFILING)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")

    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
    SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")
endif(USE_PROFILING)

set(CMAKE_EXPORT_COMPILE_COMMANDS OFF)

include(CMake/Molflow.cmake)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake/")
add_subdirectory(src_shared)
add_subdirectory(CMake/molflow_sub)
add_subdirectory(CMake/molflow_win)

IF (WIN32)
    add_subdirectory(CMake/png_zlib_win_library)
ELSE()
    #add_subdirectory(CMake/molflow_testsuit)
ENDIF(WIN32)

add_dependencies(molflow molflowSub compress)

option(USE_TESTS "Build the tests" ON)
if(USE_TESTS)
    enable_testing()
    include(GoogleTest)
    add_subdirectory(tests)
    add_dependencies(example molflowSub)
endif()

IF(WITH_GPU AND (MSVC AND WIN32) OR (NOT WIN32))
#IF(MSVC AND WITH_GPU)
    #set(cuda_obj_dir ${PROJECT_SOURCE_DIR}/cuda_kernel/)

    # load defines
    include(CMake/MolflowGPU.cmake)

    add_subdirectory(CMake/molflow_gpusub)

    add_subdirectory(CMake/gpu_testsuite)
    add_subdirectory(CMake/cuda_kernel #[["${cuda_obj_dir}"]]) #add_subdirectory(${cuda_obj_dir} "${cuda_obj_dir}" EXCLUDE_FROM_ALL) #already added in molflow_sub

    add_dependencies(gpuTestsuite flowGPU cudaKernel)
    add_dependencies(gpuSub flowGPU cudaKernel)
    add_dependencies(molflow gpuSub)

ENDIF()