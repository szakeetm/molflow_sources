cmake_minimum_required(VERSION 3.12.2 FATAL_ERROR)

################### Variables. ####################
# Change if you want modify path or other values. #
###################################################
set(PROJECT_NAME gpuSub)
project(${PROJECT_NAME} CUDA CXX)

# Folders files

set(HEADER_DIR_1 ../../src)
set(HEADER_DIR_2 ../../src_shared)
set(HEADER_DIR_3 ../../include)
set(HEADER_DIR_4 ../../include/windows_only)
set(HEADER_DIR_5 ../../include/common_cuda)

set(GPU_DIR_1 ../../src/GPUSim)
set(CUDA_DIR_1 ../../src/GPUSim/CUDA)

IF (WIN32)
    # set stuff for windows
    set(LINK_DIR_1 ../../lib/win/${MY_BUILD_TYPE})
    #set(LINK_DIR_2 ../../lib_external/win/${MY_BUILD_TYPE})
ELSE()
    # set stuff for other systems
ENDIF()

# Definition of Macros
add_definitions(
        -DWIN
        -DMOLFLOW
        -D_MBCS
        -DGPU_MODE
        -D_CRT_SECURE_NO_WARNINGS
)

option(DEBUG_MISS "Enable Debugging for Miss Program" OFF)
option(DEBUG_COUNT "Enable counters for intersection results" OFF)
option(DEBUG_POS "Enable generating hit positions from GPU calculations" OFF)
option(WITH_TRIANGLES "Enable calculations with triangles only" ON)
option(WITH_TEXTURES "Enable textures" ON)
option(WITH_PROFILES "Enable profiles" ON)

#option(USE_RANDOM_NUMBER_TYPE_64 "Use double instead of float for random numbers" ON)

IF (CMAKE_BUILD_TYPE STREQUAL "Debug")
    if (DEBUG_MISS)
        add_definitions(-DDEBUGMISS)
    endif (DEBUG_MISS)

    if (DEBUG_COUNT)
        add_definitions(-DDEBUGCOUNT)
    endif (DEBUG_COUNT)

    if (DEBUG_POS)
        add_definitions(-DDEBUGPOS)
    endif (DEBUG_POS)
ENDIF()
if (WITH_TRIANGLES)
    add_definitions(-DWITHTRIANGLES)
endif (WITH_TRIANGLES)
if (WITH_TEXTURES)
    add_definitions(-DWITH_TEX)
endif (WITH_TEXTURES)
if (WITH_PROFILES)
    add_definitions(-DWITH_PROF)
endif (WITH_PROFILES)
if (USE_RANDOM_NUMBER_TYPE_64)
    add_definitions(-DRNG64)
endif (USE_RANDOM_NUMBER_TYPE_64)

################# Flags ################
# Defines Flags for Windows and Linux. #
########################################

if(MSVC)

endif(MSVC)
if(NOT MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
    endif()
endif(NOT MSVC)

# create sub proc executable
include(${PROJECT_SOURCE_DIR}/configure_optix.cmake)

set(optix_LIBRARY "")
#message("Optix include " ${OptiX_INCLUDE})
include_directories(SYSTEM ${OptiX_INCLUDE} )

add_executable(${PROJECT_NAME}
        ${GPU_DIR_1}/gpuSub.cpp
        ${GPU_DIR_1}/SimulationGPU.cpp
        ${GPU_DIR_1}/SimulationControllerGPU.cpp
        ${GPU_DIR_1}/ModelReader.cpp
        ${GPU_DIR_1}/Poly2TriConverter.cpp
        ${HEADER_DIR_1}/GeometrySimu.cpp
        ${HEADER_DIR_1}/MolflowBuffer.cpp
        ${HEADER_DIR_1}/Parameter.cpp
        )

target_include_directories(${PROJECT_NAME} PRIVATE
        ${GPU_DIR_1}
        ${CUDA_DIR_1}
        ${common_cuda_dir}
        ${HEADER_DIR_1}
        ${HEADER_DIR_2}
        ${HEADER_DIR_3}
        ${HEADER_DIR_5}
        )

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)
target_link_libraries(${PROJECT_NAME}
        flowGPU
        simulator_shared
        )

if(MSVC)
    set(LINK_DIR_2 ../../lib_external/win/${MY_BUILD_TYPE})
    target_link_directories( ${PROJECT_NAME} PRIVATE ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
    target_link_directories(${PROJECT_NAME} PRIVATE ${LINK_DIR_1})
    target_link_directories(${PROJECT_NAME} PUBLIC
            ${LINK_DIR_2}
    )
endif(MSVC)
if(MSVC)
    # Add Windows Console output
    set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS "/SUBSYSTEM:CONSOLE ${LINK_FLAGS}")
endif(MSVC)