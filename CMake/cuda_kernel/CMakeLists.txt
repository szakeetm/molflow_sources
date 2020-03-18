cmake_minimum_required(VERSION 3.12.2 FATAL_ERROR)

################### Variables. ####################
# Change if you want modify path or other values. #
###################################################

set(PROJECT_NAME cudaKernel)
#set(CMAKE_VERBOSE_MAKEFILE ON)
#set(CUDA_VERBOSE_BUILD ON)

IF (WIN32)
    set(OS_NAME "win")
ELSE()
    IF(os_version_suffix STREQUAL ".el7")
        set(OS_NAME "linux_fedora")
    ELSE()
        set(OS_NAME "linux_debian")
    ENDIF()
ENDIF()

IF (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(MY_BUILD_TYPE "debug")
ELSE()
    set(MY_BUILD_TYPE "release")
ENDIF()

# Output Variables
set(OUTPUT_DEBUG ../lib/${OS_NAME}/debug/)
set(OUTPUT_REL ../lib/${OS_NAME}/release/)

# Folders files
#set(CPP_DIR_1 ../_molflow/GPUSim)

#set(HEADER_DIR_1 ../_molflow/Source_files)
#set(HEADER_DIR_1 ${CMAKE_SOURCE_DIR}/include/cuda)
#set(HEADER_DIR_2 ${CMAKE_SOURCE_DIR}/include)

set(GPU_DIR_1 ${CMAKE_SOURCE_DIR}/_molflow/Source_files/GPUSim)
set(CUDA_DIR_1 ${CMAKE_SOURCE_DIR}/_molflow/Source_files/GPUSim/CUDA)

############## CMake Project ################
#        The main options of project        #
#############################################

project(${PROJECT_NAME} CUDA CXX)

# Find at least a 5.0 version of CUDA.
#set(CUDA_ROOT "C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v10.2")
#set(CUDA_TOOLKIT_ROOT_DIR "C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v10.2")
#set(CMAKE_CUDA_COMPILER "C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v10.2")

#includes get added with the find_package command
#include_directories(${GPU_DIR_1})
#include_directories(${CUDA_DIR_1})

message("CUDA found as " ${CUDA_TOOLKIT_ROOT_DIR})
#set(CUDA_INCLUDE_DIR "${CUDA_TOOLKIT_ROOT_DIR}/include")

# Definition of Macros
add_definitions(
        -D_MBCS
        -DGPU_MODE
        -D_CRT_SECURE_NO_WARNINGS
)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    MESSAGE("[NVCC] Adding Debug flag...")
    ADD_DEFINITIONS(-DDEBUG)
    ADD_DEFINITIONS(-D_DEBUG)
endif()


############## Artefacts Output #################
# Defines outputs , depending Debug or Release. #
#################################################

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${OUTPUT_DEBUG}")
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${OUTPUT_DEBUG}")
    set(CMAKE_EXECUTABLE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${OUTPUT_DEBUG}")
else()
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${OUTPUT_REL}")
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${OUTPUT_REL}")
    set(CMAKE_EXECUTABLE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${OUTPUT_REL}")
endif()

# Messages
message("${PROJECT_NAME}: MAIN PROJECT: ${CMAKE_PROJECT_NAME}")
message("${PROJECT_NAME}: CURR PROJECT: ${CMAKE_CURRENT_SOURCE_DIR}")
message("${PROJECT_NAME}: CURR BIN DIR: ${CMAKE_CURRENT_BINARY_DIR}")
# Messages
message("${PROJECT_NAME}: CMAKE_LIBRARY_OUTPUT_DIRECTORY: ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
message("${PROJECT_NAME}: CMAKE_ARCHIVE_OUTPUT_DIRECTORY: ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}")
message("${PROJECT_NAME}: CMAKE_EXECUTABLE_OUTPUT_DIRECTORY: ${CMAKE_EXECUTABLE_OUTPUT_DIRECTORY}")

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

################ Files ################
#   --   Add files to project.   --   #
#######################################


#set(CUDA_SEPARABLE_COMPILATION ON)
#list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_30,code=sm_30 -rdc=true -std=c++11)

set(SRC_FILES ${SRC_FILES}
        ${CUDA_DIR_1}/RNGKernel.cu
        ${CUDA_DIR_1}/cudaRandom.cu
        #${GPU_DIR_1}/RandomWrapper.cpp
        )

set(HEADER_FILES ${HEADER_FILES}
        ${GPU_DIR_1}/RandomWrapper.h
        ${GPU_DIR_1}/LaunchParams.h
        ${GPU_DIR_1}/GPUDefines.h

        ${CUDA_DIR_1}/cudaRandom.cuh
        #${GPU_DIR_1}/jetbrains_indexing.h
        )

add_library(${PROJECT_NAME} STATIC ${SRC_FILES} )

#add_library(${PROJECT_NAME} STATIC ${CUDA_SRC_FILES} ${HEADER_FILES})
#set_target_properties(${PROJECT_NAME} PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
#set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)
#set_target_properties(${PROJECT_NAME} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

set(CUDA_INCLUDE_DIR "${CUDA_TOOLKIT_ROOT_DIR}/include")
target_include_directories(${PROJECT_NAME} PRIVATE
        ${GPU_DIR_1}
        ${CUDA_DIR_1}
        ${CUDA_INCLUDE_DIR}
        )

#target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)



#Flags on CUDA 10 for maximum compatibility
#set(CMAKE_CUDA_FLAGS ${CMAKE_CUDA_FLAGS} " -arch=sm_50 -gencode=arch=compute_50,code=sm_50 -gencode=arch=compute_52,code=sm_52 -gencode=arch=compute_60,code=sm_60 -gencode=arch=compute_61,code=sm_61 -gencode=arch=compute_70,code=sm_70 -gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_75,code=compute_75")
#set(CMAKE_CUDA_FLAGS "-arch=sm_52 -gencode=arch=compute_52,code=sm_52 -gencode=arch=compute_60,code=sm_60 -gencode=arch=compute_61,code=sm_61 -gencode=arch=compute_70,code=sm_70 -gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_75,code=compute_75")
set(cuda_flags -gencode=arch=compute_70,code=sm_70 -gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_75,code=compute_75)

# Enable fast math
#if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    option(ODL_CUDA_USE_FAST_MATH "Enable fast math in cuda (can decrease precision)" FALSE)
    if(ODL_CUDA_USE_FAST_MATH)
        MESSAGE("[NVCC] Using fast math...")
        set(cuda_flags ${cuda_flags} --use_fast_math)
    endif()
    option(ODL_CUDA_STREAM_PER_THREAD "Enable stream per-thread in cuda" TRUE)
    if(ODL_CUDA_STREAM_PER_THREAD)
        MESSAGE("[NVCC] Using stream per-thread...")
        set(cuda_flags ${cuda_flags} --default-stream per-thread)
    endif()
    MESSAGE("[NVCC] Flags: " ${cuda_flags})
#endif()

#set(CMAKE_CUDA_FLAGS "-arch=sm_70 -gencode=arch=compute_70,code=sm_70 -gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_75,code=compute_75 --default-stream per-thread")
#target_link_options(${PROJECT_NAME} PUBLIC -arch=sm_70 -gencode=arch=compute_70,code=sm_70 -gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_75,code=compute_75)

# https://docs.nvidia.com/cuda/pascal-compatibility-guide/index.html
target_compile_options(${PROJECT_NAME} PUBLIC
        #$<$<COMPILE_LANGUAGE:CUDA>:-arch=sm_50>
        # Maxwell Gen1
        #$<$<COMPILE_LANGUAGE:CUDA>:-gencode=arch=compute_50,code=sm_50>
        # Maxwell Gen2
        #$<$<COMPILE_LANGUAGE:CUDA>:-gencode=arch=compute_52,code=sm_52>
        # Pascal P100
        #$<$<COMPILE_LANGUAGE:CUDA>:-gencode=arch=compute_60,code=sm_60>
        # Pascal
        #$<$<COMPILE_LANGUAGE:CUDA>:-gencode=arch=compute_61,code=sm_61>
        # Volta
        #$<$<COMPILE_LANGUAGE:CUDA>:-gencode=arch=compute_70,code=sm_70>
        # Turing
        #$<$<COMPILE_LANGUAGE:CUDA>:-gencode=arch=compute_75,code=compute_75>
        $<$<COMPILE_LANGUAGE:CUDA>:${cuda_flags}>
        )

#[[message("CUDA_LIBRARIES: " ${CUDA_LIBRARIES})
message("CUDA_CUDA_LIBRARY: " ${CUDA_CUDA_LIBRARY})
message("CUDA_curand_LIBRARY: " ${CUDA_curand_LIBRARY})]]
target_link_libraries(${PROJECT_NAME}
        ${CUDA_LIBRARIES}
        ${CUDA_CUDA_LIBRARY}
        ${CUDA_curand_LIBRARY}
        )