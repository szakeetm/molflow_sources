cmake_minimum_required(VERSION 3.12.2 FATAL_ERROR)

################### Variables. ####################
# Change if you want modify path or other values. #
###################################################

set(PROJECT_NAME flowGPU)

message(${PROJECT_NAME} " generating in " ${CMAKE_CURRENT_SOURCE_DIR})


# Folders files

set(HEADER_DIR_1 ../../src)
set(HEADER_DIR_2 ../../src_shared)
set(HEADER_DIR_3 ../../include)
set(HEADER_DIR_4 ../../include/windows_only)
set(HEADER_DIR_5 ../../include/common_cuda)

set(GPU_DIR_1 ${CMAKE_SOURCE_DIR}/src/GPUSim)
set(CUDA_DIR_1 ${CMAKE_SOURCE_DIR}/src/GPUSim/CUDA)

IF (WIN32)
    # set stuff for windows
    set(LINK_DIR_1 ../../lib/win/${MY_BUILD_TYPE})
    #set(LINK_DIR_2 ../../lib_external/win/${MY_BUILD_TYPE})
ELSE()
    # set stuff for other systems
ENDIF()



############## CMake Project ################
#        The main options of project        #
#############################################

project(${PROJECT_NAME} CUDA CXX)

# Find at least a 5.0 version of CUDA.
#set(CUDA_ROOT "C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v10.2")
#set(CUDA_TOOLKIT_ROOT_DIR "C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v10.2")
#set(CMAKE_CUDA_COMPILER "C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v10.2")

#find_package(CUDA 10.2 REQUIRED)
message("CUDA found as " ${CUDA_TOOLKIT_ROOT_DIR})

# Definition of Macros
# ... via MolflowGPU.cmake

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

#find_library(PROTOBUF_LIBRARY bzip2 HINTS ${LINK_DIR_2})


#message(${EXECUTABLE_OUTPUT_PATH})

# Add paths to our CMake code to the module path, so they can be found automatically by
# CMake.
#[[set(CMAKE_MODULE_PATH
        "${CMAKE_SOURCE_DIR}"
        ${CMAKE_MODULE_PATH}
        )
message("Module Path is " ${CMAKE_MODULE_PATH})]]

#set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -g -G")

set(SRC_FILES ${SRC_FILES}
        ${GPU_DIR_1}/SimulationControllerGPU.cpp
        ${GPU_DIR_1}/SimulationGPU.cpp

        ${GPU_DIR_1}/SimulationOptiX.cpp
        ${GPU_DIR_1}/Model.cpp
        #${GPU_DIR_1}/RNGKernel.cu
        #${GPU_DIR_1}/RandomWrapper.cpp
        )


include_directories(${GPU_DIR_1})
include_directories(${CUDA_DIR_1})

set(common_cuda_dir ${PROJECT_SOURCE_DIR}/../../include/common_cuda/)
include_directories(SYSTEM ${common_cuda_dir})
set(common_cuda_dir ${PROJECT_SOURCE_DIR}/../../include/common_cuda/)
include_directories(SYSTEM ${HEADER_DIR_3})

include(${PROJECT_SOURCE_DIR}/configure_optix.cmake)

set(optix_LIBRARY "")
#message("Optix include " ${OptiX_INCLUDE})
include_directories(SYSTEM ${OptiX_INCLUDE})

cuda_compile_and_embed(geometry_ptx_code ${CUDA_DIR_1}/Geometry.cu)
cuda_compile_and_embed(trace_ptx_code ${CUDA_DIR_1}/TraceProcessing.cu)
cuda_compile_and_embed(ray_ptx_code ${CUDA_DIR_1}/RayGeneration.cu)

#cuda_compile_and_embed(rng2_ptx_code ${GPU_DIR_1}/cudaRandom.cu)

set(SRC_FILES ${geometry_ptx_code} ${SRC_FILES})
set(SRC_FILES ${trace_ptx_code} ${SRC_FILES})
set(SRC_FILES ${ray_ptx_code} ${SRC_FILES})

add_library(${PROJECT_NAME} ${SRC_FILES})

#set_target_properties(${PROJECT_NAME} PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)

target_include_directories(${PROJECT_NAME} PRIVATE
        ${GPU_DIR_1}
        ${CUDA_DIR_1}
        ${common_cuda_dir}
        ${HEADER_DIR_2} # shared buffer for model
        ${HEADER_DIR_3}
        )

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)

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


target_link_libraries(${PROJECT_NAME}
        cudaKernel
        ${optix_LIBRARY}
        ${CUDA_LIBRARIES}
        ${CUDA_CUDA_LIBRARY}
        )


# Now create executable
set(PROJECT_NAME gpuTestsuite)
set(INCLUDE_DIRECTORIES "")
add_executable(${PROJECT_NAME}
        ${GPU_DIR_1}/TestSimulator.cpp
        ${GPU_DIR_1}/ModelReader.cpp
        ${HEADER_DIR_1}/GeometrySimu.cpp
        ${HEADER_DIR_1}/MolflowBuffer.cpp
        ${HEADER_DIR_1}/Parameter.cpp
        ${HEADER_DIR_2}/Vector.cpp
        ${GPU_DIR_1}/Poly2TriConverter.cpp
        ${GPU_DIR_1}/SimulationControllerGPU.cpp
        )

target_include_directories(${PROJECT_NAME} PRIVATE
        ${GPU_DIR_1}
        ${CUDA_DIR_1}
        ${common_cuda_dir}
        #${HEADER_DIR_1}
        ${HEADER_DIR_2}
        #${HEADER_DIR_3}
        #${HEADER_DIR_4}
        )

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)
target_link_libraries(${PROJECT_NAME}
        flowGPU
        )

if(MSVC)
    # Add Windows Console output
    set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS "/SUBSYSTEM:CONSOLE ${LINK_FLAGS}")
endif(MSVC)