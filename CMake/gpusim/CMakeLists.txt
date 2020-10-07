cmake_minimum_required(VERSION 3.12.2 FATAL_ERROR)

################### Variables. ####################
# Change if you want modify path or other values. #
###################################################

set(PROJECT_NAME gpuTestsuite)

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

# Now create executable
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