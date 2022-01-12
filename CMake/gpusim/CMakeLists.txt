cmake_minimum_required(VERSION 3.12.2 FATAL_ERROR)

################### Variables. ####################
# Change if you want modify path or other values. #
###################################################

set(PROJECT_NAME gpuTestsuite)
project(${PROJECT_NAME} CUDA CXX)

message(${PROJECT_NAME} " generating in " ${CMAKE_CURRENT_SOURCE_DIR})


# Folders files

set(HEADER_DIR_1 ../../src)
set(HEADER_DIR_2 ../../src_shared)
set(HEADER_DIR_3 ../../include)
set(HEADER_DIR_4 ../../include/windows_only)
set(HEADER_DIR_5 ../../include/common_cuda)

set(CPP_DIR_1 ../../src/)
set(CPP_DIR_2 ../../src_shared/)
set(SIMU_DIR ../../src/Simulation)
set(IO_DIR ../../src/IO)

set(GPU_DIR_1 ../../src/GPUSim)
set(CUDA_DIR_1 ../../src/GPUSim/CUDA)

IF (WIN32)
    # set stuff for windows
    set(LINK_DIR_1 ../../lib/win/${MY_BUILD_TYPE})
    #set(LINK_DIR_2 ../../lib_external/win/${MY_BUILD_TYPE})
ELSE()
    # set stuff for other systems
ENDIF()

include(${PROJECT_SOURCE_DIR}/configure_optix.cmake)

set(optix_LIBRARY "")
# Now create executable
set(INCLUDE_DIRECTORIES "")
add_executable(${PROJECT_NAME}

        ${CPP_DIR_1}/MolflowTypes.cpp
        ${CPP_DIR_2}/FlowMPI.h #contains templates
        ${CPP_DIR_2}/FlowMPI.cpp
        ${CPP_DIR_2}/File.cpp
        ${CPP_DIR_1}/Initializer.cpp
        ${CPP_DIR_1}/ParameterParser.cpp

        ${SIMU_DIR}/Simulation.cpp
        ${SIMU_DIR}/SimulationMC.cpp
        ${SIMU_DIR}/Particle.cpp
        ${SIMU_DIR}/Physics.cpp
        ${SIMU_DIR}/AnglemapGeneration.cpp
        ${SIMU_DIR}/CDFGeneration.cpp
        ${SIMU_DIR}/IDGeneration.cpp
        ${CPP_DIR_2}/IntersectAABB_shared.cpp #isinfacet
        ${CPP_DIR_2}/Random.cpp #TruncatedGaussian::GetGaussian --> rtnorm

        ${IO_DIR}/LoaderXML.cpp
        ${IO_DIR}/WriterXML.cpp
        ${CPP_DIR_1}/TimeMoments.cpp

        ${GPU_DIR_1}/TestSimulator.cpp
        ${GPU_DIR_1}/ModelReader.cpp
        ${GPU_DIR_1}/SimulationGPU.cpp
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
        ${HEADER_DIR_1}
        ${HEADER_DIR_2}
        ${HEADER_DIR_3}
        ${HEADER_DIR_5}
        )

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)
target_link_libraries(${PROJECT_NAME}
        flowGPU
        shared_simulator
        )
target_link_libraries(${PROJECT_NAME} pugixml truncatedgaussian)

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