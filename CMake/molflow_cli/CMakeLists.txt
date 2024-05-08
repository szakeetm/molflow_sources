cmake_minimum_required(VERSION 3.12.2 FATAL_ERROR)

#Command-line interface, builds molflowCLI executable
set(PROJECT_NAME molflowCLI)

############## CMake Project ################
#        The main options of project        #
#############################################

project(${PROJECT_NAME} CXX)
set(CMAKE_CXX_STANDARD 17)

# Definition of Macros
add_definitions(
        -D_MBCS
        -D_CRT_SECURE_NO_WARNINGS
)

################# Flags ################
# Defines Flags for Windows and Linux. #
########################################

if(MSVC)

else() #not msvc
    if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
    endif()
endif()

################ Files ################
#   --   Add files to project.   --   #
#######################################

set(SRC_FILES
        ${CPP_DIR_SRC}/MolflowCLI.cpp
        ${CPP_DIR_SRC_SHARED}/SimulationManager.cpp
        ${CPP_DIR_SRC_SHARED}/ShMemory.cpp
        ${CPP_DIR_SRC_SHARED}/Process.cpp

        ${CPP_DIR_SRC_SHARED}/Vector.cpp
        ${CPP_DIR_SRC_SHARED}/Distributions.cpp

        ${IO_DIR}/LoaderXML.cpp
        ${IO_DIR}/WriterXML.cpp
        ${CPP_DIR_SRC}/Initializer.cpp
        ${CPP_DIR_SRC}/ParameterParser.cpp
        ${CPP_DIR_SRC_SHARED}/File.cpp

        ${CPP_DIR_SRC}/TimeMoments.cpp

        ${CPP_DIR_SRC}/MolflowTypes.cpp
        ${CPP_DIR_SRC_SHARED}/Buffer_shared.cpp

        #${SIMU_DIR}/MolflowSimulation.cpp
        #${SIMU_DIR}/SimulationMC.cpp
        ${CPP_DIR_SRC_SHARED}/IntersectAABB_shared.cpp
        ${CPP_DIR_SRC_SHARED}/Polygon.cpp
        ${CPP_DIR_SRC_SHARED}/Random.cpp

        ${CPP_DIR_SRC_SHARED}/FlowMPI.h #contains templates
        ${CPP_DIR_SRC_SHARED}/FlowMPI.cpp

        ${IO_DIR}/CSVExporter.cpp
        ${IO_DIR}/CSVExporter.h

        ${HELPER_DIR}/GLProgress_CLI.cpp
        ${HELPER_DIR}/GLProgress_CLI.hpp
        ${HELPER_DIR}/GLProgress_abstract.cpp
        ${HELPER_DIR}/GLProgress_abstract.hpp

        )

add_executable(${PROJECT_NAME} ${SRC_FILES})

IF(APPLE)
include(${CMAKE_HOME_DIRECTORY}/src_shared/SetOpenMP.cmake) #Help macOS find OpenMP
ENDIF()

find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    message(STATUS "Detected OpenMP version: ${OpenMP_CXX_VERSION}")
    target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
endif()

option(USE_MPI FALSE)
if(USE_MPI)
    message(STATUS "MPI use requested (USE_MPI is TRUE), finding...")
    find_package(MPI REQUIRED)
    if(MPI_FOUND)
    message(STATUS "Detected MPI version: ${MPI_CXX_VERSION}")
    set(USE_MPI TRUE)
    target_compile_definitions(${PROJECT_NAME} PUBLIC -DUSE_MPI)
    target_include_directories(${PROJECT_NAME} SYSTEM PRIVATE ${MPI_INCLUDE_PATH})
    target_link_libraries(${PROJECT_NAME} PRIVATE MPI::MPI_CXX)
    endif()
endif()

target_link_libraries(${PROJECT_NAME} PUBLIC molflow_core)

target_include_directories(${PROJECT_NAME} PRIVATE
        ${HEADER_DIR_SRC}
        ${HEADER_DIR_SRC_SHARED}
        ${HEADER_DIR_INCLUDE}
        )

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)


# set the path to the library folder
if(MSVC)
    # Add Whole Program Optimization and Link Time Code Generation
    if(NOT CMAKE_BUILD_TYPE MATCHES Debug)
        set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "/GL") #Incompatible with Edit and continue
        set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS "/LTCG") #no whole program optimization for debug
    endif()
    

    # Add Windows Console output
    message(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
    target_link_directories(${PROJECT_NAME} PRIVATE ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

endif(MSVC)
if(NOT MSVC)

endif(NOT MSVC)