cmake_minimum_required(VERSION 3.12.2 FATAL_ERROR)

#shared GUI and CLI part, built as library that GUI or CLI links against
set(PROJECT_NAME molflow_core)


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

else() #not MSVC
    if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++") #default on macOS: libstdc++
    endif()
endif()


################ Files ################
#   --   Add files to project.   --   #
#######################################

set(SRC_FILES
        ${CPP_DIR_SRC_SHARED}/SimulationManager.cpp
        ${CPP_DIR_SRC_SHARED}/ShMemory.cpp
        ${CPP_DIR_SRC_SHARED}/Process.cpp
        ${CPP_DIR_SRC_SHARED}/Vector.cpp
        ${CPP_DIR_SRC_SHARED}/Distributions.cpp
        ${CPP_DIR_SRC_SHARED}/Buffer_shared.cpp
        ${CPP_DIR_SRC_SHARED}/SimulationController.cpp
        ${CPP_DIR_SRC_SHARED}/SimulationManager.cpp
        ${CPP_DIR_SRC_SHARED}/IntersectAABB_shared.cpp
        ${CPP_DIR_SRC_SHARED}/Polygon.cpp
        ${CPP_DIR_SRC_SHARED}/Random.cpp
        
        ${GLAPP_DIR}/GLFormula.cpp

        ${IO_DIR}/LoaderXML.cpp
        ${IO_DIR}/WriterXML.cpp

        ${CPP_DIR_SRC}/ParameterParser.cpp
        ${CPP_DIR_SRC}/Initializer.cpp
        ${CPP_DIR_SRC}/TimeMoments.cpp
        ${CPP_DIR_SRC}/MolflowTypes.cpp
        
        ${SIMU_DIR}/MolflowSimulation.cpp
        ${SIMU_DIR}/SimulationMC.cpp
        ${SIMU_DIR}/Particle.cpp
        ${SIMU_DIR}/Physics.cpp
        ${SIMU_DIR}/AnglemapGeneration.cpp
        ${SIMU_DIR}/CDFGeneration.cpp
        ${SIMU_DIR}/IDGeneration.cpp
        ${SIMU_DIR}/MolflowSimFacet.cpp
        ${SIMU_DIR}/MolflowSimGeom.cpp
        
        ${HELPER_DIR}/GLProgress_CLI.cpp
        ${HELPER_DIR}/GLProgress_CLI.hpp
        ${HELPER_DIR}/GLProgress_abstract.cpp
        ${HELPER_DIR}/GLProgress_abstract.hpp
)

add_library(${PROJECT_NAME} ${SRC_FILES})

IF(APPLE)
    include(${CPP_DIR_SRC_SHARED}/cmake_modules/SetOpenMP.cmake)  #Help macOS find OpenMP
ENDIF()

find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    message(STATUS "Detected OpenMP version: ${OpenMP_CXX_VERSION}")
    target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
endif()

target_link_libraries(${PROJECT_NAME} PUBLIC src_shared_core)
target_link_libraries(${PROJECT_NAME} PRIVATE truncated_normal)

target_include_directories(${PROJECT_NAME} PRIVATE
        ${HEADER_DIR_SRC}
        ${HEADER_DIR_SRC_SHARED}
        ${HEADER_DIR_INCLUDE}
        )

if(MSVC)
    # Add Whole Program Optimization and Link Time Code Generation
    if(NOT CMAKE_BUILD_TYPE MATCHES Debug)
        set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "/GL") #Incompatible with Edit and continue
        set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS "/LTCG") #No need for whole-program optimization in debug
    endif()    

    # Add Windows Console output
    message(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
    link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
else() #Not MSVC

endif()