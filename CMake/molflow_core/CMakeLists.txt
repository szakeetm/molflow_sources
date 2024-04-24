cmake_minimum_required(VERSION 3.12.2 FATAL_ERROR)

################### Variables. ####################
# Change if you want modify path or other values. #
###################################################

set(PROJECT_NAME molflow_core) #shared GUI and CLI part

# Folders files
set(CPP_DIR_SRC ${CMAKE_HOME_DIRECTORY}/src)
set(SIMU_DIR ${CMAKE_HOME_DIRECTORY}/src/Simulation)
set(IO_DIR ${CMAKE_HOME_DIRECTORY}/src/IO)
set(CPP_DIR_SRC_SHARED ${CMAKE_HOME_DIRECTORY}/src_shared)
set(HELPER_DIR ${CMAKE_HOME_DIRECTORY}/src_shared/Helper)

set(HEADER_DIR_SRC ${CMAKE_HOME_DIRECTORY}/src)
set(HEADER_DIR_SRC_SHARED ${CMAKE_HOME_DIRECTORY}/src_shared)
set(HEADER_DIR_INCLUDE ${CMAKE_HOME_DIRECTORY}/include)


IF (WIN32)
    # set stuff for windows

    set(LINK_DIR_1 ${CMAKE_HOME_DIRECTORY}/lib_external/win/${MY_BUILD_TYPE})
ELSEIF(APPLE)
    set(LINK_DIR_1 ${CMAKE_HOME_DIRECTORY}/lib_external/mac)
ELSE()
    # set stuff for other systems

    # link to fedora libraries if EL Linux (Red Hat Enterprise Linux) has been detected
    #[[ #No more linux dependencies
    IF(os_version_suffix MATCHES "\\.el[1-9]")
        set(LINK_DIR_1 ${CMAKE_HOME_DIRECTORY}/lib_external/linux_fedora)
    ELSE()
        set(LINK_DIR_1 ${CMAKE_HOME_DIRECTORY}/lib_external/linux_debian)
    ENDIF()]]
ENDIF()



############## CMake Project ################
#        The main options of project        #
#############################################

project(${PROJECT_NAME} CXX)

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
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
    if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
    endif()
endif()

################ Files ################
#   --   Add files to project.   --   #
#######################################

set(SRC_FILES
        ${CPP_DIR_SRC_SHARED}/SimulationManager.cpp
        ${CPP_DIR_SRC_SHARED}/ShMemory.cpp
        ${CPP_DIR_SRC_SHARED}/Process.cpp
        ${CPP_DIR_SRC_SHARED}/GLApp/GLFormula.cpp

        ${CPP_DIR_SRC_SHARED}/Vector.cpp
        ${CPP_DIR_SRC_SHARED}/Distributions.cpp
        ${CPP_DIR_SRC_SHARED}/GLApp/GLFormula.h

        ${IO_DIR}/LoaderXML.cpp
        ${IO_DIR}/WriterXML.cpp

        ${CPP_DIR_SRC}/ParameterParser.cpp
        ${CPP_DIR_SRC}/Initializer.cpp

        ${CPP_DIR_SRC}/TimeMoments.cpp

        ${CPP_DIR_SRC}/MolflowTypes.cpp
        ${CPP_DIR_SRC_SHARED}/Buffer_shared.cpp

        ${SIMU_DIR}/MolflowSimulation.cpp
        ${SIMU_DIR}/SimulationMC.cpp
        ${SIMU_DIR}/Particle.cpp
        ${SIMU_DIR}/Physics.cpp
        ${SIMU_DIR}/AnglemapGeneration.cpp
        ${SIMU_DIR}/CDFGeneration.cpp
        ${SIMU_DIR}/IDGeneration.cpp
        ${SIMU_DIR}/MolflowSimFacet.cpp
        ${SIMU_DIR}/MolflowSimGeom.cpp

        ${CPP_DIR_SRC_SHARED}/SimulationController.cpp
        ${CPP_DIR_SRC_SHARED}/SimulationManager.cpp
        ${CPP_DIR_SRC_SHARED}/IntersectAABB_shared.cpp
        ${CPP_DIR_SRC_SHARED}/Polygon.cpp
        ${CPP_DIR_SRC_SHARED}/Random.cpp
        ${CMAKE_HOME_DIRECTORY}/src_shared/GLApp/GLFormula.h
        
        ${HELPER_DIR}/GLProgress_CLI.cpp
        ${HELPER_DIR}/GLProgress_CLI.hpp
        ${HELPER_DIR}/GLProgress_abstract.cpp
        ${HELPER_DIR}/GLProgress_abstract.hpp
)

# set the path to the library folder
link_directories(${LINK_DIR_1})

add_library(${PROJECT_NAME} ${SRC_FILES})

IF(APPLE)
include(${CMAKE_HOME_DIRECTORY}/src_shared/SetOpenMP.cmake)  #Help macOS find OpenMP
ENDIF()

find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    message(STATUS "Detected OpenMP version: ${OpenMP_CXX_VERSION}")
    target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
endif()

target_link_libraries(${PROJECT_NAME} PUBLIC src_shared_simulator)

find_package(pugixml CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE pugixml::pugixml)
target_link_libraries(${PROJECT_NAME} PRIVATE truncated_normal)

target_include_directories(${PROJECT_NAME} PRIVATE
        ${HEADER_DIR_SRC}
        ${HEADER_DIR_SRC_SHARED}
        ${HEADER_DIR_INCLUDE}
        )

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)

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