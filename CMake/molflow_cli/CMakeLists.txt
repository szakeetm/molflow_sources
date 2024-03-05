cmake_minimum_required(VERSION 3.12.2 FATAL_ERROR)

################### Variables. ####################
# Change if you want modify path or other values. #
###################################################

set(PROJECT_NAME molflowCLI)

# Folders files
set(CPP_DIR_1 ../../src)
set(SIMU_DIR ../../src/Simulation)
set(IO_DIR ../../src/IO)
set(CPP_DIR_2 ../../src_shared)
set(CPP_DIR_3 ../../src_shared/GLApp)
set(HELPER_DIR ../../src_shared/Helper)

set(HEADER_DIR_1 ../../src)
set(HEADER_DIR_2 ../../src_shared)
set(HEADER_DIR_3 ../../include)

set(HEADER_DIR_EXTERNAL ../../src_shared/external)

IF (WIN32)
    # set stuff for windows
    set(HEADER_DIR_4 ../../include/windows_only)

    set(LINK_DIR_1 ../../lib/win/${MY_BUILD_TYPE})
    set(LINK_DIR_2 ../../lib_external/win/${MY_BUILD_TYPE})
ELSEIF(APPLE)
    set(LINK_DIR_1 ../../lib_external/mac)
ELSE()
    # set stuff for other systems

    # link to fedora libraries if EL Linux (Red Hat Enterprise Linux) has been detected
    IF(os_version_suffix MATCHES "\\.el[1-9]")
        set(LINK_DIR_1 ../../lib_external/linux_fedora)
    ELSE()
        set(LINK_DIR_1 ../../lib_external/linux_debian)
    ENDIF()
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

else() #not msvc
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
    if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
    endif()
endif()

################ Files ################
#   --   Add files to project.   --   #
#######################################

set(SRC_FILES
        ${CPP_DIR_1}/MolflowCLI.cpp
        ${CPP_DIR_2}/SimulationManager.cpp
        ${CPP_DIR_2}/ShMemory.cpp
        ${CPP_DIR_2}/Process.cpp

        ${CPP_DIR_2}/Vector.cpp
        ${CPP_DIR_2}/Distributions.cpp

        ${IO_DIR}/LoaderXML.cpp
        ${IO_DIR}/WriterXML.cpp
        ${CPP_DIR_1}/Initializer.cpp
        ${CPP_DIR_1}/ParameterParser.cpp
        ${CPP_DIR_2}/File.cpp

        ${CPP_DIR_1}/TimeMoments.cpp

        ${CPP_DIR_1}/MolflowTypes.cpp
        ${CPP_DIR_2}/Buffer_shared.cpp

        #${SIMU_DIR}/MolflowSimulation.cpp
        #${SIMU_DIR}/SimulationMC.cpp
        ${CPP_DIR_2}/IntersectAABB_shared.cpp
        ${CPP_DIR_2}/Polygon.cpp
        ${CPP_DIR_2}/Random.cpp

        ${CPP_DIR_2}/FlowMPI.h #contains templates
        ${CPP_DIR_2}/FlowMPI.cpp

        ${IO_DIR}/CSVExporter.cpp
        ${IO_DIR}/CSVExporter.h

        ${HELPER_DIR}/GLProgress_CLI.cpp
        ${HELPER_DIR}/GLProgress_CLI.hpp
        ${HELPER_DIR}/GLProgress_abstract.cpp
        ${HELPER_DIR}/GLProgress_abstract.hpp

        )

add_executable(${PROJECT_NAME} ${SRC_FILES})

include(${CMAKE_HOME_DIRECTORY}/src_shared/SetOpenMP.cmake)
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

target_link_libraries(${PROJECT_NAME} PUBLIC molflowSim)
find_package(pugixml CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE pugixml::pugixml)

target_include_directories(${PROJECT_NAME} PRIVATE
        ${HEADER_DIR_1}
        ${HEADER_DIR_2}
        ${HEADER_DIR_3}
        ${HEADER_DIR_4}
        ${HEADER_DIR_5}
        ${HEADER_DIR_6}
        ${HEADER_DIR_7}   
        ${HEADER_DIR_EXTERNAL}
        )



if(MSVC)
        set(HEADER_DIR_6 ../include/windows_only)
        target_include_directories(${PROJECT_NAME} PUBLIC
                ${HEADER_DIR_6}
                )
endif()

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)


# set the path to the library folder
target_link_directories(${PROJECT_NAME} PRIVATE ${LINK_DIR_1} ${LINK_DIR_2})
if(MSVC)
    # Add Whole Program Optimization and Link Time Code Generation
    if(NOT CMAKE_BUILD_TYPE MATCHES Debug)
        set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "/GL") #Incompatible with Edit and continue
        set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS "/LTCG") #no whole program optimization for debug
    endif()
    

    # Add Windows Console output
    message(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
    target_link_directories(${PROJECT_NAME} PRIVATE ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

    # zip extraction
    target_link_libraries(${PROJECT_NAME} PRIVATE
            #bzip2.lib
            #lzma.lib
            #zlib.lib
            )
endif(MSVC)
if(NOT MSVC)

endif(NOT MSVC)