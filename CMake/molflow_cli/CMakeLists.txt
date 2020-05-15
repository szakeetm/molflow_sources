cmake_minimum_required(VERSION 3.12.2 FATAL_ERROR)

################### Variables. ####################
# Change if you want modify path or other values. #
###################################################

set(PROJECT_NAME molflowCLI)

# Folders files
set(CPP_DIR_1 ../../src)
set(CPP_DIR_2 ../../src_shared)

set(HEADER_DIR_1 ../../src)
set(HEADER_DIR_2 ../../src_shared)
set(HEADER_DIR_3 ../../include)

IF (WIN32)
    # set stuff for windows

    set(LINK_DIR_1 ../../lib/win/${MY_BUILD_TYPE})
    set(LINK_DIR_2 ../../lib_external/win/${MY_BUILD_TYPE})
ELSEIF(APPLE)
    set(LINK_DIR_1 ../../lib_external/mac)
ELSE()
    # set stuff for other systems

    # link to fedora libraries if EL Linux (Red Hat Enterprise Linux) has been detected
    IF(os_version_suffix STREQUAL ".el7")
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

# Definition of Macros
add_definitions(
        -DMOLFLOW
)

################# Flags ################
# Defines Flags for Windows and Linux. #
########################################

if(MSVC)

endif(MSVC)
if(NOT MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
    if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
    endif()
endif(NOT MSVC)

################ Files ################
#   --   Add files to project.   --   #
#######################################

set(SRC_FILES
        ${CPP_DIR_2}/FlowCLI.cpp
        ${CPP_DIR_2}/SimulationManager.cpp
        ${CPP_DIR_2}/ShMemory.cpp
        ${CPP_DIR_2}/Process.cpp

        ${CPP_DIR_1}/GeometrySimu.cpp
        ${CPP_DIR_1}/MolflowBuffer.cpp
        ${CPP_DIR_1}/Parameter.cpp
        ${CPP_DIR_2}/Vector.cpp
        )

# set the path to the library folder
#[[link_directories(${LINK_DIR_1}
        ${LINK_DIR_2}
        )]]

add_executable(${PROJECT_NAME} ${SRC_FILES})
#target_link_libraries(${PROJECT_NAME} simulator_shared)

target_include_directories(${PROJECT_NAME} PRIVATE
        ${HEADER_DIR_1}
        ${HEADER_DIR_2}
        ${HEADER_DIR_3}
        ${HEADER_DIR_4}
        ${HEADER_DIR_5}
        ${HEADER_DIR_6}
        ${HEADER_DIR_7}
        )

if(MSVC)
        set(HEADER_DIR_6 ../include/windows_only)
        target_include_directories(${PROJECT_NAME} PUBLIC
                ${HEADER_DIR_6}
                )
    else()
        find_package(GSL REQUIRED)
        target_include_directories(${PROJECT_NAME} PRIVATE ${GSL_INCLUDE_DIRS})
endif()

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)
target_compile_options(${PROJECT_NAME} PRIVATE
        $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:
        -Wall>
        $<$<CXX_COMPILER_ID:MSVC>:
        /W4>)

# Multi-processor compilation
if (MSVC)
    target_compile_options(${PROJECT_NAME} PRIVATE
            "$<$<CONFIG:Debug>:/MP>"
            "$<$<CONFIG:Release>:/MP>"
            )
endif ()

if(MSVC)
    # Add Whole Program Optimization and Link Time Code Generation
    set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "/GL")
    set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS "/LTCG")

    # Add Windows Console output
    message(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
    link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
endif(MSVC)
if(NOT MSVC)

endif(NOT MSVC)