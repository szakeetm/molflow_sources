cmake_minimum_required(VERSION 3.12.2 FATAL_ERROR)

################### Variables. ####################
# Change if you want modify path or other values. #
###################################################

set(PROJECT_NAME molflowSub)

# Folders files
set(CPP_DIR_1 ../../src)
set(CPP_DIR_2 ../../src_shared)
set(CPP_DIR_3 ../../src_shared/GLApp)

set(HEADER_DIR_1 ../../src)
set(HEADER_DIR_2 ../../src_shared)
set(HEADER_DIR_3 ../../src_shared/GLApp)

IF (WIN32)
    # set stuff for windows
    set(HEADER_DIR_6 ../../include)

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
IF (WIN32)
    add_definitions(
            -DWIN
    )
ENDIF()

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

set(SRC_FILES ${CPP_DIR_1}/molflowSub.cpp
        ${CPP_DIR_1}/MolflowTypes.cpp
        ${CPP_DIR_1}/Parameter.cpp
        ${CPP_DIR_1}/Simulation.cpp
        #${CPP_DIR_1}/SimulationAC.cpp
        #${CPP_DIR_1}/SimulationControl.cpp
        ${CPP_DIR_1}/SimulationMC.cpp
        ${CPP_DIR_1}/GeometrySimu.cpp

        ${CPP_DIR_2}/SimulationController.cpp
        ${CPP_DIR_3}/MathTools.cpp
        )

set(HEADER_FILES ${HEADER_DIR_1}/MolflowTypes.h
        ${HEADER_DIR_1}/Parameter.h
        ${HEADER_DIR_1}/Simulation.h
        ${HEADER_DIR_2}/SimulationController.h
        )

# set the path to the library folder
link_directories(${LINK_DIR_1}
        ${LINK_DIR_2}
        )

add_executable(${PROJECT_NAME} ${SRC_FILES} ${HEADER_FILES})
target_link_libraries(${PROJECT_NAME} simulator_shared)

target_include_directories(${PROJECT_NAME} PRIVATE
        ${HEADER_DIR_1}
        ${HEADER_DIR_2}
        ${HEADER_DIR_3}
        ${HEADER_DIR_4}
        ${HEADER_DIR_5}
        ${HEADER_DIR_6}
        ${HEADER_DIR_7}
        )

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
    target_link_libraries(${PROJECT_NAME}
            #libgsl.lib
            #libgslcblas.lib
            kernel32.lib
            user32.lib
            gdi32.lib
            winspool.lib
            comdlg32.lib
            advapi32.lib
            shell32.lib
            ole32.lib
            oleaut32.lib
            uuid.lib
            odbc32.lib
            odbccp32.lib)
    target_link_libraries(${PROJECT_NAME} truncatedgaussian)
endif(MSVC)
if(NOT MSVC)

    if(APPLE)
        target_link_libraries(${PROJECT_NAME} "-framework AppKit")
    else()
        #find_package(X11 REQUIRED)
        #target_include_directories(${PROJECT_NAME} PRIVATE ${X11_INCLUDE_DIRS})
    endif()

    find_package(GSL REQUIRED)
    target_include_directories(${PROJECT_NAME} PRIVATE ${GSL_INCLUDE_DIRS})

    find_package(CURL REQUIRED)
    target_include_directories(${PROJECT_NAME} PRIVATE ${CURL_INCLUDE_DIRS})

    set(THREADS_PREFER_PTHREAD_FLAG ON)
    find_package(Threads REQUIRED)

    target_link_libraries(${PROJECT_NAME} ${GSL_LIBRARIES})
    target_link_libraries(${PROJECT_NAME} ${CURL_LIBRARIES})
    target_link_libraries(${PROJECT_NAME} Threads::Threads)
    target_link_libraries(${PROJECT_NAME} ${X11_LIBRARIES})

    #for shared memory
    find_library(LIBRT rt)
    if(LIBRT)
        target_link_libraries(${PROJECT_NAME} ${LIBRT})
    endif()

    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        target_link_libraries(${PROJECT_NAME} c++fs)
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
        #don´t add anything for filesystem
    else()
        target_link_libraries(${PROJECT_NAME} stdc++fs)
    endif()

endif(NOT MSVC)