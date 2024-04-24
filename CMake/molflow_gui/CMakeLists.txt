cmake_minimum_required(VERSION 3.12.2 FATAL_ERROR)

################### Variables. ####################
# Change if you want modify path or other values. #
###################################################

set(PROJECT_NAME molflow)

# GUI-specific folders files
set(UI_DIR ${CMAKE_HOME_DIRECTORY}/src/Interface)

############## CMake Project ################
#        The main options of project        #
#############################################

project(${PROJECT_NAME} CXX)

# Definition of Macros
add_definitions(
        -DCURL_STATICLIB
        -D_MBCS
        -D_CRT_SECURE_NO_WARNINGS
        -D_CRT_NONSTDC_NO_DEPRECATE
)
IF (WIN32)
    add_definitions(
            #-D_WINDOWS
            -D_WITHGL
    )
ENDIF()

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
        ${UI_DIR}/FacetAdvParams.cpp
        ${UI_DIR}/FacetDetails.cpp
        ${UI_DIR}/GeometryRender.cpp
        ${UI_DIR}/GeometryViewer.cpp
        ${UI_DIR}/GlobalSettings.cpp
        ${UI_DIR}/ImportDesorption.cpp
        ${UI_DIR}/MomentsEditor.cpp
        ${UI_DIR}/Movement.cpp
        ${UI_DIR}/BackgroundGas.cpp
		${UI_DIR}/MeasureForce.cpp
        ${UI_DIR}/OutgassingMapWindow.cpp
        ${UI_DIR}/ParameterEditor.cpp
        ${UI_DIR}/PressureEvolution.cpp
        ${UI_DIR}/ProfilePlotter.cpp
        ${UI_DIR}/TexturePlotter.cpp
        ${UI_DIR}/TextureScaling.cpp
        ${UI_DIR}/TimeSettings.cpp
        ${UI_DIR}/TimewisePlotter.cpp
        ${UI_DIR}/Viewer3DSettings.cpp

        ${CPP_DIR_SRC}/FormulaEvaluator_MF.cpp
        ${CPP_DIR_SRC}/MolFlow.cpp
        ${CPP_DIR_SRC}/MolflowFacet.cpp
        ${CPP_DIR_SRC}/MolflowGeometry.cpp
        ${CPP_DIR_SRC}/MolflowTypes.cpp
        ${CPP_DIR_SRC}/MolflowWorker.cpp
        ${CPP_DIR_SRC}/TimeMoments.cpp

        ${HEADER_DIR_SRC_SHARED}/FlowMPI.h
        ${HEADER_DIR_SRC_SHARED}/FlowMPI.cpp
        ${HEADER_DIR_SRC_SHARED}/IntersectAABB_shared.cpp

        )

# Add executable to build.
if(WIN32)
    add_executable(${PROJECT_NAME} WIN32 ${SRC_FILES})
else()
    add_executable(${PROJECT_NAME} ${SRC_FILES} )
endif()

if(APPLE)
include(${CMAKE_HOME_DIRECTORY}/src_shared/SetOpenMP.cmake)  #Help macOS find OpenMP
endif()

find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    message(STATUS "Detected OpenMP version: ${OpenMP_CXX_VERSION}")
    target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
endif()

target_link_libraries(${PROJECT_NAME} PUBLIC src_shared_interface)
target_link_libraries(${PROJECT_NAME} PUBLIC molflow_core)

target_include_directories(${PROJECT_NAME} PRIVATE
        ${HEADER_DIR_SRC}
        ${UI_DIR}
        ${HEADER_DIR_SRC_SHARED}
        )

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)

if(MSVC)
    # Add Whole Program Optimization and Link Time Code Generation
    if(NOT CMAKE_BUILD_TYPE MATCHES Debug)
        set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "/GL") #Incompatible with Edit and continue
        set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS "/LTCG")
    endif()
    set_target_properties(${PROJECT_NAME} PROPERTIES
            WIN32_EXECUTABLE YES
            LINK_FLAGS "/ENTRY:mainCRTStartup"
            )

    # Add Windows Console output
    #if(CMAKE_BUILD_TYPE STREQUAL "Debug") #Always show console for CLI builds as error reporting not complete
        set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS "/SUBSYSTEM:CONSOLE ${LINK_FLAGS}")
    #endif()

    # set the path to the library folder
    target_link_directories(${PROJECT_NAME} PRIVATE ${LINK_DIR_1})

    target_link_directories(${PROJECT_NAME} PRIVATE ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
    
    # 1. link against external libs
    target_link_libraries(${PROJECT_NAME} PRIVATE
            opengl32
            user32.lib
            shell32.lib
            ole32.lib)
    # 2. link against our own libs
    
else() #not MSVC
    if(APPLE)
        #link to self-build sdl shared lib
        target_link_libraries(${PROJECT_NAME} PRIVATE "-framework AppKit")
    else()
        # Use the package PkgConfig to detect GTK+ headers/library files
        PKG_CHECK_MODULES(GTK3 REQUIRED gtk+-3.0)

        # Setup CMake to use GTK+, tell the compiler where to look for headers
        # and to the linker where to look for libraries
        INCLUDE_DIRECTORIES(${GTK3_INCLUDE_DIRS})
        LINK_DIRECTORIES(${GTK3_LIBRARY_DIRS})

        # Add other flags to the compiler
        ADD_DEFINITIONS(${GTK3_CFLAGS_OTHER})

        target_link_libraries(${PROJECT_NAME} ${GTK3_LIBRARIES})

    endif()

    #for shared memory
    find_library(LIBRT rt)
    if(LIBRT)
        target_link_libraries(${PROJECT_NAME} PRIVATE ${LIBRT})
    endif()

endif()