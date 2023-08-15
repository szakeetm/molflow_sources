cmake_minimum_required(VERSION 3.12.2 FATAL_ERROR)

################### Variables. ####################
# Change if you want modify path or other values. #
###################################################

set(PROJECT_NAME molflow)

# Folders files
set(CPP_DIR_1 ../../src)
set(UI_DIR ../../src/Interface)
set(IO_DIR ../../src/IO)

set(HEADER_DIR_1 ../../src)
set(HEADER_DIR_2 ../../src_shared)
#set(HEADER_DIR_3 ../../src_shared/GLApp)
#set(HEADER_DIR_4 ../../src_shared/GLApp/GLChart)
#set(HEADER_DIR_5 ../../_molflow/molflow_project/lib/include/)
set(HEADER_DIR_6 ../../include)

IF (WIN32)
    # set stuff for windows

    set(HEADER_DIR_7 ../../include/windows_only)
    set(HEADER_DIR_8 ../../include/windows_only/png)

    set(LINK_DIR_1 ../../lib_external/win/${MY_BUILD_TYPE})
    set(DLL_DIR ../../lib_external/win/dll)
ELSEIF(APPLE)
    # set stuff for mac os
    set(HEADER_DIR_7 "../../include/mac_only")
    set(LINK_DIR_1 ../../lib_external/mac)
ELSE()
    # set stuff for other systems

    # link to fedora libraries if EL Linux (Red Hat Enterprise Linux) has been detected
    #[[ #no more linux libraries
    IF(os_version_suffix MATCHES "\\.el[1-9]")
        set(LINK_DIR_1 ../../lib_external/linux_fedora)
    ELSE()
        set(LINK_DIR_1 ../../lib_external/linux_debian)
    ENDIF()]]
ENDIF()

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

        ${CPP_DIR_1}/FormulaEvaluator_MF.cpp
        ${CPP_DIR_1}/MolFlow.cpp
        ${CPP_DIR_1}/MolflowFacet.cpp
        ${CPP_DIR_1}/MolflowGeometry.cpp
        ${CPP_DIR_1}/MolflowTypes.cpp
        ${CPP_DIR_1}/MolflowWorker.cpp
        ${CPP_DIR_1}/TimeMoments.cpp

        ${HEADER_DIR_2}/FlowMPI.h
        ${HEADER_DIR_2}/FlowMPI.cpp
        ${HEADER_DIR_2}/IntersectAABB_shared.cpp

        )

# Add executable to build.
if(WIN32)
    add_executable(${PROJECT_NAME} WIN32 ${SRC_FILES})
else()
    add_executable(${PROJECT_NAME} ${SRC_FILES} )
endif()

include(${CMAKE_HOME_DIRECTORY}/src_shared/SetOpenMP.cmake)
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    message(STATUS "Detected OpenMP version: ${OpenMP_CXX_VERSION}")
    target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
endif()

target_link_libraries(${PROJECT_NAME} PUBLIC shared_interface)
target_link_libraries(${PROJECT_NAME} PUBLIC molflowSim)

target_include_directories(${PROJECT_NAME} PRIVATE
        ${HEADER_DIR_1}
        ${UI_DIR}
        ${HEADER_DIR_2}
        ${HEADER_DIR_3}
        ${HEADER_DIR_4}
        ${HEADER_DIR_5}
        ${HEADER_DIR_6}
        ${HEADER_DIR_7}
        ${HEADER_DIR_8}
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
    find_package(OpenGL REQUIRED)
    # 1. link against external libs
    target_link_libraries(${PROJECT_NAME} PRIVATE
            SDL2.lib
            SDL2main.lib
            opengl32
            user32.lib
            shell32.lib
            ole32.lib)
    # 2. link against our own libs
    #[[target_link_libraries(${PROJECT_NAME} PRIVATE pugixml nativefiledialog)]]
else() #not MSVC
    if(APPLE)
        #link to self-build sdl shared lib
        target_link_libraries(${PROJECT_NAME} PRIVATE "-framework AppKit")
    else()
        # Use the package PkgConfig to detect GTK+ headers/library files
        #[[FIND_PACKAGE(PkgConfig REQUIRED)
        PKG_CHECK_MODULES(GTK3 REQUIRED gtk+-3.0)

        # Setup CMake to use GTK+, tell the compiler where to look for headers
        # and to the linker where to look for libraries
        INCLUDE_DIRECTORIES(${GTK3_INCLUDE_DIRS})
        LINK_DIRECTORIES(${GTK3_LIBRARY_DIRS})

        # Add other flags to the compiler
        ADD_DEFINITIONS(${GTK3_CFLAGS_OTHER})

        target_link_libraries(${PROJECT_NAME} ${GTK3_LIBRARIES})

        find_package(X11 REQUIRED)
        target_include_directories(${PROJECT_NAME} PRIVATE ${X11_INCLUDE_DIRS})]]
    endif()

    find_package(OpenGL REQUIRED)
    target_include_directories(${PROJECT_NAME} SYSTEM PRIVATE ${OPENGL_INCLUDE_DIRS})

    find_package(SDL2 REQUIRED)
    target_include_directories(${PROJECT_NAME} SYSTEM PRIVATE ${SDL2_INCLUDE_DIRS})

    find_package(PNG REQUIRED)
    target_include_directories(${PROJECT_NAME} SYSTEM PRIVATE ${PNG_INCLUDE_DIRS})

    find_package(CURL REQUIRED)
    target_include_directories(${PROJECT_NAME} SYSTEM PRIVATE ${CURL_INCLUDE_DIRS})

    set(THREADS_PREFER_PTHREAD_FLAG ON)
    find_package(Threads REQUIRED)

    target_link_libraries(${PROJECT_NAME} PRIVATE ${OPENGL_LIBRARIES})
    target_link_libraries(${PROJECT_NAME} PRIVATE ${SDL2_LIBRARIES})
    target_link_libraries(${PROJECT_NAME} PRIVATE ${PNG_LIBRARIES})
    #target_link_libraries(${PROJECT_NAME} PRIVATE ${GSL_LIBRARIES})
    target_link_libraries(${PROJECT_NAME} PRIVATE ${CURL_LIBRARIES})
    target_link_libraries(${PROJECT_NAME} PRIVATE Threads::Threads)
    target_link_libraries(${PROJECT_NAME} PRIVATE ${X11_LIBRARIES})

    #for shared memory
    find_library(LIBRT rt)
    if(LIBRT)
        target_link_libraries(${PROJECT_NAME} PRIVATE ${LIBRT})
    endif()

    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        target_link_libraries(${PROJECT_NAME} PRIVATE stdc++)
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
        #don´t add anything for filesystem
    else()
        target_link_libraries(${PROJECT_NAME} PRIVATE stdc++fs)
    endif()

endif()