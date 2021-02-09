cmake_minimum_required(VERSION 3.12.2 FATAL_ERROR)

################### Variables. ####################
# Change if you want modify path or other values. #
###################################################

set(PROJECT_NAME molflow)

# Folders files
set(CPP_DIR_1 ../../src)
set(UI_DIR ../../src/Interface)
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
        -DCURL_STATICLIB
        -D_MBCS
        -D_CRT_SECURE_NO_WARNINGS
        -D_CRT_NONSTDC_NO_DEPRECATE
)
IF (WIN32)
    add_definitions(
            -D_WINDOWS
            -D_WITHGL
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

if (APPLE)
    set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include" CACHE INTERNAL "OpenMP flags for Xcode toolchain.")
    set(OpenMP_C_LIB_NAMES "omp" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")

    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include" CACHE INTERNAL "OpenMP flags for Xcode toolchain.")
    set(OpenMP_CXX_LIB_NAMES "omp" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")
    set(OpenMP_omp_LIBRARY "/usr/local/opt/libomp/lib/libomp.dylib" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")
endif()

set(SRC_FILES
        ${UI_DIR}/FacetAdvParams.cpp
        ${UI_DIR}/FacetDetails.cpp
        ${CPP_DIR_1}/FormulaEvaluator_MF.cpp
        ${UI_DIR}/GeometryRender.cpp
        ${UI_DIR}/GeometryViewer.cpp
        ${UI_DIR}/GlobalSettings.cpp
        ${UI_DIR}/ImportDesorption.cpp
        ${CPP_DIR_1}/MolFlow.cpp
        ${CPP_DIR_1}/MolflowBuffer.cpp
        ${CPP_DIR_1}/GeometrySimu.cpp
        ${CPP_DIR_1}/MolflowFacet.cpp
        ${CPP_DIR_1}/MolflowGeometry.cpp
        ${CPP_DIR_1}/MolflowTypes.cpp
        ${CPP_DIR_1}/MolflowWorker.cpp
        ${UI_DIR}/MomentsEditor.cpp
        ${UI_DIR}/Movement.cpp
        ${UI_DIR}/OutgassingMapWindow.cpp
        ${CPP_DIR_1}/Parameter.cpp
        ${UI_DIR}/ParameterEditor.cpp
        ${UI_DIR}/PressureEvolution.cpp
        ${UI_DIR}/ProfilePlotter.cpp
        ${UI_DIR}/TexturePlotter.cpp
        ${UI_DIR}/TextureScaling.cpp
        ${CPP_DIR_1}/TimeMoments.cpp
        ${UI_DIR}/TimeSettings.cpp
        ${UI_DIR}/TimewisePlotter.cpp
        ${UI_DIR}/Viewer3DSettings.cpp
        ${HEADER_DIR_2}/SimulationController.cpp
        ${HEADER_DIR_2}/SimulationManager.cpp
        ${HEADER_DIR_2}/IntersectAABB_shared.cpp


        )

# Add executable to build.
if(WIN32)
    add_executable(${PROJECT_NAME} WIN32 ${SRC_FILES})
else()
    add_executable(${PROJECT_NAME} ${SRC_FILES})
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

find_package(OpenMP REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)

if(MSVC)
    # Add Whole Program Optimization and Link Time Code Generation
    set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "/GL")
    #set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS "/LTCG")
    set_target_properties(${PROJECT_NAME} PROPERTIES
            WIN32_EXECUTABLE YES
            LINK_FLAGS "/ENTRY:mainCRTStartup"
            )

    # Add Windows Console output
    #if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS "/SUBSYSTEM:CONSOLE ${LINK_FLAGS}")
    #endif()

    # set the path to the library folder
    target_link_directories(${PROJECT_NAME} PRIVATE ${LINK_DIR_1})

    message(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
    target_link_directories(${PROJECT_NAME} PRIVATE ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
    find_package(OpenGL REQUIRED)
    # 1. link against external libs
    target_link_libraries(${PROJECT_NAME} PRIVATE
            bzip2.lib
            #libcurl_a_x64.lib
            SDL2.lib
            SDL2main.lib
            #libgsl.lib
            #libgslcblas.lib
            lzma.lib
            ZipLib.lib
            #zlib.lib
            opengl32#.lib
            user32.lib
            shell32.lib
            ole32.lib)
    # 2. link against our own libs
    target_link_libraries(${PROJECT_NAME} PRIVATE pugixml nativefiledialog)
elseif(NOT MSVC)
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
    target_include_directories(${PROJECT_NAME} PRIVATE ${OPENGL_INCLUDE_DIRS})

    find_package(SDL2 REQUIRED)
    target_include_directories(${PROJECT_NAME} PRIVATE ${SDL2_INCLUDE_DIRS})

    find_package(PNG REQUIRED)
    target_include_directories(${PROJECT_NAME} PRIVATE ${PNG_INCLUDE_DIRS})

    find_package(GSL REQUIRED)
    target_include_directories(${PROJECT_NAME} PRIVATE ${GSL_INCLUDE_DIRS})

    find_package(CURL REQUIRED)
    target_include_directories(${PROJECT_NAME} PRIVATE ${CURL_INCLUDE_DIRS})

    set(THREADS_PREFER_PTHREAD_FLAG ON)
    find_package(Threads REQUIRED)

    target_link_libraries(${PROJECT_NAME} PRIVATE ${OPENGL_LIBRARIES})
    target_link_libraries(${PROJECT_NAME} PRIVATE ${SDL2_LIBRARIES})
    target_link_libraries(${PROJECT_NAME} PRIVATE ${PNG_LIBRARIES})
    target_link_libraries(${PROJECT_NAME} PRIVATE ${GSL_LIBRARIES})
    target_link_libraries(${PROJECT_NAME} PRIVATE ${CURL_LIBRARIES})
    target_link_libraries(${PROJECT_NAME} PRIVATE Threads::Threads)
    target_link_libraries(${PROJECT_NAME} PRIVATE ${X11_LIBRARIES})

    #for shared memory
    find_library(LIBRT rt)
    if(LIBRT)
        target_link_libraries(${PROJECT_NAME} PRIVATE ${LIBRT})
    endif()

    # Your-external "mylib", add GLOBAL if the imported library is located in directories above the current.
    if (NOT TARGET libzip)
        add_library( libzip STATIC IMPORTED GLOBAL)
    endif()
    # You can define two import-locations: one for debug and one for release.
    get_filename_component(ABS_LINK_DIR_1 "${LINK_DIR_1}" REALPATH)

    # other static libraries need to be built with similar settings for clang (here: libc++)
    if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
        set_target_properties( libzip PROPERTIES IMPORTED_LOCATION ${ABS_LINK_DIR_1}/libzip_clang.a )
    else()
        set_target_properties( libzip PROPERTIES IMPORTED_LOCATION ${ABS_LINK_DIR_1}/libzip_gcc.a )
    endif()
    target_link_libraries(${PROJECT_NAME} PRIVATE libzip) # from ./lib/

    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        target_link_libraries(${PROJECT_NAME} PRIVATE stdc++)
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
        #don´t add anything for filesystem
    else()
        target_link_libraries(${PROJECT_NAME} PRIVATE stdc++fs)
    endif()

    target_link_libraries(${PROJECT_NAME} PRIVATE pugixml sdl_savepng truncatedgaussian)
    target_link_libraries(${PROJECT_NAME} PRIVATE nativefiledialog)

endif()