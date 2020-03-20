cmake_minimum_required(VERSION 3.12.2 FATAL_ERROR)

################### Variables. ####################
# Change if you want modify path or other values. #
###################################################

set(PROJECT_NAME molflow)

IF (WIN32)
    set(OS_NAME "win")
    set(OS_RELPATH "..")
ELSEIF(APPLE)
    set(OS_NAME "mac")
    set(OS_RELPATH "")
ELSE()
    IF(os_version_suffix STREQUAL ".el7")
        set(OS_NAME "linux_fedora")
    ELSE()
        set(OS_NAME "linux_debian")
    ENDIF()
    set(OS_RELPATH "..")
ENDIF()

IF (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(MY_BUILD_TYPE "debug")
ELSE()
    set(MY_BUILD_TYPE "release")
ENDIF()

# Output Variables
set(OUTPUT_DEBUG ${OS_RELPATH}/bin/${OS_NAME}/debug/)
set(OUTPUT_REL ${OS_RELPATH}/bin/${OS_NAME}/release/)

# Folders files
set(CPP_DIR_1 ../../src)
set(CPP_DIR_2 ../../src_shared)
set(CPP_DIR_3 ../../src_shared/GLApp)
set(CPP_DIR_4 ../../src_shared/GLApp/GLChart)
set(HEADER_DIR_1 ../../src)
set(HEADER_DIR_2 ../../src_shared)
set(HEADER_DIR_3 ../../src_shared/GLApp)
set(HEADER_DIR_4 ../../src_shared/GLApp/GLChart)
#set(HEADER_DIR_5 ../../_molflow/molflow_project/lib/include/)
set(HEADER_DIR_6 ../../include)

set(COPY_DIR ../../copy_to_build/)

IF (WIN32)
    # set stuff for windows

    set(CPP_DIR_5 ../../source/shared_code/windows_only)
    set(HEADER_DIR_7 ../../include/windows_only)
    set(HEADER_DIR_8 ../../include/windows_only/png)

    set(LINK_DIR_1 ../../lib/win/${MY_BUILD_TYPE})
    #set(LINK_DIR_2 ../../_molflow/molflow_project/lib/lib)
    set(LINK_DIR_3 ../../lib_external/win/${MY_BUILD_TYPE})
    set(DLL_DIR ../../lib_external/win/dll)
ELSEIF(APPLE)
    # set stuff for mac os
    set(HEADER_DIR_7 "../../include/mac_only")
    set(LINK_DIR_1 ../../lib_external/mac)
    #set(LINK_DIR_2 ../../lib_external/mac/SDL2)
ELSE()
    # set stuff for other systems

    set(CPP_DIR_5 "")
    set(HEADER_DIR_7 "")
    set(HEADER_DIR_8 "")

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
        -DMOLFLOW
        -D_CRT_SECURE_NO_WARNINGS
        -D_CRT_NONSTDC_NO_DEPRECATE
)
IF (WIN32)
    add_definitions(
            -D_WINDOWS
            -D_WITHGL
    )
ENDIF()
############## Artefacts Output #################
# Defines outputs , depending Debug or Release. #
#################################################

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY    "${CMAKE_BINARY_DIR}/${OUTPUT_DEBUG}")
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY    "${CMAKE_BINARY_DIR}/${OUTPUT_DEBUG}")
    set(CMAKE_EXECUTABLE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${OUTPUT_DEBUG}")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY    "${CMAKE_BINARY_DIR}/${OUTPUT_DEBUG}")
else()
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY    "${CMAKE_BINARY_DIR}/${OUTPUT_REL}")
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY    "${CMAKE_BINARY_DIR}/${OUTPUT_REL}")
    set(CMAKE_EXECUTABLE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${OUTPUT_REL}")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY    "${CMAKE_BINARY_DIR}/${OUTPUT_REL}")
endif()
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_EXECUTABLE_OUTPUT_DIRECTORY}) #to build executable in main folder
#message(${EXECUTABLE_OUTPUT_PATH})

# Messages
message("${PROJECT_NAME}: MAIN PROJECT: ${CMAKE_PROJECT_NAME}")
message("${PROJECT_NAME}: CURR PROJECT: ${CMAKE_CURRENT_SOURCE_DIR}")
message("${PROJECT_NAME}: CURR BIN DIR: ${CMAKE_CURRENT_BINARY_DIR}")

################# Flags ################
# Defines Flags for Windows and Linux. #
########################################

#[[if(MSVC)
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++17")
    set(CMAKE_GENERATOR_PLATFORM x64)
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /mwindows")
    #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mwindows")
    # main vs WinMain
    # https://stackoverflow.com/questions/2752792/whats-the-equivalent-of-gccs-mwindows-option-in-cmake
    #set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /subsystem:windows /ENTRY:mainCRTStartup")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /W3 /EHsc")
   set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /W3 /EHsc")
endif(MSVC)
if(NOT MSVC)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
   if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
       set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
   endif()
endif(NOT MSVC)]]

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

file(GLOB GLOB_FILES
        ${CPP_DIR_3}/*.cpp
        ${CPP_DIR_4}/*.cpp
        ${CPP_DIR_5}/*.cpp
        ${HEADER_DIR_3}/*.h
        ${HEADER_DIR_4}/*.h
        )

set(SRC_FILES ${CPP_DIR_1}/FacetAdvParams.cpp
        ${CPP_DIR_1}/FacetDetails.cpp
        ${CPP_DIR_1}/GeometryRender.cpp
        ${CPP_DIR_1}/GeometryViewer.cpp
        ${CPP_DIR_1}/GlobalSettings.cpp
        ${CPP_DIR_1}/ImportDesorption.cpp
        ${CPP_DIR_1}/MolFlow.cpp
        ${CPP_DIR_1}/MolflowFacet.cpp
        ${CPP_DIR_1}/MolflowGeometry.cpp
        ${CPP_DIR_1}/MolflowTypes.cpp
        ${CPP_DIR_1}/MolflowWorker.cpp
        ${CPP_DIR_1}/MomentsEditor.cpp
        ${CPP_DIR_1}/Movement.cpp
        ${CPP_DIR_1}/OutgassingMap.cpp
        ${CPP_DIR_1}/Parameter.cpp
        ${CPP_DIR_1}/ParameterEditor.cpp
        ${CPP_DIR_1}/PressureEvolution.cpp
        ${CPP_DIR_1}/ProfilePlotter.cpp
        ${CPP_DIR_1}/TexturePlotter.cpp
        ${CPP_DIR_1}/TextureScaling.cpp
        ${CPP_DIR_1}/TimeSettings.cpp
        ${CPP_DIR_1}/TimewisePlotter.cpp
        ${CPP_DIR_1}/Viewer3DSettings.cpp

        ${CPP_DIR_2}/ASELoader.cpp
        ${CPP_DIR_2}/AddVertex.cpp
        ${CPP_DIR_2}/AlignFacet.cpp
        ${CPP_DIR_2}/AppUpdater.cpp
        ${CPP_DIR_2}/BuildIntersection.cpp
        ${CPP_DIR_2}/CollapseSettings.cpp
        ${CPP_DIR_2}/CreateShape.cpp
        ${CPP_DIR_2}/Distributions.cpp
        ${CPP_DIR_2}/ExtrudeFacet.cpp
        ${CPP_DIR_2}/FacetCoordinates.cpp
        ${CPP_DIR_2}/Facet_shared.cpp
        ${CPP_DIR_2}/File.cpp
        ${CPP_DIR_2}/FormulaEditor.cpp
        ${CPP_DIR_2}/GeometryRender_shared.cpp
        ${CPP_DIR_2}/GeometryViewer_shared.cpp
        ${CPP_DIR_2}/Geometry_shared.cpp
        ${CPP_DIR_2}/GrahamScan.cpp
        ${CPP_DIR_2}/HistogramPlotter.cpp
        ${CPP_DIR_2}/HistogramSettings.cpp
        ${CPP_DIR_2}/Interface.cpp
        ${CPP_DIR_2}/LoadStatus.cpp
        ${CPP_DIR_2}/MirrorFacet.cpp
        ${CPP_DIR_2}/MirrorVertex.cpp
        ${CPP_DIR_2}/MoveFacet.cpp
        ${CPP_DIR_2}/MoveVertex.cpp
        ${CPP_DIR_2}/ParticleLogger.cpp
        ${CPP_DIR_2}/Polygon.cpp
        ${CPP_DIR_2}/Process.cpp
        ${CPP_DIR_2}/Random.cpp
        ${CPP_DIR_2}/RecoveryDialog.cpp
        ${CPP_DIR_2}/RotateFacet.cpp
        ${CPP_DIR_2}/RotateVertex.cpp
        ${CPP_DIR_2}/ScaleFacet.cpp
        ${CPP_DIR_2}/ScaleVertex.cpp
        ${CPP_DIR_2}/SelectDialog.cpp
        ${CPP_DIR_2}/SelectTextureType.cpp
        ${CPP_DIR_2}/ShMemory.cpp
        ${CPP_DIR_2}/SmartSelection.cpp
        ${CPP_DIR_2}/SplitFacet.cpp
        ${CPP_DIR_2}/Vector.cpp
        ${CPP_DIR_2}/VertexCoordinates.cpp
        ${CPP_DIR_2}/Web.cpp
        ${CPP_DIR_2}/Worker_shared.cpp
        ${CPP_DIR_2}/GeometryConverter.cpp
        )

set(HEADER_FILES ${HEADER_DIR_1}/FacetAdvParams.h
        ${HEADER_DIR_1}/FacetDetails.h
        ${HEADER_DIR_1}/GlobalSettings.h
        ${HEADER_DIR_1}/ImportDesorption.h
        ${HEADER_DIR_1}/MolFlow.h
        ${HEADER_DIR_1}/MolflowFacet.h
        ${HEADER_DIR_1}/MolflowGeometry.h
        ${HEADER_DIR_1}/MolflowTypes.h
        ${HEADER_DIR_1}/MomentsEditor.h
        ${HEADER_DIR_1}/Movement.h
        ${HEADER_DIR_1}/OutgassingMap.h
        ${HEADER_DIR_1}/Parameter.h
        ${HEADER_DIR_1}/ParameterEditor.h
        ${HEADER_DIR_1}/PressureEvolution.h
        ${HEADER_DIR_1}/ProfilePlotter.h
        ${HEADER_DIR_1}/TexturePlotter.h
        ${HEADER_DIR_1}/TextureScaling.h
        ${HEADER_DIR_1}/TimeSettings.h
        ${HEADER_DIR_1}/TimewisePlotter.h
        ${HEADER_DIR_1}/Viewer3DSettings.h

        ${HEADER_DIR_2}/ASELoader.h
        ${HEADER_DIR_2}/AddVertex.h
        ${HEADER_DIR_2}/AlignFacet.h
        ${HEADER_DIR_2}/AppUpdater.h
        ${HEADER_DIR_2}/Buffer_shared.h
        ${HEADER_DIR_2}/BuildIntersection.h
        ${HEADER_DIR_2}/CollapseSettings.h
        ${HEADER_DIR_2}/CreateShape.h
        ${HEADER_DIR_2}/Distributions.h
        ${HEADER_DIR_2}/ExtrudeFacet.h
        ${HEADER_DIR_2}/FacetCoordinates.h
        ${HEADER_DIR_2}/Facet_shared.h
        ${HEADER_DIR_2}/File.h
        ${HEADER_DIR_2}/FormulaEditor.h
        ${HEADER_DIR_2}/GeometryViewer.h
        ${HEADER_DIR_2}/Geometry_shared.h
        ${HEADER_DIR_2}/GrahamScan.h
        ${HEADER_DIR_2}/HistogramPlotter.h
        ${HEADER_DIR_2}/HistogramSettings.h
        ${HEADER_DIR_2}/Interface.h
        ${HEADER_DIR_2}/LoadStatus.h
        ${HEADER_DIR_2}/MirrorFacet.h
        ${HEADER_DIR_2}/MirrorVertex.h
        ${HEADER_DIR_2}/MoveFacet.h
        ${HEADER_DIR_2}/MoveVertex.h
        ${HEADER_DIR_2}/ParticleLogger.h
        ${HEADER_DIR_2}/Polygon.h
        ${HEADER_DIR_2}/RecoveryDialog.h
        ${HEADER_DIR_2}/RotateFacet.h
        ${HEADER_DIR_2}/RotateVertex.h
        ${HEADER_DIR_2}/ScaleFacet.h
        ${HEADER_DIR_2}/ScaleVertex.h
        ${HEADER_DIR_2}/SelectDialog.h
        ${HEADER_DIR_2}/SelectTextureType.h
        ${HEADER_DIR_2}/SmartSelection.h
        ${HEADER_DIR_2}/SMP.h
        ${HEADER_DIR_2}/SplitFacet.h
        ${HEADER_DIR_2}/Vector.h
        ${HEADER_DIR_2}/VertexCoordinates.h
        ${HEADER_DIR_1}/versionId.h
        ${HEADER_DIR_2}/Worker.h
        ${HEADER_DIR_2}/GeometryConverter.h
        )

if(WITH_GPU)
    add_definitions(-DGPUCOMPABILITY)
    set(SRC_FILES ${SRC_FILES} ${CPP_DIR_1}/GPUSim/MolflowModelParser.cpp)
    set(HEADER_FILES ${HEADER_FILES} ${HEADER_DIR_1}/GPUSim/MolflowModelParser.h)
endif(WITH_GPU)

# set the path to the library folder
link_directories(
        ${LINK_DIR_2}
        ${LINK_DIR_3}
        )
if(MSVC)
    link_directories(${LINK_DIR_1})
endif(MSVC)

#find_library(PROTOBUF_LIBRARY bzip2 HINTS ${LINK_DIR_2})

# Add executable to build.
add_executable(${PROJECT_NAME} WIN32 ${SRC_FILES} ${HEADER_FILES} ${GLOB_FILES})

target_include_directories(${PROJECT_NAME} PRIVATE
        ${HEADER_DIR_1}
        ${HEADER_DIR_2}
        ${HEADER_DIR_3}
        ${HEADER_DIR_4}
        ${HEADER_DIR_5}
        ${HEADER_DIR_6}
        ${HEADER_DIR_7}
        ${HEADER_DIR_8}
        )

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)


if (MSVC)
    set_target_properties(${PROJECT_NAME} PROPERTIES
            WIN32_EXECUTABLE YES
            LINK_FLAGS "/ENTRY:mainCRTStartup"
            )
    # Multi-processor compilation
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
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS "/SUBSYSTEM:CONSOLE ${LINK_FLAGS}")
    endif()

    message(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
    link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
    #target_link_libraries(${PROJECT_NAME} clipper.lib pugixml.lib nativefiledialog.lib sdl_savepng.lib truncatedgaussian.lib png_zlib_win_library.lib)
    #target_link_libraries(${PROJECT_NAME} libgsl_dll.lib libgslcblas_dll.lib SDL2.lib SDL2main.lib libcurl_a_x64.lib opengl32.lib Advapi32.lib user32.lib shell32.lib kernel32.lib ZipLib.lib zlib.lib lzma.lib)
    # 1. link against external libs
    target_link_libraries(${PROJECT_NAME} bzip2.lib
            libcurl_a_x64.lib
            SDL2.lib
            SDL2main.lib
            libgsl.lib
            libgslcblas.lib
            lzma.lib
            ZipLib.lib
            #zlib.lib
            opengl32.lib
            user32.lib
            shell32.lib
            ole32.lib)
    # 2. link against our own libs
    target_link_libraries(${PROJECT_NAME} pugixml clipper nativefiledialog sdl_savepng truncatedgaussian png_zlib_win_library)
endif(MSVC)
if(NOT MSVC)

    if(APPLE)
        #link to self-build sdl shared lib
        target_link_libraries(${PROJECT_NAME} "-framework AppKit")

        #add_library( libSDL2 SHARED IMPORTED GLOBAL)
        #get_filename_component(ABS_LINK_DIR_2 "${LINK_DIR_2}" REALPATH)
        #set_target_properties( libSDL2 PROPERTIES IMPORTED_LOCATION ${ABS_LINK_DIR_2}/libSDL2-2.0.dylib)
        #target_link_libraries(${PROJECT_NAME} libSDL2) # from ./lib/

    else()
        # Use the package PkgConfig to detect GTK+ headers/library files
        FIND_PACKAGE(PkgConfig REQUIRED)
        PKG_CHECK_MODULES(GTK3 REQUIRED gtk+-3.0)

        # Setup CMake to use GTK+, tell the compiler where to look for headers
        # and to the linker where to look for libraries
        INCLUDE_DIRECTORIES(${GTK3_INCLUDE_DIRS})
        LINK_DIRECTORIES(${GTK3_LIBRARY_DIRS})

        # Add other flags to the compiler
        ADD_DEFINITIONS(${GTK3_CFLAGS_OTHER})

        target_link_libraries(${PROJECT_NAME} ${GTK3_LIBRARIES})

        find_package(X11 REQUIRED)
        target_include_directories(${PROJECT_NAME} PRIVATE ${X11_INCLUDE_DIRS})
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

    #GSLCBLAS or BLAS?
    #find_package(GSLCBLAS REQUIRED)
    #target_include_directories(${PROJECT_NAME} PRIVATE ${GSLCBLAS_INCLUDE_DIRS})

    #find_package(stdc++fs REQUIRED)
    #target_include_directories(${PROJECT_NAME} PRIVATE ${stdc++fs_INCLUDE_DIRS})

    target_link_libraries(${PROJECT_NAME} ${OPENGL_LIBRARIES})
    target_link_libraries(${PROJECT_NAME} ${SDL2_LIBRARIES})
    target_link_libraries(${PROJECT_NAME} ${PNG_LIBRARIES})
    target_link_libraries(${PROJECT_NAME} ${GSL_LIBRARIES})
    target_link_libraries(${PROJECT_NAME} ${CURL_LIBRARIES})
    target_link_libraries(${PROJECT_NAME} Threads::Threads)
    target_link_libraries(${PROJECT_NAME} ${X11_LIBRARIES})
    #target_link_libraries(${PROJECT_NAME} ${GSLCBLAS_LIBRARIES})

    #for shared memory
    find_library(LIBRT rt)
    if(LIBRT)
        target_link_libraries(${PROJECT_NAME} ${LIBRT})
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
    target_link_libraries(${PROJECT_NAME} libzip) # from ./lib/

    #use prebuild nfd for GCC on MacOS
    if(APPLE AND "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        add_library( libnfd STATIC IMPORTED GLOBAL)
        set_target_properties( libnfd PROPERTIES IMPORTED_LOCATION ${ABS_LINK_DIR_1}/libnfd.a )
        target_link_libraries(${PROJECT_NAME} libnfd) # from ./lib/
    endif()

    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        target_link_libraries(${PROJECT_NAME} c++fs)
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
        #don´t add anything for filesystem
    else()
        target_link_libraries(${PROJECT_NAME} stdc++fs)
    endif()

    target_link_libraries(${PROJECT_NAME} pugixml clipper sdl_savepng truncatedgaussian)
    target_link_libraries(${PROJECT_NAME} nativefiledialog)
    endif(NOT MSVC)

# Clear previous build to prevent remaining (old) files
# file(REMOVE_RECURSE ${CMAKE_EXECUTABLE_OUTPUT_DIRECTORY})

IF (WIN32)
    file(GLOB DLL_FILFES
            ${DLL_DIR}/*.dll
            )
    file(COPY ${DLL_FILES}
            DESTINATION ${CMAKE_EXECUTABLE_OUTPUT_DIRECTORY})
ENDIF()

set(COPY_FILES ${COPY_DIR}/desorption_yields
        ${COPY_DIR}/images
        ${COPY_DIR}/parameter_catalog)

IF (WIN32)
    set(COPY_FILES ${COPY_FILES}
            ${COPY_DIR}/7za.exe
            )
ELSEIF(NOT APPLE)
    set(COPY_FILES ${COPY_FILES}
            ${COPY_DIR}/7za
            )
ENDIF()

IF(${OS_NAME} STREQUAL "linux_fedora")
    set(COPY_FILES ${COPY_FILES}
            ${COPY_DIR}/launch_molflow.sh
            )
ENDIF()

#Copy our own libSDL2 to /lib
IF(APPLE)
    #    file(COPY
    #	 ${LINK_DIR_2}/libSDL2-2.0.dylib
    #	DESTINATION
    #	${CMAKE_EXECUTABLE_OUTPUT_DIRECTORY}/lib
    #            )
ENDIF()

set(COPY_FILES ${COPY_FILES}
        ${COPY_DIR}/updater_config_default_${OS_NAME}.xml
        )

file(COPY ${COPY_FILES}
        DESTINATION ${CMAKE_EXECUTABLE_OUTPUT_DIRECTORY})

file(RENAME
        ${CMAKE_EXECUTABLE_OUTPUT_DIRECTORY}/updater_config_default_${OS_NAME}.xml
        ${CMAKE_EXECUTABLE_OUTPUT_DIRECTORY}/updater_config.xml)

message("COPIED: " ${COPY_FILES})