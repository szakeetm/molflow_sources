cmake_minimum_required(VERSION 3.12.2 FATAL_ERROR)

# GUI-specific folders files, builds molflow executable
set(PROJECT_NAME molflow)

set(UI_DIR ${CPP_DIR_SRC}/Interface)

project(${PROJECT_NAME} CXX)
set(CMAKE_CXX_STANDARD 17)

# Definition of Macros
add_definitions(
        -DCURL_STATICLIB
        -D_MBCS
        -D_CRT_SECURE_NO_WARNINGS
        -D_CRT_NONSTDC_NO_DEPRECATE
)
IF (WIN32)
    add_definitions(
            -D_WITHGL
    )
ENDIF()

if(MSVC)

else() #not MSVC
    if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
    endif()
endif()

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
        ${CPP_DIR_SRC}/MolflowWorker.cpp
)

# Add executable to build.
if(WIN32)
    add_executable(${PROJECT_NAME} WIN32 ${SRC_FILES})
else()
    add_executable(${PROJECT_NAME} ${SRC_FILES} )
endif()

target_link_libraries(${PROJECT_NAME} PUBLIC src_shared_gui)
target_link_libraries(${PROJECT_NAME} PUBLIC molflow_core)

target_include_directories(${PROJECT_NAME} PRIVATE
        ${HEADER_DIR_SRC}
        ${UI_DIR}
        ${HEADER_DIR_SRC_SHARED}
        ${IMGUI_DIR}
)

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
        target_link_libraries(${PROJECT_NAME} PRIVATE "-framework AppKit")
    else()

    endif()

    #for shared memory
    find_library(LIBRT rt)
    if(LIBRT)
        target_link_libraries(${PROJECT_NAME} PRIVATE ${LIBRT})
    endif()

endif()