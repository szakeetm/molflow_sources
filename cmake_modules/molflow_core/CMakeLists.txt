cmake_minimum_required(VERSION 3.12.2 FATAL_ERROR)

#shared GUI and CLI part, built as library that GUI or CLI links against
set(PROJECT_NAME molflow_core)

project(${PROJECT_NAME} CXX)
set(CMAKE_CXX_STANDARD 17)

# Definition of Macros
add_definitions(
        -D_MBCS
        -D_CRT_SECURE_NO_WARNINGS
)
if(MSVC)

else() #not MSVC
    if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++") #default on macOS: libstdc++
    endif()
endif()

set(SRC_FILES
        
        ${GLAPP_DIR}/GLFormula.cpp

        ${IO_DIR}/LoaderXML.cpp
        ${IO_DIR}/WriterXML.cpp

        ${CPP_DIR_SRC}/ParameterParser.cpp
        ${CPP_DIR_SRC}/Initializer.cpp
        ${CPP_DIR_SRC}/TimeMoments.cpp
        ${CPP_DIR_SRC}/MolflowTypes.cpp
        

       
        
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

target_link_libraries(${PROJECT_NAME} PUBLIC src_shared_core)

target_link_libraries(${PROJECT_NAME} PUBLIC truncated_normal)

target_include_directories(${PROJECT_NAME} PRIVATE
        ${HEADER_DIR_SRC}
        ${HEADER_DIR_SRC_SHARED}
        ${HEADER_DIR_INCLUDE}
        ${INTERFACE_DIR} #geometrytypes.h for example
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