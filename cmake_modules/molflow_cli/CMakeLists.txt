cmake_minimum_required(VERSION 3.12.2 FATAL_ERROR)

#Command-line interface, builds molflowCLI executable
set(PROJECT_NAME molflowCLI)

############## CMake Project ################
#        The main options of project        #
#############################################

project(${PROJECT_NAME} CXX)
set(CMAKE_CXX_STANDARD 17)

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
    if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
    endif()
endif()

################ Files ################
#   --   Add files to project.   --   #
#######################################

set(SRC_FILES
        ${CPP_DIR_SRC}/MolflowCLI.cpp
        ${CPP_DIR_SRC}/Initializer.cpp

        ${IO_DIR}/CSVExporter.cpp
        ${IO_DIR}/CSVExporter.h

        ${HELPER_DIR}/GLProgress_CLI.cpp
        ${HELPER_DIR}/GLProgress_CLI.hpp
        ${HELPER_DIR}/GLProgress_abstract.cpp
        ${HELPER_DIR}/GLProgress_abstract.hpp

)

add_executable(${PROJECT_NAME} ${SRC_FILES})

target_include_directories(${PROJECT_NAME} PRIVATE
        ${HEADER_DIR_SRC}
        ${HEADER_DIR_SRC_SHARED}
        ${HEADER_DIR_INCLUDE}
        ${INTERFACE_DIR} #geometrytypes.h for example
        ${HEADER_DIR_EXTERNAL} #Ziplib etc.
)

target_link_libraries(${PROJECT_NAME} PUBLIC molflow_core)
# set the path to the library folder
if(MSVC)
    # Add Whole Program Optimization and Link Time Code Generation
    if(NOT CMAKE_BUILD_TYPE MATCHES Debug)
        set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "/GL") #Incompatible with Edit and continue
        set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS "/LTCG") #no whole program optimization for debug
    endif()
    
    # Add Windows Console output
    message(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
    target_link_directories(${PROJECT_NAME} PRIVATE ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

endif(MSVC)
if(NOT MSVC)

endif(NOT MSVC)