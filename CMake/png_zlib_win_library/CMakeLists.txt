cmake_minimum_required(VERSION 3.12.2 FATAL_ERROR)

##################### Variables ############################
# Change if you want modify path or other values           #
############################################################

# Project
#get_filename_component(PROJECT_DIR "${CMAKE_CURRENT_SOURCE_DIR}" ABSOLUTE)
#set(DEPENDENCIES_DIR ${PROJECT_DIR}/dependencies)
set(PROJECT_NAME png_zlib_win_library)

# Folders files
set(C_DIR_1 ../../include/windows_only/zlib/)
set(HEADER_DIR_1 ../../include/windows_only/zlib/)
set(C_DIR_2 ../../include/windows_only/png/)
set(HEADER_DIR_2 ../../include/windows_only/png/)

################# CMake Project ############################
# The main options of project                              #
############################################################

project(${PROJECT_NAME} C)

############### Files & Targets ############################
# Files of project and target to build                     #
############################################################

# Source Files
file(GLOB SRC_FILES
        ${C_DIR_1}/*.c
        ${C_DIR_2}/*.c
        ${HEADER_DIR_1}/*.h
        ${HEADER_DIR_2}/*.h
        )

set(SRC_FILES
        ${C_DIR_1}/adler32.c
        ${C_DIR_1}/compress.c
        ${C_DIR_1}/crc32.c
        ${C_DIR_1}/deflate.c
        ${C_DIR_1}/gzclose.c
        ${C_DIR_1}/gzlib.c
        ${C_DIR_1}/gzread.c
        ${C_DIR_1}/gzwrite.c
        ${C_DIR_1}/infback.c
        ${C_DIR_1}/inffast.c
        ${C_DIR_1}/inflate.c
        ${C_DIR_1}/inftrees.c
        ${C_DIR_1}/trees.c
        ${C_DIR_1}/uncompr.c
        ${C_DIR_1}/zutil.c
        ${HEADER_DIR_1}/crc32.h
        ${HEADER_DIR_1}/deflate.h
        ${HEADER_DIR_1}/gzguts.h
        ${HEADER_DIR_1}/inffast.h
        ${HEADER_DIR_1}/inffixed.h
        ${HEADER_DIR_1}/inflate.h
        ${HEADER_DIR_1}/inftrees.h
        ${HEADER_DIR_1}/trees.h
        ${HEADER_DIR_1}/zconf.h
        ${HEADER_DIR_1}/zlib.h
        ${HEADER_DIR_1}/zutil.h

        ${C_DIR_2}/png.c
        ${C_DIR_2}/pngerror.c
        ${C_DIR_2}/pngget.c
        ${C_DIR_2}/pngmem.c
        ${C_DIR_2}/pngpread.c
        ${C_DIR_2}/pngread.c
        ${C_DIR_2}/pngrio.c
        ${C_DIR_2}/pngrtran.c
        ${C_DIR_2}/pngrutil.c
        ${C_DIR_2}/pngset.c
        ${C_DIR_2}/pngtest.c
        ${C_DIR_2}/pngtrans.c
        ${C_DIR_2}/pngwio.c
        ${C_DIR_2}/pngwrite.c
        ${C_DIR_2}/pngwtran.c
        ${C_DIR_2}/pngwutil.c
        ${HEADER_DIR_2}/png.h
        ${HEADER_DIR_2}/pngconf.h
        ${HEADER_DIR_2}/pngdebug.h
        ${HEADER_DIR_2}/pnginfo.h
        ${HEADER_DIR_2}/pnglibconf.h
        ${HEADER_DIR_2}/pngpriv.h
        ${HEADER_DIR_2}/pngstruct.h
        )

# Add library to build.
add_library(${PROJECT_NAME} STATIC
        ${SRC_FILES}
)

target_include_directories(${PROJECT_NAME} PRIVATE
        ${HEADER_DIR_1}
        ${HEADER_DIR_2}
        )

######################### Flags ############################
# Defines Flags for Windows and Linux                      #
############################################################

# Preprocessor definitions
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_definitions(${PROJECT_NAME} PRIVATE)
    if(MSVC)
        target_compile_options(${PROJECT_NAME} PRIVATE  /W3 /MDd /Od /EHsc)
        #for .pdb debugging files
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /Zi")
        set(CMAKE_SHARED_LINKER_FLAGS_DEBUG "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} /DEBUG:FULL /OPT:REF /OPT:ICF")
    endif()
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_definitions(${PROJECT_NAME} PRIVATE)
    if(MSVC)
        target_compile_options(${PROJECT_NAME} PRIVATE  /W3 /Oi /Gy /EHsc)
    endif()
endif()

########### Link & Dependencies ############################
# Add project dependencies and Link to project             #
############################################################
set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE C)

