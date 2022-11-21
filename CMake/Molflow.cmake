IF (WIN32)
    set(OS_NAME "win")
    set(OS_RELPATH "")
ELSEIF(APPLE)
    set(OS_NAME "mac")
    set(OS_RELPATH "")
ELSE()
    IF(os_version_suffix MATCHES "\\.el[1-9]")
        set(OS_NAME "linux_fedora")
    ELSE()
        set(OS_NAME "linux_debian")
    ENDIF()
    set(OS_RELPATH "")
ENDIF()

IF (CMAKE_BUILD_TYPE MATCHES Debug|RelWithDebInfo)
    set(MY_BUILD_TYPE "debug")
ELSE()
    set(MY_BUILD_TYPE "release")
ENDIF()

############## Artefacts Output #################
# Defines outputs , depending Debug or Release. #
#################################################


set(CMAKE_LIBRARY_OUTPUT_DIRECTORY    "${CMAKE_BINARY_DIR}/${OS_RELPATH}/lib/")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY    "$${CMAKE_BINARY_DIR}/{OS_RELPATH}/lib/")
set(CMAKE_EXECUTABLE_OUTPUT_DIRECTORY "$${CMAKE_BINARY_DIR}/{OS_RELPATH}/bin/")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY    "$${CMAKE_BINARY_DIR}/{OS_RELPATH}/bin/")

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_EXECUTABLE_OUTPUT_DIRECTORY}) #to build executable in main folder

# Messages
message("${PROJECT_NAME}: MAIN PROJECT: ${CMAKE_PROJECT_NAME}")
message("${PROJECT_NAME}: CURR PROJECT: ${CMAKE_CURRENT_SOURCE_DIR}")
message("${PROJECT_NAME}: CURR BIN DIR: ${CMAKE_CURRENT_BINARY_DIR}")
# Messages
message("${PROJECT_NAME}: CMAKE_LIBRARY_OUTPUT_DIRECTORY: ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
message("${PROJECT_NAME}: CMAKE_ARCHIVE_OUTPUT_DIRECTORY: ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}")
message("${PROJECT_NAME}: CMAKE_EXECUTABLE_OUTPUT_DIRECTORY: ${CMAKE_EXECUTABLE_OUTPUT_DIRECTORY}")


############## CMake Project ################
#        The main options of project        #
#############################################

# Definition of Macros
add_definitions(
        -DMOLFLOW
)

if (OS_NAME STREQUAL "linux_fedora" )
    add_definitions(
            -D__LINUX_FEDORA
    )
elseif(OS_NAME STREQUAL "linux_debian")
    add_definitions(
            -D__LINUX_DEBIAN
    )
endif ()

if(NOT MSVC)
    add_compile_options(
            "$<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>:-Wno-tautological-undefined-compare;-Wno-inconsistent-missing-override>"
            #is valid for C++/ObjC++ but not for C
            $<$<COMPILE_LANGUAGE:CXX>:-Wno-reorder>
            #-w
            #-Wextra
            -Wconversion
            -Wno-write-strings
            -Wno-unused
            -pedantic
            #-Werror -Wno-error=uninitialized
        "$<$<CONFIG:RELEASE>:-O3>"
        "$<$<CONFIG:DEBUG>:-O0>"
        "$<$<CONFIG:DEBUG>:-ggdb3>"
        "$<$<CONFIG:RELWITHDEBINFO>:-O2>"
        "$<$<CONFIG:RELWITHDEBINFO>:-ggdb3>"
        "$<$<CONFIG:RELWITHDEBINFO>:-g>"
    )
else() #MSVC
    add_compile_options(
        /W3
    )
    add_compile_options(
        "$<$<CONFIG:Release>:/GL;/O2;/EHsc>"
        "$<$<CONFIG:Debug>:/MDd;/Od;/EHsc>"
    )
    # Multi-processor compilation
    add_compile_options(
        "$<$<CONFIG:Debug>:/MP>"
        "$<$<CONFIG:Release>:/MP>"
    )
endif()
#disable generation of appname.manifest file
#alternative: use /MANIFEST:EMBED
if(MSVC)
    set(CMAKE_EXE_LINKER_FLAGS    "${CMAKE_EXE_LINKER_FLAGS} /MANIFEST:NO")
    set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} /MANIFEST:NO")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /MANIFEST:NO")
endif(MSVC)

set(COPY_DIR ./copy_to_build/)

# Windows DLL files
IF (WIN32)
    set(DLL_DIR ./lib_external/win/dll)
    file(GLOB DLL_FILES
            ${DLL_DIR}/*.dll
            )
    file(COPY ${DLL_FILES}
            DESTINATION ${CMAKE_EXECUTABLE_OUTPUT_DIRECTORY})

    message("COPIED: " ${DLL_FILES})
    message("    TO: " ${CMAKE_EXECUTABLE_OUTPUT_DIRECTORY})
ENDIF()

set(COPY_FILES ${COPY_DIR}/desorption_yields
        ${COPY_DIR}/images
        ${COPY_DIR}/parameter_catalog
        ${COPY_DIR}/Roboto-Medium.ttf
        ${COPY_DIR}/DroidSans.ttf
        ${COPY_DIR}/FreeMono.ttf
        ${COPY_DIR}/fa-regular-400.ttf
        ${COPY_DIR}/fa-solid-900.ttf
        )

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

file(COPY ${COPY_FILES}
        DESTINATION ${CMAKE_EXECUTABLE_OUTPUT_DIRECTORY})

message("COPIED: " ${COPY_FILES})