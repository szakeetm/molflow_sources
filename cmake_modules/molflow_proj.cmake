############## Artefacts Output #################
# Defines outputs , depending Debug or Release. #
#################################################

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY    "${CMAKE_BINARY_DIR}/lib/")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY    "${CMAKE_BINARY_DIR}/lib/")
set(CMAKE_EXECUTABLE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY    "${CMAKE_BINARY_DIR}/bin/")

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_EXECUTABLE_OUTPUT_DIRECTORY})

# Messages
message("${PROJECT_NAME}: Project dir: ${CMAKE_PROJECT_NAME}")
message("${PROJECT_NAME}: Source dir: ${CMAKE_CURRENT_SOURCE_DIR}")
message("${PROJECT_NAME}: Binary dir: ${CMAKE_CURRENT_BINARY_DIR}")
message("${PROJECT_NAME}: Library output dir: ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
message("${PROJECT_NAME}: Archive output dir: ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}")
message("${PROJECT_NAME}: Executable output dir: ${CMAKE_EXECUTABLE_OUTPUT_DIRECTORY}")


# Project-wide folders
set(CPP_DIR_SRC ${CMAKE_HOME_DIRECTORY}/src)
set(SIMU_DIR ${CPP_DIR_SRC}/Simulation)
set(IO_DIR ${CPP_DIR_SRC}/IO)
set(CPP_DIR_SRC_SHARED ${CMAKE_HOME_DIRECTORY}/src_shared)
set(GLAPP_DIR ${CPP_DIR_SRC_SHARED}/GLApp)
set(GLCHART_DIR ${GLAPP_DIR}/GLChart)
set(HELPER_DIR ${CPP_DIR_SRC_SHARED}/Helper)

set(HEADER_DIR_SRC ${CPP_DIR_SRC}) #same for .cpp and .h
set(HEADER_DIR_SRC_SHARED ${CPP_DIR_SRC_SHARED}) #same for .cpp and .h
set(HEADER_DIR_INCLUDE ${CMAKE_HOME_DIRECTORY}/include)

# Definition of Macros
add_definitions(
        -DMOLFLOW #to distinguish from SYNRAD in the source files
)

#required preprocessor definitions by the auto-updater
IF (WIN32)
    #WIN32 defined by default
ELSEIF(APPLE)
    #__MACOSX__ defined by default
ELSE()
    string (REGEX MATCH "\\.el[1-9]" os_version_suffix ${CMAKE_SYSTEM})
    message("-- os_version_suffix:      ${os_version_suffix}")
    IF(os_version_suffix MATCHES "\\.el[1-9]")
        add_definitions(
            -D__LINUX_FEDORA
        )
        message("Assuming Fedora-based Linux, defining __LINUX_FEDORA")
    ELSE()
        add_definitions(
            -D__LINUX_DEBIAN
        )
        message("Assuming Debian-based Linux, defining __LINUX_DEBIAN")
    ENDIF()
ENDIF()

if(USE_PROFILING)
    MESSAGE("Setting profiling flags.")
    string(APPEND CMAKE_C_FLAGS " -O0 -pg")
    string(APPEND CMAKE_CXX_FLAGS " -O0 -pg")
    string(APPEND CMAKE_EXE_LINKER_FLAGS " -pg")
    string(APPEND CMAKE_SHARED_LINKER_FLAGS " -pg")
endif(USE_PROFILING)

set(CMAKE_EXPORT_COMPILE_COMMANDS OFF) # Don't generate .json files

set(CMAKE_CXX_STANDARD 17)

if(MSVC) #MSVC
    if(CMAKE_BUILD_TYPE MATCHES Debug)
        string(REGEX REPLACE "/Zi" "/ZI" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}") #Edit and continue (/ZI) instead of default /Zi
    endif()
    add_compile_options(
        /W3
        /MP # Multi-processor compilation
        "$<$<CONFIG:Release>:/GL;/O2;/EHsc>"
        "$<$<CONFIG:Debug>:/MDd;/Od;/EHsc>"
    )

    #disable generation of appname.manifest file. Alternative: use /MANIFEST:EMBED
    set(CMAKE_EXE_LINKER_FLAGS    "${CMAKE_EXE_LINKER_FLAGS} /MANIFEST:NO")
    set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} /MANIFEST:NO")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /MANIFEST:NO")
else() #not MSVC
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
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        string(APPEND CMAKE_CXX_FLAGS " -stdlib=libc++") #instead of default macOS stdlibc++
    endif()
endif()

set(COPY_DIR copy_to_build)

# Windows DLL files (on other OS libraries are linked statically)
IF (WIN32)
    set(DLL_DIR lib_external/win/dll)
    file(GLOB DLL_FILES
            ${DLL_DIR}/*.dll
            )
    file(COPY ${DLL_FILES}
            DESTINATION ${CMAKE_EXECUTABLE_OUTPUT_DIRECTORY})

    message("COPIED: " ${DLL_FILES})
    message("    TO: " ${CMAKE_EXECUTABLE_OUTPUT_DIRECTORY})
ENDIF()

# Other files to include in the bin directory
set(COPY_FILES ${COPY_DIR}/desorption_yields
        ${COPY_DIR}/images
        ${COPY_DIR}/parameter_catalog
        ${COPY_DIR}/fonts
        )

IF (WIN32)
    set(COPY_FILES ${COPY_FILES}
            ${COPY_DIR}/7za.exe #no system-wide 7-zip installation
            )
ENDIF()

file(COPY ${COPY_FILES}
        DESTINATION ${CMAKE_EXECUTABLE_OUTPUT_DIRECTORY})

message("COPIED: " ${COPY_FILES})