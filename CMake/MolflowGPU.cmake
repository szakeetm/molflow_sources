############## CMake Project ################
#        The main options of project        #
#############################################

if(CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    option(DEBUG_MISS "Enable Debugging for Miss Program" ON)
    option(DEBUG_COUNT "Enable counters for intersection results" OFF)
    option(DEBUG_POS "Enable generating hit positions from GPU calculations" ON)
    option(DEBUG_LEAKPOS "Enable generating hit positions for leaks from GPU calculations" ON)
    option(DEBUG_BOUNDS "Enable bound checks for CUDA kernels" ON)
endif()
option(WITH_DESORPEXIT "Enable exit on desorption limit" ON)
option(WITH_TRIANGLES "Enable calculations with triangles only" ON)
option(WITH_TEXTURES "Enable textures" ON)
option(WITH_PROFILES "Enable profiles" ON)
option(WITH_TRANS "Enable transparent SBT" ON)
option(WITH_NBOUNCE "Enable NBBOUNCE Counter" ON)

option(USE_RANDOM_NUMBER_TYPE_64 "Use double instead of float for random numbers" ON)
option(USE_COUNTER_TYPE_64 "Use 64bit instead of 32bit precision for the counter structure" ON)
MESSAGE("[GPU_BUILD_OPTION] USE_RANDOM_NUMBER_TYPE_64: ${USE_RANDOM_NUMBER_TYPE_64}")
MESSAGE("[GPU_BUILD_OPTION] USE_COUNTER_TYPE_64: ${USE_COUNTER_TYPE_64}")

if(CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    MESSAGE("[GPU_BUILD_OPTION] DEBUG_MISS: ${DEBUG_MISS}")
    MESSAGE("[GPU_BUILD_OPTION] DEBUG_COUNT: ${DEBUG_COUNT}")
    MESSAGE("[GPU_BUILD_OPTION] DEBUG_POS: ${DEBUG_POS}")
    MESSAGE("[GPU_BUILD_OPTION] DEBUG_LEAKPOS: ${DEBUG_LEAKPOS}")
    MESSAGE("[GPU_BUILD_OPTION] DEBUG_BOUNDS: ${DEBUG_BOUNDS}")
endif()
MESSAGE("[GPU_BUILD_OPTION] WITH_DESORPEXIT: ${WITH_DESORPEXIT}")
MESSAGE("[GPU_BUILD_OPTION] WITH_TRIANGLES: ${WITH_TRIANGLES}")
MESSAGE("[GPU_BUILD_OPTION] WITH_TEXTURES: ${WITH_TEXTURES}")
MESSAGE("[GPU_BUILD_OPTION] WITH_PROFILES: ${WITH_PROFILES}")
MESSAGE("[GPU_BUILD_OPTION] WITH_TRANS: ${WITH_TRANS}")
MESSAGE("[GPU_BUILD_OPTION] WITH_NBOUNCE: ${WITH_NBOUNCE}")

# Definition of Macros
if(WIN32)
    add_definitions(
            -DWIN
    )
endif(WIN32)

add_definitions(
        -D_MBCS
        -D_CRT_SECURE_NO_WARNINGS

        -DMOLFLOW
        -DGPU_MODE
        -DGPUCOMPABILITY
)

if(CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    MESSAGE("[NVCC] Adding Debug flag...")
    ADD_DEFINITIONS(-DDEBUG)
    ADD_DEFINITIONS(-D_DEBUG)
    set(nvcc_flags ${nvcc_flags} -DDEBUG -D_DEBUG)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")


    if (DEBUG_MISS)
        add_definitions(-DDEBUGMISS)
        set(nvcc_flags ${nvcc_flags} -DDEBUGMISS)
    endif (DEBUG_MISS)

    if (DEBUG_COUNT)
        add_definitions(-DDEBUGCOUNT)
        set(nvcc_flags ${nvcc_flags} -DDEBUGCOUNT)
    endif (DEBUG_COUNT)

    if (DEBUG_POS)
        add_definitions(-DDEBUGPOS)
        set(nvcc_flags ${nvcc_flags} -DDEBUGPOS)
    endif (DEBUG_POS)

    if (DEBUG_LEAKPOS)
        add_definitions(-DDEBUGLEAKPOS)
        set(nvcc_flags ${nvcc_flags} -DDEBUGLEAKPOS)
    endif (DEBUG_LEAKPOS)
    if (DEBUG_BOUNDS)
        add_definitions(-DBOUND_CHECK)
        set(nvcc_flags ${nvcc_flags} -DBOUND_CHECK)
    endif (DEBUG_BOUNDS)
ENDIF()

if (WITH_DESORPEXIT)
    add_definitions(-DWITHDESORPEXIT)
    set(nvcc_flags ${nvcc_flags} -DWITHDESORPEXIT)
endif (WITH_DESORPEXIT)
if (WITH_NBOUNCE)
    add_definitions(-DGPUNBOUNCE)
    set(nvcc_flags ${nvcc_flags} -DGPUNBOUNCE)
endif (WITH_NBOUNCE)
if (WITH_TRIANGLES)
    add_definitions(-DWITHTRIANGLES)
    set(nvcc_flags ${nvcc_flags} -DWITHTRIANGLES)
endif (WITH_TRIANGLES)
if (WITH_TEXTURES)
    add_definitions(-DWITH_TEX)
    set(nvcc_flags ${nvcc_flags} -DWITH_TEX)
endif (WITH_TEXTURES)
if (WITH_PROFILES)
    add_definitions(-DWITH_PROF)
    set(nvcc_flags ${nvcc_flags} -DWITH_PROF)
endif (WITH_PROFILES)
if (WITH_TRANS)
    add_definitions(-DWITH_TRANS)
    set(nvcc_flags ${nvcc_flags} -DWITH_TRANS)
endif (WITH_TRANS)
if (USE_RANDOM_NUMBER_TYPE_64)
    add_definitions(-DRNG64)
    set(nvcc_flags ${nvcc_flags} -DRNG64)
endif (USE_RANDOM_NUMBER_TYPE_64)
if (USE_COUNTER_TYPE_64)
    add_definitions(-DHIT64)
    set(nvcc_flags ${nvcc_flags} -DHIT64)
endif (USE_COUNTER_TYPE_64)

#Flags on CUDA 10 for maximum compatibility
#set(CMAKE_CUDA_FLAGS ${CMAKE_CUDA_FLAGS} " -arch=sm_50 -gencode=arch=compute_50,code=sm_50 -gencode=arch=compute_52,code=sm_52 -gencode=arch=compute_60,code=sm_60 -gencode=arch=compute_61,code=sm_61 -gencode=arch=compute_70,code=sm_70 -gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_75,code=compute_75")
#set(CMAKE_CUDA_FLAGS "-arch=sm_52 -gencode=arch=compute_52,code=sm_52 -gencode=arch=compute_60,code=sm_60 -gencode=arch=compute_61,code=sm_61 -gencode=arch=compute_70,code=sm_70 -gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_75,code=compute_75")
#set(cuda_flags -gencode=arch=compute_70,code=sm_70 -gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_75,code=compute_75)
set(cuda_flags -gencode=arch=compute_75,code=compute_75)

# Enable fast math
#if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
option(ODL_CUDA_USE_FAST_MATH "Enable fast math in cuda (can decrease precision)" OFF)
if(ODL_CUDA_USE_FAST_MATH)
    MESSAGE("[NVCC] Using fast math...")
    set(cuda_flags ${cuda_flags} --use_fast_math)
endif()
option(ODL_CUDA_STREAM_PER_THREAD "Enable stream per-thread in cuda" ON)
if(ODL_CUDA_STREAM_PER_THREAD)
    MESSAGE("[NVCC] Using stream per-thread...")
    set(cuda_flags ${cuda_flags} --default-stream per-thread)
endif()
#endif()