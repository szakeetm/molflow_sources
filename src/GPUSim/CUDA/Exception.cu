// Created by pbahr

#include "CommonFunctions.cuh"
#include "PerRayData.h"
#include "LaunchParams.h"
#include "GPUDefines.h" // for NB_RAND
#include "helper_math.h"
#include <optix.h>
//using namespace flowgpu;

namespace flowgpu {

    /*! launch parameters in constant memory, filled in by optix upon
        optixLaunch (this gets filled in from the buffer we pass to
        optixLaunch) */
    extern "C" __constant__ LaunchParams optixLaunchParams;


    extern "C" __global__ void __exception__all()
    {
        // This assumes that the launch dimensions are matching the size of the output buffer.
        const uint3 theLaunchDim   = optixGetLaunchDimensions();
        const uint3 theLaunchIndex = optixGetLaunchIndex();

        const int theExceptionCode = optixGetExceptionCode();
        printf("Exception %d [%#08x] at (%u, %u)\n", theExceptionCode, theExceptionCode, theLaunchIndex.x, theLaunchIndex.y);


        const unsigned int fbIndex = getWorkIndex();
        atomicAdd(optixLaunchParams.sharedData.missCounter, 1);

        // Reset particle
        MolPRD& prd = optixLaunchParams.perThreadData.currentMoleculeData[fbIndex];
        initParticle(prd);

        //const unsigned int index = theLaunchIndex.y * theLaunchDim.x + theLaunchIndex.x;
        //sysParameter.outputBuffer[index] = make_float4(1000000.0f, 0.0f, 1000000.0f, 1.0f); // super magenta
    }

} // ::flowgpu
