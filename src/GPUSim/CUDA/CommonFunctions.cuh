//
// Created by pbahr on 09/03/2020.
// TODO: Find a better name
//

#ifndef MOLFLOW_PROJ_COMMONFUNCTIONS_CUH
#define MOLFLOW_PROJ_COMMONFUNCTIONS_CUH

#include <math_constants.h>
#include "cuda_runtime.h"
#include "helper_math.h"

#include "LaunchParams.h"
#include "GPUDefines.h"

static __forceinline__ __device__
float3 getNewDirection(flowgpu::MolPRD& hitData, const flowgeom::Polygon& poly,
                       const float* randFloat, unsigned int& randInd, unsigned int& randOffset)
{

    // generate ray direction
/*#ifdef BOUND_CHECK
    if(randInd + randOffset < 0 || randInd + randOffset >= optixLaunchParams.simConstants.nbRandNumbersPerThread*optixLaunchParams.simConstants.size.x*optixLaunchParams.simConstants.size.y){
            printf("randInd %u is out of bounds\n", randInd + randOffset);
        }
#endif*/
    const float theta = acosf(sqrtf((double)randFloat[(unsigned int)(randInd + randOffset++)]));


/*#ifdef BOUND_CHECK
    if(randInd + randOffset < 0 || randInd + randOffset >= optixLaunchParams.simConstants.nbRandNumbersPerThread*optixLaunchParams.simConstants.size.x*optixLaunchParams.simConstants.size.y){
            printf("randInd %u is out of bounds\n", randInd + randOffset);
        }
#endif*/
    const float phi = randFloat[(unsigned int)(randInd + randOffset++)] * 2.0f * CUDART_PI_F;


    const float u = sinf(theta)*cosf(phi);
    const float v = sinf(theta)*sinf(phi);
    const float n = cosf(theta);

    /*float3 nU = rayGenData->poly[facIndex].nU;
    float3 nV = rayGenData->poly[facIndex].nV;
    float3 N = rayGenData->poly[facIndex].N;*/
    const float3 nU = poly.nU;
    const float3 nV = poly.nV;
    const float3 N = poly.N;

    float3 rayDir = u*nU + v*nV + n*N;
    /*if ((randFloat[(unsigned int)(randInd + randOffset-2)] == 1.0f ||
            randFloat[(unsigned int)(randInd + randOffset - 1)] == 1.0f) ||
            (randFloat[(unsigned int)(randInd + randOffset-2)] == 0.0f ||
             randFloat[(unsigned int)(randInd + randOffset - 1)] == 0.0f) ||
            (rayDir.x==0.0f && rayDir.y==0.0f && rayDir.z==1.0f)){
        printf("[%d] Ray dir like normal: [%u , %u] {%lf, %lf}{0x%0x, 0x%0x} ( %lf , %lf ) - %lf , %lf , %lf -> %lf , %lf , %lf\n",
               hitData.inSystem,randInd + randOffset-2,randInd + randOffset-1,
               randFloat[(unsigned int)(randInd + randOffset-2)] ,
               randFloat[(unsigned int)(randInd + randOffset - 1)] ,
               randFloat[(unsigned int)(randInd + randOffset-2)] ,
               randFloat[(unsigned int)(randInd + randOffset - 1)] ,phi,theta,u,v,n, rayDir.x,rayDir.y,rayDir.z);
    }*/

    return u*nU + v*nV + n*N;
    /*if (rayDir.x != 0.0) rayDir.x = 1.0 / rayDir.x;
    if (rayDir.y != 0.0) rayDir.y = 1.0 / rayDir.y;
    if (rayDir.z != 0.0) rayDir.z = 1.0 / rayDir.z;
    rayDir = float3(-1.0f,-1.0f,-1.0f) * rayDir;*/
}

static __forceinline__ __device__
float3 getNewDirection(flowgpu::MolPRD& hitData, const flowgeom::Polygon& poly,
                       const double* randFloat, unsigned int& randInd, unsigned int& randOffset)
{

    // generate ray direction
/*#ifdef BOUND_CHECK
    if(randInd + randOffset < 0 || randInd + randOffset >= optixLaunchParams.simConstants.nbRandNumbersPerThread*optixLaunchParams.simConstants.size.x*optixLaunchParams.simConstants.size.y){
            printf("randInd %u is out of bounds\n", randInd + randOffset);
        }
#endif*/
    const float theta = acos(sqrt(randFloat[(unsigned int)(randInd + randOffset++)]));

/*#ifdef BOUND_CHECK
    if(randInd + randOffset < 0 || randInd + randOffset >= optixLaunchParams.simConstants.nbRandNumbersPerThread*optixLaunchParams.simConstants.size.x*optixLaunchParams.simConstants.size.y){
            printf("randInd %u is out of bounds\n", randInd + randOffset);
        }
#endif*/
    const float phi = randFloat[(unsigned int)(randInd + randOffset++)] * 2.0 * CUDART_PI;


    const float u = sinf(theta)*cosf(phi);
    const float v = sinf(theta)*sinf(phi);
    const float n = cosf(theta);

    /*float3 nU = rayGenData->poly[facIndex].nU;
    float3 nV = rayGenData->poly[facIndex].nV;
    float3 N = rayGenData->poly[facIndex].N;*/
    const float3 nU = poly.nU;
    const float3 nV = poly.nV;
    const float3 N = poly.N;

    /*if ((randFloat[(unsigned int)(randInd + randOffset-2)] == 1.0f ||
         randFloat[(unsigned int)(randInd + randOffset - 1)] == 1.0f) ||
        (randFloat[(unsigned int)(randInd + randOffset-2)] == 0.0f ||
         randFloat[(unsigned int)(randInd + randOffset - 1)] == 0.0f) ||
        (rayDir.x==0.0f && rayDir.y==0.0f && rayDir.z==1.0f)){
        printf("[%d] Ray dir like normal: [%u , %u] {%lf, %lf}{0x%0x, 0x%0x} ( %lf , %lf ) - %lf , %lf , %lf -> %lf , %lf , %lf\n",
               hitData.inSystem,randInd + randOffset-2,randInd + randOffset-1,
               randFloat[(unsigned int)(randInd + randOffset-2)] ,
               randFloat[(unsigned int)(randInd + randOffset - 1)] ,
               randFloat[(unsigned int)(randInd + randOffset-2)] ,
               randFloat[(unsigned int)(randInd + randOffset - 1)] ,phi,theta,u,v,n, rayDir.x,rayDir.y,rayDir.z);
    }*/

    return u*nU + v*nV + n*N;
    /*if (rayDir.x != 0.0) rayDir.x = 1.0 / rayDir.x;
    if (rayDir.y != 0.0) rayDir.y = 1.0 / rayDir.y;
    if (rayDir.z != 0.0) rayDir.z = 1.0 / rayDir.z;
    rayDir = float3(-1.0f,-1.0f,-1.0f) * rayDir;*/
}

#endif //MOLFLOW_PROJ_COMMONFUNCTIONS_CUH
