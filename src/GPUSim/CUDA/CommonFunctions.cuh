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

#include <cooperative_groups.h>
namespace cg = cooperative_groups;

constexpr __device__ float origin()      { return 1.0f / 32.0f; }
constexpr __device__ float float_scale() { return 1.0f / 65536.0f; }
constexpr __device__ float int_scale()   { return 256.0f; }
// Normal points outward for rays exiting the surface, else is flipped.
static __forceinline__ __device__ float3 offset_ray(const float3 p, const float3 n){
    int3 of_i(make_int3(int_scale() * n.x, int_scale() * n.y, int_scale() * n.z));
    float3 p_i(make_float3(
            int_as_float(float_as_int(p.x)+((p.x < 0) ? -of_i.x : of_i.x)),
            int_as_float(float_as_int(p.y)+((p.y < 0) ? -of_i.y : of_i.y)),
            int_as_float(float_as_int(p.z)+((p.z < 0) ? -of_i.z : of_i.z))));
    return float3(make_float3(
            fabsf(p.x) < origin() ? p.x+float_scale()*n.x : p_i.x,
            fabsf(p.y) < origin() ? p.y+float_scale()*n.y : p_i.y,
            fabsf(p.z) < origin() ? p.z+float_scale()*n.z : p_i.z));
}

static __forceinline__ __device__ float3 offset_ray_new(const float3 p, const float3 n){
    int3 of_i(make_int3(int_scale() * n.x, int_scale() * n.y, int_scale() * n.z));
    float3 p_i(make_float3(
            int_as_float(float_as_int(p.x)+((p.x < 0) ? -of_i.x : of_i.x)),
            int_as_float(float_as_int(p.y)+((p.y < 0) ? -of_i.y : of_i.y)),
            int_as_float(float_as_int(p.z)+((p.z < 0) ? -of_i.z : of_i.z))));
    return float3(make_float3(p.x+float_scale()*n.x,
                              p.y+float_scale()*n.y,
                              p.z+float_scale()*n.z));
}

// Normal points outward for rays exiting the surface, else is flipped.
static __forceinline__ __device__ float3 offset_ray_none(const float3 p, const float3 n){
    int3 of_i(make_int3(int_scale() * n.x, int_scale() * n.y, int_scale() * n.z));
    float3 p_i(make_float3(
            int_as_float(float_as_int(p.x)+((p.x < 0) ? -of_i.x : of_i.x)),
            int_as_float(float_as_int(p.y)+((p.y < 0) ? -of_i.y : of_i.y)),
            int_as_float(float_as_int(p.z)+((p.z < 0) ? -of_i.z : of_i.z))));
    return p;
}

__device__
int atomicAggInc(int *ptr, int incVal)
{
    cg::coalesced_group g = cg::coalesced_threads();
    int prev;

    // elect the first active thread to perform atomic add
    if (g.thread_rank() == 0) {
        prev = atomicAdd(ptr, g.size() * incVal);
    }

    // broadcast previous value within the warp
    // and add each active thread’s rank to it
    prev = g.shfl(prev, 0);
    return g.thread_rank() + prev;
}

__device__
uint32_t atomicAggInc(uint32_t *ptr, uint32_t incVal)
{
    cg::coalesced_group g = cg::coalesced_threads();
    int prev;

    // elect the first active thread to perform atomic add
    if (g.thread_rank() == 0) {
        prev = atomicAdd(ptr, g.size() * incVal);
    }

    // broadcast previous value within the warp
    // and add each active thread’s rank to it
    prev = g.shfl(prev, 0);
    return g.thread_rank() + prev;
}

__device__
float atomicAggInc(float *ptr, float incVal)
{
    cg::coalesced_group g = cg::coalesced_threads();
    int prev;

    // elect the first active thread to perform atomic add
    if (g.thread_rank() == 0) {
        prev = atomicAdd(ptr, g.size() * incVal);
    }

    // broadcast previous value within the warp
    // and add each active thread’s rank to it
    prev = g.thread_rank() + g.shfl(prev, 0);
    return prev;
}

//TODO: Only non maxwell for now
static __forceinline__ __device__
float getNewVelocity(const flowgpu::Polygon& poly, const float& gasMass)
{
    return 145.469*sqrt(poly.facProps.temperature / gasMass);
}

//TODO: Only cosine for now
static __forceinline__ __device__
float3 getNewDirection(flowgpu::MolPRD& hitData, const flowgpu::Polygon& poly,
                       const float* randFloat, unsigned int& randInd, unsigned int& randOffset)
{

    // generate ray direction
/*#ifdef BOUND_CHECK
    if(randInd + randOffset < 0 || randInd + randOffset >= optixLaunchParams.simConstants.nbRandNumbersPerThread*optixLaunchParams.simConstants.size.x*optixLaunchParams.simConstants.size.y){
            printf("randInd %u is out of bounds\n", randInd + randOffset);
        }
#endif*/

    float theta = 0.0f;
    if(poly.desProps.desorbType == 2 || poly.desProps.desorbType == 0) {
        theta = acosf(sqrtf((double) randFloat[(unsigned int) (randInd + randOffset++)]));
    }
    else if (poly.desProps.desorbType == 3) {
        theta = acosf(powf((double) randFloat[(unsigned int) (randInd + randOffset++)], (1.0f / (poly.desProps.cosineExponent + 1.0f))));
    }
    else if (poly.desProps.desorbType == 1) {
        theta = acosf((double) randFloat[(unsigned int) (randInd + randOffset++)]);
    }
    else{
        printf("Unsupported desorption type! %u on %d\n", poly.desProps.desorbType, poly.parentIndex);
        return;
    }
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
float3 getNewDirection(flowgpu::MolPRD& hitData, const flowgpu::Polygon& poly,
                       const double* randFloat, unsigned int& randInd, unsigned int& randOffset)
{

    // generate ray direction
/*#ifdef BOUND_CHECK
    if(randInd + randOffset < 0 || randInd + randOffset >= optixLaunchParams.simConstants.nbRandNumbersPerThread*optixLaunchParams.simConstants.size.x*optixLaunchParams.simConstants.size.y){
            printf("randInd %u is out of bounds\n", randInd + randOffset);
        }
#endif*/

    float theta = 0.0f;
    if(poly.desProps.desorbType == 2 || poly.desProps.desorbType == 0) {
        theta = acos(sqrt(randFloat[(unsigned int) (randInd + randOffset++)]));
    }
    else if (poly.desProps.desorbType == 3) {
        theta = acos(pow(randFloat[(unsigned int) (randInd + randOffset++)], (1.0 / (poly.desProps.cosineExponent + 1.0))));
    }
    else if (poly.desProps.desorbType == 1) {
        theta = acos(randFloat[(unsigned int) (randInd + randOffset++)]);
    }
    else{
        printf("Unsupported desorption type! %u on %d\n", poly.desProps.desorbType, poly.parentIndex);
        return;
    }
    const float phi = randFloat[(unsigned int)(randInd + randOffset++)] * 2.0 * CUDART_PI_F;

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

//TODO: Only cosine for now
static __forceinline__ __device__
float3 getNewReverseDirection(flowgpu::MolPRD& hitData, const flowgpu::Polygon& poly,
                       const float* randFloat, unsigned int& randInd, unsigned int& randOffset)
{
    // generate ray direction
    float theta = 0.0f;
    if(poly.desProps.desorbType == 2 || poly.desProps.desorbType == 0) {
        theta = acosf(sqrtf((double) randFloat[(unsigned int) (randInd + randOffset++)]));
    }
    else if (poly.desProps.desorbType == 3) {
        theta = acosf(powf((double) randFloat[(unsigned int) (randInd + randOffset++)], (1.0f / (poly.desProps.cosineExponent + 1.0f))));
    }
    else if (poly.desProps.desorbType == 1) {
        theta = acosf((double) randFloat[(unsigned int) (randInd + randOffset++)]);
    }
    else{
        printf("Unsupported desorption type! %u on %d\n", poly.desProps.desorbType, poly.parentIndex);
        return make_float3(0.0f);
    }
    const float phi = randFloat[(unsigned int)(randInd + randOffset++)] * 2.0f * CUDART_PI_F;

    const float u = sinf(theta)*cosf(phi);
    const float v = sinf(theta)*sinf(phi);
    const float n = cosf(theta);
    const float3 nU = poly.nU;
    const float3 nV = poly.nV;
    const float3 N = poly.N;

    return u*nU + v*nV - n*N;
}

static __forceinline__ __device__
float3 getNewReverseDirection(flowgpu::MolPRD& hitData, const flowgpu::Polygon& poly,
                       const double* randFloat, unsigned int& randInd, unsigned int& randOffset)
{

    // generate ray direction
    float theta = 0.0f;
    if(poly.desProps.desorbType == 2 || poly.desProps.desorbType == 0) {
        theta = acos(sqrt(randFloat[(unsigned int) (randInd + randOffset++)]));
    }
    else if (poly.desProps.desorbType == 3) {
        theta = acos(pow(randFloat[(unsigned int) (randInd + randOffset++)], (1.0 / (poly.desProps.cosineExponent + 1.0))));
    }
    else if (poly.desProps.desorbType == 1) {
        theta = acos(randFloat[(unsigned int) (randInd + randOffset++)]);
    }
    else{
        printf("Unsupported desorption type! %u on %d\n", poly.desProps.desorbType, poly.parentIndex);
        return make_float3(0.0f);
    }
    const float phi = randFloat[(unsigned int)(randInd + randOffset++)] * 2.0 * CUDART_PI_F;


    const float u = sinf(theta)*cosf(phi);
    const float v = sinf(theta)*sinf(phi);
    const float n = cosf(theta);

    const float3 nU = poly.nU;
    const float3 nV = poly.nV;
    const float3 N = poly.N; // reverse normal


    return u*nU + v*nV - n*N;
}
#endif //MOLFLOW_PROJ_COMMONFUNCTIONS_CUH
