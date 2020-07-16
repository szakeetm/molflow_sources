// Created by pbahr

#include <optix_device.h>
#include <math_constants.h>
#include "helper_math.h"
#include <cooperative_groups.h>

#include "jetbrains_indexing.h"
#include "LaunchParams.h"
#include "GPUDefines.h" // for NB_RAND
#include <LaunchParams.h>
#include "CommonFunctions.cuh"

namespace cg = cooperative_groups;

#define DET33(_11,_12,_13,_21,_22,_23,_31,_32,_33)  \
  ((_11)*( (_22)*(_33) - (_32)*(_23) ) +            \
   (_12)*( (_23)*(_31) - (_33)*(_21) ) +            \
   (_13)*( (_21)*(_32) - (_31)*(_22) ))

#define DOT(v1, v2)  \
  ((v1.x)*(v2.x) + (v1.y)*(v2.y) + (v1.z)*(v2.z))

#define CROSS(a, b)   \
  (a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x)

#define float3_as_args(u) \
    reinterpret_cast<unsigned int&>((u).x), \
    reinterpret_cast<unsigned int&>((u).y), \
    reinterpret_cast<unsigned int&>((u).z)


using namespace flowgpu;

namespace flowgpu {

    /*! launch parameters in constant memory, filled in by optix upon
        optixLaunch (this gets filled in from the buffer we pass to
        optixLaunch) */
    extern "C" __constant__ LaunchParams optixLaunchParams;

    static __forceinline__ __device__
    void *unpackPointer( unsigned int i0, unsigned int i1 )
    {
        const uint64_t uptr = static_cast<uint64_t>( i0 ) << 32 | i1;
        void*           ptr = reinterpret_cast<void*>( uptr );
        return ptr;
    }

    static __forceinline__ __device__
    void  packPointer( void* ptr, unsigned int& i0, unsigned int& i1 )
    {
        const uint64_t uptr = reinterpret_cast<uint64_t>( ptr );
        i0 = uptr >> 32;
        i1 = uptr & 0x00000000ffffffff;
    }

    template<typename T>
    static __forceinline__ __device__ T *getPRD()
    {
        const unsigned int u0 = optixGetPayload_0();
        const unsigned int u1 = optixGetPayload_1();
        return reinterpret_cast<T*>( unpackPointer( u0, u1 ) );
    }

    static __device__ __inline__ MolPRD getMolPRD()
    {
        MolPRD prd;
        prd.velocity = int_as_float( optixGetPayload_0() );
        prd.currentDepth = optixGetPayload_1();
        prd.inSystem = optixGetPayload_2();
        prd.hitFacetId = optixGetPayload_3();
        prd.hitT = int_as_float( optixGetPayload_4() );
#ifdef GPUNBOUNCE
        prd.nbBounces = optixGetPayload_5();
#endif
        prd.facetHitSide = optixGetPayload_6();
        return prd;
    }

    static __device__ __inline__ void setMolPRD( const MolPRD &prd )
    {
        optixSetPayload_0( float_as_int(prd.velocity) );
        optixSetPayload_1( prd.currentDepth );
        optixSetPayload_2( prd.inSystem );
        optixSetPayload_3( prd.hitFacetId );
        optixSetPayload_4( float_as_int(prd.hitT) );
#ifdef GPUNBOUNCE
        optixSetPayload_5( prd.nbBounces );
#endif
        optixSetPayload_6( prd.facetHitSide );

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

    //TODO: Account for oriRatio in counter increments
    // hitData.orientationRatio

    // --------------------------------------
    // increase facet counters for absorption
    // --------------------------------------
    static __forceinline__ __device__
    void increaseHitCounterAbsorp(CuFacetHitCounter& hitCounter, const float hitEquiv, const float absEquiv,  const float ortVelocity, const float velFactor, const float velocity)
    {
        //printf("-- %d -> Absorp ----- \n",prd.currentDepth);
        // 	IncreaseFacetCounter(iFacet, sHandle->currentParticle.flightTime,
        // 	1, 0, 1, 2.0 / ortVelocity, (sHandle->wp.useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity);

        atomicAdd(&hitCounter.nbMCHit,static_cast<uint32_t>(1));
        atomicAdd(&hitCounter.nbHitEquiv, hitEquiv);
        //atomicAdd(&hitCounter.nbDesorbed, static_cast<uint32_t>(0));
        atomicAdd(&hitCounter.nbAbsEquiv, absEquiv);

        atomicAdd(&hitCounter.sum_1_per_ort_velocity, (1.0f / ortVelocity));//prd.oriRatio * sum_1_per_v;
        atomicAdd(&hitCounter.sum_v_ort, velFactor*ortVelocity);//(sHandle->wp.useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity; //prd.oriRatio * sum_v_ort;
        atomicAdd(&hitCounter.sum_1_per_velocity, (hitEquiv) / velocity);//(hitEquiv + static_cast<double>(desorb)) / prd.velocity;
    }

    // --------------------------------------
    // increase facet counters for bounce
    // --------------------------------------
    static __forceinline__ __device__
    void increaseHitCounterBounce(CuFacetHitCounter &hitCounter, const float hitEquiv, const float ortVelocity, const float velFactor, const float velocity)
    {
        // 	IncreaseFacetCounter(iFacet, sHandle->currentParticle.flightTime, 1, 0, 0, 1.0 / ortVelocity, (sHandle->wp.useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity);

        const float oriRatio = 1.0f;

        atomicAdd(&hitCounter.nbMCHit, static_cast<uint32_t>(1));
        atomicAdd(&hitCounter.nbHitEquiv, hitEquiv);

        atomicAdd(&hitCounter.sum_1_per_ort_velocity, oriRatio * (1.0f / ortVelocity));//prd.oriRatio * sum_1_per_v;
        atomicAdd(&hitCounter.sum_v_ort, oriRatio *velFactor*ortVelocity);//(sHandle->wp.useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity; //prd.oriRatio * sum_v_ort;
        atomicAdd(&hitCounter.sum_1_per_velocity, (hitEquiv) / velocity);//(hitEquiv + static_cast<double>(desorb)) / prd.velocity;
    }

    // Same but without increased hit counter
    static __forceinline__ __device__
    void increaseHitCounterPostBounce(CuFacetHitCounter &hitCounter, const float hitEquiv, const float ortVelocity, const float velFactor, const float velocity)
    {
        const float oriRatio = 1.0f;
        // 		IncreaseFacetCounter(iFacet, sHandle->currentParticle.flightTime, 0, 0, 0, 1.0 / ortVelocity, (sHandle->wp.useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity);
        atomicAdd(&hitCounter.sum_1_per_ort_velocity, oriRatio * (1.0f / ortVelocity));//prd.oriRatio * sum_1_per_v;
        atomicAdd(&hitCounter.sum_v_ort, oriRatio * velFactor*ortVelocity);//(sHandle->wp.useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity; //prd.oriRatio * sum_v_ort;
        atomicAdd(&hitCounter.sum_1_per_velocity, (hitEquiv) / velocity);//(hitEquiv + static_cast<double>(desorb)) / prd.velocity;
    }

    // --------------------------------------
    // increase facet counters for bounce
    // --------------------------------------
    static __forceinline__ __device__
    void increaseTransparentCounter(CuFacetHitCounter &hitCounter, const float hitEquiv, const float ortVelocity, const float velFactor, const float velocity)
    {
        /*IncreaseFacetCounter(facet, currentParticle.flightTime +
                                    facet->colDist / 100.0 / currentParticle.velocity,
                                    1, 0, 0,
                             2.0 / (currentParticle.velocity * directionFactor),
                             2.0 * (wp.useMaxwellDistribution ? 1.0 : 1.1781) * currentParticle.velocity *
                             directionFactor);*/
        const float oriRatio = 1.0f;

        atomicAdd(&hitCounter.nbMCHit, static_cast<uint32_t>(1));
        atomicAdd(&hitCounter.nbHitEquiv, hitEquiv);

        atomicAdd(&hitCounter.sum_1_per_ort_velocity, oriRatio * (2.0f / ortVelocity));//prd.oriRatio * sum_1_per_v;
        atomicAdd(&hitCounter.sum_v_ort, oriRatio * 2.0f * velFactor * ortVelocity);//(sHandle->wp.useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity; //prd.oriRatio * sum_v_ort;
        atomicAdd(&hitCounter.sum_1_per_velocity, (hitEquiv) / velocity);//(hitEquiv + static_cast<double>(desorb)) / prd.velocity;
    }

#define EPS32 1e-6f
    // --------------------------------------
    // increase texture counters for absorption (same as desorp)
    // --------------------------------------
    static __forceinline__ __device__
    void RecordAbsorpTexture(const flowgpu::Polygon &poly, MolPRD &hitData, float3 rayOrigin, float3 rayDir){

        const float velocity_factor = 2.0f;
        const float ortSpeedFactor = 1.0f;

        const float3 b = rayOrigin - poly.O;

        float det = poly.U.x * poly.V.y - poly.U.y * poly.V.x; // TODO: Pre calculate
        float detU = b.x * poly.V.y - b.y * poly.V.x;
        float detV = poly.U.x * b.y - poly.U.y * b.x;

        if(fabsf(det)<=EPS32){
            det = poly.U.y * poly.V.z - poly.U.z * poly.V.y; // TODO: Pre calculate
            detU = b.y * poly.V.z - b.z * poly.V.y;
            detV = poly.U.y * b.z - poly.U.z * b.y;
            if(fabsf(det)<=EPS32){
                det = poly.U.z * poly.V.x - poly.U.x * poly.V.z; // TODO: Pre calculate
                detU = b.z * poly.V.x - b.x * poly.V.z;
                detV = poly.U.z * b.x - poly.U.x * b.z;
                if(fabsf(det)<=EPS32){
                    printf("[TP][Abs] Dangerous determinant calculated for texture hit: %lf : %lf : %lf -> %lf : %lf\n",det,detU,detV,detU/det,detV/det);
                }
            }
        }

        float hitLocationU = detU/det;
        float hitLocationV = detV/det;

        flowgpu::FacetTexture& facetTex = optixLaunchParams.sharedData.facetTextures[poly.texProps.textureOffset];

        unsigned int tu = (unsigned int)(hitLocationU * facetTex.texWidthD);
        unsigned int tv = (unsigned int)(hitLocationV * facetTex.texHeightD);
        unsigned int add = tu + tv * (facetTex.texWidth);

        float ortVelocity = (optixLaunchParams.simConstants.useMaxwell ? 1.0f : 1.1781f) * hitData.velocity*fabsf(dot(rayDir, poly.N)); //surface-orthogonal velocity component
        //printf("Hit at %lf , %lf for tex %lf , %lf\n", hitLocationU, hitLocationV, facetTex.texWidthD, facetTex.texHeightD);
#ifdef BOUND_CHECK
        if(facetTex.texelOffset + add < 0 || facetTex.texelOffset + add >= optixLaunchParams.simConstants.nbTexel){
            printf("facetTex.texelOffset + add %u + %u>= %u is out of bounds [%u , %u] - [%d , %d] (%lf * %lf)x(%lf * %lf)\n", facetTex.texelOffset, add, optixLaunchParams.simConstants.nbTexel,tu,tv,facetTex.texWidth,facetTex.texHeight, hitLocationU, facetTex.texWidthD, hitLocationV, facetTex.texHeightD);}
#endif
        flowgpu::Texel& tex = optixLaunchParams.sharedData.texels[facetTex.texelOffset + add];
        atomicAdd(&tex.countEquiv, static_cast<uint32_t>(1));
        atomicAdd(&tex.sum_1_per_ort_velocity, 1.0 * velocity_factor / ortVelocity);
        atomicAdd(&tex.sum_v_ort_per_area, 1.0 * ortSpeedFactor * ortVelocity * optixLaunchParams.sharedData.texelInc[facetTex.texelOffset + add]); // sum ortho_velocity[m/s] / cell_area[cm2]
    }


    // --------------------------------------
    // increase texture counters for reflections
    // --------------------------------------
    static __forceinline__ __device__
    unsigned int RecordBounceTexture(const flowgpu::Polygon& poly, MolPRD& hitData, const float3& rayOrigin, const float3& rayDir){

        const float3 b = rayOrigin - poly.O;

        float det = poly.U.x * poly.V.y - poly.U.y * poly.V.x; // TODO: Pre calculate
        float detU = b.x * poly.V.y - b.y * poly.V.x;
        float detV = poly.U.x * b.y - poly.U.y * b.x;

        if(fabsf(det)<=EPS32){
            det = poly.U.y * poly.V.z - poly.U.z * poly.V.y; // TODO: Pre calculate
            detU = b.y * poly.V.z - b.z * poly.V.y;
            detV = poly.U.y * b.z - poly.U.z * b.y;
            if(fabsf(det)<=EPS32){
                det = poly.U.z * poly.V.x - poly.U.x * poly.V.z; // TODO: Pre calculate
                detU = b.z * poly.V.x - b.x * poly.V.z;
                detV = poly.U.z * b.x - poly.U.x * b.z;
                if(fabsf(det)<=EPS32){
                    printf("[TP][Hit] Dangerous determinant calculated for texture hit: %lf : %lf : %lf -> %lf : %lf\n",det,detU,detV,detU/det,detV/det);
                }
            }
        }

        float hitLocationU = detU/det;
        float hitLocationV = detV/det;

        flowgpu::FacetTexture& facetTex = optixLaunchParams.sharedData.facetTextures[poly.texProps.textureOffset];

        unsigned int tu = (unsigned int)(hitLocationU * facetTex.texWidthD);
        unsigned int tv = (unsigned int)(hitLocationV * facetTex.texHeightD);
        unsigned int add = tu + tv * (facetTex.texWidth);

        //printf("Pre Bounce Tex: %f = %f * %f * %f\n",ortVelocity,(optixLaunchParams.simConstants.useMaxwell ? 1.0f : 1.1781f) ,hitData.velocity,fabsf(dot(rayDir, poly.N)));
#ifdef BOUND_CHECK
        if(facetTex.texelOffset + add < 0 || facetTex.texelOffset + add >= optixLaunchParams.simConstants.nbTexel){
            printf("facetTex.texelOffset + add %u + %u>= %u is out of bounds [%u , %u] - [%d , %d] (%lf * %lf)x(%lf * %lf)\n", facetTex.texelOffset, add, optixLaunchParams.simConstants.nbTexel,tu,tv,facetTex.texWidth,facetTex.texHeight, hitLocationU, facetTex.texWidthD, hitLocationV, facetTex.texHeightD);}
#endif
        flowgpu::Texel& tex = optixLaunchParams.sharedData.texels[facetTex.texelOffset + add];

        const float velocity_factor = 1.0f;
        const float ortSpeedFactor = 1.0f;
        const float oriRatio = 1.0f; // TODO: Part of particle
        const float ortVelocity = (optixLaunchParams.simConstants.useMaxwell ? 1.0f : 1.1781f) * hitData.velocity * fabsf(dot(rayDir, poly.N)); //surface-orthogonal velocity component

        atomicAdd(&tex.countEquiv, static_cast<uint32_t>(oriRatio));
        atomicAdd(&tex.sum_1_per_ort_velocity, (float)(oriRatio * velocity_factor / ortVelocity));
        atomicAdd(&tex.sum_v_ort_per_area, (float)(oriRatio * ortSpeedFactor * ortVelocity * optixLaunchParams.sharedData.texelInc[facetTex.texelOffset + add])); // sum ortho_velocity[m/s] / cell_area[cm2]
#if DEBUG
        if(facetTex.texelOffset + add == 4){
            printf("Pre Bounce Texoffset4[%d]: %f\n", poly.parentIndex,optixLaunchParams.sharedData.texelInc[facetTex.texelOffset + add]);
        }
        if(optixLaunchParams.sharedData.texelInc[facetTex.texelOffset + add] >= 1e5)
            printf("Pre Bounce Tex[%d]: %f , %f - %f (%f) -> [%u + %u]\n", poly.parentIndex,
                oriRatio * velocity_factor / ortVelocity,oriRatio * ortSpeedFactor * ortVelocity * optixLaunchParams.sharedData.texelInc[facetTex.texelOffset + add],
                ortVelocity,
                optixLaunchParams.sharedData.texelInc[facetTex.texelOffset + add], facetTex.texelOffset, add);
#endif
        // Return Texel index for reuse in post bounce processing
        return facetTex.texelOffset + add;
    }

    // --------------------------------------
    // increase texture counters for reflections
    // --------------------------------------
    static __forceinline__ __device__
    void RecordTransparentTexture(const flowgpu::Polygon& poly, MolPRD& hitData, const float3& rayOrigin, const float3& rayDir){

        const float3 b = rayOrigin - poly.O;

        float det = poly.U.x * poly.V.y - poly.U.y * poly.V.x; // TODO: Pre calculate
        float detU = b.x * poly.V.y - b.y * poly.V.x;
        float detV = poly.U.x * b.y - poly.U.y * b.x;

        if(fabsf(det)<=EPS32/* || detU/det > 1.0 || detV/det > 1.0*/){
            //printf("[1][TP][Hit] Dangerous determinant calculated for texture hit: %lf : %lf : %lf -> %lf : %lf\n",det,detU,detV,detU/det,detV/det);

            det = poly.U.y * poly.V.z - poly.U.z * poly.V.y; // TODO: Pre calculate
            detU = b.y * poly.V.z - b.z * poly.V.y;
            detV = poly.U.y * b.z - poly.U.z * b.y;
            if(fabsf(det)<=EPS32/* || detU/det > 1.0 || detV/det > 1.0*/){
                //printf("[2][TP][Hit] Dangerous determinant calculated for texture hit: %lf : %lf : %lf -> %lf : %lf\n",det,detU,detV,detU/det,detV/det);
                det = poly.U.z * poly.V.x - poly.U.x * poly.V.z; // TODO: Pre calculate
                detU = b.z * poly.V.x - b.x * poly.V.z;
                detV = poly.U.z * b.x - poly.U.x * b.z;
                if(fabsf(det)<=EPS32/* || detU/det > 1.0 || detV/det > 1.0*/){
                    printf("[3][TP][Hit] Dangerous determinant calculated for texture hit: %lf : %lf : %lf -> %lf : %lf\n",det,detU,detV,detU/det,detV/det);
                }
            }
        }

/*#ifdef DEBUG
        printf("[%u][%u] trans tex on poly -> %lf , %lf , %lf . %lf , %lf , %lf\n", poly.parentIndex,poly.parentIndex,
                poly.U.x,poly.U.y,poly.U.z,poly.V.x,poly.V.y,poly.V.z);
#endif
#ifdef DEBUG
        printf("[%u][%u] trans tex on poly pos -> %lf , %lf , %lf\n", poly.parentIndex,poly.parentIndex,
                rayOrigin.x,rayOrigin.y,rayOrigin.z);
#endif*/
        const float hitLocationU = detU/det;
        const float hitLocationV = detV/det;

/*#ifdef DEBUG
        printf("[%u][%u] trans tex offset -> %lf , %lf / %lf = %u!!\n", poly.parentIndex,poly.parentIndex,
                detU,detV,det,poly.texProps.textureOffset);
#endif*/
        flowgpu::FacetTexture& facetTex = optixLaunchParams.sharedData.facetTextures[poly.texProps.textureOffset];
/*#ifdef DEBUG
        printf("[%u][%u] trans tex width -> %lf , %lf!!\n", poly.parentIndex,poly.parentIndex, facetTex.texWidthD,
                facetTex.texHeightD);
#endif*/
        const unsigned int tu = (unsigned int)(hitLocationU * facetTex.texWidthD);
        const unsigned int tv = (unsigned int)(hitLocationV * facetTex.texHeightD);
        const unsigned int add = tu + tv * (facetTex.texWidth);

/*#ifdef DEBUG
        printf("[%u][%u] trans tex hitloc -> %lf , %lf -> %u , %u = %u ===> %u!!\n", poly.parentIndex,poly.parentIndex,
                hitLocationU, hitLocationV,
                tu,tv,add, facetTex.texelOffset + add);
#endif*/
#ifdef BOUND_CHECK
        if(facetTex.texelOffset + add < 0 || facetTex.texelOffset + add >= optixLaunchParams.simConstants.nbTexel){
            printf("facetTex.texelOffset + add %u + %u>= %u is out of bounds [%u , %u] - [%d , %d] (%lf * %lf)x(%lf * %lf)\n", facetTex.texelOffset, add, optixLaunchParams.simConstants.nbTexel,tu,tv,facetTex.texWidth,facetTex.texHeight, hitLocationU, facetTex.texWidthD, hitLocationV, facetTex.texHeightD);}
#endif
        flowgpu::Texel& tex = optixLaunchParams.sharedData.texels[facetTex.texelOffset + add];

        const float velocity_factor = 2.0f;
        const float ortSpeedFactor = 2.0f;
        const float oriRatio = 1.0f; // TODO: Part of particle
        const float ortVelocity = (optixLaunchParams.simConstants.useMaxwell ? 1.0f : 1.1781f) * hitData.velocity * fabsf(dot(rayDir, poly.N)); //surface-orthogonal velocity component

/*#ifdef DEBUG
        printf("[%u][%u] trans tex add!!\n", poly.parentIndex,poly.parentIndex);
#endif*/
        atomicAdd(&tex.countEquiv, static_cast<uint32_t>(oriRatio));
        atomicAdd(&tex.sum_1_per_ort_velocity, (float)(oriRatio * velocity_factor / ortVelocity));
/*#ifdef DEBUG
        printf("[%u][%u] trans tex add sum v!!\n", poly.parentIndex,poly.parentIndex);
#endif*/
        atomicAdd(&tex.sum_v_ort_per_area, (float)(oriRatio * ortSpeedFactor * ortVelocity * optixLaunchParams.sharedData.texelInc[facetTex.texelOffset + add])); // sum ortho_velocity[m/s] / cell_area[cm2]
/*#ifdef DEBUG
        printf("[%u][%u] trans tex added!!\n", poly.parentIndex,poly.parentIndex);
#endif*/
        return;
    }

    // Same but without increased hit counter
    static __forceinline__ __device__
    void RecordPostBounceTexture(const flowgpu::Polygon& poly, MolPRD& hitData, float3 rayDir, unsigned int texelIndex){

        const float velocity_factor = 1.0f;
        const float ortSpeedFactor = 1.0f;
        const float oriRatio = 1.0f; // TODO: Part of particle
        const float ortVelocity = (optixLaunchParams.simConstants.useMaxwell ? 1.0f : 1.1781f) * hitData.velocity * fabsf(dot(rayDir, poly.N)); //surface-orthogonal velocity component
        //printf("PostBounce Tex: %f = %f * %f * %f\n",ortVelocity,(optixLaunchParams.simConstants.useMaxwell ? 1.0f : 1.1781f) ,hitData.velocity,fabsf(dot(rayDir, poly.N)));

        flowgpu::Texel& tex = optixLaunchParams.sharedData.texels[texelIndex];
        atomicAdd(&tex.sum_1_per_ort_velocity, (float)(oriRatio * velocity_factor / ortVelocity));
        atomicAdd(&tex.sum_v_ort_per_area, (float)(oriRatio * ortSpeedFactor * ortVelocity * optixLaunchParams.sharedData.texelInc[texelIndex])); // sum ortho_velocity[m/s] / cell_area[cm2]
        //printf("PostBounce Tex[%d]: %f , %f - %f (%f)\n", poly.parentIndex, oriRatio * velocity_factor / ortVelocity,oriRatio * ortSpeedFactor * ortVelocity * optixLaunchParams.sharedData.texelInc[texelIndex],ortVelocity,optixLaunchParams.sharedData.texelInc[texelIndex]);
#if DEBUG
        if(optixLaunchParams.sharedData.texelInc[texelIndex] >= 1e5)
        printf("PostBounce Tex[%d]: %f , %f - %f (%f)\n", poly.parentIndex,
                oriRatio * velocity_factor / ortVelocity,
                oriRatio * ortSpeedFactor * ortVelocity * optixLaunchParams.sharedData.texelInc[texelIndex],
                ortVelocity,
                optixLaunchParams.sharedData.texelInc[texelIndex]);
#endif
    }

    // --------------------------------------
    // increase texture counters for absorption (same as desorp)
    // --------------------------------------
    static __forceinline__ __device__
    void RecordAbsorpProfile(const flowgpu::Polygon &poly, MolPRD &hitData, float3 rayOrigin, float3 rayDir){

        const float velocity_factor = 2.0f;
        const float ortSpeedFactor = 1.0f;

        const float3 b = rayOrigin - poly.O;

        float det = poly.U.x * poly.V.y - poly.U.y * poly.V.x; // TODO: Pre calculate
        float detU = b.x * poly.V.y - b.y * poly.V.x;
        float detV = poly.U.x * b.y - poly.U.y * b.x;

        if(fabsf(det)<=EPS32){
            det = poly.U.y * poly.V.z - poly.U.z * poly.V.y; // TODO: Pre calculate
            detU = b.y * poly.V.z - b.z * poly.V.y;
            detV = poly.U.y * b.z - poly.U.z * b.y;
            if(fabsf(det)<=EPS32){
                det = poly.U.z * poly.V.x - poly.U.x * poly.V.z; // TODO: Pre calculate
                detU = b.z * poly.V.x - b.x * poly.V.z;
                detV = poly.U.z * b.x - poly.U.x * b.z;
                if(fabsf(det)<=EPS32){
                    printf("[TP][Abs] Dangerous determinant calculated for profile hit: %lf : %lf : %lf -> %lf : %lf\n",det,detU,detV,detU/det,detV/det);
                }
            }
        }

        unsigned int add = 0;
        if(poly.profProps.profileType == flowgpu::PROFILE_FLAGS::profileU){
            float hitLocationU = detU/det;
            add = (unsigned int)(hitLocationU * PROFILE_SIZE);
        }
        else if(poly.profProps.profileType == flowgpu::PROFILE_FLAGS::profileV){
            float hitLocationV = detV/det;
            add = (unsigned int)(hitLocationV * PROFILE_SIZE);
        }

        float ortVelocity = hitData.velocity*fabsf(dot(rayDir, poly.N)); //surface-orthogonal velocity component
        //printf("Hit at %lf , %lf for tex %lf , %lf\n", hitLocationU, hitLocationV, facetTex.texWidthD, facetTex.texHeightD);

#ifdef BOUND_CHECK
        if(poly.profProps.profileOffset + add < 0 || poly.profProps.profileOffset + add >= optixLaunchParams.simConstants.nbProfSlices){printf("poly.profProps.profileOffset + add %u >= %u is out of bounds\n", poly.profProps.profileOffset + add, optixLaunchParams.simConstants.nbProfSlices);}
#endif
        flowgpu::Texel& tex = optixLaunchParams.sharedData.profileSlices[poly.profProps.profileOffset + add];
        atomicAdd(&tex.countEquiv, static_cast<uint32_t>(1));
        atomicAdd(&tex.sum_1_per_ort_velocity, 1.0 * velocity_factor / ortVelocity);
        atomicAdd(&tex.sum_v_ort_per_area, 1.0 * ortSpeedFactor * ortVelocity * (optixLaunchParams.simConstants.useMaxwell ? 1.0f : 1.1781f)); // sum ortho_velocity[m/s] / cell_area[cm2]
    }

    // --------------------------------------
    // increase profile counters for reflections
    // --------------------------------------
    static __forceinline__ __device__
    unsigned int RecordBounceProfile(const flowgpu::Polygon& poly, MolPRD& hitData, const float3& rayOrigin, const float3& rayDir){

        const float3 b = rayOrigin - poly.O;

        float det = poly.U.x * poly.V.y - poly.U.y * poly.V.x; // TODO: Pre calculate
        float detU = b.x * poly.V.y - b.y * poly.V.x;
        float detV = poly.U.x * b.y - poly.U.y * b.x;

        if(fabsf(det)<=EPS32){
            det = poly.U.y * poly.V.z - poly.U.z * poly.V.y; // TODO: Pre calculate
            detU = b.y * poly.V.z - b.z * poly.V.y;
            detV = poly.U.y * b.z - poly.U.z * b.y;
            if(fabsf(det)<=EPS32){
                det = poly.U.z * poly.V.x - poly.U.x * poly.V.z; // TODO: Pre calculate
                detU = b.z * poly.V.x - b.x * poly.V.z;
                detV = poly.U.z * b.x - poly.U.x * b.z;
                if(fabsf(det)<=EPS32){
                    printf("[TP][Hit][->%u] Dangerous determinant calculated for profile hit: %lf : %lf : %lf -> %lf : %lf\n",poly.parentIndex,det,detU,detV,detU/det,detV/det);
                }
            }
        }

        unsigned int add = 0;
        if(poly.profProps.profileType == flowgpu::PROFILE_FLAGS::profileU){
            float hitLocationU = detU/det;
            add = (unsigned int)(hitLocationU * PROFILE_SIZE);
        }
        else if(poly.profProps.profileType == flowgpu::PROFILE_FLAGS::profileV){
            float hitLocationV = detV/det;
            add = (unsigned int)(hitLocationV * PROFILE_SIZE);
        }

        //printf("Pre Bounce Tex: %f = %f * %f * %f\n",ortVelocity,(optixLaunchParams.simConstants.useMaxwell ? 1.0f : 1.1781f) ,hitData.velocity,fabsf(dot(rayDir, poly.N)));
#ifdef BOUND_CHECK
        if(poly.profProps.profileOffset + add < 0 || poly.profProps.profileOffset + add >= optixLaunchParams.simConstants.nbProfSlices){printf("poly.profProps.profileOffset + add %u >= %u is out of bounds\n", poly.profProps.profileOffset + add, optixLaunchParams.simConstants.nbProfSlices);}
#endif
        flowgpu::Texel& tex = optixLaunchParams.sharedData.profileSlices[poly.profProps.profileOffset + add];

        const float velocity_factor = 1.0f;
        const float ortSpeedFactor = 1.0f;
        const float oriRatio = 1.0f; // TODO: Part of particle
        const float ortVelocity = hitData.velocity * fabsf(dot(rayDir, poly.N)); //surface-orthogonal velocity component

        atomicAdd(&tex.countEquiv, static_cast<uint32_t>(oriRatio));
        atomicAdd(&tex.sum_1_per_ort_velocity, (float)(oriRatio * velocity_factor / ortVelocity));
        atomicAdd(&tex.sum_v_ort_per_area, (float)(oriRatio * ortSpeedFactor * ortVelocity * (optixLaunchParams.simConstants.useMaxwell ? 1.0f : 1.1781f))); // sum ortho_velocity[m/s] / cell_area[cm2]
        //printf("Pre Bounce Tex[%d]: %f , %f - %f (%f)\n", poly.parentIndex, oriRatio * velocity_factor / ortVelocity,oriRatio * ortSpeedFactor * ortVelocity * optixLaunchParams.sharedData.texelInc[facetTex.texelOffset + add],ortVelocity,optixLaunchParams.sharedData.texelInc[facetTex.texelOffset + add]);

        // Return Texel index for reuse in post bounce processing
        return poly.profProps.profileOffset + add;
    }

    // Same but without increased hit counter
    static __forceinline__ __device__
    void RecordPostBounceProfile(const flowgpu::Polygon& poly, MolPRD& hitData, float3 rayDir, unsigned int texelIndex){

        const float velocity_factor = 1.0f;
        const float ortSpeedFactor = 1.0f;
        const float oriRatio = 1.0f; // TODO: Part of particle
        const float ortVelocity =  hitData.velocity * fabsf(dot(rayDir, poly.N)); //surface-orthogonal velocity component
        //printf("PostBounce Tex: %f = %f * %f * %f\n",ortVelocity,(optixLaunchParams.simConstants.useMaxwell ? 1.0f : 1.1781f) ,hitData.velocity,fabsf(dot(rayDir, poly.N)));

        flowgpu::Texel& tex = optixLaunchParams.sharedData.profileSlices[texelIndex];
        atomicAdd(&tex.sum_1_per_ort_velocity, (float)(oriRatio * velocity_factor / ortVelocity));
        atomicAdd(&tex.sum_v_ort_per_area, (float)(oriRatio * ortSpeedFactor * ortVelocity * (optixLaunchParams.simConstants.useMaxwell ? 1.0f : 1.1781f))); // sum ortho_velocity[m/s] / cell_area[cm2]
        //printf("PostBounce Tex[%d]: %f , %f - %f (%f)\n", poly.parentIndex, oriRatio * velocity_factor / ortVelocity,oriRatio * ortSpeedFactor * ortVelocity * optixLaunchParams.sharedData.texelInc[texelIndex],ortVelocity,optixLaunchParams.sharedData.texelInc[texelIndex]);
    }

    // --------------------------------------
    // increase profile counters for reflections
    // --------------------------------------
    static __forceinline__ __device__
    unsigned int RecordTransparentProfile(const flowgpu::Polygon& poly, MolPRD& hitData, const float3& rayOrigin, const float3& rayDir){

        const float3 b = rayOrigin - poly.O;

        float det = poly.U.x * poly.V.y - poly.U.y * poly.V.x; // TODO: Pre calculate
        float detU = b.x * poly.V.y - b.y * poly.V.x;
        float detV = poly.U.x * b.y - poly.U.y * b.x;

        if(fabsf(det)<=EPS32){
            det = poly.U.y * poly.V.z - poly.U.z * poly.V.y; // TODO: Pre calculate
            detU = b.y * poly.V.z - b.z * poly.V.y;
            detV = poly.U.y * b.z - poly.U.z * b.y;
            if(fabsf(det)<=EPS32){
                det = poly.U.z * poly.V.x - poly.U.x * poly.V.z; // TODO: Pre calculate
                detU = b.z * poly.V.x - b.x * poly.V.z;
                detV = poly.U.z * b.x - poly.U.x * b.z;
                if(fabsf(det)<=EPS32){
                    printf("[TP][Hit][->%u] Dangerous determinant calculated for profile hit: %lf : %lf : %lf -> %lf : %lf\n",poly.parentIndex,det,detU,detV,detU/det,detV/det);
                }
            }
        }

        unsigned int add = 0;
        if(poly.profProps.profileType == flowgpu::PROFILE_FLAGS::profileU){
            float hitLocationU = detU/det;
            add = (unsigned int)(hitLocationU * PROFILE_SIZE);
        }
        else if(poly.profProps.profileType == flowgpu::PROFILE_FLAGS::profileV){
            float hitLocationV = detV/det;
            add = (unsigned int)(hitLocationV * PROFILE_SIZE);
        }

        //printf("Pre Bounce Tex: %f = %f * %f * %f\n",ortVelocity,(optixLaunchParams.simConstants.useMaxwell ? 1.0f : 1.1781f) ,hitData.velocity,fabsf(dot(rayDir, poly.N)));
#ifdef BOUND_CHECK
        if(poly.profProps.profileOffset + add < 0 || poly.profProps.profileOffset + add >= optixLaunchParams.simConstants.nbProfSlices){printf("poly.profProps.profileOffset + add %u >= %u is out of bounds\n", poly.profProps.profileOffset + add, optixLaunchParams.simConstants.nbProfSlices);}
#endif
        flowgpu::Texel& tex = optixLaunchParams.sharedData.profileSlices[poly.profProps.profileOffset + add];

        const float velocity_factor = 2.0f;
        const float ortSpeedFactor = 2.0f;
        const float oriRatio = 1.0f; // TODO: Part of particle
        const float ortVelocity = hitData.velocity * fabsf(dot(rayDir, poly.N)); //surface-orthogonal velocity component

        atomicAdd(&tex.countEquiv, static_cast<uint32_t>(oriRatio));
        atomicAdd(&tex.sum_1_per_ort_velocity, (float)(oriRatio * velocity_factor / ortVelocity));
        atomicAdd(&tex.sum_v_ort_per_area, (float)(oriRatio * ortSpeedFactor * ortVelocity * (optixLaunchParams.simConstants.useMaxwell ? 1.0f : 1.1781f))); // sum ortho_velocity[m/s] / cell_area[cm2]
        //printf("Pre Bounce Tex[%d]: %f , %f - %f (%f)\n", poly.parentIndex, oriRatio * velocity_factor / ortVelocity,oriRatio * ortSpeedFactor * ortVelocity * optixLaunchParams.sharedData.texelInc[facetTex.texelOffset + add],ortVelocity,optixLaunchParams.sharedData.texelInc[facetTex.texelOffset + add]);

        // Return Texel index for reuse in post bounce processing
        return poly.profProps.profileOffset + add;
    }

    //------------------------------------------------------------------------------
    //
    // closest hit and miss programs for the molflow ray tracer standard molecules
    //
    // Note eventually we will have to create one pair of those for each
    // ray type and each geometry type we want to render; but this
    // simple example doesn't use any actual geometries yet, so we only
    // create a single, dummy, set of them (we do have to have at least
    // one group of them to set up the SBT)
    // ------------------------------------------------------------------------------

    extern "C" __global__ void __anyhit__molecule()
    {
        // we don't need it
    }

    extern "C" __global__ void __closesthit__molecule_triangle()
    {

        const TriangleMeshSBTData &sbtData = *(const TriangleMeshSBTData*)optixGetSbtDataPointer();

        //const int   primID = optixGetPrimitiveIndex();
        const float3 ray_orig = optixGetWorldRayOrigin();
        const float3 ray_dir  = optixGetWorldRayDirection();
        const float  ray_t    = optixGetRayTmax();

        const uint2 launchIndex = make_uint2(optixGetLaunchIndex());
        const unsigned int bufferIndex = launchIndex.x + launchIndex.y * optixLaunchParams.simConstants.size.x;

#ifdef DEBUGMISS
        const unsigned int missIndex = bufferIndex*NMISSES;
        optixLaunchParams.perThreadData.missBuffer[missIndex] = 0;
#endif
        const int facetHitKind = optixGetHitKind();

        /*if(fbIndex==0)
            printf("--- start post processing ---\n");*/

        MolPRD prd = getMolPRD();
        //printf("[%d] -- %d Bounces loaded----- \n",ix, prd.currentDepth);

/*#ifdef DEBUG
        printf("[%u] normal hit detected -> %u\n", bufferIndex, prd.inSystem);
#endif*/

        // self intersection
        if(prd.hitFacetId == optixGetPrimitiveIndex()){
#ifdef DEBUG
            //printf("[%d] source and goal facet equal %d : %8.6f,%8.6f,%8.6f -> %8.6f,%8.6f,%8.6f : [%.5e .. %.5e] (tri)\n",(blockDim.x * blockIdx.x + threadIdx.x), prd.hitFacetId, ray_orig.x,ray_orig.y,ray_orig.z,ray_dir.x,ray_dir.y,ray_dir.z, optixGetRayTmin(),ray_t);
/*
if(prd.inSystem == 4)
                printf("[%d] [backface] source and goal facet equal %d : %8.6f,%8.6f,%8.6f -> %8.6f,%8.6f,%8.6f : [%.5e .. %.5e] (tri)\n",(blockDim.x * blockIdx.x + threadIdx.x), prd.hitFacetId, ray_orig.x,ray_orig.y,ray_orig.z,ray_dir.x,ray_dir.y,ray_dir.z, optixGetRayTmin(),ray_t);
*/

#endif
            prd.facetHitSide = facetHitKind;
            if(facetHitKind == OPTIX_HIT_KIND_TRIANGLE_BACK_FACE) {
                //prd.inSystem = SELF_BACK_HIT;
                /*if(bufferIndex==0)
                    printf("[%u] back self hit detected -> %u -> %u for n=%d\n", bufferIndex, prd.inSystem, sbtData.poly[prd.hitFacetId].parentIndex,prd.nbBounces);
            */}

            //prd.hitPos = ray_orig;
            optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].hitPos = ray_orig;
            //optixLaunchParams.perThreadData.currentMoleculeData[fbIndex].postHitDir = ray_dir;

            //prd.currentDepth = 0;
            setMolPRD(prd);
            return;
        }


        prd.hitT = ray_t;
        prd.hitPos = ray_orig + ray_t * ray_dir;
        prd.hitFacetId = optixGetPrimitiveIndex();
        prd.facetHitSide = facetHitKind;

        // first add facet hits
        const unsigned int counterIdx = prd.hitFacetId + (bufferIndex%(CORESPERSM * WARPSCHEDULERS)) * optixLaunchParams.simConstants.nbFacets;

        //Register (orthogonal) velocity
        const flowgpu::Polygon& poly  = sbtData.poly[prd.hitFacetId];

        // replace like with self intersection, but keep position etc.
        if(facetHitKind == OPTIX_HIT_KIND_TRIANGLE_FRONT_FACE) {
            //prd.inSystem = TRANSPARENT_HIT;
        }
        else
        if(bufferIndex==0)
            printf("[%u] back hit detected -> Mol: %u on F#%u\n", bufferIndex, prd.inSystem, poly.parentIndex);

        // TODO: Consider maxwelldistribution or not and oriRatio for lowFlux and desorbtion/absorbtion
        //IncreaseFacetCounter(
        // iFacet, sHandle->currentParticle.flightTime, 1, 0, 0, 1.0 / ortVelocity,
        // (sHandle->wp.useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity);

        // TODO: Save somewhere as a shared constant instead of repetively evaluating
        float velFactor = optixLaunchParams.simConstants.useMaxwell ? 1.0f : 1.1781f;

        // TODO: if bounce
        // TODO: queue bounces back in pipeline instead of recursively following them
        // stop for some point due to recursion limits

        float ortVelocity = prd.velocity*fabsf(dot(ray_dir, poly.N));
        const float hitEquiv = 1.0f; //1.0*prd.orientationRatio; // hit=1.0 (only changed for lowflux mode)
        const float absEquiv = 1.0f; //1.0*prd.orientationRatio; // hit=1.0 (only changed for lowflux mode)

#ifndef WITH_TRANS
        if(poly.facProps.is2sided && poly.facProps.opacity == 0){
        //if(poly.facProps.is2sided){
#ifdef DEBUG
            printf("[%d] 2sided facet hit %d / %d : %8.6f,%8.6f,%8.6f -> %8.6f,%8.6f,%8.6f : [%.5e .. %.5e] (tri)\n",(blockDim.x * blockIdx.x + threadIdx.x), prd.hitFacetId, poly.parentIndex, ray_orig.x,ray_orig.y,ray_orig.z,ray_dir.x,ray_dir.y,ray_dir.z, optixGetRayTmin(),ray_t);
#endif
            // replace like with self intersection, but keep position etc.
            prd.inSystem = TRANSPARENT_HIT;

            optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].hitPos = prd.hitPos;
            optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].postHitDir = ray_dir;

            increaseTransparentCounter(optixLaunchParams.hitCounter[counterIdx], hitEquiv, ortVelocity, velFactor, prd.velocity);
            if (poly.texProps.textureFlags & flowgpu::TEXTURE_FLAGS::countTrans) {
                RecordTransparentTexture(poly, prd, prd.hitPos, ray_dir);
            }
            setMolPRD(prd);
            return;
        }
#endif // WITH_TRANS

        const RN_T* randFloat = optixLaunchParams.randomNumbers;
        unsigned int randInd = NB_RAND*(bufferIndex);
        unsigned int randOffset = optixLaunchParams.perThreadData.randBufferOffset[bufferIndex];


        //TODO: AtomicAdd for smaller hitCounter structure, just per threadIdx
        if(poly.facProps.stickingFactor>=0.99999f || ((poly.facProps.stickingFactor > 0.0f) && (randFloat[(unsigned int)(randInd + randOffset++)] < (poly.facProps.stickingFactor)))){
            //--------------------------
            //-------- ABSORPTION ------
            //--------------------------
#ifdef BOUND_CHECK
            if(counterIdx < 0 || counterIdx >= optixLaunchParams.simConstants.nbFacets * CORESPERSM * WARPSCHEDULERS){printf("facIndex %u >= %u is out of bounds\n", counterIdx, optixLaunchParams.simConstants.nbFacets * CORESPERSM * WARPSCHEDULERS);}
#endif
            increaseHitCounterAbsorp(optixLaunchParams.hitCounter[counterIdx], hitEquiv, absEquiv, ortVelocity, velFactor, prd.velocity);

#ifdef WITH_TEX
            if (poly.texProps.textureFlags & flowgpu::TEXTURE_FLAGS::countAbs){
                //printf("[%u] Absorb texture hit: %lf , %lf , %lf -> %lf , %lf , %lf => %lf , %lf , %lf\n",poly.parentIndex,ray_orig.x,ray_orig.y,ray_orig.z,ray_dir.x,ray_dir.y,ray_dir.z,ray_orig.x+ray_dir.x,ray_orig.y+ray_dir.y,ray_orig.z+ray_dir.z);
                RecordAbsorpTexture(poly, prd, prd.hitPos, ray_dir);
            }
#endif
#ifdef WITH_PROF
            if (poly.profProps.profileType != flowgpu::PROFILE_FLAGS::noProfile){
                RecordAbsorpProfile(poly, prd, prd.hitPos, ray_dir);
            }
#endif
            optixLaunchParams.perThreadData.randBufferOffset[bufferIndex] = randOffset;

            prd.inSystem = NEW_PARTICLE;
            prd.currentDepth = 0;
#ifdef GPUNBOUNCE
            /*if(poly.parentIndex == 1)
                printf("[%u][%u == %u][%u] Absorb hit: %lf , %lf , %lf -> %lf , %lf , %lf => %lf , %lf , %lf\n",
                        bufferIndex,poly.parentIndex,prd.hitFacetId,prd.nbBounces,ray_orig.x,ray_orig.y,ray_orig.z,ray_dir.x,ray_dir.y,ray_dir.z,prd.hitPos.x,prd.hitPos.y,prd.hitPos.z);
*/
            optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].nbBounces = prd.nbBounces = 0;
#endif
            setMolPRD(prd);
            return;
        }
        // Process a memory error
        else if(poly.facProps.stickingFactor<-0.00001f){
            //

            printf("-- %lf Err ----- \n",poly.facProps.stickingFactor);
            return;
        }
        // Process a bounce/reflection
        else{

            //--------------------------
            //-------- REFLECTION ------
            //--------------------------
#ifdef BOUND_CHECK
            if(counterIdx < 0 || counterIdx >= optixLaunchParams.simConstants.nbFacets * CORESPERSM * WARPSCHEDULERS){printf("facIndex %u >= %u is out of bounds\n", counterIdx, optixLaunchParams.simConstants.nbFacets * CORESPERSM * WARPSCHEDULERS);}
#endif
            // 1. Increment counters
            increaseHitCounterBounce(optixLaunchParams.hitCounter[counterIdx], hitEquiv, ortVelocity, velFactor, prd.velocity);

#ifdef WITH_TEX
            unsigned int texelIndex = 1e8;
            if (poly.texProps.textureFlags & flowgpu::TEXTURE_FLAGS::countRefl){
                texelIndex = RecordBounceTexture(poly, prd, prd.hitPos, ray_dir);
            }
#endif // WITH_TEX
#ifdef WITH_PROF
            unsigned int profileIndex = 1e8;
            //printf("BOUNCE? %d != %d == %d\n",(int)poly.profProps.profileType, (int)flowgpu::PROFILE_FLAGS::noProfile,(int)poly.profProps.profileType != flowgpu::PROFILE_FLAGS::noProfile);
            if (poly.profProps.profileType != flowgpu::PROFILE_FLAGS::noProfile){
                profileIndex = RecordBounceProfile(poly, prd, prd.hitPos, ray_dir);
            }
#endif


            // also increase a nbBounce counter if there is need (e.g. recursion)
            prd.inSystem = ACTIVE_PARTICLE;

            // 2. Relaunch particle

            //-----------
            // new ray velocity (for now only perfect thermalization (needs accomodationFactor)
            //-----------
            //TODO: Calculate new Velocity from CFD or directly
            //if (optixLaunchParams.simConstants.useMaxwell) prd.velocity = GenerateRandomVelocity(collidedFacet->sh.CDFid);
            //else
            prd.velocity = getNewVelocity(poly, optixLaunchParams.simConstants.gasMass);
            if(facetHitKind == OPTIX_HIT_KIND_TRIANGLE_FRONT_FACE || prd.inSystem != TRANSPARENT_HIT) {
                //if (bufferIndex == 0 && poly.parentIndex == 102) printf("[%d] Front direction!\n", poly.parentIndex);
                prd.postHitDir = getNewDirection(prd, poly, randFloat, randInd, randOffset);
            }
            else {
                if (bufferIndex == 0)
                    printf("[%d] Back direction!\n", poly.parentIndex);
                //else if (bufferIndex == 0 && poly.parentIndex != 102) printf("[%d] Back direction!\n", poly.parentIndex);
                prd.postHitDir = getNewReverseDirection(prd, poly, randFloat, randInd, randOffset);
                //prd.inSystem = ACTIVE_PARTIC;
            }

/*#if defined(DEBUG)
        if(bufferIndex == optixLaunchParams.simConstants.size.x - 1)
                        printf("[%d] hit pos %lf , %lf , %lf -> %lf , %lf , %lf -> %d\n",bufferIndex,prd.hitPos.x,prd.hitPos.y,prd.hitPos.z,prd.postHitDir.x,prd.postHitDir.y,prd.postHitDir.z,poly.parentIndex);
//printf("[%d] hit pos %lf , %lf , %lf -> %lf , %lf , %lf -> %d\n",bufferIndex,prd.hitPos.x,prd.hitPos.y,prd.hitPos.z,prd.hitPos.x+prd.postHitDir.x,prd.hitPos.y+prd.postHitDir.y,prd.hitPos.z+prd.postHitDir.z,poly.parentIndex);
#endif*/

            // 3. Increment counters for post bounce / outgoing particles
            ortVelocity = prd.velocity*fabsf(dot(prd.postHitDir, poly.N));
            increaseHitCounterPostBounce(optixLaunchParams.hitCounter[counterIdx], hitEquiv, ortVelocity, velFactor, prd.velocity);

#ifdef WITH_TEX
            if (poly.texProps.textureFlags & flowgpu::TEXTURE_FLAGS::countRefl)
                RecordPostBounceTexture(poly, prd, prd.postHitDir, texelIndex);
#endif
#ifdef WITH_PROF
            if (poly.profProps.profileType != flowgpu::PROFILE_FLAGS::noProfile){
                RecordPostBounceProfile(poly, prd, prd.postHitDir, profileIndex);
            }
#endif

#if defined(GPUNBOUNCE)
           /* #ifdef DEBUG
        if(bufferIndex == optixLaunchParams.simConstants.size.x - 1)
            printf("[%d] reflection and bounces %d / %d\n",bufferIndex,optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].nbBounces,prd.nbBounces);
#endif*/
                optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].nbBounces++;
                prd.nbBounces++;
            //if(bufferIndex == optixLaunchParams.simConstants.size.x - 1) printf("[%d] hit pos %lf , %lf , %lf -> %d [%d]\n",bufferIndex,prd.hitPos.x,prd.hitPos.y,prd.hitPos.z,poly.parentIndex,texelIndex);
#endif

            optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].facetHitSide = prd.facetHitSide;

            optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].hitPos = prd.hitPos;
            optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].postHitDir = prd.postHitDir;

            // Write temporary local variables back to shared memory
            optixLaunchParams.perThreadData.randBufferOffset[bufferIndex] = randOffset;

            // No recursion or end of recursion
            if(prd.currentDepth >= optixLaunchParams.simConstants.maxDepth){
                //printf("-- [%d] Max Bounces reached on depth %d, resetting for next call ... ----- \n",fbIndex, prd.currentDepth);
                prd.currentDepth = 0;
            }
            // WITH RECURSION
            else {
                ++prd.currentDepth;

#ifdef DEBUGPOS
    if(bufferIndex==0){
        const unsigned int posIndexOff = optixLaunchParams.perThreadData.posOffsetBuffer_debug[(unsigned int)(bufferIndex)]++;
        if(posIndexOff<NBPOSCOUNTS){
            const unsigned int posIndex = bufferIndex*NBPOSCOUNTS+posIndexOff;
            //printf("[%d] my pos is %d\n", (unsigned int)(ix+iy*optixLaunchParams.simConstants.size.x), posIndex);
            optixLaunchParams.perThreadData.positionsBuffer_debug[posIndex] = prd.hitPos;
        }
    }
#endif
                optixTrace(optixLaunchParams.traversable,
                           prd.hitPos,
                           prd.postHitDir,
                           0.f,//1e-4f,//0.f,    // tmin
                           1e20f,  // tmax
                           0.0f,   // rayTime
                           OptixVisibilityMask(255),
                           OPTIX_RAY_FLAG_DISABLE_ANYHIT
                           /*| OPTIX_RAY_FLAG_CULL_BACK_FACING_TRIANGLES*/,//OPTIX_RAY_FLAG_NONE,
                           RayType::RAY_TYPE_MOLECULE,             // SBT offset
                           RayType::RAY_TYPE_COUNT,               // SBT stride
                           RayType::RAY_TYPE_MOLECULE,             // missSBTIndex
                        //u0, u1 , u2, u3);
                        //float3_as_args(hitData.hitPos),
                           reinterpret_cast<unsigned int &>(prd.velocity),
                           reinterpret_cast<unsigned int &>(prd.currentDepth),
                           reinterpret_cast<unsigned int &>(prd.inSystem),
                        /* Can't use float_as_int() because it returns rvalue but payload requires a lvalue */
                           reinterpret_cast<unsigned int &>(prd.hitFacetId),
                           reinterpret_cast<unsigned int &>(prd.hitT),
                           reinterpret_cast<unsigned int &>(prd.nbBounces),
                           reinterpret_cast<unsigned int &>(prd.facetHitSide)

                );
            }

            setMolPRD(prd);
        } // end bounce
    }

#ifdef WITH_TRANS
    extern "C" __global__ void __closesthit__transparent_triangle()
    {

        const TriangleMeshSBTData &sbtData = *(const TriangleMeshSBTData*)optixGetSbtDataPointer();

        //const int   primID = optixGetPrimitiveIndex();
        const float3 ray_orig = optixGetWorldRayOrigin();
        const float3 ray_dir  = optixGetWorldRayDirection();
        const float  ray_t    = optixGetRayTmax();

        const uint2 launchIndex = make_uint2(optixGetLaunchIndex());
        const unsigned int bufferIndex = launchIndex.x + launchIndex.y * optixLaunchParams.simConstants.size.x;

        /*#ifdef DEBUG
            printf("[%u] transparent hit detected\n", bufferIndex);
        #endif*/

        MolPRD prd = getMolPRD();
        const int facetHitKind = optixGetHitKind();
        const int hitFacetId = optixGetPrimitiveIndex();

        // self intersection
        // just try to offset continously
        if(prd.hitFacetId == hitFacetId){
#ifdef DEBUG
            /*if(bufferIndex == 0)
                printf("[%u->%u] transparent hit self inter\n", bufferIndex, prd.nbBounces);
*/#endif
            prd.hitFacetId = hitFacetId;
            prd.hitT = ray_t;
            prd.postHitDir = ray_dir;
            prd.hitPos = ray_orig + ray_t * ray_dir;
            //prd.inSystem = SELF_INTERSECTION;
            prd.inSystem = TRANSPARENT_HIT;
            optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].hitPos = ray_orig;

            setMolPRD(prd);
            return;
        }

        prd.hitT = ray_t;
        prd.hitPos = ray_orig + ray_t * ray_dir;

        // first add facet hits
        const unsigned int counterIdx = hitFacetId + (bufferIndex%(CORESPERSM * WARPSCHEDULERS)) * optixLaunchParams.simConstants.nbFacets;

        //Register (orthogonal) velocity
        const flowgpu::Polygon& poly  = sbtData.poly[hitFacetId];

        /*#ifdef DEBUG
            printf("[%u] trans 2side check\n", bufferIndex);
        #endif*/

#ifdef DEBUG
        /*if(facetHitKind != OPTIX_HIT_KIND_TRIANGLE_FRONT_FACE && !poly.facProps.is2sided) {
            //__miss__molecule();

            printf("[%u][%u] Back face hit on 1-sided facet should never happen!\n",bufferIndex,poly.parentIndex);
            prd.hitFacetId = hitFacetId;
            prd.facetHitSide = facetHitKind;
            prd.inSystem = NEW_PARTICLE;
            setMolPRD(prd);
            return;
        }*/
        if(facetHitKind != OPTIX_HIT_KIND_TRIANGLE_FRONT_FACE && !poly.facProps.is2sided) {
            //__miss__molecule();

            //printf("[%u][%u] Back face hit on 1-sided facet should never happen!\n",bufferIndex,poly.parentIndex);
            /*prd.hitFacetId = hitFacetId;
            prd.facetHitSide = facetHitKind;
            prd.inSystem = NEW_PARTICLE;
            setMolPRD(prd);
            return;*/
        }
#endif
        /*#ifdef DEBUG
            printf("[%u] trans hit!!\n", bufferIndex);
        #endif*/

        const float velFactor = optixLaunchParams.simConstants.useMaxwell ? 1.0f : 1.1781f;
        float ortVelocity = prd.velocity*fabsf(dot(ray_dir, poly.N));
        const float hitEquiv = 1.0f; //1.0*prd.orientationRatio; // hit=1.0 (only changed for lowflux mode)

        // replace like with self intersection, but keep position etc.
        prd.inSystem = TRANSPARENT_HIT;

        /*if(optixGetHitKind() == OPTIX_HIT_KIND_TRIANGLE_FRONT_FACE)
            prd.inSystem = TRANSPARENT_HIT;
        else
            prd.inSystem = TRANSPARENT_BACK_HIT;*/

        /*#ifdef DEBUG
            printf("[%u] trans count!!\n", bufferIndex);
        #endif*/

        // Only increase counters if not dealing with self intersection
        if(prd.hitFacetId != hitFacetId){

            //printf("Registering trans hit on status %u\n",prd.inSystem);

            increaseTransparentCounter(optixLaunchParams.hitCounter[counterIdx], hitEquiv, ortVelocity, velFactor, prd.velocity);
#ifdef WITH_TEX
            if (poly.texProps.textureFlags & flowgpu::TEXTURE_FLAGS::countTrans) {
                RecordTransparentTexture(poly, prd, prd.hitPos, ray_dir);
            }
#endif // WITH_TEX
#ifdef WITH_PROF
            //printf("BOUNCE? %d != %d == %d\n",(int)poly.profProps.profileType, (int)flowgpu::PROFILE_FLAGS::noProfile,(int)poly.profProps.profileType != flowgpu::PROFILE_FLAGS::noProfile);
            if (poly.profProps.profileType != flowgpu::PROFILE_FLAGS::noProfile){
                RecordTransparentProfile(poly, prd, prd.hitPos, ray_dir);
            }
#endif
        }
        prd.hitFacetId = hitFacetId;
        prd.postHitDir = ray_dir;
        prd.facetHitSide = facetHitKind;

        optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].hitPos = prd.hitPos;
        optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].postHitDir = ray_dir;

        /*#ifdef DEBUG
            printf("[%u] trans end!!\n", bufferIndex);
        #endif*/
        setMolPRD(prd);
    }
#endif // WITH_TRANS

#if defined(DEBUGMISS)
    // Parallelogram intersection based on the SDK optixWhitted example
    extern "C" __device__ void is_parallelogram_miss(uint32_t primID)
    {
        const PolygonMeshSBTData &sbtData = *(const PolygonMeshSBTData*)optixGetSbtDataPointer();

        //const int   primID = optixGetPrimitiveIndex();

        float3 ray_orig = optixGetWorldRayOrigin();
        float3 ray_dir = optixGetWorldRayDirection();

        /*if(blockDim.x * blockIdx.x + threadIdx.x == 3 && ray_dir.x > -0.26 && ray_dir.x < -0.24)
            primID=162;*/

        ray_dir = make_float3(-1.0,-1.0,-1.0) * ray_dir;

        const float ray_tmin = optixGetRayTmin(), ray_tmax = optixGetRayTmax();

        const flowgpu::Polygon& poly  = sbtData.poly[primID];
        const float det = dot(poly.Nuv, ray_dir);
        printf("[%d] intersection det = %12.10f\n", primID, det);
        if(det > 0.0) {
            const float iDet = 1.0 / det;
            float3 intZ = ray_orig - poly.O;

            const float u = iDet * DET33(intZ.x, poly.V.x, ray_dir.x,
                                         intZ.y, poly.V.y, ray_dir.y,
                                         intZ.z, poly.V.z, ray_dir.z);
            printf("[%d] intersection   u = %12.10f\n",primID, u);

            const float v = iDet * DET33(poly.U.x, intZ.x, ray_dir.x,
                                         poly.U.y, intZ.y, ray_dir.y,
                                         poly.U.z, intZ.z, ray_dir.z);
            printf("[%d] intersection   v = %12.10f\n",primID, v);

            const float d = iDet * dot(poly.Nuv, intZ);
            printf("[%d] intersection   d = %12.10f\n",primID, d);

        }
    }
#endif //DEBUGMISS

    //------------------------------------------------------------------------------
    // miss program that gets called for any ray that did not have a
    // valid intersection
    //
    // reset molecule position to start from source again and count the leak
    // ------------------------------------------------------------------------------

    extern "C" __global__ void __miss__molecule()
    {
        //int3 &prd = *(int3*)getPRD<int3>();

        MolPRD prd = getMolPRD();

#if defined(DEBUGMISS)
        const float3 ray_orig = optixGetWorldRayOrigin();
        const float3 ray_dir  = optixGetWorldRayDirection();
        const float  ray_t    = optixGetRayTmax();


        const unsigned int fbIndex = blockDim.x * blockIdx.x + threadIdx.x;
        const unsigned int missIndex = fbIndex*NMISSES;
#if not defined(GPUNBOUNCE)
        printf("--------------------(%d) miss[%d -> %d][%d] (%12.10f , %12.10f , %12.10f) --> (%12.10f , %12.10f , %12.10f) = %e \n",
               prd.inSystem, fbIndex, prd.hitFacetId, missIndex,
               ray_orig.x, ray_orig.y , ray_orig.z , ray_dir.x, ray_dir.y , ray_dir.z, ray_t);
#else
        printf("--------------------(%d , %d) miss[%d -> %d][%d] (%12.10f , %12.10f , %12.10f) --> (%12.10f , %12.10f , %12.10f) = %e \n",
               prd.inSystem, prd.nbBounces, fbIndex, prd.hitFacetId, missIndex,
               ray_orig.x, ray_orig.y , ray_orig.z , ray_dir.x, ray_dir.y , ray_dir.z, ray_t);
#endif

        for(int i=missIndex+1; i <= missIndex+optixLaunchParams.perThreadData.missBuffer[missIndex];i++){
            printf("-------------------- miss[%d -> %d -> %d] at %d\n",
                    blockDim.x * blockIdx.x + threadIdx.x,
                    missIndex,
                    optixLaunchParams.perThreadData.missBuffer[missIndex],
                   optixLaunchParams.perThreadData.missBuffer[i]);

            is_parallelogram_miss(optixLaunchParams.perThreadData.missBuffer[i]);

        }


        optixLaunchParams.perThreadData.missBuffer[missIndex] = 0;
#endif //DEBUGMISS

/*optixLaunchParams.perThreadData.currentMoleculeData[fbIndex].hitPos = prd.hitPos;
optixLaunchParams.perThreadData.currentMoleculeData[fbIndex].postHitDir = prd.postHitDir;*/

#ifdef DEBUGLEAKPOS
        {
            const float3 ray_orig = optixGetWorldRayOrigin();
            const float3 ray_dir  = optixGetWorldRayDirection();
            const float  ray_t    = optixGetRayTmax();


            const unsigned int fbIndex = blockDim.x * blockIdx.x + threadIdx.x;



#if not defined(GPUNBOUNCE)
        printf("--------------------(%d) miss[%d -> %d] (%12.10f , %12.10f , %12.10f) --> (%12.10f , %12.10f , %12.10f) = %e \n",
               prd.inSystem, fbIndex, prd.hitFacetId,
               ray_orig.x, ray_orig.y , ray_orig.z , ray_dir.x, ray_dir.y , ray_dir.z, ray_t);
#else
        printf("--------------------(%d , %d) miss[%d -> %d] (%12.10f , %12.10f , %12.10f) --> (%12.10f , %12.10f , %12.10f) = %e \n",
               prd.inSystem, prd.nbBounces, fbIndex, prd.hitFacetId,
               ray_orig.x, ray_orig.y , ray_orig.z , ray_dir.x, ray_dir.y , ray_dir.z, ray_t);
#endif

        const unsigned int posIndexOff = optixLaunchParams.perThreadData.leakPosOffsetBuffer_debug[(unsigned int)(blockDim.x * blockIdx.x + threadIdx.x)]++;
            if(posIndexOff<NBCOUNTS){
                const unsigned int posIndex = (blockDim.x * blockIdx.x + threadIdx.x)*NBCOUNTS+posIndexOff;
                //printf("[%d] my pos is %d\n", (unsigned int)(fbIndex), posIndex);
                optixLaunchParams.perThreadData.leakPositionsBuffer_debug[posIndex] = ray_orig;
                optixLaunchParams.perThreadData.leakDirectionsBuffer_debug[posIndex] = ray_dir;

            }
            }
#endif

        // Try to relaunch a molecule for some time until removing it from the system
        // Differentiate between physical and single precision leaks here

        //if(prd.inSystem > 20) {
            //optixLaunchParams.sharedData.missCounter += 1;
            atomicAdd(optixLaunchParams.sharedData.missCounter, 1);

            prd.velocity = -999.0;
            prd.hitPos = make_float3(-999.0);
            prd.postHitDir = make_float3(-999.0);
            prd.hitFacetId = -1;
            prd.hitT = -999.0f;
            prd.inSystem = NEW_PARTICLE;
#if defined(GPUNBOUNCE)
        prd.nbBounces = 0;
#endif
        /*}
        else{
            prd.inSystem = max(3,prd.inSystem+1);
        }*/
        setMolPRD(prd);
    }


} // ::flowgpu
