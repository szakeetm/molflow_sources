// Created by pbahr

#include <optix_device.h>
#include <math_constants.h>

#include "LaunchParams.h"
#include "GPUDefines.h" // for NB_RAND
#include "jetbrains_indexing.h"
#include "helper_math.h"
#include <cooperative_groups.h>
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
    reinterpret_cast<uint32_t&>((u).x), \
    reinterpret_cast<uint32_t&>((u).y), \
    reinterpret_cast<uint32_t&>((u).z)

/*// Helper routines
thrust::device_vector<float> rand_vec();
//void initialize(thrust::device_vector<float>& v){
void initialize(){
    rand_vec().resize(1000);
    thrust::default_random_engine rng(123456);
    thrust::uniform_real_distribution<float> dist(0,1); // [0,1]
    for(size_t i = 0; i < rand_vec.size(); i++)
        rand_vec[i] = dist(rng);
}*/

namespace flowgpu {

    /*! launch parameters in constant memory, filled in by optix upon
        optixLaunch (this gets filled in from the buffer we pass to
        optixLaunch) */
    extern "C" __constant__ LaunchParams optixLaunchParams;

    static __forceinline__ __device__
    void  packPointer( void* ptr, uint32_t& i0, uint32_t& i1 )
    {
        const uint64_t uptr = reinterpret_cast<uint64_t>( ptr );
        i0 = uptr >> 32;
        i1 = uptr & 0x00000000ffffffff;
    }

    static __forceinline__ __device__
    void *getRayOriginPolygon( unsigned int i0, unsigned int i1 )
    {
        const uint64_t uptr = static_cast<uint64_t>( i0 ) << 32 | i1;
        void*           ptr = reinterpret_cast<void*>( uptr );
        return ptr;
    }

    static __forceinline__ __device__
    void *getRayOriginTriangle( unsigned int i0, unsigned int i1 )
    {
        const uint64_t uptr = static_cast<uint64_t>( i0 ) << 32 | i1;
        void*           ptr = reinterpret_cast<void*>( uptr );
        return ptr;
    }

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

    extern "C" __device__ int point_in_polygon(float u, float v, const flowgeom::Polygon& poly) {
        // Fast method to check if a point is inside a polygon or not.
        // Works with convex and concave polys, orientation independent
        const PolygonRayGenData* rayGenData = (PolygonRayGenData*) optixGetSbtDataPointer();
        //const int   primID = optixGetPrimitiveIndex();

        //const Polygon& poly  = sbtData.poly[primID];

        const int nbSizeMinusOne = poly.nbVertices - 1;
        const float2* polyPoints = rayGenData->vertex2;

        int n_updown = 0;
        int n_found = 0;

        float2 p;
        p.x = u;
        p.y = v;

        //printf("[%d] poly check with %d vert at poly offset %d for p %4.2f %4.2f\n",threadIdx.x,poly.nbVertices,poly.yertOffset,p.x,p.y);
        //printf("[%d] poly first vec2 %4.2f , %4.2f\n",threadIdx.x,rayGenData->vertex2[0].x,rayGenData->vertex2[0].y);
        //printf("[%d] poly second vec2 %4.2f , %4.2f\n",threadIdx.x,rayGenData->vertex2[1].x,rayGenData->vertex2[1].y);

#ifdef BOUND_CHECK
        if(poly.indexOffset < 0 || poly.indexOffset+poly.nbVertices >= optixLaunchParams.simConstants.nbVertices){
            printf("[%d] indexOffset %u -- %u >= %u is out of bounds\n", poly.indexOffset, poly.indexOffset+poly.nbVertices, optixLaunchParams.simConstants.nbVertices);
        }
#endif
        for (int j = 0; j < nbSizeMinusOne; j++) {
            const float2& p1 = polyPoints[poly.indexOffset + j];
            const float2& p2 = polyPoints[poly.indexOffset + j + 1];

            //printf("[%d / %u] p1 %4.2f %4.2f / p2 %4.2f %4.2f \n",threadIdx.x,j,p1.x,p1.y,p2.x,p2.y);

            if (p.x<p1.x != p.x<p2.x) {
                float slope = (p2.y - p1.y) / (p2.x - p1.x);
                //printf("[%d / %u] slope %4.2f --> %4.2f < %4.2f \n",threadIdx.x,j,slope, (slope * p.x - p.y), (slope * p1.x - p1.y));

                if ((slope * p.x - p.y) < (slope * p1.x - p1.y)) {
                    n_updown++;
                }
                else {
                    n_updown--;
                }
                n_found++;
            }
        }

        //Last point. Repeating code because it's the fastest and this function is heavily used
        const float2& p1 = polyPoints[poly.indexOffset + nbSizeMinusOne];
        const float2& p2 = polyPoints[poly.indexOffset + 0];
        if (p.x<p1.x != p.x<p2.x) {
            //printf("[%d / %u] p1 %4.2f %4.2f / p2 %4.2f %4.2f \n",threadIdx.x,0,p1.x,p1.y,p2.x,p2.y);

            float slope = (p2.y - p1.y) / (p2.x - p1.x);
            //printf("[%d / %u] slope %4.2f --> %4.2f < %4.2f \n",threadIdx.x,0,slope, (slope * p.x - p.y), (slope * p1.x - p1.y));

            if ((slope * p.x - p.y) < (slope * p1.x - p1.y)) {
                n_updown++;
            }
            else {
                n_updown--;
            }
            n_found++;
        }

        //printf("[%d] found %d / updown %d = %d\n",threadIdx.x,n_found,n_updown,((n_found / 2) & 1) ^ ((n_updown / 2) & 1));
        if(((n_found / 2) & 1) ^ ((n_updown / 2) & 1)){
            return 1;
        }
        else{
            return 0;
        }
    }
    static __forceinline__ __device__
    int getSourceFacet(MolPRD& hitData,
#ifdef WITHTRIANGLES
            const TriangleRayGenData* rayGenData,
#else
                         const PolygonRayGenData* rayGenData,
#endif
                         const float* randFloat, unsigned int& randInd, unsigned int& randOffset)
    {
#ifdef BOUND_CHECK
        if(randInd + randOffset < 0 || randInd + randOffset >= optixLaunchParams.simConstants.nbRandNumbersPerThread*optixLaunchParams.simConstants.size.x*optixLaunchParams.simConstants.size.y){
            printf("randInd %u is out of bounds\n", randInd + randOffset);
        }
#endif
        float facetRnd = randFloat[(unsigned int)(randInd + randOffset++)];

        int facIndex = 0;
        bool found = false;
        do{
#ifdef BOUND_CHECK
            if(facIndex < 0 || facIndex >= optixLaunchParams.simConstants.nbFacets){
                printf("facIndex %u >= %u is out of bounds\n", facIndex, optixLaunchParams.simConstants.nbFacets);
                printf("found %u with last prob %4.2f > facetRnd %4.2f \n", (rayGenData->facetProbabilities[facIndex].y >= facetRnd), optixLaunchParams.simConstants.nbFacets,
                       rayGenData->facetProbabilities[facIndex].y,facetRnd);
            }
#endif

            found = (rayGenData->facetProbabilities[facIndex].y >= facetRnd); // find probability interval rand lies in
            if(!found)
                facIndex++;

            /*if(found){
                if(facIndex==0)
                    printf("[%d] %u/%u with last prob %8.6f > facetRnd %8.6f and rand index %u+%u\n", fbIndex, facIndex, optixLaunchParams.simConstants.nbFacets,
                           rayGenData->facetProbabilities[facIndex].y,facetRnd, randInd,randOffset-1);
                else
                    printf("found %u/%u with last prob %8.6f > facetRnd %8.6f \n", facIndex, optixLaunchParams.simConstants.nbFacets,
                       rayGenData->facetProbabilities[facIndex].y,facetRnd);
            }*/
        }while(!found && facIndex<optixLaunchParams.simConstants.nbFacets);

#ifdef BOUND_CHECK
        if(facIndex < 0 || facIndex >= optixLaunchParams.simConstants.nbFacets){
            printf("Post facIndex %u >= %u is out of bounds\n", facIndex, optixLaunchParams.simConstants.nbFacets);
            printf("Post found %u with last prob %4.2f > facetRnd %4.2f \n", found, optixLaunchParams.simConstants.nbFacets,
                   rayGenData->facetProbabilities[facIndex].y,facetRnd);

        }
#endif

        hitData.hitFacetId = facIndex;
        return facIndex;

    }

    static __forceinline__ __device__
    float3 getNewOrigin(MolPRD& hitData,
#ifdef WITHTRIANGLES
        const TriangleRayGenData* rayGenData,
#else
        const PolygonRayGenData* rayGenData,
#endif
        const int& facIndex,
        const float* randFloat, unsigned int& randInd, unsigned int& randOffset)
    {

#ifdef WITHTRIANGLES

        #ifdef BOUND_CHECK
        if(randInd + randOffset < 0 || randInd + randOffset >= optixLaunchParams.simConstants.nbRandNumbersPerThread*optixLaunchParams.simConstants.size.x*optixLaunchParams.simConstants.size.y){
            printf("randInd %u is out of bounds\n", randInd + randOffset);
        }
        #endif

        float r1_sqrt = sqrtf(randFloat[(unsigned int)(randInd + randOffset++)]);

        #ifdef BOUND_CHECK
        if(randInd + randOffset < 0 || randInd + randOffset >= optixLaunchParams.simConstants.nbRandNumbersPerThread*optixLaunchParams.simConstants.size.x*optixLaunchParams.simConstants.size.y){
            printf("randInd %u is out of bounds\n", randInd + randOffset);
        }
        #endif

            float r2 = randFloat[(unsigned int)(randInd + randOffset++)];


        #ifdef BOUND_CHECK
        if(rayGenData->index[facIndex].x < 0 || rayGenData->index[facIndex].x >= optixLaunchParams.simConstants.nbVertices
        || rayGenData->index[facIndex].y < 0 || rayGenData->index[facIndex].y >= optixLaunchParams.simConstants.nbVertices
        || rayGenData->index[facIndex].z < 0 || rayGenData->index[facIndex].z >= optixLaunchParams.simConstants.nbVertices){
            printf("rayGenData->index[facIndex] %u,%u,%u >= %u is out of bounds\n", rayGenData->index[facIndex].x, rayGenData->index[facIndex].y, rayGenData->index[facIndex].z, optixLaunchParams.simConstants.nbVertices);
        }
        #endif

            float3 vertA = rayGenData->vertex[rayGenData->index[facIndex].x];
            float3 vertB = rayGenData->vertex[rayGenData->index[facIndex].y];
            float3 vertC = rayGenData->vertex[rayGenData->index[facIndex].z];

            //rayOrigin =
            return (1-r1_sqrt) * vertA + r1_sqrt * (1 - r2) * vertB + r1_sqrt * r2 * vertC; //rayGenData
#else
        // start position of particle (U,V) -> (x,y,z)
        float uDir = 0.0f, vDir = 0.0f;

        short isInPoly = 0;
        unsigned int loopNb = 0;

        while(!isInPoly && loopNb < NB_INPOLYCHECKS){
#ifdef BOUND_CHECK
            if(randInd + randOffset < 0 || randInd + randOffset >= optixLaunchParams.simConstants.nbRandNumbersPerThread*optixLaunchParams.simConstants.size.x*optixLaunchParams.simConstants.size.y){
                printf("randInd %u is out of bounds\n", randInd + randOffset);
            }
#endif
            uDir = randFloat[(unsigned int)(randInd + randOffset++)]; // float3(randFloat[randInd++],randFloat[randInd++],randFloat[randInd++])

#ifdef BOUND_CHECK
            if(randInd + randOffset < 0 || randInd + randOffset >= optixLaunchParams.simConstants.nbRandNumbersPerThread*optixLaunchParams.simConstants.size.x*optixLaunchParams.simConstants.size.y){
                printf("randInd %u is out of bounds\n", randInd + randOffset);
            }
#endif
            vDir = randFloat[(unsigned int)(randInd + randOffset++)];

            isInPoly = point_in_polygon(uDir,vDir, rayGenData->poly[facIndex]);

            ++loopNb;
        }

        // if no suitable starting point was found
        if(!isInPoly){
            uDir = 0.5;
            vDir = 0.5;
        }

        //rayOrigin =
        return rayGenData->poly[facIndex].O + uDir * rayGenData->poly[facIndex].U + vDir * rayGenData->poly[facIndex].V;
#endif //WITHTRIANGLES
    }


    static __forceinline__ __device__
    float3 getNewDirection(MolPRD& hitData, flowgeom::Polygon& poly,
            const float* randFloat, unsigned int& randInd, unsigned int& randOffset)
    {

        // generate ray direction
#ifdef BOUND_CHECK
        if(randInd + randOffset < 0 || randInd + randOffset >= optixLaunchParams.simConstants.nbRandNumbersPerThread*optixLaunchParams.simConstants.size.x*optixLaunchParams.simConstants.size.y){
            printf("randInd %u is out of bounds\n", randInd + randOffset);
        }
#endif
        const float theta = acosf(sqrtf(randFloat[(unsigned int)(randInd + randOffset++)]));


#ifdef BOUND_CHECK
        if(randInd + randOffset < 0 || randInd + randOffset >= optixLaunchParams.simConstants.nbRandNumbersPerThread*optixLaunchParams.simConstants.size.x*optixLaunchParams.simConstants.size.y){
            printf("randInd %u is out of bounds\n", randInd + randOffset);
        }
#endif
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

        //rayDir = u*nU + v*nV + n*N;
        return u*nU + v*nV + n*N;
        /*if (rayDir.x != 0.0) rayDir.x = 1.0 / rayDir.x;
        if (rayDir.y != 0.0) rayDir.y = 1.0 / rayDir.y;
        if (rayDir.z != 0.0) rayDir.z = 1.0 / rayDir.z;
        rayDir = float3(-1.0f,-1.0f,-1.0f) * rayDir;*/
    }

    static __forceinline__ __device__
    void initMoleculeInSystem(const unsigned int bufferIndex, MolPRD& hitData, float3& rayDir , float3& rayOrigin)
    {
        hitData = optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex];
        rayDir = hitData.postHitDir;
        rayOrigin = hitData.hitPos;

#ifdef DEBUG
        //printf("--- in[%d] %4.2f , %4.2f , %4.2f - %4.2f , %4.2f , %4.2f ---\n", blockDim.x * blockIdx.x + threadIdx.x,optixLaunchParams.perThreadData.currentMoleculeData[fbIndex].hitPos.x, optixLaunchParams.perThreadData.currentMoleculeData[fbIndex].hitPos.y, optixLaunchParams.perThreadData.currentMoleculeData[fbIndex].hitPos.z,optixLaunchParams.perThreadData.currentMoleculeData[fbIndex].postHitDir.x, optixLaunchParams.perThreadData.currentMoleculeData[fbIndex].postHitDir.y, optixLaunchParams.perThreadData.currentMoleculeData[fbIndex].postHitDir.z);
#endif
    }


    // Calculate new direction, but keep old position for molecule that was stuck due to self-intersection
    static __forceinline__ __device__
    void initMoleculePostSelfIntersection(const unsigned int bufferIndex, MolPRD& hitData, float3& rayDir , float3& rayOrigin)
    {
        hitData = optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex];

        float* randFloat = optixLaunchParams.randomNumbers;
        unsigned int randInd = optixLaunchParams.simConstants.nbRandNumbersPerThread*(bufferIndex);
        unsigned int randOffset = optixLaunchParams.perThreadData.randBufferOffset[bufferIndex];

#ifdef BOUND_CHECK
        if(hitData.hitFacetId < 0 || hitData.hitFacetId >= optixLaunchParams.simConstants.nbFacets){
            printf("hitData.hitFacetId %u >= %u is out of bounds\n", hitData.hitFacetId, optixLaunchParams.simConstants.nbFacets);
        }
#endif

#ifdef WITHTRIANGLES
        const TriangleRayGenData* rayGenData = (TriangleRayGenData*) optixGetSbtDataPointer();
#else
        const PolygonRayGenData* rayGenData = (PolygonRayGenData*) optixGetSbtDataPointer();
#endif

        //rayDir = hitData.postHitDir;
        rayDir = getNewDirection(hitData,rayGenData->poly[hitData.hitFacetId],randFloat,randInd,randOffset);
        rayOrigin = hitData.hitPos;
        optixLaunchParams.perThreadData.randBufferOffset[bufferIndex] = randOffset; // set value back to buffer

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

    // --------------------------------------
    // increase facet counters for desorption
    // --------------------------------------
    static __forceinline__ __device__
    void increaseHitCounterDesorption(CuFacetHitCounter& hitCounter, const MolPRD hitData, const float3 rayDir, const float3 polyNormal)
    {
        double velFactor = optixLaunchParams.simConstants.useMaxwell ? 1.0f : 1.1781f; // TODO: Save somewhere as a shared constant instead of repetively evaluating

        // TODO: fix velocity etc
        const float velocity = hitData.velocity;
        const float ortVelocity = velocity*abs(dot(rayDir, polyNormal));
        //const float hitEquiv = 0.0f; //1.0*prd.orientationRatio; // hit=1.0 (only changed for lowflux mode)
        //const float absEquiv = 0.0f; //1.0*prd.orientationRatio; // hit=1.0 (only changed for lowflux mode)
        //atomicAdd(&optixLaunchParams.hitCounter[counterIdx].nbMCHit,1);
        //++optixLaunchParams.hitCounter[counterIdx].nbMCHit;
        //const unsigned int counterIdx = facIndex + ((bufferIndex)%(CORESPERMP*WARPSCHEDULERS))*optixLaunchParams.simConstants.nbFacets;

        //atomicAdd(&hitCounter.nbMCHit, static_cast<uint32_t>(0));
        //atomicAdd(&hitCounter.nbHitEquiv, 0.0f);
        atomicAdd(&hitCounter.nbDesorbed, static_cast<uint32_t>(1));
        //atomicAdd(&hitCounter.nbAbsEquiv, 0.0f); //static_cast<double>(absorb)*sHandle->currentParticle.oriRatio;


        atomicAdd(&hitCounter.sum_1_per_ort_velocity, (2.0f / ortVelocity));//prd.oriRatio * sum_1_per_v;
        atomicAdd(&hitCounter.sum_v_ort, velFactor*ortVelocity);//(sHandle->wp.useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity; //prd.oriRatio * sum_v_ort;
        atomicAdd(&hitCounter.sum_1_per_velocity, 1.0f / velocity);//(hitEquiv + static_cast<double>(desorb)) / prd.velocity;
    }

    // --------------------------------------
    // increase texture counters for desorption
    // --------------------------------------
    static __forceinline__ __device__
    void RecordDesorptionTexture(flowgpu::Polygon *poly, float3 rayOrigin, float velocity_factor, float ortSpeedFactor){
/*
        size_t tu = (size_t)(f->colU * poly->texWidthD);
        size_t tv = (size_t)(f->colV * poly->texHeightD);
        size_t add = tu + tv * (f->sh.texWidth);
        double ortVelocity = (sHandle->wp.useMaxwellDistribution ? 1.0 : 1.1781)*sHandle->currentParticle.velocity*abs(Dot(sHandle->currentParticle.direction, f->sh.N)); //surface-orthogonal velocity component

        if (countHit) f->texture[m][add].countEquiv += sHandle->currentParticle.oriRatio;
        f->texture[m][add].sum_1_per_ort_velocity += sHandle->currentParticle.oriRatio * velocity_factor / ortVelocity;
        f->texture[m][add].sum_v_ort_per_area += sHandle->currentParticle.oriRatio * ortSpeedFactor*ortVelocity*f->textureCellIncrements[add]; // sum ortho_velocity[m/s] / cell_area[cm2]
*/

    }

    // Calculate new direction, but keep old position
    static __forceinline__ __device__
    void initMoleculeFromStart(const unsigned int bufferIndex, MolPRD& hitData, float3& rayDir , float3& rayOrigin)
    {
        hitData.hitPos = make_float3(-999.9f,-999.9f,-999.9f);
        hitData.postHitDir = make_float3(-999.9f,-999.9f,-999.9f);
        hitData.hitT = -999.0f;
        hitData.velocity = -999.0f;
        hitData.currentDepth = 0;
        hitData.inSystem = 0;

        const float* randFloat = optixLaunchParams.randomNumbers;
        unsigned int randInd = optixLaunchParams.simConstants.nbRandNumbersPerThread*(bufferIndex);
        unsigned int randOffset = optixLaunchParams.perThreadData.randBufferOffset[bufferIndex];

#ifdef WITHTRIANGLES
        const TriangleRayGenData* rayGenData = (TriangleRayGenData*) optixGetSbtDataPointer();
#else
        const PolygonRayGenData* rayGenData = (PolygonRayGenData*) optixGetSbtDataPointer();
#endif

        const int facIndex = getSourceFacet(hitData,rayGenData,randFloat,randInd,randOffset);
#ifdef BOUND_CHECK
        if(facIndex < 0 || facIndex >= optixLaunchParams.simConstants.nbFacets){
            printf("facIndex %u >= %u is out of bounds\n", facIndex, optixLaunchParams.simConstants.nbFacets);
        }
#endif
        rayOrigin = getNewOrigin(hitData,rayGenData,facIndex,randFloat,randInd,randOffset);
        rayDir = getNewDirection(hitData,rayGenData->poly[facIndex],randFloat,randInd,randOffset);


        // Write back cached local variables to shared memory
        optixLaunchParams.perThreadData.randBufferOffset[bufferIndex] = randOffset;

        // --------------------------------------
        // increase facet counters for desorption
        // --------------------------------------
        const unsigned int counterIdx = facIndex + ((bufferIndex)%(CORESPERSM * WARPSCHEDULERS)) * optixLaunchParams.simConstants.nbFacets;
#ifdef BOUND_CHECK
        if(counterIdx < 0 || counterIdx >= optixLaunchParams.simConstants.nbFacets * CORESPERSM * WARPSCHEDULERS){
            printf("facIndex %u >= %u is out of bounds\n", counterIdx, optixLaunchParams.simConstants.nbFacets * CORESPERSM * WARPSCHEDULERS);
        }
#endif
        increaseHitCounterDesorption(optixLaunchParams.hitCounter[counterIdx],hitData,rayDir, rayGenData->poly[facIndex].N);
        /*if (rayGenData->poly[facIndex].textureFlags & TextureCounter::countDes)
            RecordHitOnTexture(rayGenData->poly[facIndex], true, 2.0, 1.0); //was 2.0, 1.0
*/

    }

    //------------------------------------------------------------------------------
    // ray gen program - the actual rendering happens in here
    //------------------------------------------------------------------------------
    extern "C" __global__ void __raygen__startFromSource()
    {
        //TODO: use thrust::random or curand

        const uint2 launchIndex = make_uint2(optixGetLaunchIndex());
        const unsigned int bufferIndex = launchIndex.x + launchIndex.y * optixLaunchParams.simConstants.size.x;

#ifdef BOUND_CHECK
        if(bufferIndex < 0 || bufferIndex >= optixLaunchParams.simConstants.nbRandNumbersPerThread * optixLaunchParams.simConstants.size.x * optixLaunchParams.simConstants.size.y){
            printf("bufferIndex %u is out of bounds\n", bufferIndex);
        }
#endif
        float3 rayOrigin;
        float3 rayDir;

        MolPRD hitData;

        switch (optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].inSystem) {
            case 1:
            {
                /* if molecule is still in the system (bounce etc.)
                 * load data from thread buffer
                 */
                initMoleculeInSystem(bufferIndex, hitData, rayDir, rayOrigin);
                break;
            }
            case 0:
            {
                /*
                 * start from a source
                 */
                initMoleculeFromStart(bufferIndex, hitData, rayDir, rayOrigin);
                break;
            }
            default:
            {
                /*
                 * if molecule suffered from self intersection with starting facet, try again with new direction
                 */
                initMoleculePostSelfIntersection(bufferIndex, hitData, rayDir, rayOrigin);
            }
        }

#ifdef DEBUGPOS
        const unsigned int posIndexOffset = optixLaunchParams.perThreadData.posOffsetBuffer_debug[bufferIndex]++;
        if(posIndexOffset<NBCOUNTS){
            const unsigned int posIndex = bufferIndex*NBCOUNTS+posIndexOffset;
            //printf("[%d] my pos is %d\n", bufferIndex, posIndex);
            optixLaunchParams.perThreadData.positionsBuffer_debug[posIndex] = rayOrigin;
        }
#endif


        {
#ifdef WITHTRIANGLES
        const TriangleRayGenData* rayGenData = (TriangleRayGenData*) optixGetSbtDataPointer();
#else
        const PolygonRayGenData* rayGenData = (PolygonRayGenData*) optixGetSbtDataPointer();
#endif

            uint32_t facIndex = hitData.hitFacetId;
#ifdef BOUND_CHECK
            if(facIndex < 0 || facIndex >= optixLaunchParams.simConstants.nbFacets){
            printf("[RayOffset] facIndex %u >= %u is out of bounds (%u)\n", facIndex, optixLaunchParams.simConstants.nbFacets, optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].inSystem);
        }
#endif
            rayOrigin = offset_ray(rayOrigin,rayGenData->poly[facIndex].N);

        }

        optixTrace(optixLaunchParams.traversable,
               rayOrigin,
               rayDir,
               0.f,//optixLaunchParams.simConstants.scene_epsilon,//1.e-6f,//0.f,    // tmin // TODO: Choose scene_epsilon wisely
               1e20f,  // tmax
               0.0f,   // rayTime
               OptixVisibilityMask( 255 ),
               OPTIX_RAY_FLAG_DISABLE_ANYHIT
               | OPTIX_RAY_FLAG_CULL_BACK_FACING_TRIANGLES,//OPTIX_RAY_FLAG_NONE,
                   RayType::RAY_TYPE_MOLECULE,             // SBT offset
                   RayType::RAY_TYPE_COUNT,               // SBT stride
                   RayType::RAY_TYPE_MOLECULE,             // missSBTIndex
               //u0, u1 , u2, u3);
               //float3_as_args(hitData.hitPos),
               reinterpret_cast<unsigned int&>(hitData.velocity),
               reinterpret_cast<unsigned int&>(hitData.currentDepth),
               reinterpret_cast<unsigned int&>(hitData.inSystem),
            /* Can't use float_as_int() because it returns rvalue but payload requires a lvalue */
               reinterpret_cast<unsigned int&>(hitData.hitFacetId),
               reinterpret_cast<unsigned int&>(hitData.hitT));

        //hitData.hitPos = hitData.hitOri + hitData.hitT * hitData.hitDir;

        // and write to thread buffer ...
        optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].velocity = hitData.velocity;
        optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].currentDepth = hitData.currentDepth;
        optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].inSystem = hitData.inSystem;
        optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].hitFacetId = hitData.hitFacetId;
        optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].hitT = hitData.hitT;
    }
} // ::flowgpu
