// Created by pbahr

#include <math_constants.h>

#include "LaunchParams.h"
#include "GPUDefines.h" // for NB_RAND
#include "helper_math.h"
#include <LaunchParams.h>
#include "CommonFunctions.cuh"




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

    extern "C" __device__ int point_in_polygon(float u, float v, const flowgpu::Polygon& poly) {
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
/*
        if(poly.indexOffset < 0 || poly.indexOffset+poly.nbVertices >= optixLaunchParams.simConstants.nbVertices){
            printf("[%d] indexOffset %u -- %u >= %u is out of bounds\n", poly.parentIndex, poly.indexOffset, poly.indexOffset+poly.nbVertices, optixLaunchParams.simConstants.nbVertices);
        }*/
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
                       const RN_T* randFloat, unsigned int& randInd, unsigned int& randOffset)
    {
#ifdef BOUND_CHECK
        if(randInd + randOffset < 0 || randInd + randOffset >= optixLaunchParams.simConstants.nbRandNumbersPerThread*optixLaunchParams.simConstants.size.x*optixLaunchParams.simConstants.size.y){
            printf("randInd %u is out of bounds\n", randInd + randOffset);
        }
        else if(randFloat[(unsigned int)(randInd + randOffset)] < 0.0 || randFloat[(unsigned int)(randInd + randOffset)] > 1.0){
            printf("randFloat %lf is out of bounds [%u + %u]\n", randFloat[(unsigned int)(randInd + randOffset)], randInd, randOffset);
        }
#endif
        float facetRnd = randFloat[(unsigned int)(randInd + randOffset++)];

        int facIndex = 0;
        bool found = false;
        do{
#ifdef BOUND_CHECK
            if(facIndex < 0 || facIndex >= optixLaunchParams.simConstants.nbFacets){
                printf("facIndex %u >= %u is out of bounds\n", facIndex, optixLaunchParams.simConstants.nbFacets);
                printf("[%d] found %u with last prob %4.2f > facetRnd %4.2f \n", (rayGenData->facetProbabilities[facIndex].y >= facetRnd), optixLaunchParams.simConstants.nbFacets,
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
            printf("[%d] Post found %u with last prob %4.2f > facetRnd %4.2f \n", found, optixLaunchParams.simConstants.nbFacets,
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
                        const RN_T* randFloat, unsigned int& randInd, unsigned int& randOffset)
    {

#ifdef WITHTRIANGLES

#ifdef BOUND_CHECK
        if(randInd + randOffset < 0 || randInd + randOffset >= optixLaunchParams.simConstants.nbRandNumbersPerThread*optixLaunchParams.simConstants.size.x*optixLaunchParams.simConstants.size.y){
            printf("randInd %u is out of bounds\n", randInd + randOffset);
        }
#endif
        hitData.rndOrigin[0] = randFloat[(unsigned int)(randInd + randOffset)];
        hitData.rndOrigin[1] = randFloat[(unsigned int)(randInd + randOffset+1)];

        RN_T r1_sqrt =
                #ifdef RNG64
                    sqrt((double)randFloat[(unsigned int)(randInd + randOffset++)]);
                #else
                    sqrtf(randFloat[(unsigned int)(randInd + randOffset++)]);
                #endif


#ifdef BOUND_CHECK
        if(randInd + randOffset < 0 || randInd + randOffset >= optixLaunchParams.simConstants.nbRandNumbersPerThread*optixLaunchParams.simConstants.size.x*optixLaunchParams.simConstants.size.y){
            printf("randInd %u is out of bounds\n", randInd + randOffset);
        }
#endif

        RN_T r2 = randFloat[(unsigned int)(randInd + randOffset++)];

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

#ifdef RNG64
        return make_float3((1.0 - r1_sqrt) * vertA + r1_sqrt * (1.0 - r2) * vertB + r1_sqrt * r2 * vertC); //rayGenData
#else
        return (1.0f-r1_sqrt) * vertA + r1_sqrt * (1.0f - r2) * vertB + r1_sqrt * r2 * vertC; //rayGenData
#endif

#else //WITHTRIANGLES

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

        hitData.rndOrigin[0] = uDir;
        hitData.rndOrigin[1] = vDir;

        //rayOrigin =
        return rayGenData->poly[facIndex].O + uDir * rayGenData->poly[facIndex].U + vDir * rayGenData->poly[facIndex].V;
#endif //WITHTRIANGLES
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

        RN_T* randFloat = optixLaunchParams.randomNumbers;
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
        /*if(rayGenData->poly[hitData.hitFacetId].facProps.is2sided && rayGenData->poly[hitData.hitFacetId].facProps.opacity == 0.0f){
            // todo: just offset later, no direction change
            if(bufferIndex == 0) printf(" Self intersection back facet hit %d\n",rayGenData->poly[hitData.hitFacetId].parentIndex);
            rayDir = getNewReverseDirection(hitData, rayGenData->poly[hitData.hitFacetId], randFloat, randInd,
                                            randOffset);
        }
        else *//*if(
                (hitData.facetHitSide == OPTIX_HIT_KIND_TRIANGLE_FRONT_FACE)
                )*//*{
            rayDir = getNewDirection(hitData, rayGenData->poly[hitData.hitFacetId], randFloat, randInd, randOffset);
            if(hitData.facetHitSide != OPTIX_HIT_KIND_TRIANGLE_FRONT_FACE)
                if(bufferIndex == 0)
                    printf(" Self intersection back facet hit %d\n",rayGenData->poly[hitData.hitFacetId].parentIndex);

        }*/
        /*else {
            if(bufferIndex == 0) printf(" Self intersection back facet hit %d\n",rayGenData->poly[hitData.hitFacetId].parentIndex);
            rayDir = getNewReverseDirection(hitData, rayGenData->poly[hitData.hitFacetId], randFloat, randInd,
                                            randOffset);
        }*/
        //rayDir = getNewDirection(hitData,rayGenData->poly[hitData.hitFacetId],randFloat,randInd,randOffset);


        /*if(hitData.facetHitSide == OPTIX_HIT_KIND_TRIANGLE_BACK_FACE) {
            hitData.hitPos = offset_ray(hitData.hitPos, -1.0f*rayGenData->poly[hitData.hitFacetId].N);
        }
        else{
            hitData.hitPos = offset_ray_old(hitData.hitPos,1.0f* rayGenData->poly[hitData.hitFacetId].N);
        }*/

        rayOrigin = hitData.hitPos;
        optixLaunchParams.perThreadData.randBufferOffset[bufferIndex] = randOffset; // set value back to buffer

    }

// Calculate new position for molecule with an offset
static __forceinline__ __device__
void initMoleculeTransparentHit(const unsigned int bufferIndex, MolPRD& hitData, float3& rayDir , float3& rayOrigin)
{
    hitData = optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex];

    //rayDir = hitData.postHitDir;
    rayDir = hitData.postHitDir;
    rayOrigin = hitData.hitPos;

}

    // --------------------------------------
    // increase facet counters for desorption
    // --------------------------------------
    static __forceinline__ __device__
    void increaseHitCounterDesorption(CuFacetHitCounter& hitCounter, const MolPRD &hitData, const float3 &rayDir, const float3 &polyNormal)
    {

        //const float hitEquiv = 0.0f; //1.0*prd.orientationRatio; // hit=1.0 (only changed for lowflux mode)
        //const float absEquiv = 0.0f; //1.0*prd.orientationRatio; // hit=1.0 (only changed for lowflux mode)
        //atomicAdd(&optixLaunchParams.hitCounter[counterIdx].nbMCHit,1);
        //++optixLaunchParams.hitCounter[counterIdx].nbMCHit;
        //const unsigned int counterIdx = facIndex + ((bufferIndex)%(CORESPERMP*WARPSCHEDULERS))*optixLaunchParams.simConstants.nbFacets;

        //atomicAdd(&hitCounter.nbMCHit, static_cast<uint32_t>(0));
        //atomicAdd(&hitCounter.nbHitEquiv, 0.0f);

#ifdef HIT64
        const FLOAT_T velFactor = optixLaunchParams.simConstants.useMaxwell ? 1.0 : 1.1781;
        const FLOAT_T velocity = hitData.velocity;
        const FLOAT_T ortVelocity = velocity*fabsf(dot(rayDir, polyNormal));

        atomicAdd(&hitCounter.nbDesorbed, static_cast<uint64_t>(1));
        atomicAdd(&hitCounter.sum_1_per_ort_velocity, (2.0 / ortVelocity));//prd.oriRatio * sum_1_per_v;
        atomicAdd(&hitCounter.sum_v_ort, velFactor*ortVelocity);//(sHandle->wp.useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity; //prd.oriRatio * sum_v_ort;
        atomicAdd(&hitCounter.sum_1_per_velocity, 1.0 / velocity);//(hitEquiv + static_cast<double>(desorb)) / prd.velocity;
#else
        const FLOAT_T velFactor = optixLaunchParams.simConstants.useMaxwell ? 1.0f : 1.1781f;
        const FLOAT_T velocity = hitData.velocity;
        const FLOAT_T ortVelocity = velocity*fabsf(dot(rayDir, polyNormal));

        atomicAdd(&hitCounter.nbDesorbed, static_cast<uint32_t>(1));
        atomicAdd(&hitCounter.sum_1_per_ort_velocity, (2.0f / ortVelocity));//prd.oriRatio * sum_1_per_v;
        atomicAdd(&hitCounter.sum_v_ort, velFactor*ortVelocity);//(sHandle->wp.useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity; //prd.oriRatio * sum_v_ort;
        atomicAdd(&hitCounter.sum_1_per_velocity, 1.0f / velocity);//(hitEquiv + static_cast<double>(desorb)) / prd.velocity;
#endif
    }

    // --------------------------------------
    // increase texture counters for desorption
    // --------------------------------------
    static __forceinline__ __device__
    void RecordDesorptionTexture(const flowgpu::Polygon& poly, const MolPRD& hitData, const float3& rayOrigin, const float3& rayDir){

        const float2 hitLocation = getHitLocation(poly,rayOrigin);

        flowgpu::FacetTexture& facetTex = optixLaunchParams.sharedData.facetTextures[poly.texProps.textureOffset];

        const unsigned int tu = (unsigned int)(__saturatef(hitLocation.x) * facetTex.texWidthD); // saturate for rounding errors that can result for values larger than 1.0
        const unsigned int tv = (unsigned int)(__saturatef(hitLocation.y) * facetTex.texHeightD);
        unsigned int add = tu + tv * (facetTex.texWidth);
        //printf("Hit at %lf , %lf for tex %lf , %lf\n", hitLocation.x, hitLocation.y, facetTex.texWidthD, facetTex.texHeightD);

#ifdef BOUND_CHECK
        if(facetTex.texelOffset + add < 0 || facetTex.texelOffset + add >= optixLaunchParams.simConstants.nbTexel){printf("facetTex.texelOffset + add %u >= %u is out of bounds\n", facetTex.texelOffset + add, optixLaunchParams.simConstants.nbTexel);}
#endif

        flowgpu::Texel& tex = optixLaunchParams.sharedData.texels[facetTex.texelOffset + add];

#ifdef HIT64
        const FLOAT_T velocity_factor = 2.0;
        const FLOAT_T ortSpeedFactor = 1.0;
        const FLOAT_T ortVelocity = (optixLaunchParams.simConstants.useMaxwell ? 1.0 : 1.1781) * hitData.velocity*fabs(dot(rayDir, poly.N)); //surface-orthogonal velocity component
        atomicAdd(&tex.countEquiv, static_cast<uint64_t>(1));
        atomicAdd(&tex.sum_1_per_ort_velocity, 1.0 * velocity_factor / ortVelocity);
        atomicAdd(&tex.sum_v_ort_per_area, 1.0 * ortSpeedFactor * ortVelocity * optixLaunchParams.sharedData.texelInc[facetTex.texelOffset + add]); // sum ortho_velocity[m/s] / cell_area[cm2]
#else
        const FLOAT_T velocity_factor = 2.0f;
        const FLOAT_T ortSpeedFactor = 1.0f;
        FLOAT_T ortVelocity = (optixLaunchParams.simConstants.useMaxwell ? 1.0f : 1.1781f) * hitData.velocity*fabsf(dot(rayDir, poly.N)); //surface-orthogonal velocity component
        atomicAdd(&tex.countEquiv, static_cast<uint32_t>(1));
        atomicAdd(&tex.sum_1_per_ort_velocity, 1.0f * velocity_factor / ortVelocity);
        atomicAdd(&tex.sum_v_ort_per_area, 1.0f * ortSpeedFactor * ortVelocity * optixLaunchParams.sharedData.texelInc[facetTex.texelOffset + add]); // sum ortho_velocity[m/s] / cell_area[cm2]
#endif


        //printf("Desorbing on facet#%u %u + %u = %u\n", poly.parentIndex, tu, tv, add);

    }

    // --------------------------------------
    // increase texture counters for desorption
    // --------------------------------------
    static __forceinline__ __device__
    void RecordDesorptionProfile(const flowgpu::Polygon& poly, const MolPRD& hitData, const float3& rayOrigin, const float3& rayDir){

        const float2 hitLocation = getHitLocation(poly,rayOrigin);

        unsigned int add = 0;
        if(poly.profProps.profileType == flowgpu::PROFILE_FLAGS::profileU){
            add = (unsigned int)(__saturatef(hitLocation.x) * PROFILE_SIZE);
        }
        else if(poly.profProps.profileType == flowgpu::PROFILE_FLAGS::profileV){
            add = (unsigned int)(__saturatef(hitLocation.y) * PROFILE_SIZE);
        }


#ifdef BOUND_CHECK
        if(poly.profProps.profileOffset + add < 0 || poly.profProps.profileOffset + add >= optixLaunchParams.simConstants.nbProfSlices){
            printf("[DES] poly.profProps.profileOffset + add %u >= %u is out of bounds\n", poly.profProps.profileOffset + add, optixLaunchParams.simConstants.nbProfSlices);}
#endif
        flowgpu::Texel& tex = optixLaunchParams.sharedData.profileSlices[poly.profProps.profileOffset + add];

#ifdef HIT64
        const FLOAT_T velocity_factor = 2.0;
        const FLOAT_T ortSpeedFactor = 1.0;
        const FLOAT_T ortVelocity = hitData.velocity*fabsf(dot(rayDir, poly.N)); //surface-orthogonal velocity component

        atomicAdd(&tex.countEquiv, static_cast<uint64_t>(1));
        atomicAdd(&tex.sum_1_per_ort_velocity, velocity_factor / ortVelocity);
        atomicAdd(&tex.sum_v_ort_per_area, ortSpeedFactor * ortVelocity * (optixLaunchParams.simConstants.useMaxwell ? 1.0 : 1.1781)); // sum ortho_velocity[m/s] / cell_area[cm2]
#else
        const FLOAT_T velocity_factor = 2.0f;
        const FLOAT_T ortSpeedFactor = 1.0f;
        const FLOAT_T ortVelocity = hitData.velocity*fabsf(dot(rayDir, poly.N)); //surface-orthogonal velocity component

        atomicAdd(&tex.countEquiv, static_cast<uint32_t>(1));
        atomicAdd(&tex.sum_1_per_ort_velocity, velocity_factor / ortVelocity);
        atomicAdd(&tex.sum_v_ort_per_area, ortSpeedFactor * ortVelocity * (optixLaunchParams.simConstants.useMaxwell ? 1.0f : 1.1781f)); // sum ortho_velocity[m/s] / cell_area[cm2]
#endif
    }

static __forceinline__ __device__
void recordDesorption(const unsigned int& counterIdx, const flowgpu::Polygon& poly, const MolPRD& hitData, const float3& rayDir, const float3& rayOrigin)
{
#ifdef BOUND_CHECK
    if(counterIdx < 0 || counterIdx >= optixLaunchParams.simConstants.nbFacets * CORESPERSM * WARPSCHEDULERS){
            printf("facIndex %u >= %u is out of bounds\n", counterIdx, optixLaunchParams.simConstants.nbFacets * CORESPERSM * WARPSCHEDULERS);
        }
#endif
    increaseHitCounterDesorption(optixLaunchParams.hitCounter[counterIdx],hitData,rayDir, poly.N);
#ifdef WITH_TEX
    if (poly.texProps.textureFlags & flowgpu::TEXTURE_FLAGS::countDes)
            RecordDesorptionTexture(poly, hitData, rayOrigin, rayDir);
#endif
#ifdef WITH_PROF
    if (poly.profProps.profileType != flowgpu::PROFILE_FLAGS::noProfile)
            RecordDesorptionProfile(poly, hitData, rayOrigin, rayDir);
#endif
}

    // Calculate new direction, but keep old position
    static __forceinline__ __device__
    void initMoleculeFromStart(const unsigned int bufferIndex, MolPRD& hitData, float3& rayDir , float3& rayOrigin) {
        hitData.hitPos = make_float3(-999.9f, -999.9f, -999.9f);
        hitData.postHitDir = make_float3(-999.9f, -999.9f, -999.9f);
        hitData.hitT = -999.0f;
        hitData.velocity = -999.0f;
        hitData.currentDepth = 0;
#ifdef GPUNBOUNCE
        hitData.nbBounces = 0;
#endif
        hitData.inSystem = NEW_PARTICLE;
        hitData.orientationRatio = 1.0f;
        hitData.rndOrigin[0] = -999999.9;
        hitData.rndOrigin[1] = -999999.9;
        hitData.rndDirection[0] = -999999.9;
        hitData.rndDirection[1] = -999999.9;

        const RN_T *randFloat = optixLaunchParams.randomNumbers;
        unsigned int randInd = optixLaunchParams.simConstants.nbRandNumbersPerThread * (bufferIndex);
        unsigned int randOffset = optixLaunchParams.perThreadData.randBufferOffset[bufferIndex];

#ifdef WITHTRIANGLES
        const TriangleRayGenData* rayGenData = (TriangleRayGenData*) optixGetSbtDataPointer();
#else
        const PolygonRayGenData *rayGenData = (PolygonRayGenData *) optixGetSbtDataPointer();
#endif

        const int facIndex = getSourceFacet(hitData, rayGenData, randFloat, randInd, randOffset);
#ifdef BOUND_CHECK
        if(facIndex < 0 || facIndex >= optixLaunchParams.simConstants.nbFacets){
            printf("facIndex %u >= %u is out of bounds\n", facIndex, optixLaunchParams.simConstants.nbFacets);
        }
#endif
/*#ifdef DEBUG
        printf("[%u] generating direction\n", bufferIndex);
#endif*/

        rayOrigin = getNewOrigin(hitData, rayGenData, facIndex, randFloat, randInd, randOffset);
        rayDir = getNewDirection(hitData, rayGenData->poly[facIndex], randFloat, randInd, randOffset);

        hitData.velocity = getNewVelocity(rayGenData->poly[facIndex], optixLaunchParams.simConstants.gasMass);

        // Write back cached local variables to shared memory
        optixLaunchParams.perThreadData.randBufferOffset[bufferIndex] = randOffset;

        // --------------------------------------
        // increase facet counters for desorption
        // --------------------------------------
        const unsigned int counterIdx =
                facIndex + ((bufferIndex) % (CORESPERSM * WARPSCHEDULERS)) * optixLaunchParams.simConstants.nbFacets;
        recordDesorption(counterIdx, rayGenData->poly[facIndex], hitData, rayDir, rayOrigin);
    }


    //const __device__ float offset_val = 1.0f/64.0f;
    //const __device__ float offset_val_n = -1.0f/64.0f;
    const __device__ float offset_val = 1.0f/1.0f;
    const __device__ float offset_val_n = -1.0f/1.0f;
    //------------------------------------------------------------------------------
    // ray gen program - the actual rendering happens in here
    //------------------------------------------------------------------------------
    extern "C" __global__ void __raygen__startFromSource()
    {
        //TODO: use thrust::random or curand

        const unsigned int bufferIndex = getWorkIndex();

#ifdef BOUND_CHECK
        if(bufferIndex < 0 || bufferIndex >= optixLaunchParams.simConstants.nbRandNumbersPerThread * optixLaunchParams.simConstants.size.x * optixLaunchParams.simConstants.size.y){
            printf("bufferIndex %u > %u is out of bounds\n", bufferIndex,optixLaunchParams.simConstants.nbRandNumbersPerThread * optixLaunchParams.simConstants.size.x * optixLaunchParams.simConstants.size.y);
        }
#endif
        float3 rayOrigin;
        float3 rayDir;

        MolPRD& hitData = optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex];

#ifdef WITHDESORPEXIT
        if(hitData.hasToTerminate==2){
            return;
        }
#endif
        switch (hitData.inSystem) {
            case TRANSPARENT_HIT: {
                initMoleculeTransparentHit(bufferIndex, hitData, rayDir, rayOrigin);
                break;
            }
            case ACTIVE_PARTICLE: {
                /* if molecule is still in the system (bounce etc.)
                 * load data from thread buffer
                 */
                initMoleculeInSystem(bufferIndex, hitData, rayDir, rayOrigin);
                break;
            }
            case NEW_PARTICLE: {
#ifdef WITHDESORPEXIT
                if(hitData.hasToTerminate==1){
                    hitData.hasToTerminate=2;
                    return;
                }
#endif
                /*
                 * start from a source
                 */
                initMoleculeFromStart(bufferIndex, hitData, rayDir, rayOrigin);
                break;
            }
            case SELF_INTERSECTION: {
                /*
                 * if molecule suffered from self intersection with starting facet, try again with new direction
                 */
                initMoleculePostSelfIntersection(bufferIndex, hitData, rayDir, rayOrigin);
                break;
            }
            default: {
                printf("Unknown launch status!: %u\n",hitData.inSystem);
                break;
            }
        }

#ifdef DEBUGPOS
    if(bufferIndex==0){
        const unsigned int posIndexOffset = optixLaunchParams.perThreadData.posOffsetBuffer_debug[bufferIndex]++;
        if(posIndexOffset<NBPOSCOUNTS){
            const unsigned int posIndex = bufferIndex*NBPOSCOUNTS+posIndexOffset;
            //printf("[%d] my pos is %d\n", bufferIndex, posIndex);
            optixLaunchParams.perThreadData.positionsBuffer_debug[posIndex] = rayOrigin;
        }
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
                printf("[RayOffset] facIndex %u >= %u is out of bounds (%u)\n", facIndex, optixLaunchParams.simConstants.nbFacets, hitData.inSystem);
            }
#endif
            //do not offset a transparent hit
            //if(optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].inSystem != TRANSPARENT_HIT){
            float3 facNormal = rayGenData->poly[facIndex].N;
            /*if((hitData.facetHitSide == OPTIX_HIT_KIND_TRIANGLE_BACK_FACE && hitData.inSystem != TRANSPARENT_HIT)
               || (hitData.facetHitSide == OPTIX_HIT_KIND_TRIANGLE_FRONT_FACE && hitData.inSystem == TRANSPARENT_HIT))*/

            // Don't flip normal on a normal backface hit (should only be allowed for 2sided transparent facets)
            if((/*rayGenData->poly[hitData.hitFacetId].facProps.is2sided && */rayGenData->poly[hitData.hitFacetId].facProps.opacity == 0.0f)
                    && (hitData.facetHitSide == OPTIX_HIT_KIND_TRIANGLE_FRONT_FACE) )

            // previous backface hit
            {
                /*if(bufferIndex == 0)
                    printf("[%d] reverting offset -> %d -> %d\n",bufferIndex,hitData.inSystem, hitData.nbBounces);
                */facNormal *= (offset_val_n);
                //rayOrigin = offset_ray(rayOrigin, (-1.0f) * rayGenData->poly[facIndex].N);
            }
            else{
                facNormal *= (offset_val);
            }
            //rayOrigin = offset_ray(rayOrigin,facNormal);

            //}

        }

        int hi_vel = __double2hiint(hitData.velocity);
        int lo_vel = __double2loint(hitData.velocity);
        if(isnan(rayDir.x)) {
            unsigned int randInd = optixLaunchParams.simConstants.nbRandNumbersPerThread*(bufferIndex);
            unsigned int randOffset = optixLaunchParams.perThreadData.randBufferOffset[bufferIndex];

            DEBUG_PRINT("[%d] ray dir NaN -> %d / %d + %d\n",bufferIndex,hitData.inSystem, randInd, randOffset);
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
               reinterpret_cast<unsigned int&>(hitData.currentDepth),
               reinterpret_cast<unsigned int&>(hitData.inSystem),
            /* Can't use float_as_int() because it returns rvalue but payload requires a lvalue */
               reinterpret_cast<unsigned int&>(hitData.hitFacetId),
               reinterpret_cast<unsigned int&>(hitData.hitT),
               reinterpret_cast<unsigned int&>( hi_vel ),
               reinterpret_cast<unsigned int&>( lo_vel ),
               reinterpret_cast<unsigned int&>(hitData.nbBounces),
               reinterpret_cast<unsigned int&>(hitData.facetHitSide)

        );

        //hitData.hitPos = hitData.hitOri + hitData.hitT * hitData.hitDir;
        hitData.velocity = __hiloint2double(hi_vel,lo_vel );
#if defined(DEBUG) && defined(GPUNBOUNCE)
        if(bufferIndex == optixLaunchParams.simConstants.size.x - 1 && (hitData.nbBounces % (int)1e4 == 1e4 - 1))
            printf("[%d] has new launch status -> %d and bounces %d / %d\n",bufferIndex,hitData.inSystem,optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].nbBounces,hitData.nbBounces);
#endif
/*
        // and write to thread buffer ...
        optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].velocity = hitData.velocity;
        hitData.currentDepth = hitData.currentDepth;
#ifdef GPUNBOUNCE
#ifdef DEBUG
        if(bufferIndex == optixLaunchParams.simConstants.size.x - 1 && (hitData.nbBounces % (int)1e4 == 1e4 - 1))
            printf("[%d] has new launch status -> %d and bounces %d / %d\n",bufferIndex,hitData.inSystem,optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].nbBounces,hitData.nbBounces);
#endif
        optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].nbBounces = hitData.nbBounces;
#endif
        optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].inSystem = hitData.inSystem;
        optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].hitFacetId = hitData.hitFacetId;
        optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].hitT = hitData.hitT;

        optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].facetHitSide = hitData.facetHitSide;
*/

    }
} // ::flowgpu
