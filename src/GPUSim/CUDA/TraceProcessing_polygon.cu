// Created by pbahr

#include <optix_device.h>
#include <math_constants.h>
#include "helper_math.h"
#include <cooperative_groups.h>

#include "LaunchParams.h"
#include "GPUDefines.h" // for NB_RAND
#include <LaunchParams.h>
#include "CommonFunctions.cuh"

using namespace flowgpu;

namespace flowgpu {

    /*! launch parameters in constant memory, filled in by optix upon
        optixLaunch (this gets filled in from the buffer we pass to
        optixLaunch) */
    extern "C" __constant__ LaunchParams optixLaunchParams;

#if defined(PAYLOAD_DIRECT)
    static __device__ __inline__ MolPRD getMolPRD()
    {
        MolPRD prd;
        prd.currentDepth = optixGetPayload_0();
        prd.inSystem = optixGetPayload_1();
        prd.hitFacetId = optixGetPayload_2();
        prd.hitT = int_as_float( optixGetPayload_3() );
        prd.velocity = __hiloint2double( optixGetPayload_4(),optixGetPayload_5() );

#ifdef GPUNBOUNCE
        prd.nbBounces = optixGetPayload_6();
#endif
        prd.facetHitSide = optixGetPayload_7();
        return prd;
    }

    static __device__ __inline__ void setMolPRD( const MolPRD &prd )
    {
        optixSetPayload_0( prd.currentDepth );
        optixSetPayload_1( prd.inSystem );
        optixSetPayload_2( prd.hitFacetId );
        optixSetPayload_3( float_as_int(prd.hitT) );
        optixSetPayload_4( __double2hiint(prd.velocity) );
        optixSetPayload_5( __double2loint(prd.velocity) );

#ifdef GPUNBOUNCE
        optixSetPayload_6( prd.nbBounces );
#endif
        optixSetPayload_7( prd.facetHitSide );

    }
#endif

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
    // increase texture counters for absorption (same as desorp)
    // --------------------------------------
    static __forceinline__ __device__
    void RecordAbsorpTexture(const flowgpu::Polygon &poly, MolPRD &hitData, float3 rayOrigin, float3 rayDir){



        const float hitLocationU = uint_as_float(optixGetAttribute_3());
        const float hitLocationV = uint_as_float(optixGetAttribute_4());

        flowgpu::FacetTexture& facetTex = optixLaunchParams.sharedData.facetTextures[poly.texProps.textureOffset];
        const auto tu = (unsigned int)(hitLocationU * facetTex.texWidthD);
        const auto tv = (unsigned int)(hitLocationV * facetTex.texHeightD);
        unsigned int add = tu + tv * (facetTex.texWidth);

        //printf("Hit at %lf , %lf for tex %lf , %lf : %d , %d\n", hitLocationU, hitLocationV, facetTex.texWidthD, facetTex.texHeightD, tu, tv);
        const float velocity_factor = 2.0f;
        const float ortSpeedFactor = 1.0f;
        const float ortVelocity = (optixLaunchParams.simConstants.useMaxwell ? 1.0f : 1.1781f) * hitData.velocity*fabsf(dot(rayDir, poly.N)); //surface-orthogonal velocity component

#ifdef BOUND_CHECK
        if(facetTex.texelOffset + add < 0 || facetTex.texelOffset + add >= optixLaunchParams.simConstants.nbTexel){printf("facetTex.texelOffset + add %u >= %u is out of bounds\n", facetTex.texelOffset + add, optixLaunchParams.simConstants.nbTexel);}
#endif
        flowgpu::Texel& tex = optixLaunchParams.sharedData.texels[facetTex.texelOffset + add];
        atomicAdd(&tex.countEquiv, static_cast<uint32_t>(1u));
        atomicAdd(&tex.sum_1_per_ort_velocity, 1.0 * velocity_factor / ortVelocity);
        atomicAdd(&tex.sum_v_ort_per_area, 1.0 * ortSpeedFactor * ortVelocity * optixLaunchParams.sharedData.texelInc[facetTex.texelOffset + add]); // sum ortho_velocity[m/s] / cell_area[cm2]
    }


    // --------------------------------------
    // increase texture counters for reflections
    // --------------------------------------
    static __forceinline__ __device__
    unsigned int RecordBounceTexture(const flowgpu::Polygon& poly, MolPRD& hitData, const float3& rayOrigin, const float3& rayDir){

        flowgpu::FacetTexture& facetTex = optixLaunchParams.sharedData.facetTextures[poly.texProps.textureOffset];

        const float hitLocationU = uint_as_float(optixGetAttribute_3());
        const float hitLocationV = uint_as_float(optixGetAttribute_4());
        const auto tu = (unsigned int)(hitLocationU * facetTex.texWidthD);
        const auto tv = (unsigned int)(hitLocationV * facetTex.texHeightD);
        unsigned int add = tu + tv * (facetTex.texWidth);

        //printf("Pre Bounce Tex: %f = %f * %f * %f\n",ortVelocity,(optixLaunchParams.simConstants.useMaxwell ? 1.0f : 1.1781f) ,hitData.velocity,fabsf(dot(rayDir, poly.N)));
#ifdef BOUND_CHECK
        if(facetTex.texelOffset + add < 0 || facetTex.texelOffset + add >= optixLaunchParams.simConstants.nbTexel){printf("facetTex.texelOffset + add %u >= %u is out of bounds\n", facetTex.texelOffset + add, optixLaunchParams.simConstants.nbTexel);}
#endif
        flowgpu::Texel& tex = optixLaunchParams.sharedData.texels[facetTex.texelOffset + add];

        const float velocity_factor = 1.0f;
        const float ortSpeedFactor = 1.0f;
        const float oriRatio = 1.0f; // TODO: Part of particle
        const float ortVelocity = (optixLaunchParams.simConstants.useMaxwell ? 1.0f : 1.1781f) * hitData.velocity * fabsf(dot(rayDir, poly.N)); //surface-orthogonal velocity component

        atomicAdd(&tex.countEquiv, static_cast<uint32_t>(oriRatio));
        atomicAdd(&tex.sum_1_per_ort_velocity, (float)(oriRatio * velocity_factor / ortVelocity));
        atomicAdd(&tex.sum_v_ort_per_area, (float)(oriRatio * ortSpeedFactor * ortVelocity * optixLaunchParams.sharedData.texelInc[facetTex.texelOffset + add])); // sum ortho_velocity[m/s] / cell_area[cm2]
        //printf("Pre Bounce Tex[%d]: %f , %f - %f (%f)\n", poly.parentIndex, oriRatio * velocity_factor / ortVelocity,oriRatio * ortSpeedFactor * ortVelocity * optixLaunchParams.sharedData.texelInc[facetTex.texelOffset + add],ortVelocity,optixLaunchParams.sharedData.texelInc[facetTex.texelOffset + add]);

        // Return Texel index for reuse in post bounce processing
        return facetTex.texelOffset + add;
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
    }

    // --------------------------------------
    // increase texture counters for absorption (same as desorp)
    // --------------------------------------
    static __forceinline__ __device__
    void RecordAbsorpProfile(const flowgpu::Polygon &poly, MolPRD &hitData, float3 rayOrigin, float3 rayDir){

        const float velocity_factor = 2.0f;
        const float ortSpeedFactor = 1.0f;

        float hitLocation;
        if(poly.profProps.profileType == flowgpu::PROFILE_FLAGS::profileU){
            hitLocation = uint_as_float(optixGetAttribute_3());
        }
        else if(poly.profProps.profileType == flowgpu::PROFILE_FLAGS::profileV){
            hitLocation = uint_as_float(optixGetAttribute_4());
        }
        auto add = (unsigned int)(hitLocation * PROFILE_SIZE);

        float ortVelocity = hitData.velocity*fabsf(dot(rayDir, poly.N)); //surface-orthogonal velocity component
        //printf("Hit at %lf , %lf for tex %lf , %lf\n", hitLocationU, hitLocationV, facetTex.texWidthD, facetTex.texHeightD);
#ifdef BOUND_CHECK
        if(poly.profProps.profileOffset + add < 0 || poly.profProps.profileOffset + add >= optixLaunchParams.simConstants.nbProfSlices){
            printf("[ABS] poly.profProps.profileOffset + add %u >= %u is out of bounds\n", poly.profProps.profileOffset + add, optixLaunchParams.simConstants.nbProfSlices);}
#endif
        flowgpu::Texel& tex = optixLaunchParams.sharedData.profileSlices[poly.profProps.profileOffset + add];
        atomicAdd(&tex.countEquiv, static_cast<uint32_t>(1u));
        atomicAdd(&tex.sum_1_per_ort_velocity, 1.0 * velocity_factor / ortVelocity);
        atomicAdd(&tex.sum_v_ort_per_area, 1.0 * ortSpeedFactor * ortVelocity * (optixLaunchParams.simConstants.useMaxwell ? 1.0f : 1.1781f)); // sum ortho_velocity[m/s] / cell_area[cm2]
    }

    // --------------------------------------
    // increase texture counters for reflections
    // --------------------------------------
    static __forceinline__ __device__
    unsigned int RecordBounceProfile(const flowgpu::Polygon& poly, MolPRD& hitData, const float3& rayOrigin, const float3& rayDir){

        float hitLocation;
        if(poly.profProps.profileType == flowgpu::PROFILE_FLAGS::profileU){
            hitLocation = uint_as_float(optixGetAttribute_3());
        }
        else if(poly.profProps.profileType == flowgpu::PROFILE_FLAGS::profileV){
            hitLocation = uint_as_float(optixGetAttribute_4());
        }
        auto add = (unsigned int)(hitLocation * PROFILE_SIZE);

        //printf("Pre Bounce Tex: %f = %f * %f * %f\n",ortVelocity,(optixLaunchParams.simConstants.useMaxwell ? 1.0f : 1.1781f) ,hitData.velocity,fabsf(dot(rayDir, poly.N)));
#ifdef BOUND_CHECK
        if(poly.profProps.profileOffset + add < 0 || poly.profProps.profileOffset + add >= optixLaunchParams.simConstants.nbProfSlices){
            printf("[HIT] poly.profProps.profileOffset + add %u >= %u is out of bounds\n", poly.profProps.profileOffset + add, optixLaunchParams.simConstants.nbProfSlices);}
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

    static __forceinline__ __device__
    void recordAbsorption(const unsigned int& counterIdx, const flowgpu::Polygon& poly, MolPRD& prd, const float3& rayDir, const float3 rayOrigin)
    {
#ifdef HIT64
        const FLOAT_T velFactor = optixLaunchParams.simConstants.useMaxwell ? 1.0 : 1.1781;
        const FLOAT_T ortVelocity = prd.velocity*fabsf(dot(rayDir, poly.N));
        const FLOAT_T hitEquiv = 1.0; //1.0*prd.orientationRatio; // hit=1.0 (only changed for lowflux mode)
        const FLOAT_T absEquiv = 1.0; //1.0*prd.orientationRatio; // hit=1.0 (only changed for lowflux mode)
#else
        const FLOAT_T velFactor = optixLaunchParams.simConstants.useMaxwell ? 1.0f : 1.1781f;
        const FLOAT_T ortVelocity = prd.velocity*fabsf(dot(rayDir, poly.N));
        const FLOAT_T hitEquiv = 1.0f; //1.0*prd.orientationRatio; // hit=1.0 (only changed for lowflux mode)
        const FLOAT_T absEquiv = 1.0f; //1.0*prd.orientationRatio; // hit=1.0 (only changed for lowflux mode)

#endif

#ifdef BOUND_CHECK
        if(counterIdx < 0 || counterIdx >= optixLaunchParams.simConstants.nbFacets * EXTRAFACETCOUNTERS){
            printf("facIndex %u >= %u is out of bounds\n", counterIdx, optixLaunchParams.simConstants.nbFacets * EXTRAFACETCOUNTERS);
        }
#endif
        increaseHitCounterAbsorp(optixLaunchParams.hitCounter[counterIdx], hitEquiv, absEquiv, ortVelocity, velFactor, prd.velocity);

#ifdef WITH_TEX
        if (poly.texProps.textureFlags & flowgpu::TEXTURE_FLAGS::countAbs){
            RecordAbsorpTexture(poly, prd, rayOrigin, rayDir);
        }
#endif
#ifdef WITH_PROF
        if (poly.profProps.profileType != flowgpu::PROFILE_FLAGS::noProfile){
            RecordAbsorpProfile(poly, prd, rayOrigin, rayDir);
        }
#endif

// terminate if exit has been called
#ifdef WITHDESORPEXIT
#ifdef PAYLOAD_DIRECT
        if(optixLaunchParams.perThreadData.currentMoleculeData[getWorkIndex()].hasToTerminate==1){
            optixLaunchParams.perThreadData.currentMoleculeData[getWorkIndex()].hasToTerminate=2;
        }
#else
        if(prd.hasToTerminate==1){
            prd.hasToTerminate=2;
        }
#endif //PAYLOAD_DIRECT
#endif
    }

    static __forceinline__ __device__
    void recordBounce(const unsigned int& bufferIndex, const unsigned int& counterIdx, const flowgpu::Polygon& poly, MolPRD& prd, const float3& rayDir, const float3 rayOrigin,
#ifdef RNG_BULKED
            const RN_T* randFloat, unsigned int& randInd, unsigned int& randOffset)
#else
            curandState_t* states)
#endif
    {
        // TODO: Save somewhere as a shared constant instead of repetively evaluating
#ifdef HIT64
        const FLOAT_T velFactor = optixLaunchParams.simConstants.useMaxwell ? 1.0 : 1.1781;
        FLOAT_T ortVelocity = prd.velocity*fabsf(dot(rayDir, poly.N));
        const FLOAT_T hitEquiv = 1.0; //1.0*prd.orientationRatio; // hit=1.0 (only changed for lowflux mode)
        const FLOAT_T absEquiv = 1.0; //1.0*prd.orientationRatio; // hit=1.0 (only changed for lowflux mode)
#else
        const FLOAT_T velFactor = optixLaunchParams.simConstants.useMaxwell ? 1.0f : 1.1781f;
        FLOAT_T ortVelocity = prd.velocity*fabsf(dot(rayDir, poly.N));
        const FLOAT_T hitEquiv = 1.0f; //1.0*prd.orientationRatio; // hit=1.0 (only changed for lowflux mode)
        const FLOAT_T absEquiv = 1.0f; //1.0*prd.orientationRatio; // hit=1.0 (only changed for lowflux mode)
#endif

#ifdef BOUND_CHECK
        if(counterIdx < 0 || counterIdx >= optixLaunchParams.simConstants.nbFacets * EXTRAFACETCOUNTERS){
            printf("facIndex %u >= %u is out of bounds\n", counterIdx, optixLaunchParams.simConstants.nbFacets * EXTRAFACETCOUNTERS);
        }
#endif
        // 1. Increment counters
        increaseHitCounterBounce(optixLaunchParams.hitCounter[counterIdx], hitEquiv, ortVelocity, velFactor, prd.velocity);

#ifdef WITH_TEX
        unsigned int texelIndex = 1e8;
        if (poly.texProps.textureFlags & flowgpu::TEXTURE_FLAGS::countRefl){
            texelIndex = RecordBounceTexture(poly, prd, rayOrigin, rayDir);
        }
#endif // WITH_TEX
#ifdef WITH_PROF
        unsigned int profileIndex = 1e8;
        //printf("BOUNCE? %d != %d == %d\n",(int)poly.profProps.profileType, (int)flowgpu::PROFILE_FLAGS::noProfile,(int)poly.profProps.profileType != flowgpu::PROFILE_FLAGS::noProfile);
        if (poly.profProps.profileType != flowgpu::PROFILE_FLAGS::noProfile){
            profileIndex = RecordBounceProfile(poly, prd, rayOrigin, rayDir);
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
        if(prd.facetHitSide == OPTIX_HIT_KIND_TRIANGLE_FRONT_FACE || prd.inSystem != TRANSPARENT_HIT) {
            //if (bufferIndex == 0 && poly.parentIndex == 102) printf("[%d] Front direction!\n", poly.parentIndex);
#ifdef RNG_BULKED
            prd.postHitDir = getNewDirection(prd, poly, randFloat, randInd, randOffset);
#else
            prd.postHitDir = getNewDirection(prd, poly, states);
#endif
        }
        else {
            /*if (bufferIndex == 0)
                printf("[%d] Back direction!\n", poly.parentIndex);*/
            //else if (bufferIndex == 0 && poly.parentIndex != 102) printf("[%d] Back direction!\n", poly.parentIndex);
#ifdef RNG_BULKED
            prd.postHitDir = getNewReverseDirection(prd, poly, randFloat, randInd, randOffset);
#else
            prd.postHitDir = getNewReverseDirection(prd, poly, states);
#endif
            //prd.inSystem = ACTIVE_PARTIC;
        }

        // 3. Increment counters for post bounce / outgoing particles
#ifdef HIT64
        ortVelocity = prd.velocity*fabs(dot(prd.postHitDir, poly.N));
#else
        ortVelocity = prd.velocity*fabsf(dot(prd.postHitDir, poly.N));
#endif

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
        ++optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].nbBounces;
        ++prd.nbBounces;
#endif
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

    extern "C" __global__ void __closesthit__molecule_polygon()
    {
#ifdef DEBUG
        //printf("--- pre %d data ---\n", threadIdx.x);
#endif

        const PolygonMeshSBTData &sbtData = *(const PolygonMeshSBTData*)optixGetSbtDataPointer();

        //const int   primID = optixGetPrimitiveIndex();
        const float3 ray_orig = optixGetWorldRayOrigin();
        const float3 ray_dir  = optixGetWorldRayDirection();
        const float  ray_t    = optixGetRayTmax();

        //const uint2 launchIndex = make_uint2(optixGetLaunchIndex());
        const unsigned int bufferIndex = getWorkIndex();

#ifdef DEBUGMISS
        const unsigned int missIndex = bufferIndex*NMISSES;
        optixLaunchParams.perThreadData.missBuffer[missIndex] = 0;
#endif
        const unsigned int facetHitKind = optixGetHitKind();

#ifdef PAYLOAD_DIRECT
        MolPRD myPrd = getMolPRD();
        MolPRD* prd = &myPrd;
#else
        MolPRD* prd = mergePointer(optixGetPayload_0(), optixGetPayload_1());
#endif

        if(optixGetPrimitiveIndex()==79 && ray_t > -0.00000032 && ray_t < -0.00000029)
            DEBUG_PRINT("[%d <- %d] / /got hit after all \n", optixGetPrimitiveIndex(), bufferIndex);


// self intersection
        if(prd->hitFacetId == optixGetPrimitiveIndex()){
#ifdef DEBUG
            //printf("[%d] source and goal facet equal %d : %8.6f,%8.6f,%8.6f -> %8.6f,%8.6f,%8.6f : [%.5e .. %.5e] (tri)\n",(blockDim.x * blockIdx.x + threadIdx.x), prd->hitFacetId, ray_orig.x,ray_orig.y,ray_orig.z,ray_dir.x,ray_dir.y,ray_dir.z, optixGetRayTmin(),ray_t);
/*
if(prd->inSystem == 4)
                printf("[%d] [backface] source and goal facet equal %d : %8.6f,%8.6f,%8.6f -> %8.6f,%8.6f,%8.6f : [%.5e .. %.5e] (tri)\n",(blockDim.x * blockIdx.x + threadIdx.x), prd->hitFacetId, ray_orig.x,ray_orig.y,ray_orig.z,ray_dir.x,ray_dir.y,ray_dir.z, optixGetRayTmin(),ray_t);
*/

#endif
            prd->facetHitSide = facetHitKind;
            if(facetHitKind == OPTIX_HIT_KIND_TRIANGLE_BACK_FACE) {
                prd->inSystem = SELF_INTERSECTION;
                /*if(bufferIndex==0)
                    printf("[%u] back self hit detected -> %u -> %u for n=%d\n", bufferIndex, prd->inSystem, sbtData.poly[prd->hitFacetId].parentIndex,prd->nbBounces);
            */}
            else{
                prd->inSystem = SELF_INTERSECTION;
            }
            const flowgpu::Polygon& poly  = sbtData.poly[prd->hitFacetId];


            prd->hitPos = ray_orig;
            prd->postHitDir = ray_dir;

            //prd->currentDepth = 0;
#ifdef PAYLOAD_DIRECT
            optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].hitPos = prd->hitPos;
            optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].postHitDir = prd->postHitDir;
            setMolPRD(myPrd);
#endif
            return;
        }

        prd->hitT = ray_t;
        prd->hitPos = ray_orig + ray_t * ray_dir;
        prd->hitFacetId = optixGetPrimitiveIndex();
        prd->facetHitSide = facetHitKind;
        //setMolPRD(prd);

        //TODO: only assume bounce for now
        //TODO: Need to account for post-bounce effects on same facet

        // first add facet hits
        //TODO: Check if counters can be used on threads instead of launch id
        //const unsigned int counterIdx = prd->hitFacetId+ix*optixLaunchParams.simConstants.nbFacets+iy*optixLaunchParams.simConstants.nbFacets*optixLaunchParams.simConstants.size.x;
        //const unsigned int counterIdx = prd->hitFacetId + (blockDim.x * blockIdx.x + threadIdx.x)*optixLaunchParams.simConstants.nbFacets;
        const unsigned int counterIdx = prd->hitFacetId + (bufferIndex%(EXTRAFACETCOUNTERS)) * optixLaunchParams.simConstants.nbFacets;

        //Register (orthogonal) velocity
        const flowgpu::Polygon& poly  = sbtData.poly[prd->hitFacetId];

#ifdef HIT64
        const FLOAT_T velFactor = optixLaunchParams.simConstants.useMaxwell ? 1.0 : 1.1781;
        const FLOAT_T ortVelocity = prd->velocity*fabsf(dot(ray_dir, poly.N));
        const FLOAT_T hitEquiv = 1.0; //1.0*prd->orientationRatio; // hit=1.0 (only changed for lowflux mode)
        const FLOAT_T absEquiv = 1.0; //1.0*prd->orientationRatio; // hit=1.0 (only changed for lowflux mode)
#else
        const float velFactor = optixLaunchParams.simConstants.useMaxwell ? 1.0f : 1.1781f;
        const float ortVelocity = prd->velocity*fabsf(dot(ray_dir, poly.N));
        const float hitEquiv = 1.0f; //1.0*prd->orientationRatio; // hit=1.0 (only changed for lowflux mode)
        const float absEquiv = 1.0f; //1.0*prd->orientationRatio; // hit=1.0 (only changed for lowflux mode)

#endif

#ifdef RNG_BULKED
        const RN_T* randFloat = optixLaunchParams.randomNumbers;
        unsigned int randInd = NB_RAND(optixLaunchParams.simConstants.maxDepth)*(bufferIndex);
        unsigned int randOffset = optixLaunchParams.perThreadData.randBufferOffset[bufferIndex];
#endif


        //TODO: AtomicAdd for smaller hitCounter structure, just per threadIdx
        if(poly.facProps.stickingFactor>=0.99999f || ((poly.facProps.stickingFactor > 0.0f) && (
#ifdef RNG_BULKED
        randFloat[(unsigned int)(randInd + randOffset++)]
#else
        generate_rand(optixLaunchParams.randomNumbers, bufferIndex)
#endif

        < (poly.facProps.stickingFactor)))){
            //--------------------------
            //-------- ABSORPTION ------
            //--------------------------

            recordAbsorption(counterIdx,poly,*prd, ray_dir, prd->hitPos);
#ifdef RNG_BULKED
            optixLaunchParams.perThreadData.randBufferOffset[bufferIndex] = randOffset;
#endif

            prd->inSystem = NEW_PARTICLE;
            prd->currentDepth = 0;
#ifdef GPUNBOUNCE
            /*optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].nbBounces =*/ prd->nbBounces = 0;
#ifdef PAYLOAD_DIRECT
            optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].nbBounces = 0;
#endif // PAYLOAD_DIRECT
#endif // GPUNBOUNCE

#ifdef PAYLOAD_DIRECT
            setMolPRD(myPrd);
#endif
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
#ifdef RNG_BULKED
            recordBounce(bufferIndex, counterIdx, poly, *prd, ray_dir, prd->hitPos, randFloat, randInd, randOffset);
            // Write temporary local variables back to shared memory
            optixLaunchParams.perThreadData.randBufferOffset[bufferIndex] = randOffset;
#else
            recordBounce(bufferIndex, counterIdx, poly, *prd, ray_dir, prd->hitPos, optixLaunchParams.randomNumbers);
#endif

            optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].rndDirection[0] = prd->rndDirection[0];
            optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].rndDirection[1] = prd->rndDirection[1];

#ifdef PAYLOAD_DIRECT
            optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].facetHitSide = prd->facetHitSide;
            optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].hitPos = prd->hitPos;
            optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].postHitDir = prd->postHitDir;
#endif

            // No recursion or end of recursion
            if(prd->currentDepth >= optixLaunchParams.simConstants.maxDepth){
                //printf("-- [%d] Max Bounces reached on depth %d, resetting for next call ... ----- \n",fbIndex, prd->currentDepth);
                prd->currentDepth = 0;
            }
            // WITH RECURSION
            else {
                ++prd->currentDepth;

#ifdef DEBUGPOS
                if(bufferIndex==0){
                    const unsigned int posIndexOff = optixLaunchParams.perThreadData.posOffsetBuffer_debug[(unsigned int)(bufferIndex)]++;
                    if(posIndexOff<NBPOSCOUNTS){
                        const unsigned int posIndex = bufferIndex*NBPOSCOUNTS+posIndexOff;
                        //printf("[%d] my pos is %d\n", (unsigned int)(ix+iy*optixLaunchParams.simConstants.size.x), posIndex);
                        optixLaunchParams.perThreadData.positionsBuffer_debug[posIndex] = prd->hitPos;
                    }
                }
#endif

                apply_offset(*prd, optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].hitPos);

#ifdef PAYLOAD_DIRECT
                int hi_vel = __double2hiint(prd->velocity);
                int lo_vel = __double2loint(prd->velocity);
#else
                // Pass the current payload registers through to the shadow ray.
                unsigned int p0 = optixGetPayload_0();
                unsigned int p1 = optixGetPayload_1();
#endif
                optixTrace(optixLaunchParams.traversable,
                           optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].hitPos,
                           optixLaunchParams.perThreadData.currentMoleculeData[bufferIndex].postHitDir,
                           0.f,//1e-4f,//0.f,    // tmin
                           1e20f,  // tmax
                           0.0f,   // rayTime
                           OptixVisibilityMask(255),
                           OPTIX_RAY_FLAG_DISABLE_ANYHIT
                           | OPTIX_RAY_FLAG_CULL_BACK_FACING_TRIANGLES,//OPTIX_RAY_FLAG_NONE,
                           RayType::RAY_TYPE_MOLECULE,             // SBT offset
                           RayType::RAY_TYPE_COUNT,               // SBT stride
                           RayType::RAY_TYPE_MOLECULE,             // missSBTIndex
#ifdef PAYLOAD_DIRECT
                        // direct
                           reinterpret_cast<unsigned int&>(prd->currentDepth),
                           reinterpret_cast<unsigned int&>(prd->inSystem),
                           reinterpret_cast<unsigned int&>(prd->hitFacetId),
                           reinterpret_cast<unsigned int&>(prd->hitT),
                           reinterpret_cast<unsigned int&>( hi_vel ),
                           reinterpret_cast<unsigned int&>( lo_vel ),
                           reinterpret_cast<unsigned int&>(prd->nbBounces),
                           reinterpret_cast<unsigned int&>(prd->facetHitSide)
#else
                        // ptr
                           p0, p1
#endif
                );

#ifdef PAYLOAD_DIRECT
                prd->velocity = __hiloint2double(hi_vel,lo_vel );
#endif
            }

#ifdef PAYLOAD_DIRECT
            setMolPRD(myPrd);
#endif
        } // end bounce
    }

    extern "C" __global__ void __anyhit__molecule()
    { /*! for this simple example, this will remain empty */ }

    //------------------------------------------------------------------------------
    // miss program that gets called for any ray that did not have a
    // valid intersection
    //
    // as with the anyhit/closest hit programs, in this example we only
    // need to have _some_ dummy function to set up a valid SBT
    // ------------------------------------------------------------------------------

#if defined(DEBUGMISS)

    extern "C" __device__ void is_pnp_miss(const int primID, float d, float u, float v) {
        // Fast method to check if a point is inside a polygon or not.
        // Works with convex and concave polys, orientation independent

        const PolygonMeshSBTData &sbtData = *(const PolygonMeshSBTData*)optixGetSbtDataPointer();
        const unsigned int fbIndex = getWorkIndex();

#ifdef BOUND_CHECK
        if(primID < 0 || primID >= optixLaunchParams.simConstants.nbFacets){
            printf("primID %u >= %u is out of bounds\n", primID, optixLaunchParams.simConstants.nbFacets);
#if defined(DEBUG)
            optixThrowException(7);
#endif
        }
#endif
        const flowgpu::Polygon& poly  = sbtData.poly[primID];

        const int nbSizeMinusOne = (int)poly.nbVertices - 1;
        const float2* polyPoints = sbtData.vertex2;

        float2 p;
        p.x = u;
        p.y = v;

        int i, j, c = 0;
        DEBUG_PRINT("[%d][%d] pnp! [%d : %d] -- %lf , %lf , %lf\n", primID,fbIndex, nbSizeMinusOne,poly.nbVertices, u, v, d);
#ifdef BOUND_CHECK
        if(poly.nbVertices < 0 || poly.nbVertices >= optixLaunchParams.simConstants.nbVertices){
            printf("poly.nbVertices %u >= %u is out of bounds\n", poly.nbVertices, optixLaunchParams.simConstants.nbVertices);
#if defined(DEBUG)
            optixThrowException(8);
#endif
        }
#endif
#ifdef BOUND_CHECK
        if(nbSizeMinusOne < 0 || poly.indexOffset + nbSizeMinusOne >= optixLaunchParams.simConstants.nbIndices * optixLaunchParams.simConstants.nbFacets){
            printf("poly.indexOffset %u >= %u [%d < 0]is out of bounds\n", poly.indexOffset, poly.indexOffset + nbSizeMinusOne, nbSizeMinusOne);
#if defined(DEBUG)
            optixThrowException(9);
#endif
        }
#endif
        for (i = 0, j = nbSizeMinusOne; i < poly.nbVertices; j = i++) {
            //DEBUG_PRINT("[%d] pnp! [%d : %d] --> %d : %d -- %lf , %lf , %lf\n", primID, nbSizeMinusOne,poly.nbVertices, poly.indexOffset + i, poly.indexOffset + j, u, v, d);
            const float2& p1 = polyPoints[poly.indexOffset + i];
            const float2& p2 = polyPoints[poly.indexOffset + j];
            if ((p1.y>v) != (p2.y>v)) {
#if defined(DEBUG)
            if(p2.y - p1.y == 0.0) optixThrowException(300);
#endif

                float slope = (v - p1.y) / (p2.y - p1.y);
                if(u < (p2.x - p1.x) * slope + p1.x)
                    c = !c;
            }
        }

        if(!c){
            DEBUG_PRINT("point in poly![%d] -- %lf , %lf , %lf\n", primID, u, v, d);
        }
    }

    extern "C" __device__ void is_polygon_miss_double(int primID, double d, double u, double v) {
        // Fast method to check if a point is inside a polygon or not.
        // Works with convex and concave polys, orientation independent

        const PolygonMeshSBTData &sbtData = *(const PolygonMeshSBTData*)optixGetSbtDataPointer();


        const flowgpu::Polygon& poly  = sbtData.poly[primID];

        const int nbSizeMinusOne = poly.nbVertices - 1;
        const double2* polyPoints = sbtData.vertex2x64;

        int n_updown = 0;
        int n_found = 0;

        double2 p;
        p.x = u;
        p.y = v;

        for (size_t j = 0; j < nbSizeMinusOne; j++) {
            const double2& p1 = polyPoints[poly.indexOffset + j];
            const double2& p2 = polyPoints[poly.indexOffset + j + 1];
/*
            if(primID==2 && optixGetLaunchIndex().x+optixGetLaunchIndex().y*optixLaunchParams.frame.size.x % 500 == 0)
                */
/*printf("[%d] -- %10.4f / %10.4f / %10.4f / %10.4f > %10.4f for Ray from %10.4f , %10.4f , %10.4f to %10.4f , %10.4f , %10.4f \n",
                       primID, det, u, v, d, ray_tmin, ray_orig.x,ray_orig.y,ray_orig.z,ray_dir.x,ray_dir.y,ray_dir.z);*//*

                printf("[%d] -- %10.4f , %10.4f -- %10.4f , %10.4f \n", sbtData.index[poly.indexOffset+j], p1.x, p1.y, p2.x, p2.y);
*/

            if (p.x<p1.x != p.x<p2.x) {
                double slope = (p2.y - p1.y) / (p2.x - p1.x);
                if ((slope * p.x - p.y) < (slope * p1.x - p1.y)) {
                    n_updown++;
                }
                else {
                    n_updown--;
                }
                n_found++;
            }
        }

/*
        if(primID==2 && optixGetLaunchIndex().x+optixGetLaunchIndex().y*optixLaunchParams.frame.size.x % 500 == 0)
            */
/*printf("[%d] -- %10.4f / %10.4f / %10.4f / %10.4f > %10.4f for Ray from %10.4f , %10.4f , %10.4f to %10.4f , %10.4f , %10.4f \n",
                   primID, det, u, v, d, ray_tmin, ray_orig.x,ray_orig.y,ray_orig.z,ray_dir.x,ray_dir.y,ray_dir.z);*//*

            printf("[%d]half -- found would be %d with %d and %d \n",nbSizeMinusOne, ((n_found / 2) & 1) ^ ((n_updown / 2) & 1), n_found, n_updown);
*/

        //Last point. Repeating code because it's the fastest and this function is heavily used
        const double2& p1 = polyPoints[poly.indexOffset + nbSizeMinusOne];
        const double2& p2 = polyPoints[poly.indexOffset + 0];
        if (p.x<p1.x != p.x<p2.x) {
            double slope = (p2.y - p1.y) / (p2.x - p1.x);
            if ((slope * p.x - p.y) < (slope * p1.x - p1.y)) {
                n_updown++;
            }
            else {
                n_updown--;
            }
            n_found++;
        }

        /*if(((n_found / 2u) & 1u) ^ ((n_updown / 2u) & 1u)){

        }
        else{

        }*/
        DEBUG_PRINT("inpoly_d[%d] -- [%d] %u ^ %u (%d / %d) -- %10.4f / %10.4f / %10.4f \n",
                    primID, (((n_found / 2u) & 1u) ^ ((n_updown / 2u) & 1u)), ((n_found / 2u) & 1u) , ((n_updown / 2u) & 1u), n_found, n_updown, u, v, d);

    }

    extern "C" __device__ void is_polygon_miss(int primID, float d, float u, float v) {
        // Fast method to check if a point is inside a polygon or not.
        // Works with convex and concave polys, orientation independent

        const PolygonMeshSBTData &sbtData = *(const PolygonMeshSBTData*)optixGetSbtDataPointer();

        const flowgpu::Polygon& poly  = sbtData.poly[primID];

        const int nbSizeMinusOne = poly.nbVertices - 1;
        const float2* polyPoints = sbtData.vertex2;

        int n_updown = 0;
        int n_found = 0;

        float2 p;
        p.x = u;
        p.y = v;

        for (size_t j = 0; j < nbSizeMinusOne; j++) {
            const float2& p1 = polyPoints[poly.indexOffset + j];
            const float2& p2 = polyPoints[poly.indexOffset + j + 1];
/*
            if(primID==2 && optixGetLaunchIndex().x+optixGetLaunchIndex().y*optixLaunchParams.frame.size.x % 500 == 0)
                */
/*printf("[%d] -- %10.4f / %10.4f / %10.4f / %10.4f > %10.4f for Ray from %10.4f , %10.4f , %10.4f to %10.4f , %10.4f , %10.4f \n",
                       primID, det, u, v, d, ray_tmin, ray_orig.x,ray_orig.y,ray_orig.z,ray_dir.x,ray_dir.y,ray_dir.z);*//*

                printf("[%d] -- %10.4f , %10.4f -- %10.4f , %10.4f \n", sbtData.index[poly.indexOffset+j], p1.x, p1.y, p2.x, p2.y);
*/

            if (p.x<p1.x != p.x<p2.x) {
                float slope = (p2.y - p1.y) / (p2.x - p1.x);
                if ((slope * p.x - p.y) < (slope * p1.x - p1.y)) {
                    n_updown++;
                }
                else {
                    n_updown--;
                }
                n_found++;
            }
        }

/*
        if(primID==2 && optixGetLaunchIndex().x+optixGetLaunchIndex().y*optixLaunchParams.frame.size.x % 500 == 0)
            */
/*printf("[%d] -- %10.4f / %10.4f / %10.4f / %10.4f > %10.4f for Ray from %10.4f , %10.4f , %10.4f to %10.4f , %10.4f , %10.4f \n",
                   primID, det, u, v, d, ray_tmin, ray_orig.x,ray_orig.y,ray_orig.z,ray_dir.x,ray_dir.y,ray_dir.z);*//*

            printf("[%d]half -- found would be %d with %d and %d \n",nbSizeMinusOne, ((n_found / 2) & 1) ^ ((n_updown / 2) & 1), n_found, n_updown);
*/

        //Last point. Repeating code because it's the fastest and this function is heavily used
        const float2& p1 = polyPoints[poly.indexOffset + nbSizeMinusOne];
        const float2& p2 = polyPoints[poly.indexOffset + 0];
        if (p.x<p1.x != p.x<p2.x) {
            float slope = (p2.y - p1.y) / (p2.x - p1.x);
            if ((slope * p.x - p.y) < (slope * p1.x - p1.y)) {
                n_updown++;
            }
            else {
                n_updown--;
            }
            n_found++;
        }

        /*if(((n_found / 2u) & 1u) ^ ((n_updown / 2u) & 1u)){

        }
        else{

        }*/
        DEBUG_PRINT("inpoly[%d] -- [%d] %u ^ %u (%d / %d) -- %10.4f / %10.4f / %10.4f \n",
                    primID, (((n_found / 2u) & 1u) ^ ((n_updown / 2u) & 1u)), ((n_found / 2u) & 1u) , ((n_updown / 2u) & 1u), n_found, n_updown, u, v, d);

    }

    // Parallelogram intersection based on the SDK optixWhitted example
    extern "C" __device__ void is_parallelogram_miss_double(uint32_t primID)
    {
        const PolygonMeshSBTData &sbtData = *(const PolygonMeshSBTData*)optixGetSbtDataPointer();

        //const int   primID = optixGetPrimitiveIndex();

        const double3 ray_orig = make_double3(optixGetWorldRayOrigin());
        //const double3 ray_dir = make_double3(optixGetWorldRayDirection().x,optixGetWorldRayDirection().y,optixGetWorldRayDirection().z);
        //ray_dir = make_double3(-1.0,-1.0,-1.0) * ray_dir;
        const double3 ray_dir = make_double3(-1.0f * optixGetWorldRayDirection());

        /*if(blockDim.x * blockIdx.x + threadIdx.x == 3 && ray_dir.x > -0.26 && ray_dir.x < -0.24)
            primID=162;*/

        const float ray_tmin = optixGetRayTmin(), ray_tmax = optixGetRayTmax();

        const flowgpu::Polygon& poly  = sbtData.poly[primID];
        const double det = dot(poly.Nuvx64, ray_dir);
        printf("missd[%d] intersection det = %12.10f\n", primID, det);
        if(det > 0.0f) {
            const double iDet = 1.0f / det;
            double3 intZ = ray_orig - poly.Ox64;

            const double u = iDet * DET33_ROW(intZ, poly.Vx64, ray_dir);

            const double v = iDet * DET33_ROW(poly.Ux64, intZ, ray_dir);

            const double d = iDet * dot(poly.Nuvx64, intZ);

            printf("missdist[%d]--> %lf * (%lf , %lf, %lf) = %lf\n", primID, iDet, intZ.x, intZ.y, intZ.z, d);
            printf("miss[float] ray %lf , %lf , %lf ---> %lf , %lf , %lf\n", ray_orig.x,
                   ray_orig.y, ray_orig.z, ray_dir.x,
                   ray_dir.y, ray_dir.z);

            DEBUG_PRINT("d[%d] intersection   u = %lf\n",primID, u);
            DEBUG_PRINT("d[%d] intersection   v = %lf\n",primID, v);
            DEBUG_PRINT("d[%d] intersection   d = %lf\n",primID, d);
            //is_pnp_miss(primID, static_cast<double>(d),static_cast<double>(u),static_cast<double>(v));

            double doubleEps = 1e-5;
            double newU = u;
            double newV = v;
            if (u >= 0.0 - doubleEps && u <= 1.0 + doubleEps && v >= 0.0 - doubleEps && v <= 1.0 + doubleEps) {
                DEBUG_PRINT("d[%d && %d] pre clamp\n",u >= 0.0 - doubleEps && u <= 1.0 + doubleEps,v >= 0.0 - doubleEps && v <= 1.0 + doubleEps);
            }
            if(newU >= 0.0 - doubleEps && newU <= 1.0 + doubleEps)
                newU = clamp(newU, 0.0, 1.0);
            if(newV >= 0.0 - doubleEps && newV <= 1.0 + doubleEps)
                newV = clamp(newV, 0.0, 1.0);
            if(newU >= 0.0 - doubleEps && newU <= 1.0 + doubleEps && (newV >= 0.0 - doubleEps && newV <= 1.0 + doubleEps))
                is_pnp_miss(primID, static_cast<double>(d),static_cast<double>(newU),static_cast<double>(newV));
            DEBUG_PRINT("d[%d] %e > %e\n",d > -doubleEps,d, doubleEps);
            doubleEps = 1e-4;
            if(newU >= 0.0 - doubleEps && newU <= 1.0 + doubleEps)
                newU = clamp(newU, 0.0, 1.0);
            if(newV >= 0.0 - doubleEps && newV <= 1.0 + doubleEps)
                newV = clamp(newV, 0.0, 1.0);
            if(newU >= 0.0 - doubleEps && newU <= 1.0 + doubleEps && (newV >= 0.0 - doubleEps && newV <= 1.0 + doubleEps))
                is_pnp_miss(primID, static_cast<double>(d),static_cast<double>(newU),static_cast<double>(newV));
            DEBUG_PRINT("d[%d] %e > %e\n",d > -doubleEps,d, doubleEps);
            doubleEps = 1e-3;
            if(newU >= 0.0 - doubleEps && newU <= 1.0 + doubleEps)
                newU = clamp(newU, 0.0, 1.0);
            if(newV >= 0.0 - doubleEps && newV <= 1.0 + doubleEps)
                newV = clamp(newV, 0.0, 1.0);
            if(newU >= 0.0 - doubleEps && newU <= 1.0 + doubleEps && (newV >= 0.0 - doubleEps && newV <= 1.0 + doubleEps))
                is_pnp_miss(primID, static_cast<double>(d),static_cast<double>(newU),static_cast<double>(newV));
            DEBUG_PRINT("d[%d] %e > %e\n",d > -doubleEps,d, doubleEps);

            newU = u > 0.5 ? u - doubleEps : u + doubleEps;
            newV = v > 0.5 ? v - doubleEps : v + doubleEps;

            /*if(newU >= 0.0 - doubleEps && newU <= 1.0 + doubleEps && (newV >= 0.0 - doubleEps && newV <= 1.0 + doubleEps))
                is_pnp_miss(primID, static_cast<double>(d),static_cast<double>(newU),static_cast<double>(newV));
*/
            DEBUG_PRINT("d[%d] %e > %e\n",d > -doubleEps,d, doubleEps);

            doubleEps = 1e-2;
            newU = u > 0.5 ? u - doubleEps : u + doubleEps;
            newV = v > 0.5 ? v - doubleEps : v + doubleEps;

            /*if(newU >= 0.0 - doubleEps && newU <= 1.0 + doubleEps && (newV >= 0.0 - doubleEps && newV <= 1.0 + doubleEps))
                is_pnp_miss(primID, static_cast<double>(d),static_cast<double>(newU),static_cast<double>(newV));
*/
            DEBUG_PRINT("d[%d] %e > %e\n",d > -doubleEps,d, doubleEps);

            //is_pnp_miss(primID, static_cast<float>(d),static_cast<float>(u),static_cast<float>(v));
        }
    }

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

#ifdef BOUND_CHECK
        if(primID < 0 || primID >= optixLaunchParams.simConstants.nbFacets){
            printf("primID %u >= %u is out of bounds\n", primID, optixLaunchParams.simConstants.nbFacets);
#if defined(DEBUG)
            optixThrowException(10);
#endif
        }
#endif
        const flowgpu::Polygon& poly  = sbtData.poly[primID];
        const float det = dot(poly.Nuv, ray_dir);
        DEBUG_PRINT("f[%d] intersection det = %12.10f\n", primID, det);
        if(det > 0.0f) {
            const float iDet = 1.0f / det;
            float3 intZ = ray_orig - poly.O;

            float u = iDet * DET33(intZ.x, poly.V.x, ray_dir.x,
                                         intZ.y, poly.V.y, ray_dir.y,
                                         intZ.z, poly.V.z, ray_dir.z);

            float v = iDet * DET33(poly.U.x, intZ.x, ray_dir.x,
                                         poly.U.y, intZ.y, ray_dir.y,
                                         poly.U.z, intZ.z, ray_dir.z);

            const float d = iDet * dot(poly.Nuv, intZ);

            DEBUG_PRINT("f[%d] intersection   u = %12.10f\n",primID, u);
            DEBUG_PRINT("f[%d] intersection   v = %12.10f\n",primID, v);
            DEBUG_PRINT("f[%d] intersection   d = %12.10f\n",primID, d);

            //is_polygon_miss(primID, d,u,v);
            is_pnp_miss(primID, d,u,v);

            float floatEps = 1e-5;
            float newU = u;
            float newV = v;
            if (u >= 0.0f - floatEps && u <= 1.0f + floatEps && v >= 0.0f - floatEps && v <= 1.0f + floatEps) {
                DEBUG_PRINT("[%d && %d] pre clamp\n",u >= 0.0f - floatEps && u <= 1.0f + floatEps,(v >= 0.0f - floatEps && v <= 1.0f + floatEps));
            }

            if(newU >= 0.0f - floatEps && newU <= 1.0f + floatEps)
                newU = clamp(newU, 0.0f, 1.0f);
            if(newV >= 0.0f - floatEps && newV <= 1.0f + floatEps)
                newV = clamp(newV, 0.0f, 1.0f);
            if(newU >= 0.0f - floatEps && newU <= 1.0f + floatEps && (newV >= 0.0f - floatEps && newV <= 1.0f + floatEps))
                is_pnp_miss(primID, static_cast<float>(d),static_cast<float>(newU),static_cast<float>(newV));
            DEBUG_PRINT("[%d] %e > %e\n",d > -floatEps,d, floatEps);
            floatEps = 1e-4;
            if(newU >= 0.0f - floatEps && newU <= 1.0f + floatEps)
                newU = clamp(newU, 0.0f, 1.0f);
            if(newV >= 0.0f - floatEps && newV <= 1.0f + floatEps)
                newV = clamp(newV, 0.0f, 1.0f);
            if(newU >= 0.0f - floatEps && newU <= 1.0f + floatEps && (newV >= 0.0f - floatEps && newV <= 1.0f + floatEps))
                is_pnp_miss(primID, static_cast<float>(d),static_cast<float>(newU),static_cast<float>(newV));
            DEBUG_PRINT("[%d] %e > %e\n",d > -floatEps,d, floatEps);
            floatEps = 1e-3;
            if(newU >= 0.0f - floatEps && newU <= 1.0f + floatEps)
                newU = clamp(newU, 0.0f, 1.0f);
            if(newV >= 0.0f - floatEps && newV <= 1.0f + floatEps)
                newV = clamp(newV, 0.0f, 1.0f);
            if(newU >= 0.0f - floatEps && newU <= 1.0f + floatEps && (newV >= 0.0f - floatEps && newV <= 1.0f + floatEps))
                is_pnp_miss(primID, static_cast<float>(d),static_cast<float>(newU),static_cast<float>(newV));
            DEBUG_PRINT("[%d] %e > %e\n",d > -floatEps,d, floatEps);

            if(newU >= 0.0f - floatEps && newU <= 1.0f + floatEps)
                newU = clamp(newU, u + floatEps, u - floatEps);
            if(newV >= 0.0f - floatEps && newV <= 1.0f + floatEps)
                newV = clamp(newV, v + floatEps, v - floatEps);
            if(newU >= 0.0f - floatEps && newU <= 1.0f + floatEps && (newV >= 0.0f - floatEps && newV <= 1.0f + floatEps))
                is_pnp_miss(primID, static_cast<float>(d),static_cast<float>(newU),static_cast<float>(newV));
            DEBUG_PRINT("[%d] %e > %e\n",d > -floatEps,d, floatEps);
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

#ifdef PAYLOAD_DIRECT
        MolPRD myPrd = getMolPRD();
        MolPRD* prd = &myPrd;
#else
        MolPRD* prd = mergePointer(optixGetPayload_0(), optixGetPayload_1());
#endif

#if defined(DEBUGMISS)
        const float3 ray_orig = optixGetWorldRayOrigin();
        const float3 ray_dir  = optixGetWorldRayDirection();
        const float  ray_t    = optixGetRayTmax();


        const unsigned int fbIndex = getWorkIndex();
        const unsigned int missIndex = fbIndex*NMISSES;
        const PolygonMeshSBTData &sbtData = *(const PolygonMeshSBTData*)optixGetSbtDataPointer();

#ifdef BOUND_CHECK
        if(prd->hitFacetId < 0 || prd->hitFacetId >= optixLaunchParams.simConstants.nbFacets){
            printf("prd->hitFacetId %u >= %u is out of bounds\n", prd->hitFacetId, optixLaunchParams.simConstants.nbFacets);
#if defined(DEBUG)
            optixThrowException(11);
#endif
        }
#endif
#if not defined(GPUNBOUNCE)
        DEBUG_PRINT("DEBMISS(%d) miss[%d -> %u -> %u][%d] "
               "(%12.10f , %12.10f , %12.10f) --> (%12.10f , %12.10f , %12.10f) = %e \n",
               prd->inSystem, fbIndex, prd->hitFacetId, sbtData.poly[prd->hitFacetId].parentIndex, missIndex,
               ray_orig.x, ray_orig.y , ray_orig.z , ray_dir.x, ray_dir.y , ray_dir.z, ray_t);
#else
        DEBUG_PRINT("[DEBMISS](%d , %d) miss[%d -> %u -> %u][%d] "
               "(%12.10f , %12.10f , %12.10f) --> (%12.10f , %12.10f , %12.10f) = %e \n",
               prd->inSystem, prd->nbBounces, fbIndex, prd->hitFacetId, sbtData.poly[prd->hitFacetId].parentIndex, missIndex,
               ray_orig.x, ray_orig.y , ray_orig.z , ray_dir.x, ray_dir.y , ray_dir.z, ray_t);
#endif

#ifdef BOUND_CHECK
        if(missIndex < 0 || missIndex >= NMISSES*optixLaunchParams.simConstants.size.x*optixLaunchParams.simConstants.size.y){
            printf("missIndex %u >= %u is out of bounds\n", missIndex, NMISSES*optixLaunchParams.simConstants.size.x*optixLaunchParams.simConstants.size.y);
#if defined(DEBUG)
            optixThrowException(12);
#endif
        }
#endif
#ifdef BOUND_CHECK
        if(missIndex + optixLaunchParams.perThreadData.missBuffer[missIndex] >= NMISSES*optixLaunchParams.simConstants.size.x*optixLaunchParams.simConstants.size.y){
            printf("missIndex+n %u >= %u is out of bounds\n", missIndex + optixLaunchParams.perThreadData.missBuffer[missIndex], NMISSES*optixLaunchParams.simConstants.size.x*optixLaunchParams.simConstants.size.y);
#if defined(DEBUG)
            optixThrowException(13);
#endif
        }
#endif

#ifdef WITHTRIANGLES
        const flowgpu::TriangleRayGenData* rayGenData = (flowgpu::TriangleRayGenData*) optixGetSbtDataPointer();
#else
        const flowgpu::PolygonRayGenData* rayGenData = (flowgpu::PolygonRayGenData*) optixGetSbtDataPointer();
#endif
        const flowgpu::Polygon& poly  = sbtData.poly[prd->hitFacetId];

        // recalculate origin in double precision only if new particle failed
        prd->rndOrigin[0] = optixLaunchParams.perThreadData.currentMoleculeData[fbIndex].rndOrigin[0];
        prd->rndOrigin[1] = optixLaunchParams.perThreadData.currentMoleculeData[fbIndex].rndOrigin[1];
        double3 dray_orig = (prd->inSystem == NEW_PARTICLE) ? (getOrigin_double(rayGenData, prd->hitFacetId, prd->rndOrigin[0],prd->rndOrigin[1])) : make_double3(ray_orig);
        //const double3 ray_dir = make_double3(optixGetWorldRayDirection().x,optixGetWorldRayDirection().y,optixGetWorldRayDirection().z);
        //ray_dir = make_double3(-1.0,-1.0,-1.0) * ray_dir;
        prd->rndDirection[0] = optixLaunchParams.perThreadData.currentMoleculeData[fbIndex].rndDirection[0];
        prd->rndDirection[1] = optixLaunchParams.perThreadData.currentMoleculeData[fbIndex].rndDirection[1];

        double3 dray_dir = /*-1.0 * */getDirection_double(*prd, poly, prd->rndDirection[0],prd->rndDirection[1]);
        double3 nU = poly.nUx64;
        double3 nV = poly.nVx64;
        double3 N = poly.Nx64;
        DEBUG_PRINT("Double Dir: (%lf , %lf) -> [%lf,%lf,%lf] + [%lf,%lf,%lf] + [%lf,%lf,%lf]\n", prd->rndDirection[0], prd->rndDirection[1] ,nU.x , nU.y , nU.z ,nV.x , nV.y , nV.z ,N.x , N.y , N.z);

        // if double ray is too similar, skip tracing new ray and finally declare a miss
        float orig_sim = (prd->inSystem == NEW_PARTICLE) ? (cosine_sim(make_float3(dray_orig), ray_orig)) : 1.0f;
        float dir_sim = cosine_sim(make_float3(dray_dir), ray_dir);

        if(orig_sim < 0.99f || dir_sim < 0.99f){

#if not defined(GPUNBOUNCE)
            DEBUG_PRINT("[DEBMISS-RELAUNCH](%d) miss[%d -> %u -> %u][%d] "
               "(%12.10f , %12.10f , %12.10f) --> (%12.10f , %12.10f , %12.10f) = %e \n",
               prd->inSystem, fbIndex, prd->hitFacetId, sbtData.poly[prd->hitFacetId].parentIndex, missIndex,
               ray_orig.x, ray_orig.y , ray_orig.z , ray_dir.x, ray_dir.y , ray_dir.z, ray_t);
#else
            DEBUG_PRINT("[DEBMISS-RELAUNCH](%d , %d) miss[%d -> %u -> %u][%d] "
                        "(%12.10f , %12.10f , %12.10f) --> (%12.10f , %12.10f , %12.10f) = %lf , %lf \n",
                        prd->inSystem, prd->nbBounces, fbIndex, prd->hitFacetId, sbtData.poly[prd->hitFacetId].parentIndex, missIndex,
                        dray_orig.x, dray_orig.y , dray_orig.z , dray_dir.x, dray_dir.y , dray_dir.z, orig_sim, dir_sim);
#endif

#ifdef PAYLOAD_DIRECT
            int hi_vel = __double2hiint(prd->velocity);
            int lo_vel = __double2loint(prd->velocity);
#else
            // Pass the current payload registers through to the shadow ray.
            unsigned int p0 = optixGetPayload_0();
            unsigned int p1 = optixGetPayload_1();
#endif
            prd->hitPos = make_float3(dray_orig);
            prd->postHitDir = make_float3(dray_dir);
            // reset misses and launch again

            optixLaunchParams.perThreadData.missBuffer[missIndex] = 0;
            optixLaunchParams.perThreadData.currentMoleculeData[fbIndex].hitPos = prd->hitPos;
            optixLaunchParams.perThreadData.currentMoleculeData[fbIndex].postHitDir = prd->postHitDir;

#if defined(PAYLOAD_DIRECT)
        setMolPRD(myPrd);
#endif

            DEBUG_PRINT("[START-DOUBLE-RELAUNCH](%d , %d) miss[%d -> %u -> %u][%d] "
                        "(%12.10f , %12.10f , %12.10f) --> (%12.10f , %12.10f , %12.10f) = %lf , %lf \n",
                        prd->inSystem, prd->nbBounces, fbIndex, prd->hitFacetId, sbtData.poly[prd->hitFacetId].parentIndex, missIndex,
                        prd->hitPos.x, prd->hitPos.y , prd->hitPos.z , prd->postHitDir.x, prd->postHitDir.y , prd->postHitDir.z, orig_sim, dir_sim);

            optixTrace(optixLaunchParams.traversable,
                       prd->hitPos,
                       prd->postHitDir,
                       0.f,//1e-4f,//0.f,    // tmin
                       1e20f,  // tmax
                       0.0f,   // rayTime
                       OptixVisibilityMask(255),
                       OPTIX_RAY_FLAG_DISABLE_ANYHIT
                       | OPTIX_RAY_FLAG_CULL_BACK_FACING_TRIANGLES,//OPTIX_RAY_FLAG_NONE,
                       RayType::RAY_TYPE_MOLECULE,             // SBT offset
                       RayType::RAY_TYPE_COUNT,               // SBT stride
                       RayType::RAY_TYPE_MOLECULE,             // missSBTIndex
#ifdef PAYLOAD_DIRECT
                    // direct
               reinterpret_cast<unsigned int&>(prd->currentDepth),
               reinterpret_cast<unsigned int&>(prd->inSystem),
               reinterpret_cast<unsigned int&>(prd->hitFacetId),
               reinterpret_cast<unsigned int&>(prd->hitT),
               reinterpret_cast<unsigned int&>( hi_vel ),
               reinterpret_cast<unsigned int&>( lo_vel ),
               reinterpret_cast<unsigned int&>(prd->nbBounces),
               reinterpret_cast<unsigned int&>(prd->facetHitSide)
#else
                    // ptr
                       p0, p1
#endif
            );
            DEBUG_PRINT("Finished double relaunch\n");
            return;
            DEBUG_PRINT("... but went further\n");
        }

/*for(int i=missIndex+1; i <= missIndex+optixLaunchParams.perThreadData.missBuffer[missIndex]; i++){
            DEBUG_PRINT("[MISSINDEX] miss[%d -> %d -> %d] at %d\n",
                   blockDim.x * blockIdx.x + threadIdx.x,
                   missIndex,
                   optixLaunchParams.perThreadData.missBuffer[missIndex],
                   optixLaunchParams.perThreadData.missBuffer[i]);

            is_parallelogram_miss(optixLaunchParams.perThreadData.missBuffer[i]);
            //is_parallelogram_miss_double(optixLaunchParams.perThreadData.missBuffer[i]);
        }*/


    /*if(optixLaunchParams.perThreadData.missBuffer[missIndex] == 0)
        return;*/

        optixLaunchParams.perThreadData.missBuffer[missIndex] = 0;
#endif //DEBUGMISS

/*optixLaunchParams.perThreadData.currentMoleculeData[fbIndex].hitPos = prd.hitPos;
optixLaunchParams.perThreadData.currentMoleculeData[fbIndex].postHitDir = prd.postHitDir;*/

#ifdef DEBUGLEAKPOS
        {
            const float3 ray_orig = optixGetWorldRayOrigin();
            const float3 ray_dir  = optixGetWorldRayDirection();
            const float  ray_t    = optixGetRayTmax();


            const unsigned int fbIndex = getWorkIndex();



/*#if not defined(GPUNBOUNCE)
            DEBUG_PRINT("[LEAKPOS](%d) miss[%d -> %d] (%12.10f , %12.10f , %12.10f) --> (%12.10f , %12.10f , %12.10f) = %e \n",
               prd.inSystem, fbIndex, prd.hitFacetId,
               ray_orig.x, ray_orig.y , ray_orig.z , ray_dir.x, ray_dir.y , ray_dir.z, ray_t);
#else
            DEBUG_PRINT("[LEAKPOS](%d , %d) miss[%d -> %d] (%12.10f , %12.10f , %12.10f) --> (%12.10f , %12.10f , %12.10f) = %e \n",
                   prd.inSystem, prd.nbBounces, fbIndex, prd.hitFacetId,
                   ray_orig.x, ray_orig.y , ray_orig.z , ray_dir.x, ray_dir.y , ray_dir.z, ray_t);
#endif*/

            const unsigned int posIndexOff = optixLaunchParams.perThreadData.leakPosOffsetBuffer_debug[fbIndex]++;
            if(posIndexOff<NBCOUNTS){
                const unsigned int posIndex = (fbIndex)*NBCOUNTS+posIndexOff;
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

        // Reset particle
        initParticle(myPrd);

        // terminate if exit has been called
#ifdef WITHDESORPEXIT
#ifdef PAYLOAD_DIRECT
        if(optixLaunchParams.perThreadData.currentMoleculeData[getWorkIndex()].hasToTerminate==1){
            optixLaunchParams.perThreadData.currentMoleculeData[getWorkIndex()].hasToTerminate=2;
        }
#else
        if(prd->hasToTerminate==1){
            prd->hasToTerminate=2;
        }
#endif //PAYLOAD_DIRECT
#endif
        /*}
        else{
            prd.inSystem = max(3,prd.inSystem+1);
        }*/
#if defined(PAYLOAD_DIRECT)
        setMolPRD(myPrd);
#endif
    }

} // ::flowgpu
