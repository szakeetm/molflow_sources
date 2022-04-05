// ======================================================================== //
// Copyright 2018-2019 Ingo Wald                                            //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

#pragma once

#include "OptixPolygon.h"
#include "GPUDefines.h"
#include "PerRayData.h"

#if !defined(RNG_BULKED)
#include <curand_kernel.h>
#endif

#ifdef DEBUGCOUNT
#define NCOUNTBINS 100
#define NBCOUNTS 100
#define DETLOW (-0.2f)
#define DETHIGH 0.2f
#define ULOW (-0.2f)
#define UHIGH 1.2f
#define VLOW (-0.2f)
#define VHIGH 1.2f
#endif
#ifdef DEBUGPOS
#define NBPOSCOUNTS 128
#endif
#ifdef DEBUGLEAKPOS
#define NBCOUNTS 3
#endif
namespace flowgpu {

    struct CuFacetHitCounter64;
    struct CuFacetHitCounter32;

    // Texel typedef
#ifdef HIT64
    using CuFacetHitCounter = CuFacetHitCounter64;
#else
    using CuFacetHitCounter = CuFacetHitCounter32;
#endif

    struct TriangleMeshSBTData {
        float3 *vertex;
        int3 *index;
        float2* texcoord;
        flowgpu::Polygon *poly;
    };

    struct TriangleRayGenData {
        float3 *vertex;
        double3 *vertex64;
        int3 *index;
        flowgpu::Polygon *poly;

        // -- data for launch parameters --
        // --------------------------------
        // probability for facet selection
        float2* facetProbabilities;
        //double* prob_lower;
        //double* prob_upper;
    };

    struct PolygonRayGenData {
        float3 *vertex;
        float2 *vertex2;
        double2 *vertex2x64;
        uint32_t *index;
        flowgpu::Polygon *poly;

        // -- data for launch parameters --
        // --------------------------------
        // probability for facet selection
        float2* facetProbabilities;
        //double* prob_lower;
        //double* prob_upper;
        // CFD for velocity calculation (temperature, v-bin)
        float* cdfs;
        //double *cdf_x; //index of map
        //double *cdf_val; //value of map
    };

    struct PolygonMeshSBTData {
        float3 *vertex;
        float2 *vertex2;
        double2 *vertex2x64;
        uint32_t *index;
        flowgpu::Polygon *poly;
    };

    enum MOLSTATUS : uint8_t {
        NEW_PARTICLE = 0u,
        ACTIVE_PARTICLE,
        SELF_INTERSECTION,
        TRANSPARENT_HIT,
        ACTIVE_BACK_HIT,
        SELF_BACK_HIT,
        TRANSPARENT_BACK_HIT,
        END_STATUS
    };

    // Globalhitbuffer
    /*FacetHitBuffer globalHits;               // Global counts (as if the whole geometry was one extra facet)
    size_t hitCacheSize;              // Number of valid hits in cache
    size_t lastHitIndex;					//Index of last recorded hit in gHits (turns over when reaches HITCACHESIZE)
    HIT hitCache[HITCACHESIZE];       // Hit history
    LEAK leakCache[LEAKCACHESIZE];      // Leak history
    size_t  lastLeakIndex;		  //Index of last recorded leak in gHits (turns over when reaches LEAKCACHESIZE)
    size_t  leakCacheSize;        //Number of valid leaks in the cache
    size_t  nbLeakTotal;*/

    // Facethitbuffer

    struct CuFacetHitCounter64{
        CuFacetHitCounter64() :
                nbDesorbed(0) , nbMCHit(0) , nbHitEquiv(0.0) ,
                nbAbsEquiv(0.0) , sum_1_per_ort_velocity(0.0) ,
                sum_1_per_velocity(0.0) , sum_v_ort(0.0)
        {};
        // Counts
        uint64_cu nbMCHit;                   // Number of hits
        uint64_cu nbDesorbed;                // Number of desorbed molec
        double nbAbsEquiv;                   // Equivalent number of absorbed molecules
        double nbHitEquiv;                   //Equivalent number of hits, used for low-flux impingement rate and density calculation
        double sum_1_per_ort_velocity;       // sum of reciprocials of orthogonal velocity components, used to determine the density, regardless of facet orientation
        double sum_1_per_velocity;           //For average molecule speed calculation
        double sum_v_ort;                    // sum of orthogonal speeds of incident velocities, used to determine the pressure
    };

    //TODO: Molflow counter
    struct CuFacetHitCounter32{
        CuFacetHitCounter32() :
        nbDesorbed(0) , nbMCHit(0) , nbHitEquiv(0.0f) ,
        nbAbsEquiv(0.0f) , sum_1_per_ort_velocity(0.0f) ,
        sum_1_per_velocity(0.0f) , sum_v_ort(0.0f)
        {};
        // Counts
        uint32_t nbMCHit;                   // Number of hits
        uint32_t nbDesorbed;                // Number of desorbed molec
        float nbAbsEquiv;                   // Equivalent number of absorbed molecules
        float nbHitEquiv;                   //Equivalent number of hits, used for low-flux impingement rate and density calculation
        float sum_1_per_ort_velocity;       // sum of reciprocials of orthogonal velocity components, used to determine the density, regardless of facet orientation
        float sum_1_per_velocity;           //For average molecule speed calculation
        float sum_v_ort;                    // sum of orthogonal speeds of incident velocities, used to determine the pressure
    };

    struct LaunchParams
    {

        struct {
            uint32_t *missCounter;
            flowgpu::FacetTexture *facetTextures;
            float *texelInc;
            flowgpu::Texel *texels;
            flowgpu::Texel *profileSlices;

            // CFD for velocity calculation (temperature, v-bin)
            float* cdfs1;
            float* cdfs2;
            //double *cdf_x; //index of map
            //double *cdf_val; //value of map
        } sharedData;

        struct {
            MolPRD *currentMoleculeData;
            uint32_t* randBufferOffset; /*! offset to remember which random number is next */
#ifdef DEBUGPOS
            uint32_t* posOffsetBuffer_debug;
            float3* positionsBuffer_debug;
            uint16_t* positionsType_debug;
#endif
#ifdef DEBUGLEAKPOS
            uint32_t* leakPosOffsetBuffer_debug;
            float3* leakPositionsBuffer_debug;
            float3* leakDirectionsBuffer_debug;
#endif
#ifdef DEBUGMISS
            uint32_t* missBuffer; //first value is the amount of primitives N, followed by N primitive IDs
#endif
        } perThreadData;

        struct {
            //uint32_t *colorBuffer;
            uint2     size;
            uint32_t nbFacets;
            uint32_t nbIndices;
            uint32_t nbVertices;

            // globalSettings
            bool useMaxwell;
            float gasMass;

            uint32_t maxDepth; // for recursion
            float scene_epsilon; // to prevent self intersection (currently unused due to other techniques)
            uint32_t nbRandNumbersPerThread;
#ifdef BOUND_CHECK
            uint32_t nbTexel;
            uint32_t nbProfSlices;
#endif
            float offset_center_magnitude;
            float offset_normal_magnitude;
        } simConstants;

#if defined(RNG_BULKED)
        RN_T* randomNumbers;
#else
        curandState_t* randomNumbers;
#endif
        CuFacetHitCounter* hitCounter;

#ifdef DEBUGCOUNT
        struct {
            uint32_t* detCount;
            uint32_t* uCount;
            uint32_t* vCount;
        } debugCounter;
#endif
        OptixTraversableHandle traversable;
    };



} // ::flowgpu
