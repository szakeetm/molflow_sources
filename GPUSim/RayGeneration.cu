// Created by pbahr

#include <optix_device.h>
#include <math_constants.h>

#include "LaunchParams.h"
#include "GPUDefines.h" // for NB_RAND
#include "jetbrains_indexing.h"

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

    // for this simple example, we have a single ray type
    enum { SURFACE_RAY_TYPE=0, RAY_TYPE_COUNT };

    static __forceinline__ __device__
    void  packPointer( void* ptr, uint32_t& i0, uint32_t& i1 )
    {
        const uint64_t uptr = reinterpret_cast<uint64_t>( ptr );
        i0 = uptr >> 32;
        i1 = uptr & 0x00000000ffffffff;
    }

    /*//------------------------------------------------------------------------------
    // ray gen program - the actual rendering happens in here
    //------------------------------------------------------------------------------
    extern "C" __global__ void __raygen__renderFrame()
    {
        // compute a test pattern based on pixel ID
        const int ix = optixGetLaunchIndex().x;
        const int iy = optixGetLaunchIndex().y;

        const auto &camera = optixLaunchParams.camera;

        // normalized screen plane position, in [0,1]^2
        const vec2f screen(vec2f(ix+.5f,iy+.5f)
                           / vec2f(optixLaunchParams.simConstants.size));

        // generate ray direction
        vec3f rayDir = normalize(camera.direction
                                 + (screen.x - 0.5f) * camera.horizontal
                                 + (screen.y - 0.5f) * camera.vertical);

        // our per-ray data for this example. what we initialize it to
        // won't matter, since this value will be overwritten by either
        // the miss or hit program, anyway
        vec3f pixelColorPRD = camera.position;
        vec3f pixelRay = rayDir;

        // the values we store the PRD pointer in:
        uint32_t u0, u1, u2, u3;
        packPointer( &pixelColorPRD, u0, u1 );
        packPointer( &pixelRay, u2, u3 );

        optixTrace(optixLaunchParams.traversable,
                   camera.position,
                   rayDir,
                   0.f,    // tmin
                   1e20f,  // tmax
                   0.0f,   // rayTime
                   OptixVisibilityMask( 255 ),
                   OPTIX_RAY_FLAG_DISABLE_ANYHIT,//OPTIX_RAY_FLAG_NONE,
                   SURFACE_RAY_TYPE,             // SBT offset
                   RAY_TYPE_COUNT,               // SBT stride
                   SURFACE_RAY_TYPE,             // missSBTIndex
                   u0, u1 , u2, u3);

        const int r = int(255.99f*pixelColorPRD.x);
        const int g = int(255.99f*pixelColorPRD.y);
        const int b = int(255.99f*pixelColorPRD.z);

        // convert to 32-bit rgba value (we explicitly set alpha to 0xff
        // to make stb_image_write happy ...
        const uint32_t rgba = 0xff000000
                              | (r<<0) | (g<<8) | (b<<16);

        // and write to frame buffer ...
        const uint32_t fbIndex = ix+iy*optixLaunchParams.simConstants.size.x;
        //optixLaunchParams.frame.colorBuffer[fbIndex] = rgba;

        //pixelColorPRD = camera.position;
        //pixelRay = rayDir;

        //optixLaunchParams.simConstants.dir[fbIndex] = pixelRay;
        //optixLaunchParams.simConstants.origin[fbIndex] = pixelColorPRD;
    }*/

    static __forceinline__ __device__
    void *getRayOriginPolygon( cuuint32_t i0, cuuint32_t i1 )
    {
        const uint64_t uptr = static_cast<uint64_t>( i0 ) << 32 | i1;
        void*           ptr = reinterpret_cast<void*>( uptr );
        return ptr;
    }

    static __forceinline__ __device__
    void *getRayOriginTriangle( cuuint32_t i0, cuuint32_t i1 )
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

    extern "C" __device__ int point_in_polygon(float u, float v, const Polygon& poly) {
        // Fast method to check if a point is inside a polygon or not.
        // Works with convex and concave polys, orientation independent
        const PolygonRayGenData* rayGenData = (PolygonRayGenData*) optixGetSbtDataPointer();
        //const int   primID = optixGetPrimitiveIndex();

        //const Polygon& poly  = sbtData.poly[primID];

        const int nbSizeMinusOne = poly.nbVertices - 1;
        const vec2f* polyPoints = rayGenData->vertex2;

        int n_updown = 0;
        int n_found = 0;

        vec2f p;
        p.u = u;
        p.v = v;

        //printf("[%d] poly check with %d vert at poly offset %d for p %4.2f %4.2f\n",threadIdx.x,poly.nbVertices,poly.vertOffset,p.u,p.v);
        //printf("[%d] poly first vec2 %4.2f , %4.2f\n",threadIdx.x,rayGenData->vertex2[0].u,rayGenData->vertex2[0].v);
        //printf("[%d] poly second vec2 %4.2f , %4.2f\n",threadIdx.x,rayGenData->vertex2[1].u,rayGenData->vertex2[1].v);

#ifdef BOUND_CHECK
        if(poly.vertOffset < 0 || poly.vertOffset+poly.nbVertices >= optixLaunchParams.simConstants.nbVertices){
            printf("[%d] vertOffset %u -- %u >= %u is out of bounds\n", poly.vertOffset, poly.vertOffset+poly.nbVertices, optixLaunchParams.simConstants.nbVertices);
        }
#endif
        for (int j = 0; j < nbSizeMinusOne; j++) {
            const vec2f& p1 = polyPoints[poly.vertOffset+j];
            const vec2f& p2 = polyPoints[poly.vertOffset+j+1];

            //printf("[%d / %u] p1 %4.2f %4.2f / p2 %4.2f %4.2f \n",threadIdx.x,j,p1.u,p1.v,p2.u,p2.v);

            if (p.u<p1.u != p.u<p2.u) {
                float slope = (p2.v - p1.v) / (p2.u - p1.u);
                //printf("[%d / %u] slope %4.2f --> %4.2f < %4.2f \n",threadIdx.x,j,slope, (slope * p.u - p.v), (slope * p1.u - p1.v));

                if ((slope * p.u - p.v) < (slope * p1.u - p1.v)) {
                    n_updown++;
                }
                else {
                    n_updown--;
                }
                n_found++;
            }
        }

        //Last point. Repeating code because it's the fastest and this function is heavily used
        const vec2f& p1 = polyPoints[poly.vertOffset+nbSizeMinusOne];
        const vec2f& p2 = polyPoints[poly.vertOffset+0];
        if (p.u<p1.u != p.u<p2.u) {
            //printf("[%d / %u] p1 %4.2f %4.2f / p2 %4.2f %4.2f \n",threadIdx.x,0,p1.u,p1.v,p2.u,p2.v);

            float slope = (p2.v - p1.v) / (p2.u - p1.u);
            //printf("[%d / %u] slope %4.2f --> %4.2f < %4.2f \n",threadIdx.x,0,slope, (slope * p.u - p.v), (slope * p1.u - p1.v));

            if ((slope * p.u - p.v) < (slope * p1.u - p1.v)) {
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
                         const float* randFloat, cuuint32_t& randInd, cuuint32_t& randOffset)
    {
#ifdef BOUND_CHECK
        if(randInd + randOffset < 0 || randInd + randOffset >= optixLaunchParams.simConstants.nbRandNumbersPerThread*optixLaunchParams.simConstants.size.x*optixLaunchParams.simConstants.size.y){
            printf("randInd %u is out of bounds\n", randInd + randOffset);
        }
#endif
        float facetRnd = randFloat[(cuuint32_t)(randInd + randOffset++)];

        int facIndex = 0;
        bool found = false;
        do{
#ifdef BOUND_CHECK
            if(facIndex < 0 || facIndex >= optixLaunchParams.simConstants.nbFacets){
                printf("facIndex %u >= %u is out of bounds\n", facIndex, optixLaunchParams.simConstants.nbFacets);
                printf("found %u with last prob %4.2f > facetRnd %4.2f \n", (rayGenData->facetProbabilities[facIndex].t >= facetRnd), optixLaunchParams.simConstants.nbFacets,
                       rayGenData->facetProbabilities[facIndex].t,facetRnd);
            }
#endif

            found = (rayGenData->facetProbabilities[facIndex].t >= facetRnd); // find probability interval rand lies in
            if(!found)
                facIndex++;

            /*if(found){
                if(facIndex==0)
                    printf("[%d] %u/%u with last prob %8.6f > facetRnd %8.6f and rand index %u+%u\n", fbIndex, facIndex, optixLaunchParams.simConstants.nbFacets,
                           rayGenData->facetProbabilities[facIndex].t,facetRnd, randInd,randOffset-1);
                else
                    printf("found %u/%u with last prob %8.6f > facetRnd %8.6f \n", facIndex, optixLaunchParams.simConstants.nbFacets,
                       rayGenData->facetProbabilities[facIndex].t,facetRnd);
            }*/
        }while(!found && facIndex<optixLaunchParams.simConstants.nbFacets);

#ifdef BOUND_CHECK
        if(facIndex < 0 || facIndex >= optixLaunchParams.simConstants.nbFacets){
            printf("Post facIndex %u >= %u is out of bounds\n", facIndex, optixLaunchParams.simConstants.nbFacets);
            printf("Post found %u with last prob %4.2f > facetRnd %4.2f \n", found, optixLaunchParams.simConstants.nbFacets,
                   rayGenData->facetProbabilities[facIndex].t,facetRnd);

        }
#endif

        hitData.hitFacetId = facIndex;
        return facIndex;

    }

    static __forceinline__ __device__
    vec3f getNewOrigin(MolPRD& hitData,
#ifdef WITHTRIANGLES
        const TriangleRayGenData* rayGenData,
#else
        const PolygonRayGenData* rayGenData,
#endif
        const int& facIndex,
        const float* randFloat, cuuint32_t& randInd, cuuint32_t& randOffset)
    {

#ifdef WITHTRIANGLES

        #ifdef BOUND_CHECK
        if(randInd + randOffset < 0 || randInd + randOffset >= optixLaunchParams.simConstants.nbRandNumbersPerThread*optixLaunchParams.simConstants.size.x*optixLaunchParams.simConstants.size.y){
            printf("randInd %u is out of bounds\n", randInd + randOffset);
        }
        #endif

        float r1_sqrt = sqrtf(randFloat[(cuuint32_t)(randInd + randOffset++)]);

        #ifdef BOUND_CHECK
        if(randInd + randOffset < 0 || randInd + randOffset >= optixLaunchParams.simConstants.nbRandNumbersPerThread*optixLaunchParams.simConstants.size.x*optixLaunchParams.simConstants.size.y){
            printf("randInd %u is out of bounds\n", randInd + randOffset);
        }
        #endif

            float r2 = randFloat[(cuuint32_t)(randInd + randOffset++)];


        #ifdef BOUND_CHECK
        if(rayGenData->index[facIndex].x < 0 || rayGenData->index[facIndex].x >= optixLaunchParams.simConstants.nbVertices
        || rayGenData->index[facIndex].y < 0 || rayGenData->index[facIndex].y >= optixLaunchParams.simConstants.nbVertices
        || rayGenData->index[facIndex].z < 0 || rayGenData->index[facIndex].z >= optixLaunchParams.simConstants.nbVertices){
            printf("rayGenData->index[facIndex] %u,%u,%u >= %u is out of bounds\n", rayGenData->index[facIndex].x, rayGenData->index[facIndex].y, rayGenData->index[facIndex].z, optixLaunchParams.simConstants.nbVertices);
        }
        #endif

            vec3f vertA = rayGenData->vertex[rayGenData->index[facIndex].x];
            vec3f vertB = rayGenData->vertex[rayGenData->index[facIndex].y];
            vec3f vertC = rayGenData->vertex[rayGenData->index[facIndex].z];

            //rayOrigin =
            return (1-r1_sqrt) * vertA + r1_sqrt * (1 - r2) * vertB + r1_sqrt * r2 * vertC; //rayGenData

            /*printf("--- [%d..%d] triangle origin %8.6f , %8.6f , %8.6f - %4.2f ---\n", threadIdx.x,facIndex, rayOrigin.x, rayOrigin.y, rayOrigin.z);
            printf("--- [%d..%d/0] triangle vert %8.6f , %8.6f , %8.6f - %4.2f ---\n", threadIdx.x,facIndex,vertA.x, vertA.y, vertA.z);
            printf("--- [%d..%d/1] triangle vert %8.6f , %8.6f , %8.6f - %4.2f ---\n", threadIdx.x,facIndex,vertB.x, vertB.y, vertB.z);
            printf("--- [%d..%d/2] triangle vert %8.6f , %8.6f , %8.6f - %4.2f ---\n", threadIdx.x,facIndex,vertC.x, vertC.y, vertC.z);*/
#else
        // start position of particle (U,V) -> (x,y,z)
        float uDir = 0.0f, vDir = 0.0f;

        short isInPoly = 0;
        cuuint32_t loopNb = 0;

        while(!isInPoly && loopNb < NB_INPOLYCHECKS){
#ifdef BOUND_CHECK
            if(randInd + randOffset < 0 || randInd + randOffset >= optixLaunchParams.simConstants.nbRandNumbersPerThread*optixLaunchParams.simConstants.size.x*optixLaunchParams.simConstants.size.y){
                printf("randInd %u is out of bounds\n", randInd + randOffset);
            }
#endif
            uDir = randFloat[(cuuint32_t)(randInd + randOffset++)]; // vec3f(randFloat[randInd++],randFloat[randInd++],randFloat[randInd++])

#ifdef BOUND_CHECK
            if(randInd + randOffset < 0 || randInd + randOffset >= optixLaunchParams.simConstants.nbRandNumbersPerThread*optixLaunchParams.simConstants.size.x*optixLaunchParams.simConstants.size.y){
                printf("randInd %u is out of bounds\n", randInd + randOffset);
            }
#endif
            vDir = randFloat[(cuuint32_t)(randInd + randOffset++)];

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
#endif WITHTRIANGLES
    }


    static __forceinline__ __device__
    vec3f getNewDirection(MolPRD& hitData, Polygon& poly,
            const float* randFloat, cuuint32_t& randInd, cuuint32_t& randOffset)
    {

        // generate ray direction
#ifdef BOUND_CHECK
        if(randInd + randOffset < 0 || randInd + randOffset >= optixLaunchParams.simConstants.nbRandNumbersPerThread*optixLaunchParams.simConstants.size.x*optixLaunchParams.simConstants.size.y){
            printf("randInd %u is out of bounds\n", randInd + randOffset);
        }
#endif
        const float theta = acosf(sqrtf(randFloat[(cuuint32_t)(randInd + randOffset++)]));


#ifdef BOUND_CHECK
        if(randInd + randOffset < 0 || randInd + randOffset >= optixLaunchParams.simConstants.nbRandNumbersPerThread*optixLaunchParams.simConstants.size.x*optixLaunchParams.simConstants.size.y){
            printf("randInd %u is out of bounds\n", randInd + randOffset);
        }
#endif
        const float phi = randFloat[(cuuint32_t)(randInd + randOffset++)] * 2.0 * CUDART_PI_F;


        const float u = sinf(theta)*cosf(phi);
        const float v = sinf(theta)*sinf(phi);
        const float n = cosf(theta);

        /*vec3f nU = rayGenData->poly[facIndex].nU;
        vec3f nV = rayGenData->poly[facIndex].nV;
        vec3f N = rayGenData->poly[facIndex].N;*/
        const vec3f nU = poly.nU;
        const vec3f nV = poly.nV;
        const vec3f N = poly.N;

        //rayDir = u*nU + v*nV + n*N;
        return u*nU + v*nV + n*N;
        /*if (rayDir.x != 0.0) rayDir.x = 1.0 / rayDir.x;
        if (rayDir.y != 0.0) rayDir.y = 1.0 / rayDir.y;
        if (rayDir.z != 0.0) rayDir.z = 1.0 / rayDir.z;
        rayDir = vec3f(-1.0f,-1.0f,-1.0f) * rayDir;*/
    }

    static __forceinline__ __device__
    void initMoleculeInSystem(const cuuint32_t thrInd, MolPRD& hitData, vec3f& rayDir , vec3f& rayOrigin)
    {
        hitData = optixLaunchParams.perThreadData.currentMoleculeData[thrInd];
        rayDir = hitData.postHitDir;
        rayOrigin = hitData.hitPos;

#ifdef DEBUG
        //printf("--- in[%d] %4.2f , %4.2f , %4.2f - %4.2f , %4.2f , %4.2f ---\n", blockDim.x * blockIdx.x + threadIdx.x,optixLaunchParams.perThreadData.currentMoleculeData[fbIndex].hitPos.x, optixLaunchParams.perThreadData.currentMoleculeData[fbIndex].hitPos.y, optixLaunchParams.perThreadData.currentMoleculeData[fbIndex].hitPos.z,optixLaunchParams.perThreadData.currentMoleculeData[fbIndex].postHitDir.x, optixLaunchParams.perThreadData.currentMoleculeData[fbIndex].postHitDir.y, optixLaunchParams.perThreadData.currentMoleculeData[fbIndex].postHitDir.z);
#endif
    }


    // Calculate new direction, but keep old position for molecule that was stuck due to self-intersection
    static __forceinline__ __device__
    void initMoleculePostSelfIntersection(const cuuint32_t thrInd, MolPRD& hitData, vec3f& rayDir , vec3f& rayOrigin)
    {
        hitData = optixLaunchParams.perThreadData.currentMoleculeData[thrInd];

        float* randFloat = optixLaunchParams.randomNumbers;
        cuuint32_t randInd = optixLaunchParams.simConstants.nbRandNumbersPerThread*(thrInd);
        cuuint32_t randOffset = optixLaunchParams.perThreadData.randBufferOffset[thrInd];

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
        optixLaunchParams.perThreadData.randBufferOffset[thrInd] = randOffset; // set value back to buffer

    }

    // --------------------------------------
    // increase facet counters for desorption
    // --------------------------------------
    static __forceinline__ __device__
    void increaseHitCounterDesorption(CuFacetHitCounter& hitCounter, MolPRD& hitData, vec3f& rayDir, vec3f& polyNormal)
    {
        double velFactor = optixLaunchParams.simConstants.useMaxwell ? 1.0f : 1.1781f; // TODO: Save somewhere as a shared constant instead of repetively evaluating

        // TODO: fix velocity etc
        const float velocity = hitData.velocity;
        const float ortVelocity = velocity*abs(dot(rayDir, polyNormal));
        //const float hitEquiv = 0.0f; //1.0*prd.orientationRatio; // hit=1.0 (only changed for lowflux mode)
        //const float absEquiv = 0.0f; //1.0*prd.orientationRatio; // hit=1.0 (only changed for lowflux mode)
        //atomicAdd(&optixLaunchParams.hitCounter[counterIdx].nbMCHit,1);
        //++optixLaunchParams.hitCounter[counterIdx].nbMCHit;
        //const cuuint32_t counterIdx = facIndex + ((thrInd)%(CORESPERMP*WARPSCHEDULERS))*optixLaunchParams.simConstants.nbFacets;

        hitCounter.nbHitEquiv += 0.0f;
        hitCounter.nbDesorbed += 1;
        hitCounter.nbAbsEquiv += 0.0f; //static_cast<double>(absorb)*sHandle->currentParticle.oriRatio;
        hitCounter.sum_1_per_ort_velocity += (2.0f / ortVelocity);//prd.oriRatio * sum_1_per_v;
        hitCounter.sum_v_ort += velFactor*ortVelocity;//(sHandle->wp.useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity; //prd.oriRatio * sum_v_ort;
        hitCounter.sum_1_per_velocity += 0.0f / velocity;//(hitEquiv + static_cast<double>(desorb)) / prd.velocity;
    }
    
    // Calculate new direction, but keep old position
    static __forceinline__ __device__
    void initMoleculeFromStart(const cuuint32_t thrInd, MolPRD& hitData, vec3f& rayDir , vec3f& rayOrigin)
    {
        hitData.hitPos = vec3f(-999.9f,-999.9f,-999.9f);
        hitData.postHitDir = vec3f(-999.9f,-999.9f,-999.9f);
        hitData.hitT = -999.0f;
        hitData.velocity = -999.0f;
        hitData.currentDepth = 0;
        hitData.inSystem = 0;

        const float* randFloat = optixLaunchParams.randomNumbers;
        cuuint32_t randInd = optixLaunchParams.simConstants.nbRandNumbersPerThread*(thrInd);
        cuuint32_t randOffset = optixLaunchParams.perThreadData.randBufferOffset[thrInd];

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
        optixLaunchParams.perThreadData.randBufferOffset[thrInd] = randOffset;

        // --------------------------------------
        // increase facet counters for desorption
        // --------------------------------------
        const cuuint32_t counterIdx = facIndex + ((thrInd)%(CORESPERMP*WARPSCHEDULERS))*optixLaunchParams.simConstants.nbFacets;
#ifdef BOUND_CHECK
        if(counterIdx < 0 || counterIdx >= optixLaunchParams.simConstants.nbFacets*CORESPERMP*WARPSCHEDULERS){
            printf("facIndex %u >= %u is out of bounds\n", counterIdx, optixLaunchParams.simConstants.nbFacets*CORESPERMP*WARPSCHEDULERS);
        }
#endif
        increaseHitCounterDesorption(optixLaunchParams.hitCounter[counterIdx],hitData,rayDir, rayGenData->poly[facIndex].N);
    }

    //------------------------------------------------------------------------------
    // ray gen program - the actual rendering happens in here
    //------------------------------------------------------------------------------
    extern "C" __global__ void __raygen__startFromSource()
    {
        //TODO: use thrust::random or curand

        const int ix = optixGetLaunchIndex().x;
        const int iy = optixGetLaunchIndex().y;
        // load from thread buffer ...
        const cuuint32_t threadIndex = ix + iy * optixLaunchParams.simConstants.size.x;

#ifdef BOUND_CHECK
        if(threadIndex < 0 || threadIndex >= optixLaunchParams.simConstants.nbRandNumbersPerThread * optixLaunchParams.simConstants.size.x * optixLaunchParams.simConstants.size.y){
            printf("threadIndex %u is out of bounds\n", threadIndex);
        }
#endif
        vec3f rayOrigin;
        vec3f rayDir;

        MolPRD hitData;



        switch (optixLaunchParams.perThreadData.currentMoleculeData[threadIndex].inSystem) {
            case 1:
            {
                /* if molecule is still in the system (bounce etc.)
                 * load data from thread buffer
                 */
                initMoleculeInSystem(threadIndex, hitData, rayDir, rayOrigin);
                break;
            }
            case 0:
            {
                /*
                 * start from a source
                 */
                initMoleculeFromStart(threadIndex, hitData, rayDir, rayOrigin);
                break;
            }
            default:
            {
                /*
                 * if molecule suffered from self intersection with starting facet, try again with new direction
                 */
                initMoleculePostSelfIntersection(threadIndex, hitData, rayDir, rayOrigin);
            }
        }

#ifdef DEBUGPOS
        const cuuint32_t posIndexOffset = optixLaunchParams.perThreadData.posOffsetBuffer_debug[threadIndex]++;
        if(posIndexOffset<NBCOUNTS){
            const cuuint32_t posIndex = threadIndex*NBCOUNTS+posIndexOffset;
            //printf("[%d] my pos is %d\n", threadIndex, posIndex);
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
            printf("[RayOffset] facIndex %u >= %u is out of bounds (%u)\n", facIndex, optixLaunchParams.simConstants.nbFacets, optixLaunchParams.perThreadData.currentMoleculeData[threadIndex].inSystem);
        }
#endif
            rayOrigin = offset_ray(make_float3(rayOrigin.x,rayOrigin.y,rayOrigin.z),make_float3(rayGenData->poly[facIndex].N.x,rayGenData->poly[facIndex].N.y,rayGenData->poly[facIndex].N.z));

        }
        optixTrace(optixLaunchParams.traversable,
               rayOrigin,
               rayDir,
                   0.f,//optixLaunchParams.simConstants.scene_epsilon,//1.e-6f,//0.f,    // tmin // TODO: Choose scene_epsilon wisely
               1e20f,  // tmax
               0.0f,   // rayTime
               OptixVisibilityMask( 255 ),
               OPTIX_RAY_FLAG_DISABLE_ANYHIT,//OPTIX_RAY_FLAG_NONE,
               SURFACE_RAY_TYPE,             // SBT offset
               RAY_TYPE_COUNT,               // SBT stride
               SURFACE_RAY_TYPE,             // missSBTIndex
               //u0, u1 , u2, u3);
               //float3_as_args(hitData.hitPos),
               reinterpret_cast<cuuint32_t&>(hitData.velocity),
               reinterpret_cast<cuuint32_t&>(hitData.currentDepth),
               reinterpret_cast<cuuint32_t&>(hitData.inSystem),
            /* Can't use float_as_int() because it returns rvalue but payload requires a lvalue */
               reinterpret_cast<cuuint32_t&>(hitData.hitFacetId),
               reinterpret_cast<cuuint32_t&>(hitData.hitT));

        //hitData.hitPos = hitData.hitOri + hitData.hitT * hitData.hitDir;

        // and write to thread buffer ...
        optixLaunchParams.perThreadData.currentMoleculeData[threadIndex].velocity = hitData.velocity;
        optixLaunchParams.perThreadData.currentMoleculeData[threadIndex].currentDepth = hitData.currentDepth;
        optixLaunchParams.perThreadData.currentMoleculeData[threadIndex].inSystem = hitData.inSystem;
        optixLaunchParams.perThreadData.currentMoleculeData[threadIndex].hitFacetId = hitData.hitFacetId;
        optixLaunchParams.perThreadData.currentMoleculeData[threadIndex].hitT = hitData.hitT;
    }
} // ::flowgpu
