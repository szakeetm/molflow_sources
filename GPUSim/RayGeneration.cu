// Created by pbahr

#include <optix_device.h>
#include <math_constants.h>

#include "LaunchParams.h"
#include "GPUDefines.h" // for NB_RAND

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
            printf("vertOffset %u -- %u >= %u is out of bounds\n", poly.vertOffset, poly.vertOffset+poly.nbVertices, optixLaunchParams.simConstants.nbVertices);
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

    //------------------------------------------------------------------------------
    // ray gen program - the actual rendering happens in here
    //------------------------------------------------------------------------------
    extern "C" __global__ void __raygen__startFromSource()
    {



        //TODO: check if u,v are in polygon
        //TODO: import nU,nV,N from facet
        //TODO: use thrust::random or curand

        const int ix = optixGetLaunchIndex().x;
        const int iy = optixGetLaunchIndex().y;
        // load from thread buffer ...
        const cuuint32_t fbIndex = ix+iy*optixLaunchParams.simConstants.size.x;

#ifdef BOUND_CHECK
        if(fbIndex<0 || fbIndex>=LAUNCH_SIZE_X*LAUNCH_SIZE_Y*LAUNCH_SIZE_Z){
            printf("fbIndex %u is out of bounds\n", fbIndex);
        }
#endif
        vec3f rayOrigin;
        vec3f rayDir;

        MolPRD hitData;
        //hitData.orientationRatio = 1.0f; // for lowflux
        //vec3i hitData(0);
        //packPointer( &hitData, u0, u1 );
        //packPointer( &rayDir, u2, u3 );

#ifdef DEBUG
        printf("--- [%d] insystem %d with bounce #%d ---\n", threadIdx.x,optixLaunchParams.perThreadData.hitBuffer[fbIndex].inSystem,optixLaunchParams.perThreadData.hitBuffer[fbIndex].currentDepth);
#endif

        if(optixLaunchParams.perThreadData.hitBuffer[fbIndex].inSystem){
            /* if molecule is still in the system (bounce etc.)
             * load data from thread buffer
             */
            hitData = optixLaunchParams.perThreadData.hitBuffer[fbIndex];
            rayDir = hitData.postHitDir;
            rayOrigin = hitData.hitPos;

#ifdef DEBUG
            //printf("--- in[%d] %4.2f , %4.2f , %4.2f - %4.2f , %4.2f , %4.2f ---\n", threadIdx.x,optixLaunchParams.perThreadData.hitBuffer[fbIndex].hitPos.x, optixLaunchParams.perThreadData.hitBuffer[fbIndex].hitPos.y, optixLaunchParams.perThreadData.hitBuffer[fbIndex].hitPos.z,optixLaunchParams.perThreadData.hitBuffer[fbIndex].postHitDir.x, optixLaunchParams.perThreadData.hitBuffer[fbIndex].postHitDir.y, optixLaunchParams.perThreadData.hitBuffer[fbIndex].postHitDir.z);
#endif
        }
        else {
            /* start from a source
             */
            hitData.hitPos = vec3f(-999.9f,-999.9f,-999.9f);
            hitData.postHitDir = vec3f(-999.9f,-999.9f,-999.9f);
            hitData.hitT = -999.0f;
            hitData.velocity = -999.0f;
            hitData.currentDepth = 0;

            /*if(ix%4==0){
            printf("Griddim %d, blockdim %d Blockidx %d threadidx %d and launch idx %d \n",
                   gridDim.x, blockDim.x , blockIdx.x , threadIdx.x, optixGetLaunchIndex().x);
        }*/

            float* randFloat = optixLaunchParams.randomNumbers;
            cuuint32_t randInd = NB_RAND*(ix+iy*optixLaunchParams.simConstants.size.x);
            cuuint32_t randOffset = optixLaunchParams.perThreadData.randBufferOffset[fbIndex];

            //printf("Optix vals: %d %d - %d %d - %d - %10.4f\n", ix, iy, optixGetLaunchDimensions().x, optixGetLaunchDimensions().y, randInd++, randFloat[randInd++]);

            const PolygonRayGenData* rayGenData = (PolygonRayGenData*) optixGetSbtDataPointer();

#ifdef BOUND_CHECK
            if(randInd + randOffset < 0 || randInd + randOffset >= NB_RAND*LAUNCH_SIZE_X*LAUNCH_SIZE_Y*LAUNCH_SIZE_Z){
                printf("randInd %u is out of bounds\n", randInd + randOffset);
            }
#endif
            float facetRnd = randFloat[(cuuint32_t)(randInd + randOffset++)];
            /*if(randOffset > NB_RAND-(NB_RAND*0.15))
                printf("[%d] Dangerous Offset for Rand: %u \n", fbIndex, randOffset);*/

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
            }while(!found);

#ifdef BOUND_CHECK
            if(facIndex < 0 || facIndex >= optixLaunchParams.simConstants.nbFacets){
                printf("Post facIndex %u >= %u is out of bounds\n", facIndex, optixLaunchParams.simConstants.nbFacets);
                printf("Post found %u with last prob %4.2f > facetRnd %4.2f \n", found, optixLaunchParams.simConstants.nbFacets,
                        rayGenData->facetProbabilities[facIndex].t,facetRnd);

            }
#endif

            // start position of particle (U,V) -> (x,y,z)
            float uDir = 0.0f, vDir = 0.0f;
            cuuint32_t loopNb = 0;
            short isInPoly = 0;
            //printf("[%d] INIT POLY CHECK: %d vert at poly offset %d\n",threadIdx.x,rayGenData->poly[facIndex].nbVertices,rayGenData->poly[facIndex].vertOffset);

            while(!isInPoly && loopNb < NB_INPOLYCHECKS){
#ifdef BOUND_CHECK
                if(randInd + randOffset < 0 || randInd + randOffset >= NB_RAND*LAUNCH_SIZE_X*LAUNCH_SIZE_Y*LAUNCH_SIZE_Z){
                    printf("randInd %u is out of bounds\n", randInd + randOffset);
                }
#endif
                uDir = randFloat[(cuuint32_t)(randInd + randOffset++)]; // vec3f(randFloat[randInd++],randFloat[randInd++],randFloat[randInd++])
                /*if(randOffset > NB_RAND-(NB_RAND*0.15))
                    printf("[%d] Dangerous Offset for Rand: %u \n", fbIndex, randOffset);*/
#ifdef BOUND_CHECK
                if(randInd + randOffset < 0 || randInd + randOffset >= NB_RAND*LAUNCH_SIZE_X*LAUNCH_SIZE_Y*LAUNCH_SIZE_Z){
                    printf("randInd %u is out of bounds\n", randInd + randOffset);
                }
#endif
                vDir = randFloat[(cuuint32_t)(randInd + randOffset++)];
                /*if(randOffset > NB_RAND-(NB_RAND*0.15))
                    printf("[%d] Dangerous Offset for Rand: %u \n", fbIndex, randOffset);*/

                isInPoly = point_in_polygon(uDir,vDir, rayGenData->poly[facIndex]);

                /*if(isInPoly){

                    vec3f rayO = rayGenData->poly[facIndex].O
                                 + uDir * rayGenData->poly[facIndex].U
                                 + vDir * rayGenData->poly[facIndex].V;
                    //printf("--- [%d..%d] point inside %4.2f , %4.2f , %4.2f - %4.2f , %4.2f , %4.2f ---\n", threadIdx.x,loopNb,rayO.x, rayO.y, rayO.z);

                    break;
                }
                else{
#ifdef DEBUG

                    vec3f rayO = rayGenData->poly[facIndex].O
                                + uDir * rayGenData->poly[facIndex].U
                                + vDir * rayGenData->poly[facIndex].V;
                    //printf("--- [%d..%d] point outside %4.2f , %4.2f , %4.2f - %4.2f , %4.2f , %4.2f ---\n", threadIdx.x,loopNb,rayO.x, rayO.y, rayO.z);
#endif
                }*/
                ++loopNb;
            }

            // if no suitable starting point was found
            if(!isInPoly){
                uDir = 0.5;
                vDir = 0.5;
            }
            /*const float uDir = randFloat[(cuuint32_t)(randInd + optixLaunchParams.perThreadData.randBufferOffset[ix+iy*optixLaunchParams.simConstants.size.x]++)]; // vec3f(randFloat[randInd++],randFloat[randInd++],randFloat[randInd++])
            const float vDir = randFloat[(cuuint32_t)(randInd + optixLaunchParams.perThreadData.randBufferOffset[ix+iy*optixLaunchParams.simConstants.size.x]++)];
            */
            rayOrigin = rayGenData->poly[facIndex].O
                              + uDir * rayGenData->poly[facIndex].U
                              + vDir * rayGenData->poly[facIndex].V;

#ifdef DEBUG
            //printf("--- point inside [%d] %4.2f , %4.2f , %4.2f - %4.2f , %4.2f , %4.2f ---\n", threadIdx.x,rayOrigin.x, rayOrigin.y, rayOrigin.z);
#endif

            //TODO: verify that u,v is inside the actual polygon
            hitData.hitFacetId = facIndex;

            // generate ray direction
#ifdef BOUND_CHECK
            if(randInd + randOffset < 0 || randInd + randOffset >= NB_RAND*LAUNCH_SIZE_X*LAUNCH_SIZE_Y*LAUNCH_SIZE_Z){
                printf("randInd %u is out of bounds\n", randInd + randOffset);
            }
#endif
            const float theta = acosf(sqrtf(randFloat[(cuuint32_t)(randInd + randOffset++)]));
            /*if(randOffset > NB_RAND-(NB_RAND*0.15))
                printf("[%d] Dangerous Offset for Rand: %u \n", fbIndex, randOffset);*/

#ifdef BOUND_CHECK
            if(randInd + randOffset < 0 || randInd + randOffset >= NB_RAND*LAUNCH_SIZE_X*LAUNCH_SIZE_Y*LAUNCH_SIZE_Z){
                printf("randInd %u is out of bounds\n", randInd + randOffset);
            }
#endif
            const float phi = randFloat[(cuuint32_t)(randInd + randOffset++)] * 2.0 * CUDART_PI_F;
            /*if(randOffset > NB_RAND-(NB_RAND*0.15))
                printf("[%d] Dangerous Offset for Rand: %u \n", fbIndex, randOffset);*/

            const float u = sinf(theta)*cosf(phi);
            const float v = sinf(theta)*sinf(phi);
            const float n = cosf(theta);

            vec3f nU = rayGenData->poly[facIndex].nU;
            vec3f nV = rayGenData->poly[facIndex].nV;
            vec3f N = rayGenData->poly[facIndex].N;

            rayDir = u*nU + v*nV + n*N;
            /*if (rayDir.x != 0.0) rayDir.x = 1.0 / rayDir.x;
            if (rayDir.y != 0.0) rayDir.y = 1.0 / rayDir.y;
            if (rayDir.z != 0.0) rayDir.z = 1.0 / rayDir.z;
            rayDir = vec3f(-1.0f,-1.0f,-1.0f) * rayDir;*/

            // our per-ray data for this example. what we initialize it to
            // won't matter, since this value will be overwritten by either
            // the miss or hit program, anyway

            // the values we store the PRD pointer in:
            //cuuint32_t u0, u1, u2, u3;

            // test ray for debugging
            //rayOrigin = vec3f(0.557846,-0.252638,0);
            //rayDir = vec3f( 0.242345,0.621272,0.745178);

            //TODO: Increase FacetCounters according to "desorption"
#ifdef DEBUG
            //printf("---out[%d] %4.2f , %4.2f , %4.2f - %4.2f , %4.2f , %4.2f ---\n", threadIdx.x, rayOrigin.x, rayOrigin.y, rayOrigin.z, rayDir.x, rayDir.y, rayDir.z);
#endif

            // Write back cached local variables to shared memory
            optixLaunchParams.perThreadData.randBufferOffset[fbIndex] = randOffset;
        }
#ifdef DEBUG
        //printf("--- [%d] start trace with insystem %d with bounce #%d ---\n", threadIdx.x,hitData.inSystem,hitData.currentDepth);
        //printf("--- [%d] %4.2f , %4.2f , %4.2f - %4.2f , %4.2f , %4.2f ---\n", threadIdx.x, rayOrigin.x, rayOrigin.y, rayOrigin.z, rayDir.x, rayDir.y, rayDir.z);
#endif

        /*if(fbIndex==0)
            printf("--- pre trace ---\n");*/

        optixTrace(optixLaunchParams.traversable,
               rayOrigin,
               rayDir,
               0.f,//1e-4f,//0.f,    // tmin
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
        optixLaunchParams.perThreadData.hitBuffer[fbIndex].velocity = hitData.velocity;
        optixLaunchParams.perThreadData.hitBuffer[fbIndex].currentDepth = hitData.currentDepth;
        optixLaunchParams.perThreadData.hitBuffer[fbIndex].inSystem = hitData.inSystem;
        optixLaunchParams.perThreadData.hitBuffer[fbIndex].hitFacetId = hitData.hitFacetId;
        optixLaunchParams.perThreadData.hitBuffer[fbIndex].hitT = hitData.hitT;


    }

} // ::flowgpu
