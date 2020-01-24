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
//printf("[%d->%d] -- 1DET33: %20.18f > 0 \n",threadIdx.x, primID, DET33(poly.U.x, intZ.x, ray_dir.x,poly.U.y, intZ.y, ray_dir.y,poly.U.z, intZ.z, ray_dir.z));
//printf("[%d->%d] -- DET33: %20.18f > 0 \n",threadIdx.x, primID, m[0]*m[4]*m[8] + m[1]*m[5]*m[6] + m[2]*m[3]*m[7] - m[0]*m[5]*m[7] - m[1]*m[3]*m[8] - m[2]*m[4]*m[6]);
//printf("[%d->%d] -- 2DET33: %20.18f > 0 \n",threadIdx.x, primID, dot(vec3d(cross(ray_dir,poly.U)), intZ));
//printf("[%d->%d] -- 3DET33: %20.18f > 0 = %d\n",threadIdx.x, primID, poly.U.x*intZ.y*ray_dir.z + intZ.x*ray_dir.y*poly.U.z + ray_dir.x*poly.U.y*intZ.z - poly.U.x*ray_dir.y*intZ.z - intZ.x*poly.U.y*ray_dir.z - ray_dir.x*intZ.y*poly.U.z,poly.U.x*intZ.y*ray_dir.z + intZ.x*ray_dir.y*poly.U.z + ray_dir.x*poly.U.y*intZ.z - poly.U.x*ray_dir.y*intZ.z - intZ.x*poly.U.y*ray_dir.z - ray_dir.x*intZ.y*poly.U.z > 0);

#define DOT(v1, v2)  \
  ((v1.x)*(v2.x) + (v1.y)*(v2.y) + (v1.z)*(v2.z))

#define CROSS(a, b)   \
  (a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x)

#define float3_as_args(u) \
    reinterpret_cast<cuuint32_t&>((u).x), \
    reinterpret_cast<cuuint32_t&>((u).y), \
    reinterpret_cast<cuuint32_t&>((u).z)


using namespace flowgpu;

namespace flowgpu {

    /*! launch parameters in constant memory, filled in by optix upon
        optixLaunch (this gets filled in from the buffer we pass to
        optixLaunch) */
    extern "C" __constant__ LaunchParams optixLaunchParams;

    // for this simple example, we have a single ray type
    enum { SURFACE_RAY_TYPE=0, RAY_TYPE_COUNT };

    static __forceinline__ __device__
    void *unpackPointer( cuuint32_t i0, cuuint32_t i1 )
    {
        const uint64_t uptr = static_cast<uint64_t>( i0 ) << 32 | i1;
        void*           ptr = reinterpret_cast<void*>( uptr );
        return ptr;
    }

    static __forceinline__ __device__
    void  packPointer( void* ptr, cuuint32_t& i0, cuuint32_t& i1 )
    {
        const uint64_t uptr = reinterpret_cast<uint64_t>( ptr );
        i0 = uptr >> 32;
        i1 = uptr & 0x00000000ffffffff;
    }

    template<typename T>
    static __forceinline__ __device__ T *getPRD()
    {
        const cuuint32_t u0 = optixGetPayload_0();
        const cuuint32_t u1 = optixGetPayload_1();
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
        return prd;
    }

    static __device__ __inline__ void setMolPRD( const MolPRD &prd )
    {
        optixSetPayload_0( float_as_int(prd.velocity) );
        optixSetPayload_1( prd.currentDepth );
        optixSetPayload_2( prd.inSystem );
        optixSetPayload_3( prd.hitFacetId );
        optixSetPayload_4( float_as_int(prd.hitT) );
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

    extern "C" __global__ void __closesthit__molecule()
    {
#ifdef DEBUG
        //printf("--- pre %d data ---\n", threadIdx.x);
#endif

        const PolygonMeshSBTData &sbtData = *(const PolygonMeshSBTData*)optixGetSbtDataPointer();

        //const int   primID = optixGetPrimitiveIndex();
        const vec3f ray_orig = optixGetWorldRayOrigin();
        const vec3f ray_dir  = optixGetWorldRayDirection();
        const float  ray_t    = optixGetRayTmax();

        const int ix = optixGetLaunchIndex().x;
        const int iy = optixGetLaunchIndex().y;
        const cuuint32_t fbIndex = ix+iy*optixLaunchParams.simConstants.size.x;

#ifdef DEBUGMISS
        const unsigned int missIndex = fbIndex*NMISSES;
        optixLaunchParams.perThreadData.missBuffer[missIndex] = 0;
#endif

        MolPRD prd = getMolPRD();

        prd.hitT = ray_t;
        prd.hitPos = ray_orig + ray_t * ray_dir;

#ifdef DEBUG
        if(prd.hitFacetId == optixGetPrimitiveIndex()){
            printf("[%d] source and goal facet equal %d : %4.2f,%4.2f,%4.2f -> %4.2f,%4.2f,%4.2f----- (poly)\n",(blockDim.x * blockIdx.x + threadIdx.x), prd.hitFacetId, ray_orig.x,ray_orig.y,ray_orig.z,prd.hitPos.x,prd.hitPos.y,prd.hitPos.z);
        }
#endif
        prd.hitFacetId = optixGetPrimitiveIndex();

        //setMolPRD(prd);

        //TODO: only assume bounce for now
        //TODO: Need to account for post-bounce effects on same facet

        // first add facet hits
        //TODO: Check if counters can be used on threads instead of launch id
        //const cuuint32_t counterIdx = prd.hitFacetId+ix*optixLaunchParams.simConstants.nbFacets+iy*optixLaunchParams.simConstants.nbFacets*optixLaunchParams.simConstants.size.x;
        //const cuuint32_t counterIdx = prd.hitFacetId + (blockDim.x * blockIdx.x + threadIdx.x)*optixLaunchParams.simConstants.nbFacets;
        const cuuint32_t counterIdx = prd.hitFacetId + ((blockDim.x * blockIdx.x + threadIdx.x)%(CORESPERMP*WARPSCHEDULERS))*optixLaunchParams.simConstants.nbFacets;
        /*if(blockIdx.x > 0){
            printf("x[%u] Here is thread [%u , %u] from dim #%u / grid #%u! \n",blockDim.x * blockIdx.x + threadIdx.x, threadIdx.x, blockIdx.x, blockDim.x, gridDim.x);
        }*/
        /*if(threadIdx.x > 1662){
            printf("XHere is thread [%u]! \n",threadIdx.x);
        }
        if(threadIdx.y > 0){
            printf("YHere is thread [%u]! \n",threadIdx.y);
        }
        if(threadIdx.z > 0){
            printf("ZHere is thread [%u]! \n",threadIdx.z);
        }*/

        //Register (orthogonal) velocity
        const Polygon& poly  = sbtData.poly[prd.hitFacetId];

        // TODO: Consider maxwelldistribution or not and oriRatio for lowFlux and desorbtion/absorbtion
        //IncreaseFacetCounter(
        // iFacet, sHandle->currentParticle.flightTime, 1, 0, 0, 1.0 / ortVelocity,
        // (sHandle->wp.useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity);

        // TODO: Save somewhere as a shared constant instead of repetively evaluating
        double velFactor = optixLaunchParams.simConstants.useMaxwell ? 1.0 : 1.1781;

        // TODO: if bounce
        // TODO: queue bounces back in pipeline instead of recursively following them
        // stop for some point due to recursion limits

        const float ortVelocity = prd.velocity*abs(dot(ray_dir, poly.N));
        const float hitEquiv = 1.0; //1.0*prd.orientationRatio; // hit=1.0 (only changed for lowflux mode)
        const float absEquiv = 1.0; //1.0*prd.orientationRatio; // hit=1.0 (only changed for lowflux mode)

        //--------------------------
        //-------- ABSORPTION ------
        //--------------------------

#ifdef DEBUG
        /*if(prd.currentDepth > 1){
            printf("[%d] -- %d Bounces ----- \n",ix, prd.currentDepth);
            return;
        }*/
#endif

        /*if(fbIndex==0)
            printf("--- pre absorb ---\n");*/

        const float* randFloat = optixLaunchParams.randomNumbers;
        cuuint32_t randInd = NB_RAND*(fbIndex);
        cuuint32_t randOffset = optixLaunchParams.perThreadData.randBufferOffset[fbIndex];
#ifdef DEBUG

        /*if(prd.hitPos.z < -12.999 && prd.hitPos.z > -13.001){
            printf("hitting origin [%d-%d] setting hitpos: %4.2f,%4.2f,%4.2f----- \n",(blockDim.x * blockIdx.x + threadIdx.x), prd.hitFacetId, prd.hitPos.x,prd.hitPos.y,prd.hitPos.z);
        }*/
#endif

        //TODO: AtomicAdd for smaller hitCounter structure, just per threadIdx
        if(poly.stickingFactor>=0.99999f || ((poly.stickingFactor > 0.0) && (randFloat[(cuuint32_t)(randInd + randOffset++)] < (poly.stickingFactor)))){
#ifdef DEBUG

           /* if(prd.hitPos.z < -12.999 && prd.hitPos.z > -13.001){
            printf("sticking origin [%d-%d] setting hitpos: %4.2f,%4.2f,%4.2f----- \n",(blockDim.x * blockIdx.x + threadIdx.x), prd.hitFacetId, prd.hitPos.x,prd.hitPos.y,prd.hitPos.z);
        }*/
#endif

            //printf("-- %d -> Absorp ----- \n",prd.currentDepth);
            // 	IncreaseFacetCounter(iFacet, sHandle->currentParticle.flightTime,
            // 	1, 0, 1, 2.0 / ortVelocity, (sHandle->wp.useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity);

            atomicAdd(&optixLaunchParams.hitCounter[counterIdx].nbMCHit,1);
            //++optixLaunchParams.hitCounter[counterIdx].nbMCHit;
            optixLaunchParams.hitCounter[counterIdx].nbHitEquiv += hitEquiv;
            optixLaunchParams.hitCounter[counterIdx].nbDesorbed += 0.0;
            optixLaunchParams.hitCounter[counterIdx].nbAbsEquiv += absEquiv; //static_cast<double>(absorb)*sHandle->currentParticle.oriRatio;
            optixLaunchParams.hitCounter[counterIdx].sum_1_per_ort_velocity += (2.0 / ortVelocity);//prd.oriRatio * sum_1_per_v;
            optixLaunchParams.hitCounter[counterIdx].sum_v_ort += velFactor*ortVelocity;//(sHandle->wp.useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity; //prd.oriRatio * sum_v_ort;
            optixLaunchParams.hitCounter[counterIdx].sum_1_per_velocity += (hitEquiv) / prd.velocity;//(hitEquiv + static_cast<double>(desorb)) / prd.velocity;

            optixLaunchParams.perThreadData.randBufferOffset[fbIndex] = randOffset;

            #ifdef DEBUG
            //printf("--- [%d] absorbed on facet #%d after bounce #%d ---\n", threadIdx.x,prd.hitFacetId,prd.currentDepth);
#endif

            prd.inSystem = 0;
            prd.currentDepth = 0;
            setMolPRD(prd);
            return;
        }
        else if(/*poly.stickingFactor<0.00001f || */poly.stickingFactor>-0.00001f){ // bounce
#ifdef DEBUG

            /*if(prd.hitPos.z < -12.999 && prd.hitPos.z > -13.001){
            printf("bouncing origin [%d-%d] setting hitpos: %4.2f,%4.2f,%4.2f----- \n",(blockDim.x * blockIdx.x + threadIdx.x), prd.hitFacetId, prd.hitPos.x,prd.hitPos.y,prd.hitPos.z);
        }*/
#endif
            atomicAdd(&optixLaunchParams.hitCounter[counterIdx].nbMCHit,1);
            //++optixLaunchParams.hitCounter[counterIdx].nbMCHit;
            optixLaunchParams.hitCounter[counterIdx].nbHitEquiv += hitEquiv;
            optixLaunchParams.hitCounter[counterIdx].nbDesorbed += 0.0;
            optixLaunchParams.hitCounter[counterIdx].nbAbsEquiv += 0.0; //static_cast<double>(absorb)*sHandle->currentParticle.oriRatio;
            optixLaunchParams.hitCounter[counterIdx].sum_1_per_ort_velocity += (1.0 / ortVelocity);//prd.oriRatio * sum_1_per_v;
            optixLaunchParams.hitCounter[counterIdx].sum_v_ort += velFactor*ortVelocity;//(sHandle->wp.useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity; //prd.oriRatio * sum_v_ort;
            optixLaunchParams.hitCounter[counterIdx].sum_1_per_velocity += (hitEquiv) / prd.velocity;//(hitEquiv + static_cast<double>(desorb)) / prd.velocity;

            prd.inSystem = 1;
        }
        else{
            //optixThrowException( 42);

            printf("-- %lf Err ----- \n",poly.stickingFactor);
            return;
        }

        /*if(fbIndex==0)
            printf("--- pre new dir ---\n");*/

        /*if(optixLaunchParams.perThreadData.randBufferOffset[fbIndex] > NB_RAND-(NB_RAND*0.15))
            printf("[%d] Dangerous Offset for Rand: %u \n", fbIndex, optixLaunchParams.perThreadData.randBufferOffset[fbIndex]);*/


#ifdef DEBUG
        //printf("--- [%d] Starting bounce #%d with randInd = %d + %d ---\n", threadIdx.x, prd.currentDepth, randInd, optixLaunchParams.perThreadData.randBufferOffset[ix + iy * optixLaunchParams.simConstants.size.x]);
#endif
        // new ray direction (for now only diffuse)

        const float theta = acosf(sqrtf(randFloat[(cuuint32_t)(randInd + randOffset++)]));
        /*if(optixLaunchParams.perThreadData.randBufferOffset[fbIndex] > NB_RAND-(NB_RAND*0.15))
            printf("[%d] Dangerous Offset for Rand: %u \n", fbIndex, optixLaunchParams.perThreadData.randBufferOffset[fbIndex]);*/
        const float phi = randFloat[(cuuint32_t)(randInd + randOffset++)] * 2.0 * CUDART_PI_F;
        /*if(optixLaunchParams.perThreadData.randBufferOffset[fbIndex] > NB_RAND-(NB_RAND*0.15))
            printf("[%d] Dangerous Offset for Rand: %u \n", fbIndex, optixLaunchParams.perThreadData.randBufferOffset[fbIndex]);*/
        const float u = sinf(theta)*cosf(phi);
        const float v = sinf(theta)*sinf(phi);
        const float n = cosf(theta);

        vec3f nU = poly.nU;
        vec3f nV = poly.nV;
        vec3f N = poly.N;

        prd.postHitDir = u*nU + v*nV + n*N;




        optixLaunchParams.perThreadData.currentMoleculeData[fbIndex].hitPos = prd.hitPos;
        optixLaunchParams.perThreadData.currentMoleculeData[fbIndex].postHitDir = prd.postHitDir;


        // Write temporary local variables back to shared memory
        optixLaunchParams.perThreadData.randBufferOffset[fbIndex] = randOffset;

        if(prd.currentDepth >= optixLaunchParams.simConstants.maxDepth){
            //printf("-- [%d] Max Bounces reached on depth %d, resetting for next call ... ----- \n",fbIndex, prd.currentDepth);
            prd.currentDepth = 0;
        }
        else {
            ++prd.currentDepth;

#ifdef DEBUGPOS
            const cuuint32_t posIndexOff = optixLaunchParams.perThreadData.posOffsetBuffer_debug[(cuuint32_t)(ix+iy*optixLaunchParams.simConstants.size.x)]++;
            if(posIndexOff<NBCOUNTS){
        const cuuint32_t posIndex = (ix+iy*optixLaunchParams.simConstants.size.x)*NBCOUNTS+posIndexOff;
        //printf("[%d] my pos is %d\n", (cuuint32_t)(ix+iy*optixLaunchParams.simConstants.size.x), posIndex);
        optixLaunchParams.perThreadData.positionsBuffer_debug[posIndex] = prd.hitPos;
            }
#endif

            optixTrace(optixLaunchParams.traversable,
                       prd.hitPos,
                       prd.postHitDir,
                       0.f,//1e-4f,//0.f,    // tmin
                       1e20f,  // tmax
                       0.0f,   // rayTime
                       OptixVisibilityMask(255),
                       OPTIX_RAY_FLAG_DISABLE_ANYHIT,//OPTIX_RAY_FLAG_NONE,
                       SURFACE_RAY_TYPE,             // SBT offset
                       RAY_TYPE_COUNT,               // SBT stride
                       SURFACE_RAY_TYPE,             // missSBTIndex
                    //u0, u1 , u2, u3);
                    //float3_as_args(hitData.hitPos),
                       reinterpret_cast<cuuint32_t &>(prd.velocity),
                       reinterpret_cast<cuuint32_t &>(prd.currentDepth),
                       reinterpret_cast<cuuint32_t &>(prd.inSystem),
                       /* Can't use float_as_int() because it returns rvalue but payload requires a lvalue */
                       reinterpret_cast<cuuint32_t &>(prd.hitFacetId),
                       reinterpret_cast<cuuint32_t &>(prd.hitT));
        }

        setMolPRD(prd);
#ifdef DEBUG
        //printf("--- %d post launch %d ---\n", threadIdx.x,prd.currentDepth);
#endif
        /*if(fbIndex==0)
            printf("--- [%d > %d] post recurs ---\n", prd.currentDepth,optixLaunchParams.simConstants.maxDepth);*/

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

// Parallelogram intersection based on the SDK optixWhitted example
    extern "C" __device__ void is_parallelogram_miss(uint32_t primID)
    {
        const PolygonMeshSBTData &sbtData = *(const PolygonMeshSBTData*)optixGetSbtDataPointer();

        //const int   primID = optixGetPrimitiveIndex();

        vec3f ray_orig = optixGetWorldRayOrigin();
        vec3f ray_dir = optixGetWorldRayDirection();

        /*if(blockDim.x * blockIdx.x + threadIdx.x == 3 && ray_dir.x > -0.26 && ray_dir.x < -0.24)
            primID=162;*/

        ray_dir = vec3f(-1.0,-1.0,-1.0) * ray_dir;

        const float ray_tmin = optixGetRayTmin(), ray_tmax = optixGetRayTmax();

        const Polygon& poly  = sbtData.poly[primID];
        const float det = dot(poly.Nuv, ray_dir);
        printf("[%d] intersection det = %12.10f\n", primID, det);
        if(det > 0.0) {
            const float iDet = 1.0 / det;
            vec3f intZ = ray_orig - poly.O;

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




    extern "C" __global__ void __closesthit__molecule_triangle()
    {

        const TriangleMeshSBTData &sbtData = *(const TriangleMeshSBTData*)optixGetSbtDataPointer();

        //const int   primID = optixGetPrimitiveIndex();
        const vec3f ray_orig = optixGetWorldRayOrigin();
        const vec3f ray_dir  = optixGetWorldRayDirection();
        const float  ray_t    = optixGetRayTmax();

        const int ix = optixGetLaunchIndex().x;
        const int iy = optixGetLaunchIndex().y;
        const cuuint32_t fbIndex = ix+iy*optixLaunchParams.simConstants.size.x;

#ifdef DEBUGMISS
        const unsigned int missIndex = fbIndex*NMISSES;
        optixLaunchParams.perThreadData.missBuffer[missIndex] = 0;
#endif

        /*if(fbIndex==0)
            printf("--- start post processing ---\n");*/

        MolPRD prd = getMolPRD();
        //printf("[%d] -- %d Bounces loaded----- \n",ix, prd.currentDepth);

#ifdef DEBUG

#endif
        // self intersection
        if(prd.hitFacetId == optixGetPrimitiveIndex()){
            //printf("[%d] source and goal facet equal %d : %8.6f,%8.6f,%8.6f -> %8.6f,%8.6f,%8.6f : [%.5e .. %.5e] (tri)\n",(blockDim.x * blockIdx.x + threadIdx.x), prd.hitFacetId, ray_orig.x,ray_orig.y,ray_orig.z,ray_dir.x,ray_dir.y,ray_dir.z, optixGetRayTmin(),ray_t);
            prd.inSystem = 2;
            //prd.hitPos = ray_orig;
            optixLaunchParams.perThreadData.currentMoleculeData[fbIndex].hitPos = ray_orig;
            //optixLaunchParams.perThreadData.currentMoleculeData[fbIndex].postHitDir = ray_dir;

            //prd.currentDepth = 0;
            setMolPRD(prd);
            return;
        }

        prd.hitT = ray_t;
        prd.hitPos = ray_orig + ray_t * ray_dir;
        prd.hitFacetId = optixGetPrimitiveIndex();

        //setMolPRD(prd);

        //TODO: only assume bounce for now
        //TODO: Need to account for post-bounce effects on same facet

        // first add facet hits
        //TODO: Check if counters can be used on threads instead of launch id
        //const cuuint32_t counterIdx = prd.hitFacetId+ix*optixLaunchParams.simConstants.nbFacets+iy*optixLaunchParams.simConstants.nbFacets*optixLaunchParams.simConstants.size.x;
        //const cuuint32_t counterIdx = prd.hitFacetId + (blockDim.x * blockIdx.x + threadIdx.x)*optixLaunchParams.simConstants.nbFacets;
        const cuuint32_t counterIdx = prd.hitFacetId + ((blockDim.x * blockIdx.x + threadIdx.x)%(CORESPERMP*WARPSCHEDULERS))*optixLaunchParams.simConstants.nbFacets;
        /*if(blockIdx.x > 0){
            printf("x[%u] Here is thread [%u , %u] from dim #%u / grid #%u! \n",blockDim.x * blockIdx.x + threadIdx.x, threadIdx.x, blockIdx.x, blockDim.x, gridDim.x);
        }*/
        /*if(threadIdx.x > 1662){
            printf("XHere is thread [%u]! \n",threadIdx.x);
        }
        if(threadIdx.y > 0){
            printf("YHere is thread [%u]! \n",threadIdx.y);
        }
        if(threadIdx.z > 0){
            printf("ZHere is thread [%u]! \n",threadIdx.z);
        }*/

        //Register (orthogonal) velocity
        const Polygon& poly  = sbtData.poly[prd.hitFacetId];

        // TODO: Consider maxwelldistribution or not and oriRatio for lowFlux and desorbtion/absorbtion
        //IncreaseFacetCounter(
        // iFacet, sHandle->currentParticle.flightTime, 1, 0, 0, 1.0 / ortVelocity,
        // (sHandle->wp.useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity);

        // TODO: Save somewhere as a shared constant instead of repetively evaluating
        double velFactor = optixLaunchParams.simConstants.useMaxwell ? 1.0 : 1.1781;

        // TODO: if bounce
        // TODO: queue bounces back in pipeline instead of recursively following them
        // stop for some point due to recursion limits

        const float ortVelocity = prd.velocity*abs(dot(ray_dir, poly.N));
        const float hitEquiv = 1.0; //1.0*prd.orientationRatio; // hit=1.0 (only changed for lowflux mode)
        const float absEquiv = 1.0; //1.0*prd.orientationRatio; // hit=1.0 (only changed for lowflux mode)

        //--------------------------
        //-------- ABSORPTION ------
        //--------------------------

#ifdef DEBUG
        /*if(prd.currentDepth > 1){
            printf("[%d] -- %d Bounces ----- \n",ix, prd.currentDepth);
            return;
        }*/
#endif

        /*if(fbIndex==0)
            printf("--- pre absorb ---\n");*/

        const float* randFloat = optixLaunchParams.randomNumbers;
        cuuint32_t randInd = NB_RAND*(fbIndex);
        cuuint32_t randOffset = optixLaunchParams.perThreadData.randBufferOffset[fbIndex];

        //TODO: AtomicAdd for smaller hitCounter structure, just per threadIdx
        if(poly.stickingFactor>=0.99999f || ((poly.stickingFactor > 0.0) && (randFloat[(cuuint32_t)(randInd + randOffset++)] < (poly.stickingFactor)))){

            //printf("-- %d -> Absorp ----- \n",prd.currentDepth);
            // 	IncreaseFacetCounter(iFacet, sHandle->currentParticle.flightTime,
            // 	1, 0, 1, 2.0 / ortVelocity, (sHandle->wp.useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity);

            atomicAdd(&optixLaunchParams.hitCounter[counterIdx].nbMCHit,1);
            //++optixLaunchParams.hitCounter[counterIdx].nbMCHit;
            optixLaunchParams.hitCounter[counterIdx].nbHitEquiv += hitEquiv;
            optixLaunchParams.hitCounter[counterIdx].nbDesorbed += 0.0f;
            optixLaunchParams.hitCounter[counterIdx].nbAbsEquiv += absEquiv; //static_cast<double>(absorb)*sHandle->currentParticle.oriRatio;
            optixLaunchParams.hitCounter[counterIdx].sum_1_per_ort_velocity += (2.0f / ortVelocity);//prd.oriRatio * sum_1_per_v;
            optixLaunchParams.hitCounter[counterIdx].sum_v_ort += velFactor*ortVelocity;//(sHandle->wp.useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity; //prd.oriRatio * sum_v_ort;
            optixLaunchParams.hitCounter[counterIdx].sum_1_per_velocity += (hitEquiv) / prd.velocity;//(hitEquiv + static_cast<double>(desorb)) / prd.velocity;

            optixLaunchParams.perThreadData.randBufferOffset[fbIndex] = randOffset;

            prd.inSystem = 0;
            prd.currentDepth = 0;
            setMolPRD(prd);
            return;
        }
        else if(poly.stickingFactor > -0.00001f){ // bounce
            atomicAdd(&optixLaunchParams.hitCounter[counterIdx].nbMCHit,1);
            //++optixLaunchParams.hitCounter[counterIdx].nbMCHit;
            optixLaunchParams.hitCounter[counterIdx].nbHitEquiv += hitEquiv;
            optixLaunchParams.hitCounter[counterIdx].nbDesorbed += 0.0f;
            optixLaunchParams.hitCounter[counterIdx].nbAbsEquiv += 0.0f; //static_cast<double>(absorb)*sHandle->currentParticle.oriRatio;
            optixLaunchParams.hitCounter[counterIdx].sum_1_per_ort_velocity += (1.0f / ortVelocity);//prd.oriRatio * sum_1_per_v;
            optixLaunchParams.hitCounter[counterIdx].sum_v_ort += velFactor*ortVelocity;//(sHandle->wp.useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity; //prd.oriRatio * sum_v_ort;
            optixLaunchParams.hitCounter[counterIdx].sum_1_per_velocity += (hitEquiv) / prd.velocity;//(hitEquiv + static_cast<double>(desorb)) / prd.velocity;

            prd.inSystem = 1;
        }
        else{
            //optixThrowException( 42);

            printf("-- %lf Err ----- \n",poly.stickingFactor);
            return;
        }

        /*if(fbIndex==0)
            printf("--- pre new dir ---\n");*/

        /*if(optixLaunchParams.perThreadData.randBufferOffset[fbIndex] > NB_RAND-(NB_RAND*0.15))
            printf("[%d] Dangerous Offset for Rand: %u \n", fbIndex, optixLaunchParams.perThreadData.randBufferOffset[fbIndex]);*/


#ifdef DEBUG
        //printf("--- [%d] Starting bounce #%d with randInd = %d + %d ---\n", threadIdx.x, prd.currentDepth, randInd, optixLaunchParams.perThreadData.randBufferOffset[ix + iy * optixLaunchParams.simConstants.size.x]);
#endif
        // new ray direction (for now only diffuse)
        const float theta = acosf(sqrtf(randFloat[(cuuint32_t)(randInd + randOffset++)]));
        /*if(optixLaunchParams.perThreadData.randBufferOffset[fbIndex] > NB_RAND-(NB_RAND*0.15))
            printf("[%d] Dangerous Offset for Rand: %u \n", fbIndex, optixLaunchParams.perThreadData.randBufferOffset[fbIndex]);*/
        const float phi = randFloat[(cuuint32_t)(randInd + randOffset++)] * 2.0 * CUDART_PI_F;
        /*if(optixLaunchParams.perThreadData.randBufferOffset[fbIndex] > NB_RAND-(NB_RAND*0.15))
            printf("[%d] Dangerous Offset for Rand: %u \n", fbIndex, optixLaunchParams.perThreadData.randBufferOffset[fbIndex]);*/
        const float u = sinf(theta)*cosf(phi);
        const float v = sinf(theta)*sinf(phi);
        const float n = cosf(theta);

        vec3f nU = poly.nU;
        vec3f nV = poly.nV;
        vec3f N = poly.N;

        prd.postHitDir = u*nU + v*nV + n*N;




        optixLaunchParams.perThreadData.currentMoleculeData[fbIndex].hitPos = prd.hitPos;
        optixLaunchParams.perThreadData.currentMoleculeData[fbIndex].postHitDir = prd.postHitDir;


        // Write temporary local variables back to shared memory
        optixLaunchParams.perThreadData.randBufferOffset[fbIndex] = randOffset;

        if(prd.currentDepth >= optixLaunchParams.simConstants.maxDepth){
            //printf("-- [%d] Max Bounces reached on depth %d, resetting for next call ... ----- \n",fbIndex, prd.currentDepth);
            prd.currentDepth = 0;
        }
        else {
            ++prd.currentDepth;

#ifdef DEBUGPOS
            const cuuint32_t posIndexOff = optixLaunchParams.perThreadData.posOffsetBuffer_debug[(cuuint32_t)(ix+iy*optixLaunchParams.simConstants.size.x)]++;
            if(posIndexOff<NBCOUNTS){
        const cuuint32_t posIndex = (ix+iy*optixLaunchParams.simConstants.size.x)*NBCOUNTS+posIndexOff;
        //printf("[%d] my pos is %d\n", (cuuint32_t)(ix+iy*optixLaunchParams.simConstants.size.x), posIndex);
        optixLaunchParams.perThreadData.positionsBuffer_debug[posIndex] = prd.hitPos;
            }
#endif

            optixTrace(optixLaunchParams.traversable,
                       prd.hitPos,
                       prd.postHitDir,
                       0.f,//1e-4f,//0.f,    // tmin
                       1e20f,  // tmax
                       0.0f,   // rayTime
                       OptixVisibilityMask(255),
                       OPTIX_RAY_FLAG_DISABLE_ANYHIT,//OPTIX_RAY_FLAG_NONE,
                       SURFACE_RAY_TYPE,             // SBT offset
                       RAY_TYPE_COUNT,               // SBT stride
                       SURFACE_RAY_TYPE,             // missSBTIndex
                    //u0, u1 , u2, u3);
                    //float3_as_args(hitData.hitPos),
                       reinterpret_cast<cuuint32_t &>(prd.velocity),
                       reinterpret_cast<cuuint32_t &>(prd.currentDepth),
                       reinterpret_cast<cuuint32_t &>(prd.inSystem),
                    /* Can't use float_as_int() because it returns rvalue but payload requires a lvalue */
                       reinterpret_cast<cuuint32_t &>(prd.hitFacetId),
                       reinterpret_cast<cuuint32_t &>(prd.hitT));
        }

        setMolPRD(prd);
#ifdef DEBUG
        //printf("--- %d post launch %d ---\n", threadIdx.x,prd.currentDepth);
#endif
        /*if(fbIndex==0)
            printf("--- [%d > %d] post recurs ---\n", prd.currentDepth,optixLaunchParams.simConstants.maxDepth);*/

    }

    extern "C" __global__ void __miss__molecule()
    {
        //vec3i &prd = *(vec3i*)getPRD<vec3i>();

        MolPRD prd = getMolPRD();

#ifdef DEBUG
        const vec3f ray_orig = optixGetWorldRayOrigin();
        const vec3f ray_dir  = optixGetWorldRayDirection();
        const float  ray_t    = optixGetRayTmax();

        printf("--------------------(%d) miss[%d -> %d] (%12.10f , %12.10f , %12.10f) --> (%12.10f , %12.10f , %12.10f) = %12.10f \n",
               prd.inSystem, blockDim.x * blockIdx.x + threadIdx.x, prd.hitFacetId,
               ray_orig.x, ray_orig.y , ray_orig.z , ray_dir.x, ray_dir.y , ray_dir.z, ray_t);
#endif

#ifdef DEBUGMISS
        const cuuint32_t fbIndex = blockDim.x * blockIdx.x + threadIdx.x;
        const unsigned int missIndex = fbIndex*NMISSES;
        for(int i=missIndex+1; i <= missIndex+optixLaunchParams.perThreadData.missBuffer[missIndex];i++){
            printf("-------------------- miss[%d -> %d -> %d] at %d\n",
                    blockDim.x * blockIdx.x + threadIdx.x,
                    missIndex,
                    optixLaunchParams.perThreadData.missBuffer[missIndex],
                   optixLaunchParams.perThreadData.missBuffer[i]);

            is_parallelogram_miss(optixLaunchParams.perThreadData.missBuffer[i]);

        }


        optixLaunchParams.perThreadData.missBuffer[missIndex] = 0;

#endif

/*optixLaunchParams.perThreadData.currentMoleculeData[fbIndex].hitPos = prd.hitPos;
optixLaunchParams.perThreadData.currentMoleculeData[fbIndex].postHitDir = prd.postHitDir;*/

        // Try to relaunch a molecule for some time until removing it from the system
        // Differentiate between physical and single precision leaks here
        //if(prd.inSystem > 20) {
            optixLaunchParams.sharedData.missCounter += 1;
            atomicAdd(optixLaunchParams.sharedData.missCounter, 1);

            prd.velocity = -999.0;
            prd.hitPos = vec3f(-999.0, -999.0, -999.0);
            prd.postHitDir = vec3f(-999.0, -999.0, -999.0);
            prd.hitFacetId = -1;
            prd.hitT = -999.0f;
            prd.inSystem = 0;
        /*}
        else{
            prd.inSystem = max(3,prd.inSystem+1);
        }*/
        setMolPRD(prd);
    }


} // ::flowgpu
