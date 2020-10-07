// Created by pbahr

#include "PerRayData.h"
#include "LaunchParams.h"
#include "GPUDefines.h" // for NB_RAND
#include "helper_math.h"

#define DET33(_11,_12,_13,_21,_22,_23,_31,_32,_33)  \
  ((_11)*( (_22)*(_33) - (_32)*(_23) ) +            \
   (_12)*( (_23)*(_31) - (_33)*(_21) ) +            \
   (_13)*( (_21)*(_32) - (_31)*(_22) ))

#define DET33_ROW(_1,_2,_3)  \
  ((_1.x)*( (_2.y)*(_3.z) - (_2.z)*(_3.y) ) +            \
   (_2.x)*( (_3.y)*(_1.z) - (_3.z)*(_1.y) ) +            \
   (_3.x)*( (_1.y)*(_2.z) - (_1.z)*(_2.y) ))

#define DOT(v1, v2)  \
  ((v1.x)*(v2.x) + (v1.y)*(v2.y) + (v1.z)*(v2.z))

#define CROSS(a, b)   \
  (a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x)

#define float3_as_args(u) \
    reinterpret_cast<unsigned int&>((u).x), \
    reinterpret_cast<unsigned int&>((u).y), \
    reinterpret_cast<unsigned int&>((u).z)


//using namespace flowgpu;

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
        prd.currentDepth = optixGetPayload_0();
        prd.inSystem = optixGetPayload_1();
        prd.hitFacetId = optixGetPayload_2();
        prd.hitT = int_as_float( optixGetPayload_3() );
        prd.velocity = int_as_float( optixGetPayload_4() );

        return prd;
    }

    static __device__ __inline__ void setMolPRD( const MolPRD &prd )
    {
        optixSetPayload_0( prd.currentDepth );
        optixSetPayload_1( prd.inSystem );
        optixSetPayload_2( prd.hitFacetId );
        optixSetPayload_3( float_as_int(prd.hitT) );
        optixSetPayload_4( float_as_int(prd.velocity) );

    }

    //------------------------------------------------------------------------------
    // closest hit and anyhit programs for radiance-type rays.
    //
    // Note eventually we will have to create one pair of those for each
    // ray type and each geometry type we want to render; but this
    // simple example doesn't use any actual geometries yet, so we only
    // create a single, dummy, set of them (we do have to have at least
    // one group of them to set up the SBT)
    //------------------------------------------------------------------------------

    // Parallelogram intersection from the SDK optixWhitted example
    extern "C" __device__ void intersection__parallelogram__camera()
    {
        //const Parallelogram* floor = reinterpret_cast<Parallelogram*>( optixGetSbtDataPointer() );
        const PolygonMeshSBTData &sbtData = *(const PolygonMeshSBTData*)optixGetSbtDataPointer();

        const int   primID = optixGetPrimitiveIndex();
        const flowgpu::Polygon& poly  = sbtData.poly[primID];
        const float3 &Aa     = sbtData.vertex[sbtData.index[poly.indexOffset + 0]];
        const float3 &Bb     = sbtData.vertex[sbtData.index[poly.indexOffset + 1]];
        const float3 &Cc     = sbtData.vertex[sbtData.index[poly.indexOffset + 2]];

        float3 v1 = Bb-Aa; // v1 = P0P1
        float3 v2 = Cc-Aa; // v2 = P1P2
        float3 n = cross(v1,v2);
        v1 *= 1.0f / dot( v1, v1 );
        v2 *= 1.0f / dot( v2, v2 );

        //printf("Normal (on device): %10.4f %10.4f %10.4f \n", n.x, n.y, n.z);
        //f->sh.N = CrossProduct(v1, v2);
        //float3 n = cross(A,B);

        /*int ind = 2;
        while (ind < poly.nbVertices) {
            int i2 = sbtData.index[poly.indexOffset+ind++];

            v1 = Bb - Aa; // v1 = P0P1
            v2 = sbtData.vertex[i2] - Bb; // v2 = P1P2
            n = cross(v1, v2);              // Cross product
        }*/

        //n = poly.Nuv;
        //v1 = poly.U;
        //v2 = poly.V;

        float3 ray_orig = optixGetWorldRayOrigin();
        float3 ray_dir  = optixGetWorldRayDirection();

        const float ray_tmin = optixGetRayTmin(), ray_tmax = optixGetRayTmax();

        //make_float3( floor->plane );
        float dt = dot(ray_dir, n );
        float t = (dot(n,Aa) - dot(n, ray_orig))/dt;
        if( t > ray_tmin && t < ray_tmax )
        {
            float3 p = ray_orig + ray_dir * t;
            float3 vi = p - Aa;
            float a1 = dot(v1, vi);
            if(a1 >= 0 && a1 <= 1)
            {
                float a2 = dot(v2, vi);
                if(a2 >= 0 && a2 <= 1)
                {
                    //intersection__polygon(t,a1,a2,n);
                    optixReportIntersection(
                            t,
                            0,
                            float3_as_args(n),
                            float_as_int( a1 ), float_as_int( a2 )
                    );
                }
            }
        }
    }

    extern "C" __device__ void intersection__polygon(float d, float u, float v, float3 n) {
        // Fast method to check if a point is inside a polygon or not.
        // Works with convex and concave polys, orientation independent

        const PolygonMeshSBTData &sbtData = *(const PolygonMeshSBTData*)optixGetSbtDataPointer();
        const int   primID = optixGetPrimitiveIndex();

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

        /*if(primID==2 && optixGetLaunchIndex().x+optixGetLaunchIndex().y*optixLaunchParams.frame.size.x % 500 == 0)
            *//*printf("[%d] -- %10.4f / %10.4f / %10.4f / %10.4f > %10.4f for Ray from %10.4f , %10.4f , %10.4f to %10.4f , %10.4f , %10.4f \n",
                   primID, det, u, v, d, ray_tmin, ray_orig.x,ray_orig.y,ray_orig.z,ray_dir.x,ray_dir.y,ray_dir.z);*//*
            printf("[%d] -- found would be %d with %d and %d \n",nbSizeMinusOne, ((n_found / 2) & 1) ^ ((n_updown / 2) & 1), n_found, n_updown);*/
        if(((n_found / 2) & 1) ^ ((n_updown / 2) & 1)){
            optixReportIntersection(
                    d,
                    0,
                    float3_as_args(n),
                    float_as_int( u ), float_as_int( v )
            );
        }
    }

    extern "C" __device__ void intersection__polygon_double(double d, double u, double v, float3 n) {
        // Fast method to check if a point is inside a polygon or not.
        // Works with convex and concave polys, orientation independent

        const PolygonMeshSBTData &sbtData = *(const PolygonMeshSBTData*)optixGetSbtDataPointer();
        const int   primID = optixGetPrimitiveIndex();

        const flowgpu::Polygon& poly  = sbtData.poly[primID];

        const int nbSizeMinusOne = poly.nbVertices - 1;
        const float2* polyPoints = sbtData.vertex2;

        int n_updown = 0;
        int n_found = 0;

        double2 p;
        p.x = u;
        p.y = v;

        for (size_t j = 0; j < nbSizeMinusOne; j++) {
            const float2& p1 = polyPoints[poly.indexOffset + j];
            const float2& p2 = polyPoints[poly.indexOffset + j + 1];

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

        //Last point. Repeating code because it's the fastest and this function is heavily used
        const double2& p1 = make_double2(polyPoints[poly.indexOffset + nbSizeMinusOne].x, polyPoints[poly.indexOffset + nbSizeMinusOne].y);
        const double2& p2 = make_double2(polyPoints[poly.indexOffset + 0].x, polyPoints[poly.indexOffset + 0].y);
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

           if(((n_found / 2) & 1) ^ ((n_updown / 2) & 1)){
            optixReportIntersection(
                    d,
                    0,
                    float3_as_args(n),
                    float_as_int( float(u) ), float_as_int( float(v) )
            );
        }
    }

    // Parallelogram intersection based on the SDK optixWhitted example
    extern "C" __device__ void intersection__parallelogram()
    {
        const PolygonMeshSBTData &sbtData = *(const PolygonMeshSBTData*)optixGetSbtDataPointer();

        const int   primID = optixGetPrimitiveIndex();

        float3 ray_orig = optixGetWorldRayOrigin();
        float3 ray_dir = optixGetWorldRayDirection();

        /*if(blockDim.x * blockIdx.x + threadIdx.x == 3 && ray_dir.x > -0.26 && ray_dir.x < -0.24)
            primID=162;*/

        ray_dir = make_float3(-1.0,-1.0,-1.0) * ray_dir;

        const float ray_tmin = optixGetRayTmin(), ray_tmax = optixGetRayTmax();

        const flowgpu::Polygon& poly  = sbtData.poly[primID];
        const float det = dot(poly.Nuv, ray_dir);
        //atomicAdd(optixLaunchParams.debugCounter.detCount,1);

#ifdef DEBUGMISS
        const int ix = optixGetLaunchIndex().x;
        const int iy = optixGetLaunchIndex().y;
        const unsigned int fbIndex = ix+iy*optixLaunchParams.simConstants.size.x;
        const unsigned int missIndex = fbIndex*NMISSES;
        optixLaunchParams.perThreadData.missBuffer[missIndex]++;
        if(optixLaunchParams.perThreadData.missBuffer[missIndex] < NMISSES)
            optixLaunchParams.perThreadData.missBuffer[missIndex+optixLaunchParams.perThreadData.missBuffer[missIndex]] = primID;
#endif

#ifdef DEBUGCOUNT
        unsigned int counterIndex = (det-DETLOW)/(DETHIGH-DETLOW)*NCOUNTBINS;
        if(counterIndex<0) counterIndex = 0;
        else if(counterIndex>=NCOUNTBINS) counterIndex = NCOUNTBINS-1;
#ifdef BOUND_CHECK
        //printf("detCount Index %u >= %u is out of bounds: %4.2f\n", counterIndex, NCOUNTBINS, det);
        if(counterIndex < 0 || counterIndex >= NCOUNTBINS){
            printf("detCount Index %u >= %u is out of bounds\n", counterIndex, NCOUNTBINS);
        }
#endif
        atomicAdd(&optixLaunchParams.debugCounter.detCount[counterIndex],1);
#endif

        /*if(fabsf(det) < 1e-6){
            MolPRD prd = getMolPRD();
            printf("det[%d->%d] insystem %d det: %12.10f\n", blockDim.x * blockIdx.x + threadIdx.x,primID, prd.inSystem, det);
        }*/

        if(det > 0.0) {
            const float iDet = 1.0 / det;
            float3 intZ = ray_orig - poly.O;
#ifdef DEBUG
            /*printf("[%d] IntZ = (%12.10f , %12.10f , %12.10f) - (%12.10f , %12.10f , %12.10f)\n",
                   primID, ray_orig.x, ray_orig.y, ray_orig.z, poly.O.x, poly.O.y, poly.O.z);*/
#endif

            const float u = iDet * DET33(intZ.x, poly.V.x, ray_dir.x,
                                   intZ.y, poly.V.y, ray_dir.y,
                                   intZ.z, poly.V.z, ray_dir.z);
#ifdef DEBUGCOUNT
            counterIndex = (int)((u-ULOW)/(UHIGH-ULOW)*NCOUNTBINS);
        if(counterIndex<0) counterIndex = 0;
        else if(counterIndex>=NCOUNTBINS) counterIndex = NCOUNTBINS-1;
            atomicAdd(&optixLaunchParams.debugCounter.uCount[counterIndex],1);

#endif

            /*if(fabsf(u) < 1e-6f || (u > 1.0 && u < 1.0f + 1e-6f)){
                printf("u[%d] u = %12.10f\n", primID, u);
                *//*printf("u[%d] PolyNorm: %12.10f , %12.10f , %12.10f\n", primID, poly.Nuv.x, poly.Nuv.y, poly.Nuv.z);
                printf("u[%d] PolyU: %12.10f , %12.10f , %12.10f\n", primID, poly.U.x, poly.U.y, poly.U.z);
                printf("u[%d] IntZ: %12.10f , %12.10f , %12.10f\n", primID, intZ.x, intZ.y, intZ.z);
                printf("u[%d] DET33: %12.10f = %12.10f * %12.10f\n", primID, u, iDet,DET33(intZ.x, poly.V.x, ray_dir.x,
                                                                                          intZ.y, poly.V.y, ray_dir.y,
                                                                                          intZ.z, poly.V.z, ray_dir.z));
            *//*
            }*/
            if (u >= 0.0 && u <= 1.0) {

                const float v = iDet * DET33(poly.U.x, intZ.x, ray_dir.x,
                                 poly.U.y, intZ.y, ray_dir.y,
                                 poly.U.z, intZ.z, ray_dir.z);

#ifdef DEBUGCOUNT
                counterIndex = (int)((v-VLOW)/(VHIGH-VLOW)*NCOUNTBINS);

        if(counterIndex<0)counterIndex = 0;
        else if(counterIndex>=NCOUNTBINS)counterIndex = NCOUNTBINS-1;
            atomicAdd(&optixLaunchParams.debugCounter.vCount[counterIndex],1);
#endif

               /* if(fabsf(v) < 1e-6f || (v > 1.0 && v < 1.0f + 1e-6f)){
                    printf("v[%d] v = %12.10f\n", primID, v);
                    *//*printf("v[%d] PolyNorm: %12.10f , %12.10f , %12.10f\n", primID, poly.Nuv.x, poly.Nuv.y, poly.Nuv.z);
                    printf("v[%d] PolyU: %12.10f , %12.10f , %12.10f\n", primID, poly.U.x, poly.U.y, poly.U.z);
                    printf("v[%d] IntZ: %12.10f , %12.10f , %12.10f\n", primID, intZ.x, intZ.y, intZ.z);
                    printf("v[%d] DET33: %12.10f = %12.10f * (%12.10f + %12.10f + %12.10f)\n", primID, v, iDet,
                           (poly.U.x)*( (intZ.y)*(ray_dir.z) - (intZ.z)*(ray_dir.y)) ,
                           (intZ.x)*( (ray_dir.y)*(poly.U.z) - (ray_dir.z)*(poly.U.y)),
                           (ray_dir.x)*( (poly.U.y)*(intZ.z) - (poly.U.z)*(intZ.y)));*//*
                    //printf("[%d] 1DET33: %12.10f * (%12.10f - %12.10f) = %12.10f\n", primID, (poly.U.x), (intZ.y)*(ray_dir.z), (intZ.z)*(ray_dir.y), (poly.U.x) * ((intZ.y)*(ray_dir.z) - (intZ.z)*(ray_dir.y)));
                    //printf("[%d] 2DET33: %12.10f * (%12.10f - %12.10f) = %12.10f\n", primID, (intZ.x), (ray_dir.y)*(poly.U.z), (ray_dir.z)*(poly.U.y), (intZ.x) * ((ray_dir.y)*(poly.U.z) - (ray_dir.z)*(poly.U.y)));
                    //printf("[%d] 3DET33: %12.10f * (%12.10f - %12.10f) = %12.10f\n", primID, (ray_dir.x), (poly.U.y)*(intZ.z) , (poly.U.z)*(intZ.y), (ray_dir.x) * ( (poly.U.y)*(intZ.z) - (poly.U.z)*(intZ.y)));
                }*/

                if (v >= 0.0 && v <= 1.0) {

                    const float d = iDet * dot(poly.Nuv, intZ);

                    /*if(primID==2)
                        printf("[%d] -- %10.4f / %10.4f / %10.4f / %10.4f > %10.4f for Ray from %10.4f , %10.4f , %10.4f to %10.4f , %10.4f , %10.4f \n",
                               primID, det, u, v, d, ray_tmin, ray_orig.x,ray_orig.y,ray_orig.z,ray_dir.x,ray_dir.y,ray_dir.z);
*/
                    if (d>ray_tmin) {
                        //intersection__polygon(u,v);
                        intersection__polygon(d,u,v,poly.Nuv);
                        /*if(inPoly > 0){
                            optixReportIntersection(
                                    d,
                                    0,
                                    //float_as_int( poly.Nuv.x ), float_as_int( poly.Nuv.y ), float_as_int( n.z ),
                                    float3_as_args(poly.Nuv),
                                    float_as_int( u ), float_as_int( v )
                            );
                        }*/
                    }
                }
            }

        }
    }

    // Parallelogram intersection based on the SDK optixWhitted example
    extern "C" __device__ void intersection__parallelogram_double()
    {
        const PolygonMeshSBTData &sbtData = *(const PolygonMeshSBTData*)optixGetSbtDataPointer();

        const int   primID = optixGetPrimitiveIndex();

        const double3 ray_orig = make_double3(optixGetWorldRayOrigin());
        //const double3 ray_dir = make_double3(optixGetWorldRayDirection().x,optixGetWorldRayDirection().y,optixGetWorldRayDirection().z);
        //ray_dir = make_double3(-1.0,-1.0,-1.0) * ray_dir;
        const double3 ray_dir = make_double3(-1.0f * optixGetWorldRayDirection());

        const double ray_tmin = optixGetRayTmin(), ray_tmax = optixGetRayTmax();

        const flowgpu::Polygon& poly  = sbtData.poly[primID];
        const double det = dot(make_double3(poly.Nuv), ray_dir);

        if(det > 0.0) {
            const double iDet = 1.0 / det;
            double3 intZ = ray_orig - make_double3(poly.O);
            const double u = iDet * DET33_ROW(intZ, make_double3(poly.V), ray_dir);

            if (u >= 0.0 && u <= 1.0) {

                const double v = iDet * DET33_ROW(make_double3(poly.U), intZ, ray_dir);

                if (v >= 0.0 && v <= 1.0) {

                    const double d = iDet * dot(make_double3(poly.Nuv), intZ);

                    /*if(primID==2)
                        printf("[%d] -- %10.4f / %10.4f / %10.4f / %10.4f > %10.4f for Ray from %10.4f , %10.4f , %10.4f to %10.4f , %10.4f , %10.4f \n",
                               primID, det, u, v, d, ray_tmin, ray_orig.x,ray_orig.y,ray_orig.z,ray_dir.x,ray_dir.y,ray_dir.z);
*/
                    if (d>ray_tmin) {
                        //intersection__polygon(u,v);
                        intersection__polygon_double(d,u,v,poly.Nuv);
                        /*if(inPoly > 0){
                            optixReportIntersection(
                                    d,
                                    0,
                                    //float_as_int( poly.Nuv.x ), float_as_int( poly.Nuv.y ), float_as_int( n.z ),
                                    float3_as_args(poly.Nuv),
                                    float_as_int( u ), float_as_int( v )
                            );
                        }*/
                    }
                }
            }

        }
    }

    extern "C" __global__ void __intersection__polygon()
    {
        intersection__parallelogram();
        //intersection__parallelogram_double();
    }

} // ::flowgpu
