// Created by pbahr

#include <optix_device.h>
#include "LaunchParams.h"

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


using namespace flowgpu;

namespace flowgpu {

    /*! launch parameters in constant memory, filled in by optix upon
        optixLaunch (this gets filled in from the buffer we pass to
        optixLaunch) */
    extern "C" __constant__ LaunchParams optixLaunchParams;

    static __forceinline__ __device__
    void *unpackPointer( uint32_t i0, uint32_t i1 )
    {
        const uint64_t uptr = static_cast<uint64_t>( i0 ) << 32 | i1;
        void*           ptr = reinterpret_cast<void*>( uptr );
        return ptr;
    }

    static __forceinline__ __device__
    void  packPointer( void* ptr, uint32_t& i0, uint32_t& i1 )
    {
        const uint64_t uptr = reinterpret_cast<uint64_t>( ptr );
        i0 = uptr >> 32;
        i1 = uptr & 0x00000000ffffffff;
    }

    template<typename T>
    static __forceinline__ __device__ T *getPRD()
    {
        const uint32_t u0 = optixGetPayload_0();
        const uint32_t u1 = optixGetPayload_1();
        return reinterpret_cast<T*>( unpackPointer( u0, u1 ) );
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

    extern "C" __global__ void __closesthit__radiance()
    {
        /*const TriangleMeshSBTData &sbtData
                = *(const TriangleMeshSBTData*)optixGetSbtDataPointer();
                // compute normal:
        const int   primID = optixGetPrimitiveIndex();
        const vec3i index  = sbtData.index[primID];
        const vec3f &A     = sbtData.vertex[index.x];
        const vec3f &B     = sbtData.vertex[index.y];
        const vec3f &C     = sbtData.vertex[index.z];
        const vec3f Ng     = normalize(cross(B-A,C-A));
                */
        const PolygonMeshSBTData &sbtData = *(const PolygonMeshSBTData*)optixGetSbtDataPointer();

        //const int   primID = optixGetPrimitiveIndex();
        //vec3f &prd = *(vec3f*)getPRD<vec3f>();
        //prd=vec3f(0.0,0.0,0.0);
        //prd = gdt::randomColor(primID);

        // compute normal:
        const int   primID = optixGetPrimitiveIndex();
        const Polygon& poly  = sbtData.poly[primID];
        vec3f &prd = *(vec3f*)getPRD<vec3f>();
        prd = gdt::randomColor(primID);
        const vec3f &A     = sbtData.vertex[sbtData.index[poly.vertOffset+0]];
        const vec3f &B     = sbtData.vertex[sbtData.index[poly.vertOffset+1]];
        const vec3f &C     = sbtData.vertex[sbtData.index[poly.vertOffset+2]];
        //const vec3f Ng     = normalize(cross(B-A,C-A));
        const vec3f Ng = normalize(
                vec3f(
                        int_as_float( optixGetAttribute_0() ),
                        int_as_float( optixGetAttribute_1() ),
                        int_as_float( optixGetAttribute_2() ))
        );
        const vec3f rayDir = optixGetWorldRayDirection();
        const float cosDN  = 0.2f + .8f*fabsf(dot(rayDir,Ng));
        /*vec3f &prd = *(vec3f*)getPRD<vec3f>();
        prd = cosDN * sbtData.color;*/
        prd = cosDN * gdt::randomColor(primID);

    }

    extern "C" __global__ void __anyhit__radiance()
    { /*! for this simple example, this will remain empty */ }



    //------------------------------------------------------------------------------
    // miss program that gets called for any ray that did not have a
    // valid intersection
    //
    // as with the anyhit/closest hit programs, in this example we only
    // need to have _some_ dummy function to set up a valid SBT
    // ------------------------------------------------------------------------------

    extern "C" __global__ void __miss__radiance()
    {
        vec3f &prd = *(vec3f*)getPRD<vec3f>();
        // set to constant white as background color
        prd = vec3f(0.9f);
    }

} // ::flowgpu
