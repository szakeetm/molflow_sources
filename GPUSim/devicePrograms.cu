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


using namespace osc;



namespace osc {

    /*! launch parameters in constant memory, filled in by optix upon
        optixLaunch (this gets filled in from the buffer we pass to
        optixLaunch) */
    extern "C" __constant__ LaunchParams optixLaunchParams;

    // for this simple example, we have a single ray type
    enum { SURFACE_RAY_TYPE=0, RAY_TYPE_COUNT };

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

    // Parallelogram intersection from the SDK optixWhitted example
    extern "C" __device__ void intersection__parallelogram()
    {
        /*const Parallelogram* floor = reinterpret_cast<Parallelogram*>( optixGetSbtDataPointer() );

        const float3 ray_orig = optixGetWorldRayOrigin();
        const float3 ray_dir  = optixGetWorldRayDirection();
        const float ray_tmin = optixGetRayTmin(), ray_tmax = optixGetRayTmax();

        float3 n = make_float3( floor->plane );
        float dt = dot(ray_dir, n );
        float t = (floor->plane.w - dot(n, ray_orig))/dt;
        if( t > ray_tmin && t < ray_tmax )
        {
            float3 p = ray_orig + ray_dir * t;
            float3 vi = p - floor->anchor;
            float a1 = dot(floor->v1, vi);
            if(a1 >= 0 && a1 <= 1)
            {
                float a2 = dot(floor->v2, vi);
                if(a2 >= 0 && a2 <= 1)
                {
                    optixReportIntersection(
                            t,
                            0,
                            float3_as_args(n),
                            float_as_int( a1 ), float_as_int( a2 )
                    );
                }
            }
        }*/
        //const Parallelogram* floor = reinterpret_cast<Parallelogram*>( optixGetSbtDataPointer() );
        //const PolygonMeshSBTData &sbtData = *(const PolygonMeshSBTData*)optixGetSbtDataPointer();
        const Polygon* poly = reinterpret_cast<Polygon*>( optixGetSbtDataPointer() );

        const vec3f ray_orig = optixGetWorldRayOrigin();
        const vec3f ray_dir  = optixGetWorldRayDirection();
        const float ray_tmin = optixGetRayTmin(), ray_tmax = optixGetRayTmax();

        vec3f n = poly->Nuv;
        float dt = dot(ray_dir, n );
        float t = (dot(n, ray_orig))/dt;
        if( t > ray_tmin && t < ray_tmax )
        {
            vec3f p = ray_orig + ray_dir * t;
            vec3f vi = p;// - floor->anchor;
            float a1 = dot(poly->U, vi);
            if(a1 >= 0 && a1 <= 1)
            {
                float a2 = dot(poly->V, vi);
                if(a2 >= 0 && a2 <= 1)
                {
                    optixReportIntersection(
                            t,
                            0,
                            float_as_int( n.x ), float_as_int( n.y ), float_as_int( n.z ),
                             //float3_as_args(n),
                            float_as_int( a1 ), float_as_int( a2 )
                    );
                }
            }
        }
    }

    extern "C" __global__ void __intersection__polygon()
    {

        /*const PolygonMeshSBTData &sbtData
                = *(const PolygonMeshSBTData*)optixGetSbtDataPointer();

        // compute normal
        const int   primID = optixGetPrimitiveIndex();
        const Polygon poly  = sbtData.poly[primID];
        const vec3f &A     = poly.vertices2d[0];
        const vec3f &B     = poly.indices[1];
        const vec3f &C     = poly.indices[2];
        const vec3f Ng     = poly.Nuv; //normalize(cross(B-A,C-A));
        //float3 Nuv = CROSS(sbtData.poly,sbtData.poly);

        const vec3f rayDirOpposite = optixGetWorldRayDirection(); // really opposite?
        const vec3f rayPos = optixGetWorldRayOrigin(); // really opposite?

        double det = DOT(Ng, rayDirOpposite);
        // Eliminate "back facet"
        if ((f->sh.is2sided) || (det > 0.0)) { //If 2-sided or if ray going opposite facet normal

            double u, v, d;
            // Ray/rectangle instersection. Find (u,v,dist) and check 0<=u<=1, 0<=v<=1, dist>=0

            if (det != 0.0) {

                double iDet = 1.0 / det;
                vec3f intZ = rayPos - poly.O;

                u = iDet * DET33(intZ.x, poly.V.x, rayDirOpposite.x,
                                 intZ.y, poly.V.y, rayDirOpposite.y,
                                 intZ.z, poly.V.z, rayDirOpposite.z);

                if (u >= 0.0 && u <= 1.0) {

                    v = iDet * DET33(poly.U.x, intZ.x, rayDirOpposite.x,
                                     poly.U.y, intZ.y, rayDirOpposite.y,
                                     poly.U.z, intZ.z, rayDirOpposite.z);

                    if (v >= 0.0 && v <= 1.0) {

                        d = iDet * DOT(poly.Nuv, intZ);

                        if (d > 0.0) {
                            // go on with point in poly
                        }
                    }
                }
            }
        }*/

        intersection__parallelogram();

        // 2d point in poly check
        /*int n_updown = 0;
        int n_found = 0;

        sbtData.poly->nbVertices;
        size_t nbSizeMinusOne = sbtData.poly->nbVertices - 1;
        for (size_t j = 0; j < nbSizeMinusOne; j++) {
            const vec2d& p1 = sbtData.poly->vertices2d[j];
            const vec2d& p2 = sbtData.poly->vertices2d[j+1];

            if (p.u<p1.u != p.u<p2.u) {
                double slope = (p2.v - p1.v) / (p2.u - p1.u);
                if ((slope * p.u - p.v) < (slope * p1.u - p1.v)) {
                    n_updown++;
                }
                else {
                    n_updown--;
                }
                n_found++;
            }
        }
s
        //Last point. Repeating code because it's the fastest and this function is heavily used
        const vec2d& p1 = polyPoints[nbSizeMinusOne];
        const vec2d& p2 = polyPoints[0];

        if (p.u<p1.u != p.u<p2.u) {
            double slope = (p2.v - p1.v) / (p2.u - p1.u);
            if ((slope * p.u - p.v) < (slope * p1.u - p1.v)) {
                n_updown++;
            }
            else {
                n_updown--;
            }
            n_found++;
        }

        if(((n_found / 2) & 1) ^ ((n_updown / 2) & 1)){
            optixReportIntersection(
                    t,
                    0,
                    float3_as_args(n),
                    float_as_int( a1 ), float_as_int( a2 )
            );
        }*/


    }

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
        //const PolygonMeshSBTData &sbtData = *(const PolygonMeshSBTData*)optixGetSbtDataPointer();

        const int   primID = optixGetPrimitiveIndex();
        vec3f &prd = *(vec3f*)getPRD<vec3f>();
        prd = gdt::randomColor(primID);
        // compute normal:
        /*const int   primID = optixGetPrimitiveIndex();
        const Polygon poly  = sbtData.poly[primID];
        const vec3f &A     = poly.indices[0];
        const vec3f &B     = poly.indices[1];
        const vec3f &C     = poly.indices[2];
        const vec3f Ng     = normalize(cross(B-A,C-A));

        const vec3f rayDir = optixGetWorldRayDirection();
        const float cosDN  = 0.2f + .8f*fabsf(dot(rayDir,Ng));
        vec3f &prd = *(vec3f*)getPRD<vec3f>();
        prd = cosDN * sbtData.color;*/
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
        prd = vec3f(1.0f);
    }

    //------------------------------------------------------------------------------
    // ray gen program - the actual rendering happens in here
    //------------------------------------------------------------------------------
    extern "C" __global__ void __raygen__renderFrame()
    {
        // compute a test pattern based on pixel ID
        const int ix = optixGetLaunchIndex().x;
        const int iy = optixGetLaunchIndex().y;

        const auto &camera = optixLaunchParams.camera;

        // our per-ray data for this example. what we initialize it to
        // won't matter, since this value will be overwritten by either
        // the miss or hit program, anyway
        vec3f pixelColorPRD = vec3f(0.f);

        // the values we store the PRD pointer in:
        uint32_t u0, u1;
        packPointer( &pixelColorPRD, u0, u1 );

        // normalized screen plane position, in [0,1]^2
        const vec2f screen(vec2f(ix+.5f,iy+.5f)
                           / vec2f(optixLaunchParams.frame.size));

        // generate ray direction
        vec3f rayDir = normalize(camera.direction
                                 + (screen.x - 0.5f) * camera.horizontal
                                 + (screen.y - 0.5f) * camera.vertical);

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
                   u0, u1 );

        const int r = int(255.99f*pixelColorPRD.x);
        const int g = int(255.99f*pixelColorPRD.y);
        const int b = int(255.99f*pixelColorPRD.z);

        // convert to 32-bit rgba value (we explicitly set alpha to 0xff
        // to make stb_image_write happy ...
        const uint32_t rgba = 0xff000000
                              | (r<<0) | (g<<8) | (b<<16);

        // and write to frame buffer ...
        const uint32_t fbIndex = ix+iy*optixLaunchParams.frame.size.x;
        optixLaunchParams.frame.colorBuffer[fbIndex] = rgba;
    }

} // ::osc
