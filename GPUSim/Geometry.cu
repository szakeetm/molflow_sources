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

    //------------------------------------------------------------------------------
    // closest hit and anyhit programs for radiance-type rays.
    //
    // Note eventually we will have to create one pair of those for each
    // ray type and each geometry type we want to render; but this
    // simple example doesn't use any actual geometries yet, so we only
    // create a single, dummy, set of them (we do have to have at least
    // one group of them to set up the SBT)
    //------------------------------------------------------------------------------

    extern "C" __device__ void intersection__polygon(float d, float u, float v, vec3f n) {
        // Fast method to check if a point is inside a polygon or not.
        // Works with convex and concave polys, orientation independent

        const PolygonMeshSBTData &sbtData = *(const PolygonMeshSBTData*)optixGetSbtDataPointer();
        const int   primID = optixGetPrimitiveIndex();

        const Polygon& poly  = sbtData.poly[primID];

        const int nbSizeMinusOne = poly.nbVertices - 1;
        const vec2f* polyPoints = sbtData.vertex2;

        int n_updown = 0;
        int n_found = 0;

        vec2f p;
        p.u = u;
        p.v = v;

        const vec3f& p0 = sbtData.vertex[sbtData.index[poly.vertOffset]];

        for (size_t j = 0; j < nbSizeMinusOne; j++) {
            const vec2f& p1 = polyPoints[sbtData.index[poly.vertOffset+j]];
            const vec2f& p2 = polyPoints[sbtData.index[poly.vertOffset+j+1]];

            const vec3f& p31 = sbtData.vertex[sbtData.index[poly.vertOffset+j]];
            const vec3f& p32 = sbtData.vertex[sbtData.index[poly.vertOffset+j+1]];

            const vec3f v2U1 = normalize(p31-p0);
            const vec3f v2V1 = cross(n,v2U1);
            const float u1 = dot(v2U1, p31-p0);
            const float v1 = dot(v2V1, p31-p0);
            const vec3f v2U2 = normalize(p32-p0);
            const vec3f v2V2 = cross(n,v2U2);
            const float u2 = dot(v2U2, p32-p0);
            const float v2 = dot(v2V2, p32-p0);

            if (p.u<u1 != p.u<u2) {
                float slope = (v2 - v1) / (u2 - u1);
                if ((slope * p.u - p.v) < (slope * u1 - v1)) {
                    n_updown++;
                }
                else {
                    n_updown--;
                }
                n_found++;
            }
        }
        //Last point. Repeating code because it's the fastest and this function is heavily used
        const vec2f& p1 = polyPoints[sbtData.index[poly.vertOffset+nbSizeMinusOne]];
        const vec2f& p2 = polyPoints[sbtData.index[poly.vertOffset+0]];

        const vec3f& p31 = sbtData.vertex[sbtData.index[poly.vertOffset+nbSizeMinusOne]];
        const vec3f v2U1 = normalize(p31-p0);
        const vec3f v2V1 = cross(n,v2U1);
        const float u1 = dot(v2U1, p31-p0);
        const float v1 = dot(v2V1, p31-p0);

        if (p.u<p1.u != p.u<0.0) {
            float slope = (0.0 - v1) / (0.0 - u1);
            if ((slope * p.u - p.v) < (slope * u1 - v1)) {
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
                    float_as_int( u ), float_as_int( v )
            );
        }
    }

    // Parallelogram intersection from the SDK optixWhitted example
    extern "C" __device__ void intersection__parallelogram()
    {
        //const Parallelogram* floor = reinterpret_cast<Parallelogram*>( optixGetSbtDataPointer() );
        const PolygonMeshSBTData &sbtData = *(const PolygonMeshSBTData*)optixGetSbtDataPointer();

        const int   primID = optixGetPrimitiveIndex();
        const Polygon& poly  = sbtData.poly[primID];
        const vec3f &Aa     = sbtData.vertex[sbtData.index[poly.vertOffset+0]];
        const vec3f &Bb     = sbtData.vertex[sbtData.index[poly.vertOffset+1]];
        const vec3f &Cc     = sbtData.vertex[sbtData.index[poly.vertOffset+2]];

        vec3f v1 = Bb-Aa; // v1 = P0P1
        vec3f v2 = Cc-Aa; // v2 = P1P2
        vec3f n = cross(v1,v2);
        v1 *= 1.0f / dot( v1, v1 );
        v2 *= 1.0f / dot( v2, v2 );
        //f->sh.N = CrossProduct(v1, v2);
        //vec3f n = cross(A,B);

        /*int ind = 2;
        while (ind < poly.nbVertices) {
            int i2 = sbtData.index[poly.vertOffset+ind++];

            v1 = Bb - Aa; // v1 = P0P1
            v2 = sbtData.vertex[i2] - Bb; // v2 = P1P2
            n = cross(v1, v2);              // Cross product
        }*/

        //n = poly.Nuv;
        //v1 = poly.U;
        //v2 = poly.V;

        vec3f ray_orig = optixGetWorldRayOrigin();
        vec3f ray_dir  = optixGetWorldRayDirection();

        const float ray_tmin = optixGetRayTmin(), ray_tmax = optixGetRayTmax();

        //make_float3( floor->plane );
        float dt = dot(ray_dir, n );
        float t = (dot(n,Aa) - dot(n, ray_orig))/dt;
        if( t > ray_tmin && t < ray_tmax )
        {
            vec3f p = ray_orig + ray_dir * t;
            vec3f vi = p - Aa;
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



        /*const int   primID = optixGetPrimitiveIndex();

        const vec3f ray_orig = optixGetWorldRayOrigin();
        const vec3f ray_dir  = optixGetWorldRayDirection();

        const float ray_tmin = optixGetRayTmin(), ray_tmax = optixGetRayTmax();

        const Polygon& poly  = sbtData.poly[primID];
        const vec3f &A     = sbtData.vertex[sbtData.index[poly.vertOffset+0]];
        const vec3f &B     = sbtData.vertex[sbtData.index[poly.vertOffset+1]];
        const vec3f &C     = sbtData.vertex[sbtData.index[poly.vertOffset+2]];
        const vec3f n     = normalize(cross(A,B));
        //const vec3f n = vec3f(-1.0) * poly->Nuv;
        const float det = DOT(n, ray_dir);

        if(det > 0.0) {
            const float iDet = 1.0 / det;
            vec3f intZ = ray_orig - C;
            const float u = iDet * DET33(intZ.x, poly.V.x, ray_dir.x,
                                   intZ.y, poly.V.y, ray_dir.y,
                                   intZ.z, poly.V.z, ray_dir.z);

            if (u >= 0.0 && u <= 1.0) {

                const float v = iDet * DET33(poly.U.x, intZ.x, ray_dir.x,
                                 poly.U.y, intZ.y, ray_dir.y,
                                 poly.U.z, intZ.z, ray_dir.z);

                if (v >= 0.0 && v <= 1.0) {

                    const float d = iDet * DOT(n, intZ);
                    if (d>0.0) {
                        //intersection__polygon(u,v);
                        intersection__polygon(d,u,v,n);
                        *//*if(inPoly > 0){
                            optixReportIntersection(
                                    d,
                                    0,
                                    float_as_int( n.x ), float_as_int( n.y ), float_as_int( n.z ),
                                    //float3_as_args(n),
                                    float_as_int( u ), float_as_int( v )
                            );
                        }*//*
                    }

                        *//*optixReportIntersection(
                                d,
                                0,
                                float_as_int( n.x ), float_as_int( n.y ), float_as_int( n.z ),
                                //float3_as_args(n),
                                float_as_int( u ), float_as_int( v )
                        );*//*
                }
            }

        }*/
    }

    extern "C" __global__ void __intersection__polygon()
    {
        intersection__parallelogram();
    }

} // ::flowgpu
