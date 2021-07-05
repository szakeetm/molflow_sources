// Created by pbahr

#include "CommonFunctions.cuh"
#include "PerRayData.h"
#include "LaunchParams.h"
#include "GPUDefines.h" // for NB_RAND
#include "helper_math.h"

//using namespace flowgpu;

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

    // Parallelogram intersection from the SDK optixWhitted example
    /*extern "C" __device__ void intersection__parallelogram__camera()
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

        *//*int ind = 2;
        while (ind < poly.nbVertices) {
            int i2 = sbtData.index[poly.indexOffset+ind++];

            v1 = Bb - Aa; // v1 = P0P1
            v2 = sbtData.vertex[i2] - Bb; // v2 = P1P2
            n = cross(v1, v2);              // Cross product
        }*//*

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
                            float_as_uint( a1 ), float_as_uint( a2 )
                    );
                }
            }
        }
    }*/
    extern "C" __device__ void intersection__pnp_double(double d, double u, double v, float3 n) {
        // Fast method to check if a point is inside a polygon or not.
        // Works with convex and concave polys, orientation independent

        const PolygonMeshSBTData &sbtData = *(const PolygonMeshSBTData*)optixGetSbtDataPointer();
        const int   primID = optixGetPrimitiveIndex();

        const flowgpu::Polygon& poly  = sbtData.poly[primID];

        const int nbSizeMinusOne = (int)poly.nbVertices - 1;
        const double2* polyPoints = sbtData.vertex2x64;

        //double2 p{u,v};
        //p.x = u;
        //p.y = v;

        int i, j, c = 0;
        for (i = 0, j = nbSizeMinusOne; i < poly.nbVertices; j = i++) {
            const double2& p1 = polyPoints[poly.indexOffset + i];
            const double2& p2 = polyPoints[poly.indexOffset + j];
            if ((p1.y>v) != (p2.y>v)) {
#if defined(DEBUG)
                if(p2.y - p1.y == 0.0) optixThrowException(100);
#endif
                double slope = (v - p1.y) / (p2.y - p1.y);
                if(u < (p2.x - p1.x) * slope + p1.x)
                    c = !c;
            }
        }

        if(c){
            /*if (fbIndex == 1601 && primID == 1)
                printf("poly in --> (%lf , %lf) -- %lf --> (%d)\n", u, v, d, (c));

            if(!(optixGetRayTmin() <= static_cast<float>(d) && static_cast<float>(d) <= optixGetRayTmax())){
                if (fbIndex == 1601 && primID == 1)
                    printf("but poly out --> (%lf <= %lf <= %lf)\n", optixGetRayTmin() , static_cast<float>(d) , optixGetRayTmax());
            }*/
            optixReportIntersection(
                    static_cast<float>(d),
                    0,
                    float3_as_args(n),
                    float_as_uint( static_cast<float>(u) ), float_as_uint( static_cast<float>(v) )
            );
        }
        else{
            /*if (fbIndex == 1601 && primID == 1)
                printf("poly out--> (%lf , %lf) -- %lf --> (%d)\n", u, v, d, (c));*/
        }
    }

    extern "C" __device__ void intersection__pnp(float d, float u, float v, float3 n) {
        // Fast method to check if a point is inside a polygon or not.
        // Works with convex and concave polys, orientation independent

        const PolygonMeshSBTData &sbtData = *(const PolygonMeshSBTData*)optixGetSbtDataPointer();
        const int   primID = optixGetPrimitiveIndex();

#ifdef BOUND_CHECK
        if(primID < 0 || primID >= optixLaunchParams.simConstants.nbFacets){
            printf("primID %u >= %u is out of bounds\n", primID, optixLaunchParams.simConstants.nbFacets);
            optixThrowException(1);
        }
#endif
        const flowgpu::Polygon& poly  = sbtData.poly[primID];

        const int nbSizeMinusOne = (int)poly.nbVertices - 1;
        const float2* polyPoints = sbtData.vertex2;

        //float2 p{u,v};
        //p.x = u;
        //p.y = v;

#ifdef BOUND_CHECK
        if(poly.nbVertices < 0 || poly.nbVertices >= optixLaunchParams.simConstants.nbVertices){
            printf("poly.nbVertices %u >= %u is out of bounds\n", poly.nbVertices, optixLaunchParams.simConstants.nbVertices);
            optixThrowException(2);
        }
#endif
#ifdef BOUND_CHECK
        if(nbSizeMinusOne < 0 || poly.indexOffset + nbSizeMinusOne >= optixLaunchParams.simConstants.nbIndices * optixLaunchParams.simConstants.nbFacets){
            printf("poly.indexOffset %u >= %u [%d < 0]is out of bounds\n", poly.indexOffset, poly.indexOffset + nbSizeMinusOne, nbSizeMinusOne);
            optixThrowException(3);
        }
#endif

        int i, j, c = 0;
        for (i = 0, j = nbSizeMinusOne; i < poly.nbVertices; j = i++) {
            const float2& p1 = polyPoints[poly.indexOffset + i];
            const float2& p2 = polyPoints[poly.indexOffset + j];
            if ((p1.y>v) != (p2.y>v)) {
#if defined(DEBUG)
                if(p2.y - p1.y == 0.0) optixThrowException(100);
#endif
                float slope = (v - p1.y) / (p2.y - p1.y);
                if(u < (p2.x - p1.x) * slope + p1.x)
                    c = !c;
            }
        }

        if(c){
            /*if (fbIndex == 1601 && primID == 1)
                printf("poly in --> (%lf , %lf) -- %lf --> (%d)\n", u, v, d, (c));*/

            /*if(!(optixGetRayTmin() <= static_cast<float>(d) && static_cast<float>(d) <= optixGetRayTmax())){
                if (fbIndex == 1601 && primID == 1)
                    printf("but poly out --> (%lf <= %lf <= %lf)\n", optixGetRayTmin() , static_cast<float>(d) , optixGetRayTmax());
            }*/
            optixReportIntersection(
                    static_cast<float>(d),
                    0,
                    float3_as_args(n),
                    float_as_uint( static_cast<float>(u) ), float_as_uint( static_cast<float>(v) )
            );
        }
        /*else{
            if (fbIndex == 1601 && primID == 1)
                printf("poly out--> (%lf , %lf) -- %lf --> (%d)\n", u, v, d, (c));
        }*/
        /*if(c){
            optixReportIntersection(
                    d,
                    0,
                    float3_as_args(n),
                    float_as_uint( u ), float_as_uint( v )
            );
        }*/
    }

    extern "C" __device__ void intersection__polygon_double(double d, double u, double v, float3 n) {
        // Fast method to check if a point is inside a polygon or not.
        // Works with convex and concave polys, orientation independent

        const PolygonMeshSBTData &sbtData = *(const PolygonMeshSBTData*)optixGetSbtDataPointer();
        const int   primID = optixGetPrimitiveIndex();

        const flowgpu::Polygon& poly  = sbtData.poly[primID];

        const int nbSizeMinusOne = (int)poly.nbVertices - 1;
        const double2* polyPoints = sbtData.vertex2x64;

        int n_updown = 0;
        int n_found = 0;

        double2 p;
        p.x = u;
        p.y = v;

        for (size_t j = 0; j < nbSizeMinusOne; j++) {
            const double2& p1 = polyPoints[poly.indexOffset + j];
            const double2& p2 = polyPoints[poly.indexOffset + j + 1];

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

        if(((n_found / 2u) & 1u) ^ ((n_updown / 2u) & 1u)){
            /*if (fbIndex == 1601 && primID == 1)
                printf("poly in --> (%lf , %lf) -- %lf --> (%d , %d)\n", u, v, d, ((n_found / 2u) & 1u), ((n_updown / 2u) & 1u));

            if(!(optixGetRayTmin() <= static_cast<float>(d) && static_cast<float>(d) <= optixGetRayTmax())){
                if (fbIndex == 1601 && primID == 1)
                    printf("but poly out --> (%lf <= %lf <= %lf)\n", optixGetRayTmin() , static_cast<float>(d) , optixGetRayTmax());
            }*/
            optixReportIntersection(
                    static_cast<float>(d),
                    0,
                    float3_as_args(n),
                    float_as_uint( static_cast<float>(u) ), float_as_uint( static_cast<float>(v) )
            );
        }
        /*else{
            if (fbIndex == 1601 && primID == 1)
                printf("poly out--> (%lf , %lf) -- %lf --> (%d , %d)\n", u, v, d, ((n_found / 2u) & 1u), ((n_updown / 2u) & 1u));

        }*/
    }

    // Parallelogram intersection based on the SDK optixWhitted example
    extern "C" __device__ void intersection__parallelogram_double() {
        const PolygonMeshSBTData &sbtData = *(const PolygonMeshSBTData*)optixGetSbtDataPointer();

        const int   primID = optixGetPrimitiveIndex();

        const flowgpu::Polygon& poly  = sbtData.poly[primID];

        const unsigned int fbIndex = getWorkIndex();


        double3 ray_orig = make_double3(optixGetWorldRayOrigin());
        //const double3 ray_dir = make_double3(optixGetWorldRayDirection().x,optixGetWorldRayDirection().y,optixGetWorldRayDirection().z);
        //ray_dir = make_double3(-1.0,-1.0,-1.0) * ray_dir;
        double3 ray_dir = make_double3(-1.0f * optixGetWorldRayDirection());

        //if(fbIndex==0) printf("[float] ray %lf , %lf , %lf ---> %lf , %lf , %lf\n", optixGetWorldRayOrigin().x , optixGetWorldRayOrigin().y, optixGetWorldRayOrigin().z , optixGetWorldRayDirection().x , optixGetWorldRayDirection().y, optixGetWorldRayDirection().z);
        //if(fbIndex==0) printf("[double] ray %lf , %lf , %lf ---> %lf , %lf , %lf\n", ray_orig.x , ray_orig.y, ray_orig.z , ray_dir.x , ray_dir.y, ray_dir.z);

        const double ray_tmin = optixGetRayTmin(), ray_tmax = optixGetRayTmax();
        /*if (fbIndex == 1601 && primID == 1) {
            printf("d[%d]--> (%lf , %lf)\n", primID, ray_tmin, ray_tmax);
            printf("d[float] ray %lf , %lf , %lf ---> %lf , %lf , %lf\n", ray_orig.x,
                   ray_orig.y, ray_orig.z, ray_dir.x,
                   ray_dir.y, ray_dir.z);
        }*/
        {
            const double det = dot(poly.Nuvx64, ray_dir);
            if (det > 0.0) {
                const double iDet = 1.0 / det;
                double3 intZ = ray_orig - poly.Ox64;

                double d = iDet * dot(poly.Nuvx64, intZ);

                /*if (fbIndex == 1601 && primID == 1)
                    printf("dist[%d]--> %lf * (%lf , %lf, %lf) = %lf\n", primID, iDet, intZ.x, intZ.y, intZ.z, d);
*/
                if (d <= 0 && d > -1e-5) {
                    d = 1e-9;
                }
                if (d > ray_tmin) {



                double u = iDet * DET33_ROW(intZ, poly.Vx64, ray_dir);
                double v = iDet * DET33_ROW(poly.Ux64, intZ, ray_dir);

                    /*if (fbIndex == 1601 && primID == 1) {
                        printf("df[%d, %d]--> (%lf , %lf)\n", fbIndex, primID, ray_tmin, ray_tmax);
                        printf("df[%d, %d]--> no interp (%lf , %lf) -- %lf\n", fbIndex, primID, u, v, d);
                        printf("df[%d, %d] ray %lf , %lf , %lf ---> %lf , %lf , %lf\n", fbIndex, primID, ray_orig.x,
                               ray_orig.y, ray_orig.z, ray_dir.x,
                               ray_dir.y, ray_dir.z);

                    }*/

                const double doubleEps = 1e-4;
                if (u >= 0.0 && u <= 1.0 && v >= 0.0 && v <= 1.0) {
                    //intersection__polygon_double(d, u, v, poly.Nuv);
                    intersection__pnp_double(d, u, v, poly.Nuv);
                }
                else if (u >= 0.0 - doubleEps && u <= 1.0 + doubleEps && v >= 0.0 - doubleEps &&
                    v <= 1.0 + doubleEps) {
                    if (u >= 0.0 - doubleEps && u <= 1.0 + doubleEps)
                        u = clamp(u, 0.0, 1.0);
                    if (v >= 0.0 - doubleEps && v <= 1.0 + doubleEps)
                        v = clamp(v, 0.0, 1.0);
                    intersection__pnp_double(d, u, v, poly.Nuv);
                    /*double u_new = clamp(u, 0.0, 1.0);
                    double v_new = clamp(v, 0.0, 1.0);
                    //intersection__polygon_double(d, u_new, v_new, poly.Nuv);
                    intersection__pnp_double(d, u_new, v_new, poly.Nuv);

                    u_new = clamp(u_new, u_new + doubleEps, u_new - doubleEps);
                    v_new = clamp(v_new, v_new + doubleEps, v_new - doubleEps);
                    intersection__pnp_double(d, u_new, v_new, poly.Nuv);

                    u_new = u > 0.5 ? u - 1e-3 : u + 1e-3;
                    v_new = v > 0.5 ? v - 1e-3 : v + 1e-3;
                    intersection__pnp_double(d, u_new, v_new, poly.Nuv);

                    u_new = u > 0.5 ? u - 1e-2 : u + 1e-2;
                    v_new = v > 0.5 ? v - 1e-2 : v + 1e-2;
                    intersection__pnp_double(d, u_new, v_new, poly.Nuv);*/
                    /*if (fbIndex == 1601 && primID == 1)
                        printf("f[%d]--> no interp (%lf , %lf) -- %lf --> (%lf , %lf)\n", primID, u, v, d, u_new, v_new);
                */
                }
                }
            }
        }

        // double ray
        {
            MolPRD& hitData = optixLaunchParams.perThreadData.currentMoleculeData[fbIndex];
#ifdef WITHTRIANGLES
            const flowgpu::TriangleRayGenData* rayGenData = (flowgpu::TriangleRayGenData*) optixGetSbtDataPointer();
#else
            const flowgpu::PolygonRayGenData* rayGenData = (flowgpu::PolygonRayGenData*) optixGetSbtDataPointer();
#endif
            ray_orig = getOrigin_double(rayGenData, hitData.hitFacetId, hitData.rndDirection[0],hitData.rndDirection[1]);
            //const double3 ray_dir = make_double3(optixGetWorldRayDirection().x,optixGetWorldRayDirection().y,optixGetWorldRayDirection().z);
            //ray_dir = make_double3(-1.0,-1.0,-1.0) * ray_dir;
            ray_dir = -1.0 * getDirection_double(hitData, poly, hitData.rndDirection[0],hitData.rndDirection[1]);

            const double det = dot(poly.Nuvx64, ray_dir);

            /*if (fbIndex == 1601 && primID == 1) {
                printf("dd[%d, %d]--> (%lf , %lf)\n", fbIndex, primID, ray_tmin, ray_tmax);
                printf("dd[%d, %d]--- %lf\n", fbIndex, primID, det);
                printf("dd[%d, %d] ray %lf , %lf , %lf ---> %lf , %lf , %lf\n", fbIndex, primID, ray_orig.x,
                       ray_orig.y, ray_orig.z, ray_dir.x,
                       ray_dir.y, ray_dir.z);
            }*/

            if (det > 0.0) {
                const double iDet = 1.0 / det;
                double3 intZ = ray_orig - poly.Ox64;

                double d = iDet * dot(poly.Nuvx64, intZ);
                if (d <= 0 && d > -1e-5) {
                    d = 1e-9;
                }
                if (d > ray_tmin) {



                double u = iDet * DET33_ROW(intZ, poly.Vx64, ray_dir);
                double v = iDet * DET33_ROW(poly.Ux64, intZ, ray_dir);

                    /*if (fbIndex == 1601 && primID == 1) {
                        printf("dd[%d, %d]--> no interp (%lf , %lf) -- %lf\n", fbIndex, primID, u, v, d);
                    }*/
                    const double doubleEps = 1e-4;
                    if (u >= 0.0 && u <= 1.0 && v >= 0.0 && v <= 1.0) {
                        //intersection__polygon_double(d, u, v, poly.Nuv);
                        intersection__pnp_double(d, u, v, poly.Nuv);
                    }
                    else if (u >= 0.0 - doubleEps && u <= 1.0 + doubleEps && v >= 0.0 - doubleEps &&
                        v <= 1.0 + doubleEps) {
                        if (u >= 0.0 - doubleEps && u <= 1.0 + doubleEps)
                            u = clamp(u, 0.0, 1.0);
                        if (v >= 0.0 - doubleEps && v <= 1.0 + doubleEps)
                            v = clamp(v, 0.0, 1.0);
                        intersection__pnp_double(d, u, v, poly.Nuv);
                        /*double u_new = clamp(u, 0.0, 1.0);
                        double v_new = clamp(v, 0.0, 1.0);
                        //intersection__polygon_double(d, u_new, v_new, poly.Nuv);
                        intersection__pnp_double(d, u_new, v_new, poly.Nuv);

                        u_new = clamp(u_new, u_new + doubleEps, u_new - doubleEps);
                        v_new = clamp(v_new, v_new + doubleEps, v_new - doubleEps);
                        intersection__pnp_double(d, u_new, v_new, poly.Nuv);

                        u_new = u > 0.5 ? u - 1e-3 : u + 1e-3;
                        v_new = v > 0.5 ? v - 1e-3 : v + 1e-3;
                        intersection__pnp_double(d, u_new, v_new, poly.Nuv);

                        u_new = u > 0.5 ? u - 1e-2 : u + 1e-2;
                        v_new = v > 0.5 ? v - 1e-2 : v + 1e-2;
                        intersection__pnp_double(d, u_new, v_new, poly.Nuv);*/
                        /*if (fbIndex == 1601 && primID == 1)
                            printf("f[%d]--> no interp (%lf , %lf) -- %lf --> (%lf , %lf)\n", primID, u, v, d, u_new, v_new);
                    */
                    }
                    /*printf("d[%d, %d]--> (%lf , %lf)\n", fbIndex, primID, ray_tmin, ray_tmax);
                    printf("d[%d, %d]--> no interp (%lf , %lf) -- %lf\n", fbIndex, primID, u, v, d);
                    printf("d[%d, %d] ray %lf , %lf , %lf ---> %lf , %lf , %lf\n", fbIndex, primID, ray_orig.x,
                           ray_orig.y, ray_orig.z, ray_dir.x,
                           ray_dir.y, ray_dir.z);*/
                }
            }
        }
    }

    extern "C" __device__ void intersection__polygon(float d, float u, float v, float3 n) {
        // Fast method to check if a point is inside a polygon or not.
        // Works with convex and concave polys, orientation independent

        const flowgpu::PolygonMeshSBTData &sbtData = *(const flowgpu::PolygonMeshSBTData*)optixGetSbtDataPointer();
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

        if(((n_found / 2u) & 1u) ^ ((n_updown / 2u) & 1u)){
            optixReportIntersection(
                    d,
                    0,
                    float3_as_args(n),
                    float_as_uint( u ), float_as_uint( v )
            );
        }
    }

    // Parallelogram intersection based on the SDK optixWhitted example
    extern "C" __device__ void intersection__parallelogram()
    {
        const flowgpu::PolygonMeshSBTData &sbtData = *(const flowgpu::PolygonMeshSBTData*)optixGetSbtDataPointer();

        const int   primID = optixGetPrimitiveIndex();

        float3 ray_orig = optixGetWorldRayOrigin();
        float3 ray_dir = optixGetWorldRayDirection();
        ray_dir = make_float3(-1.0,-1.0,-1.0) * ray_dir;

        const float ray_tmin = optixGetRayTmin(), ray_tmax = optixGetRayTmax();
#ifdef BOUND_CHECK
        if(primID < 0 || primID >= optixLaunchParams.simConstants.nbFacets){
            printf("primID %u >= %u is out of bounds\n", primID, optixLaunchParams.simConstants.nbFacets);
            #if defined(DEBUG)
                            optixThrowException(4);
            #endif
        }
#endif
        const flowgpu::Polygon& poly  = sbtData.poly[primID];

#ifdef DEBUGMISS
        const unsigned int fbIndex = getWorkIndex();
        const unsigned int missIndex = fbIndex*NMISSES;

#ifdef BOUND_CHECK
        if(missIndex < 0 || missIndex >= NMISSES*optixLaunchParams.simConstants.size.x*optixLaunchParams.simConstants.size.y){
            printf("missIndex %u >= %u is out of bounds\n", missIndex, NMISSES*optixLaunchParams.simConstants.size.x*optixLaunchParams.simConstants.size.y);
            #if defined(DEBUG)
                            optixThrowException(5);
            #endif
        }
#endif
        optixLaunchParams.perThreadData.missBuffer[missIndex]++;
        //if(optixLaunchParams.perThreadData.missBuffer[missIndex] > optixLaunchParams.simConstants.nbFacets) optixThrowException(2);
        if(optixLaunchParams.perThreadData.missBuffer[missIndex] < NMISSES) {
#ifdef BOUND_CHECK
            if(missIndex + optixLaunchParams.perThreadData.missBuffer[missIndex] >= NMISSES*optixLaunchParams.simConstants.size.x*optixLaunchParams.simConstants.size.y){
                printf("missIndex+n %u >= %u is out of bounds\n", missIndex + optixLaunchParams.perThreadData.missBuffer[missIndex], NMISSES*optixLaunchParams.simConstants.size.x*optixLaunchParams.simConstants.size.y);
                #if defined(DEBUG)
                            optixThrowException(6);
                #endif
            }
#endif
            optixLaunchParams.perThreadData.missBuffer[missIndex +
                                                       optixLaunchParams.perThreadData.missBuffer[missIndex]] = primID;
        }
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

        const float det = dot(poly.Nuv, ray_dir);
        //atomicAdd(optixLaunchParams.debugCounter.detCount,1);

        if(det > 0.0f) {
            const float iDet = 1.0f / det;
            float3 intZ = ray_orig - poly.O;

            float u = iDet * DET33_ROW(intZ, poly.V, ray_dir);
            float v = iDet * DET33_ROW(poly.U, intZ, ray_dir);

            float d = iDet * dot(poly.Nuv, intZ);

#ifdef DEBUGCOUNT
            counterIndex = (int)((u-ULOW)/(UHIGH-ULOW)*NCOUNTBINS);
            if(counterIndex<0) counterIndex = 0;
            else if(counterIndex>=NCOUNTBINS) counterIndex = NCOUNTBINS-1;
                atomicAdd(&optixLaunchParams.debugCounter.uCount[counterIndex],1);
#endif
#ifdef DEBUGCOUNT
            counterIndex = (int)((v-VLOW)/(VHIGH-VLOW)*NCOUNTBINS);
            if(counterIndex<0)counterIndex = 0;
            else if(counterIndex>=NCOUNTBINS)counterIndex = NCOUNTBINS-1;
                atomicAdd(&optixLaunchParams.debugCounter.vCount[counterIndex],1);
#endif
            /*if (fbIndex == 1601 && primID == 1) {
                printf("ff[%d, %d]--> (%lf , %lf)\n", fbIndex, primID, ray_tmin, ray_tmax);
                printf("ff[%d, %d]--> no interp (%lf , %lf) -- %lf\n", fbIndex, primID, u, v, d);
                printf("ff[%d, %d] ray %lf , %lf , %lf ---> %lf , %lf , %lf\n", fbIndex, primID, ray_orig.x,
                       ray_orig.y, ray_orig.z, ray_dir.x,
                       ray_dir.y, ray_dir.z);
            }*/

            if (d <= 0 && d > -1e-5) {
                d = 1e-9;
            }
            if (d > ray_tmin) {
                const float floatEps = 1e-4;
                if (u >= 0.0f && u <= 1.0f && v >= 0.0f && v <= 1.0f) {
                        //intersection__polygon(d, u, v, poly.Nuv);
                        intersection__pnp(d, u, v, poly.Nuv);
                }
                else if (u >= 0.0f - floatEps && u <= 1.0f + floatEps && v >= 0.0f - floatEps && v <= 1.0f + floatEps) {
                    //intersection__parallelogram_double();
                    if (u >= 0.0f - floatEps && u <= 1.0f + floatEps)
                        u = clamp(u, 0.0f, 1.0f);
                    if (v >= 0.0f - floatEps && v <= 1.0f + floatEps)
                        v = clamp(v, 0.0f, 1.0f);
                    intersection__pnp(d, u, v, poly.Nuv);
                }
                /*float floatEps = 1e-5;
                if (d>-floatEps*//*ray_tmin*//*) {
                if(d<0.0f) {
                    d *= -1.0f;
                    ray_dir = optixGetWorldRayDirection();
                }
            }
            else {
                floatEps = 1e-4;
                if (d>-floatEps*//*ray_tmin*//*) {
                    if(d<0.0f) {
                        d *= -1.0f;
                        ray_dir = optixGetWorldRayDirection();
                    }
                }
                floatEps = 1e-3;
                if (d>-floatEps*//*ray_tmin*//*) {
                    if(d<0.0f) {
                        d *= -1.0f;
                        ray_dir = optixGetWorldRayDirection();
                    }
                }
            }

            float u_orig = u;
            float v_orig = v;
            bool checkedOnce = false;
            do {
                if (u >= 0.0f - floatEps && u <= 1.0f + floatEps && v >= 0.0f - floatEps && v <= 1.0f + floatEps) {
                    checkedOnce = true;
                    //intersection__polygon(u,v);
                    intersection__polygon(d, u, v, poly.Nuv);
                    // retry with clamped values
                    if (u >= 0.0f - floatEps && u <= 1.0f + floatEps)
                        u = clamp(u, 0.0f, 1.0f);
                    if (v >= 0.0f - floatEps && v <= 1.0f + floatEps)
                        v = clamp(v, 0.0f, 1.0f);
                    intersection__polygon(d, u, v, poly.Nuv);
                    // retry with clamp harder
                    if (u_orig >= 0.0f - floatEps && u_orig <= 1.0f + floatEps)
                        u = clamp(u_orig, -u_orig, 1.0f - (u_orig - 1.0f));
                    if (v_orig >= 0.0f - floatEps && v_orig <= 1.0f + floatEps)
                        v = clamp(v_orig, -v_orig, 1.0f - (v_orig - 1.0f));
                    intersection__polygon(d, u, v, poly.Nuv);
                    // retry with clamp harder
                    if (u_orig >= 0.0f - floatEps && u_orig <= 1.0f + floatEps)
                        u = clamp(u_orig, floatEps, 1.0f - floatEps);
                    if (v_orig >= 0.0f - floatEps && v_orig <= 1.0f + floatEps)
                        v = clamp(v_orig, floatEps, 1.0f - floatEps);
                    intersection__polygon(d, u, v, poly.Nuv);
                    // retry with clamp harder
                    if (u_orig >= 0.0f - floatEps && u_orig <= 1.0f + floatEps)
                        u = clamp(u_orig, u_orig + floatEps, u_orig - floatEps);
                    if (v_orig >= 0.0f - floatEps && v_orig <= 1.0f + floatEps)
                        v = clamp(v_orig, v_orig + floatEps, v_orig - floatEps);
                    intersection__polygon(d, u, v, poly.Nuv);
                }
                floatEps *= 10.0;
            } while (floatEps <= 0.99e-3);*/
            }
        }
    }

    extern "C" __global__ void __intersection__polygon()
    {
        intersection__parallelogram();
        //intersection__parallelogram_double();
    }

} // ::flowgpu
