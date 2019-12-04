//
// Created by pbahr on 29/10/2019.
//

#pragma once

#ifndef MOLFLOW_PROJ_OPTIXPOLYGON_H
#define MOLFLOW_PROJ_OPTIXPOLYGON_H
#include "gdt/math/vec.h"

namespace flowgpu {
    using namespace gdt;

    class Polygon {
    public:
        Polygon()
        : stickingFactor(-1.0), nbVertices(0), vertOffset(0), O(), U(), V(), Nuv(), nU(), nV(), N(){
        }
        Polygon(int32_t nbOfVertices)
        : stickingFactor(-1.0), nbVertices(nbOfVertices), vertOffset(0), O(), U(), V(), Nuv(), nU(), nV(), N(){
        }
        Polygon(Polygon&& o){
            *this = std::move(o);
        }
        Polygon(const Polygon& o){
            *this = o;
        }

        ~Polygon(){

        }

        Polygon& operator=(Polygon&& o){
            if (this != &o)
            {
                this->stickingFactor = o.stickingFactor;

                this->nbVertices = o.nbVertices;
                this->vertOffset = o.vertOffset;
                this->O = o.O;
                this->U = o.U;
                this->V = o.V;
                this->Nuv = o.Nuv;
                this->nU = o.nU;
                this->nV = o.nV;
                this->N = o.N;

                o.nbVertices = 0;
                o.vertOffset = 0;
                o.O = vec3f();
                o.U = vec3f();
                o.V = vec3f();
                o.Nuv = vec3f();
                o.nU = vec3f();
                o.nV = vec3f();
                o.N = vec3f();
            }
            return *this;
        }
        Polygon& operator=(const Polygon& o){
            this->stickingFactor = o.stickingFactor;

            this->nbVertices = o.nbVertices;
            this->vertOffset = o.vertOffset;
            this->O = o.O;
            this->U = o.U;
            this->V = o.V;
            this->Nuv = o.Nuv;
            this->N = o.N;
            this->nU = o.nU;
            this->nV = o.nV;
            return *this;
        }

        // attributes that don't describe the geometry
        float stickingFactor;

        // variables for access to  global memory (indices, vertices)
        uint32_t nbVertices;
        uint32_t vertOffset;

        // variables for ray-plane (3d space) intersection
        vec3f O;
        vec3f U;
        vec3f V;
        vec3f Nuv;

        // normalized facet vectors
        vec3f nU;
        vec3f nV;
        vec3f N;

        template<class Archive>
        void serialize(Archive & archive)
        {
            archive(
                    stickingFactor,
                    nbVertices,
                    vertOffset,
                    O,
                    U,
                    V,
                    Nuv,
                    nU,
                    nV,
                    N
            );
        }
    };
}

#endif //MOLFLOW_PROJ_OPTIXPOLYGON_H
