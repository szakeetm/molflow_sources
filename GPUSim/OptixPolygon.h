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
        Polygon() : nbVertices(0), vertOffset(0), O(), U(), V(), Nuv(){

        }
        Polygon(int32_t nbOfVertices) : nbVertices(nbOfVertices), vertOffset(0), O(), U(), V(), Nuv(){

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
                this->nbVertices = o.nbVertices;
                this->vertOffset = o.vertOffset;
                this->O = o.O;
                this->U = o.U;
                this->V = o.V;

                o.nbVertices = 0;
                o.vertOffset = 0;
                o.O = vec3f();
                o.U = vec3f();
                o.V = vec3f();
            }
            return *this;
        }
        Polygon& operator=(const Polygon& o){
            this->nbVertices = o.nbVertices;
            this->vertOffset = o.vertOffset;
            this->O = o.O;
            this->U = o.U;
            this->V = o.V;
            return *this;
        }

        // variables for access to  global memory (indices, vertices)
        uint32_t nbVertices;
        uint32_t vertOffset;

        // variables for ray-plane (3d space) intersection
        vec3f O;
        vec3f U;
        vec3f V;
        vec3f Nuv;
    };
}

#endif //MOLFLOW_PROJ_OPTIXPOLYGON_H
