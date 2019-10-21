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

#pragma once

#include "gdt/math/vec.h"
#include "optix7.h"

namespace osc {
    using namespace gdt;
    struct vector2f{
        float u;
        float v;
    };

    class Polygon {
    public:
        Polygon(){

            vertices2d = nullptr;
            indices = nullptr;
            this->nbVertices = 0;
            this->nbIndices = 0;
            this->O = vec3f();
            this->U = vec3f();
            this->V = vec3f();

        }
        Polygon(int32_t nbOfVertices) : nbVertices(nbOfVertices), nbIndices(nbOfVertices),O(),U(),V(),Nuv(){
            /*if(vertices2d)
                delete[] vertices2d;
            if(indices)
                delete[] indices;*/
            vertices2d = new vector2f[nbOfVertices];
            indices = new int32_t[nbOfVertices];
            //std::cout << vertices2d << " - " << indices << " / " << &vertices2d << " - " << &indices <<std::endl;
        }
        ~Polygon(){
            if(vertices2d) delete[] vertices2d;
            if(indices) delete[] indices;
        }
        Polygon& operator=(Polygon&& o){
            if (this != &o)
            {
                if(this->vertices2d) delete[] this->vertices2d;
                if(this->indices) delete[] this->indices;

                this->vertices2d = o.vertices2d;
                this->indices = o.indices;
                this->nbVertices = o.nbVertices;
                this->nbIndices = o.nbVertices;
                this->O = o.O;
                this->U = o.U;
                this->V = o.V;

                o.vertices2d = nullptr;
                o.indices = nullptr;
                o.nbVertices = 0;
                o.nbIndices = 0;
                o.O = vec3f();
                o.U = vec3f();
                o.V = vec3f();
            }

            return *this;
        }
        Polygon(Polygon&& o){
            *this = std::move(o);
        }
        Polygon& operator=(const Polygon& o){

            this->vertices2d = o.vertices2d;
            this->indices = o.indices;
            this->nbVertices = o.nbVertices;
            this->nbIndices = o.nbVertices;
            this->O = o.O;
            this->U = o.U;
            this->V = o.V;

            return *this;
        }
        Polygon(const Polygon& o){
            *this = o;
        }



        //vec2f* vertices2d;
        vector2f* vertices2d;
        int32_t* indices;
        int32_t nbVertices;
        int32_t nbIndices;
        vec3f O;
        vec3f U;
        vec3f V;
        vec3f Nuv;
    };

    struct TriangleMeshSBTData {
        vec3f  color;
        vec3f *vertex;
        vec3i *index;
    };

    struct PolygonMeshSBTData {
        vec3f  color;
        vec3f *vertex;
        Polygon *poly;
    };

    struct LaunchParams
    {
        struct {
            uint32_t *colorBuffer;
            vec2i     size;
        } frame;

        struct {
            vec3f position;
            vec3f direction;
            vec3f horizontal;
            vec3f vertical;
        } camera;

        OptixTraversableHandle traversable;
    };

} // ::osc
