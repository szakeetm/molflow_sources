//
// Created by pbahr on 08/10/2019.
//

#ifndef MOLFLOW_PROJ_MODEL_H
#define MOLFLOW_PROJ_MODEL_H


#pragma once

#include "gdt/math/AffineSpace.h"
#include <vector>
//#include "Vector.h"
//#include "Geometry_shared.h"
#include "Simulation.h"
#include "LaunchParams.h"

/*! \namespace flowgpu - Molflow GPU code */
namespace flowgpu {
    using namespace gdt;

    /*! a simple indexed triangle mesh that our sample renderer will
        render */
    struct TriangleMesh {
        std::vector<vec3f> vertex;
        std::vector<vec3f> normal;
        std::vector<vec2f> texcoord;
        std::vector<vec3i> index;

        // material data:
        vec3f              diffuse;
    };

    struct PolygonMesh {
        std::vector<uint32_t> indices;
        std::vector<vec2f> vertices2d;
        std::vector<vec3f> vertices3d;
        std::vector<Polygon> poly;
        int32_t nbFacets;
        int32_t nbIndices;
        int32_t nbVertices;
    };

    struct Model {
        ~Model()
        {
            //for (auto mesh : meshes) delete mesh;
            for (auto mesh : poly_meshes) {
                delete mesh;
            }
        }

        //std::vector<TriangleMesh *> meshes;
        std::vector<PolygonMesh *> poly_meshes;
        //! bounding box of all vertices in the model
        box3f bounds;
    };

    //Model *loadOBJ(const std::string &objFile);
    Model *loadFromMolflow(const std::vector<Vector3d> &geomVertices, const std::vector<SuperStructure> &structures);
}



#endif //MOLFLOW_PROJ_MODEL_H
