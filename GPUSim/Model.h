//
// Created by pbahr on 08/10/2019.
//

#ifndef MOLFLOW_PROJ_MODEL_H
#define MOLFLOW_PROJ_MODEL_H


#pragma once

#include <vector>
//#include "Vector.h"
//#include "Geometry_shared.h"
//#include "Simulation.h"
#include "LaunchParams.h"

/*! \namespace flowgpu - Molflow GPU code */
namespace flowgpu {

    struct Mesh {
        Mesh() : nbFacets(0), nbIndices(0), nbVertices(0){};
        Mesh(const Mesh&) = delete;

        //std::vector<uint32_t> indices;
        //std::vector<float2> vertices2d;
        std::vector<float3> vertices3d;
        std::vector<Polygon> poly;

        std::vector<float2> facetProbabilities;
        std::vector<float> cdfs; // should not be part of each mesh, but the model itself

        uint32_t nbFacets;
        uint32_t nbIndices;
        uint32_t nbVertices;
    };

    /*! a simple indexed triangle mesh that our sample renderer will
        render */
    struct TriangleMesh : public Mesh {
        TriangleMesh() : Mesh(){};
        TriangleMesh(const TriangleMesh&) = delete;
        //std::vector<float3> vertex;
        //std::vector<float3> normal;
        //std::vector<int3> index;
        std::vector<int3> indices;

        // material data:
        //float3              diffuse;
    };

    struct PolygonMesh : public Mesh {
        PolygonMesh() : Mesh(){};
        PolygonMesh(const PolygonMesh&) = delete;

        std::vector<uint32_t> indices;
        std::vector<float2> vertices2d;
    };

    struct Model {
        Model() : nbFacets_total(), nbIndices_total(), nbVertices_total(), useMaxwell(){};
        Model(const Model&) = delete;
        ~Model()
        {
            //for (auto mesh : meshes) delete mesh;
            for (auto& mesh : poly_meshes) {
                delete mesh;
            }
            for (auto& mesh : triangle_meshes) {
                delete mesh;
            }
        }

        std::vector<TriangleMesh *> triangle_meshes;
        std::vector<PolygonMesh *> poly_meshes;
        uint32_t nbFacets_total;
        uint32_t nbIndices_total;
        uint32_t nbVertices_total;
        //! bounding box of all vertices in the model
        //float3 bounds.lower;
        //float3 bounds.upper;

        // Global Settings
        // Should they be here as well or only LaunchParams?
        bool useMaxwell;
    };

    //Model *loadOBJ(const std::string &objFile);
    //Model *loadFromMolflow(const std::vector<Vector3d> &geomVertices, const std::vector<SuperStructure> &structures, const WorkerParams& wp, const std::vector<std::vector<std::pair<double,double>>> CDFs);
}



#endif //MOLFLOW_PROJ_MODEL_H
