//
// Created by pbahr on 12/02/2020.
//

#ifndef MOLFLOW_PROJ_MODELREADER_H
#define MOLFLOW_PROJ_MODELREADER_H

#include "Model.h"
#include "cereal/types/vector.hpp"
/*! \namespace flowgeom - Molflow Geometry code */
namespace flowgpu {
    /*int loadFromMolflow(
                Geometry& geometry, const WorkerParams& wp, const std::vector<std::vector<std::pair<double,double>>> CDFs);*/
    /*flowgpu::Model *loadFromMolflowSimu(
            const std::vector<Vector3d> &geomVertices, const std::vector<SuperStructure> &structures,
            const WorkerParams& wp, const std::vector<std::vector<std::pair<double,double>>> CDFs);*/

    // for loading from a Molflow/cereal export
    struct TempFacet {
        FacetProperties facetProperties;
        std::vector<size_t> indices;
        std::vector<Vector2d> vertices2;
        std::vector<double> outgassingMap;
        std::vector<size_t> angleMapPDF;
        std::vector<double> texelInc;

        template<class Archive>
        void serialize(Archive &archive) {
            archive(
                    cereal::make_nvp("FacetProperties", facetProperties),
                    indices,
                    vertices2,
                    outgassingMap,
                    cereal::make_nvp("angleMapVector", angleMapPDF),
                    cereal::make_nvp("textIncVector", texelInc)
            );
        }
    };

    int parseGeomFromSerialization(flowgpu::Model* model, std::vector<flowgpu::TempFacet>& facets, std::vector<float3>& vertices3d);

    flowgpu::Model *loadFromSerialization(const std::string &inputString);
    int
    loadFromSimModel(std::shared_ptr<Model> &model, std::shared_ptr<flowgpu::MolflowGPUSettings> &settings,
                     const SimulationModel &simModel);
}


#endif //MOLFLOW_PROJ_MODELREADER_H
