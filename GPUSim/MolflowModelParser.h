//
// Created by pbahr on 15/11/2019.
//

#ifndef MOLFLOW_PROJ_MOLFLOWMODELPARSER_H
#define MOLFLOW_PROJ_MOLFLOWMODELPARSER_H

#include "Model.h"
#include "Simulation.h"

// debug output
/*#include <fstream>
#include <cereal/archives/xml.hpp>
#include <cereal/types/vector.hpp>*/

/*#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/sequence.h>
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/replace.h>
#include <thrust/functional.h>*/

/*! \namespace flowgpu - Molflow GPU code */
namespace flowgpu {
    Model *loadFromMolflow(
            const std::vector<Vector3d> &geomVertices, const std::vector<SuperStructure> &structures,
            const WorkerParams& wp, const std::vector<std::vector<std::pair<double,double>>> CDFs);
}

#endif //MOLFLOW_PROJ_MOLFLOWMODELPARSER_H
