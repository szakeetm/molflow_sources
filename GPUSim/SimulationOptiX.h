//
// Created by pbahr on 04/10/2019.
//

#ifndef MOLFLOW_SIMULATIONOPTIX_H
#define MOLFLOW_SIMULATIONOPTIX_H


//#include "Simulation.h"
#include "OptixController.h"

class SimulationOptiX {
protected:
    //flowgpu::SampleWindow* window;
    flowgpu::OptixController *optixHandle;
    flowgpu::Model *model;
    uint2 kernelDimensions; // blocks and threads per block
    std::vector<uint32_t> pixels;
public:
    SimulationOptiX();
    ~SimulationOptiX();
    int LoadSimulation(flowgpu::Model* loaded_model, size_t launchSize);
    //int LoadSimulation(const std::vector<Vector3d> &geomVertices, const std::vector<SuperStructure> &structures, const WorkerParams& wp, const std::vector<std::vector<std::pair<double,double>>> CDFs);
    int RunSimulation();
    unsigned long long int GetSimulationData();
    int CloseSimulation();
};


#endif //MOLFLOW_SIMULATIONOPTIX_H
