//
// Created by pbahr on 04/10/2019.
//

#ifndef MOLFLOW_SIMULATIONOPTIX_H
#define MOLFLOW_SIMULATIONOPTIX_H


#include "Simulation.h"
#include "SampleWindow.h"

class SimulationOptiX {
protected:
    flowgpu::SampleWindow* window;
    flowgpu::Model *model;
public:
    SimulationOptiX();
    ~SimulationOptiX();
    int LoadSimulation(const std::vector<Vector3d> &geomVertices, const std::vector<SuperStructure> &structures);
    int RunSimulation();
    int CloseSimulation();
};


#endif //MOLFLOW_SIMULATIONOPTIX_H
