//
// Created by pbahr on 04/10/2019.
//

#ifndef MOLFLOW_SIMULATIONOPTIX_H
#define MOLFLOW_SIMULATIONOPTIX_H


//#include "Simulation.h"
#include "OptixController.h"
#include "HostData.h"

class SimulationOptiX {
protected:
    //flowgpu::SampleWindow* window;
    flowgpu::OptixController *optixHandle;
    flowgpu::Model *model;
    uint2 kernelDimensions; // blocks and threads per block

    HostData data;
    GlobalCounter globalCounter;

    void Resize();
    unsigned long long int GetTotalHits();
public:
    SimulationOptiX();
    ~SimulationOptiX();

    int LoadSimulation(flowgpu::Model* loaded_model, size_t launchSize);
    int RunSimulation();
    int CloseSimulation();

    unsigned long long int GetSimulationData(bool silent = true);
    void IncreaseGlobalCounters(HostData* tempData);
    void PrintData();
    void PrintDataForParent();
    void PrintTotalCounters();
    void WriteDataToFile(std::string fileName);

};


#endif //MOLFLOW_SIMULATIONOPTIX_H
