//
// Created by pbahr on 04/10/2019.
//

#ifndef MOLFLOW_SIMULATIONOPTIX_H
#define MOLFLOW_SIMULATIONOPTIX_H


//#include "Simulation.h"
#include "SimulationOptiX.h"
#include "HostData.h"

class SimulationControllerGPU {
protected:
    //flowgpu::SampleWindow* window;
    flowgpu::SimulationOptiX *optixHandle;
    flowgpu::Model *model;
    uint2 kernelDimensions; // blocks and threads per block

    HostData data;
    GlobalCounter globalCounter;

    void Resize();
    unsigned long long int GetTotalHits();
public:
    SimulationControllerGPU();
    ~SimulationControllerGPU();

    int LoadSimulation(flowgpu::Model* loaded_model, size_t launchSize);
    uint64_t RunSimulation();
    int CloseSimulation();
    int ResetSimulation();
    void AllowNewParticles();
    void CheckAndBlockDesorption();
    unsigned long long int GetSimulationData(bool silent = true);
    void IncreaseGlobalCounters(HostData* tempData);
    void PrintData();
    void PrintDataForParent();
    void PrintTotalCounters();
    void WriteDataToFile(std::string fileName);
    GlobalCounter* GetGlobalCounter() ;
    double GetTransProb(size_t polyIndex);

    bool hasEnded;
};


#endif //MOLFLOW_SIMULATIONOPTIX_H
