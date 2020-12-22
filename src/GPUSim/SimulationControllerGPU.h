//
// Created by pbahr on 04/10/2019.
//

#ifndef MOLFLOW_SIMULATIONOPTIX_H
#define MOLFLOW_SIMULATIONOPTIX_H


//#include "Simulation.h"
#include "SimulationOptiX.h"
#include "HostData.h"

struct RuntimeFigures {
    uint32_t runCount {0};
    uint32_t runCountNoEnd {0};
    double desPerRun {0.0};
    unsigned long long int total_counter = 0;
    unsigned long long int total_abs = 0;
    double total_absd = 0;
    unsigned long long int total_des = 0;
};

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
    void CheckAndBlockDesorption_exact();
    unsigned long long int GetSimulationData(bool silent = true);
    void IncreaseGlobalCounters(HostData* tempData);
    void PrintData();
    void PrintDataForParent();
    void PrintTotalCounters();
    void WriteDataToFile(std::string fileName);
    GlobalCounter* GetGlobalCounter() ;
    double GetTransProb(size_t polyIndex);

    void CalcRuntimeFigures();
    int RemainingStepsUntilStop();

    RuntimeFigures figures;
    bool hasEnded{false};
    bool endCalled{false};
};


#endif //MOLFLOW_SIMULATIONOPTIX_H
