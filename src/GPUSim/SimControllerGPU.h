//
// Created by pbahr on 18/05/2020.
//

#ifndef MOLFLOW_PROJ_SIMCONTROLLERGPU_H
#define MOLFLOW_PROJ_SIMCONTROLLERGPU_H

#include "SimulationController.h"
#include "Model.h"
#include "SimulationControllerGPU.h"
#include "../Parameter.h"
#include "../GeometrySimu.h"

class SimControllerGPU : public SimulationController{
public:
    SimControllerGPU(std::string appName , std::string dpName, size_t parentPID, size_t procIdx);
    ~SimControllerGPU();

    int SanityCheckGeom() override;
    void ClearSimulation() override;
    bool LoadSimulation(Dataport *loader) override;
    bool UpdateOntheflySimuParams(Dataport *loader);

private:
    bool Load() override;
    bool UpdateParams() override;
    size_t GetHitsSize();
    void ResetTmpCounters();
    void UpdateHits(Dataport *dpHit, Dataport* dpLog,int prIdx, DWORD timeout) override;
public:
    void ResetSimulation() override;

public:
    SimulationControllerGPU gpuSim;
    flowgpu::Model* model;

    GlobalHitBuffer tmpGlobalResult; //Global results since last UpdateMCHits


};


#endif //MOLFLOW_PROJ_SIMCONTROLLERGPU_H
