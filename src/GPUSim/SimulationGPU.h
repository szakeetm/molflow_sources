//
// Created by pbahr on 18/05/2020.
//

#ifndef MOLFLOW_PROJ_SIMULATIONGPU_H
#define MOLFLOW_PROJ_SIMULATIONGPU_H

#include "SimulationController.h"
#include "Model.h"
#include "SimulationControllerGPU.h"
#include "../Parameter.h"
#include "../GeometrySimu.h"

class SimulationGPU  : public SimulationUnit {
public:
    SimulationGPU();
    ~SimulationGPU();

    int SanityCheckGeom() override;
    void ClearSimulation() override;
    bool LoadSimulation(Dataport *loader, char *loadStatus) override;
    bool UpdateOntheflySimuParams(Dataport *loader) override;
    int ReinitializeParticleLog() override;
private:
    size_t GetHitsSize() override;

    void ResetTmpCounters();
    void UpdateHits(Dataport *dpHit, Dataport* dpLog,int prIdx, DWORD timeout) override;
public:
    void ResetSimulation() override;
    bool SimulationMCStep(size_t nbStep) override;

public:
    SimulationControllerGPU gpuSim;
    flowgpu::Model* model;

    GlobalHitBuffer tmpGlobalResult; //Global results since last UpdateMCHits


};


#endif //MOLFLOW_PROJ_SIMULATIONGPU_H
