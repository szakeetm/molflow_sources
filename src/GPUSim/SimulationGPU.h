//
// Created by pbahr on 18/05/2020.
//

#ifndef MOLFLOW_PROJ_SIMULATIONGPU_H
#define MOLFLOW_PROJ_SIMULATIONGPU_H

#include "SimulationController.h"
#include "../Parameter.h"
#include "../Simulation/MolflowSimGeom.h"

class SimulationControllerGPU;
namespace flowgpu {
    class Model;
}

class SimulationGPU  : public SimulationUnit {
public:
    SimulationGPU();
    ~SimulationGPU() override;

    std::pair<int, std::optional<std::string>> SanityCheckModel(bool strictCheck) override;
    void ClearSimulation() override;
    bool LoadSimulation(Dataport *loader, char *loadStatus);
    size_t LoadSimulation(char *loadStatus) override;
    bool UpdateOntheflySimuParams(Dataport *loader);
    int ReinitializeParticleLog() override;
private:
    size_t GetHitsSize() override;

    void ResetTmpCounters();
    bool UpdateHits(Dataport *dpHit, Dataport* dpLog,int prIdx, DWORD timeout);
public:
    void ResetSimulation() override;
    bool SimulationMCStep(size_t nbStep);

    int RebuildAccelStructure() override;

    MFSim::Particle * GetParticle(size_t i) override;
    void SetNParticle(size_t n, bool fixedSeed) override;

public:
    SimulationControllerGPU* gpuSim;
    std::shared_ptr<flowgpu::Model> model;

    GlobalHitBuffer tmpGlobalResult; //Global results since last UpdateMCHits


};


#endif //MOLFLOW_PROJ_SIMULATIONGPU_H
