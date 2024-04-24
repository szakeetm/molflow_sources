
#pragma once

#include <tuple>
#include <vector>
#include <string>
#include <mutex>

#include "SimulationUnit.h"
#include "Buffer_shared.h" //Facetproperties
#include "SMP.h"
#include "Vector.h"
#include "Random.h"
#include "ProcessControl.h"
#include "SimulationController.h"
#include "MolflowSimGeom.h"
#include "RayTracing/RTHelper.h"

class ParticleTracer;

// Local simulation structure
class MolflowSimulation : public Simulation_Abstract {
public:

    MolflowSimulation()=default;
    //MolflowSimulation(MolflowSimulation&& o) noexcept;
    virtual ~MolflowSimulation() = default;

    std::vector<std::string> SanityCheckModel(bool strictCheck) override;
    //void ClearSimulation() override;
    size_t LoadSimulation(ProcCommData& procInfo, LoadStatus_abstract* loadStatus=nullptr) override;
    int RebuildAccelStructure() override;

    void ResetSimulation() override;

    size_t GetHitsSize() override;

    int ReinitializeParticleLog() override;
    std::shared_ptr<MFSim::ParticleTracer> GetParticleTracerPtr(size_t i) override;
    void ConstructParticleTracers(size_t n, bool fixedSeed) override;
    bool lastLogUpdateOK=true; // Last log update timeout
    std::vector<std::shared_ptr<MFSim::ParticleTracer>> particleTracers; //they have their own tmp counters
    mutable std::timed_mutex simuStateMutex;

};
// -- Methods ---------------------------------------------------