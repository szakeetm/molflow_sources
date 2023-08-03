#include "MolflowSimulation.h"
#include "IntersectAABB_shared.h"
#include <cstring>
#include <cereal/archives/binary.hpp>
#include <Helper/Chronometer.h>
#include <Helper/ConsoleLogger.h>
#include "Helper/StringHelper.h"
#if defined(USE_OLD_BVH)
// currently always have SuperStructure
#else
#include <RayTracing/KDTree.h>
#include <RayTracing/BVH.h>
#endif

/*
MolflowSimulation::MolflowSimulation(MolflowSimulation&& o) noexcept {

    totalDesorbed = o.totalDesorbed;

    lastLogUpdateOK = o.lastLogUpdateOK;

    model = o.model;

    particleTracers = o.particleTracers;
    for(auto& particleTracer : particleTracers) {
        particleTracer.lastHitFacet = nullptr;
        particleTracer.ray.lastIntersected = -1;
        particleTracer.model = (MolflowSimulationModel*) model.get();
    }

    hasVolatile =  o.hasVolatile;

    globalState = o.globalState;
    globParticleLog = o.globParticleLog;

}
*/

int MolflowSimulation::ReinitializeParticleLog() {

bool result = 0;
#pragma omp parallel for shared(result)
    // New GlobalSimuState structure for threads
    for (int i = 0; i < particleTracers.size(); i++)
    {
        auto& particleTracer = particleTracers[i];
        if (!particleTracer.tmpParticleLog.particleLogMutex.try_lock_for(std::chrono::seconds(10))) {
#pragma omp critical            
            result = -1;
        }
        particleTracer.tmpParticleLog.clear();
        particleTracer.tmpParticleLog.pLog.shrink_to_fit();
        if (model->otfParams.enableLogging) {
           particleTracer.tmpParticleLog.pLog.reserve(model->otfParams.logLimit);
        }
        particleTracer.tmpParticleLog.particleLogMutex.unlock();
    }
    return result;
}

MFSim::ParticleTracer * MolflowSimulation::GetParticleTracerPtr(size_t i) {
    if(i < particleTracers.size())
        return &particleTracers.at(i);
    else
        return nullptr;
}

void MolflowSimulation::ConstructParticleTracers(size_t n, bool fixedSeed) {
    particleTracers.clear();
    particleTracers.resize(n);
    size_t pid = 0;
    for(auto& particleTracer : particleTracers){
        if(fixedSeed)
         particleTracer.randomGenerator.SetSeed(42424242 + pid);
        else
         particleTracer.randomGenerator.SetSeed(GenerateSeed(pid));
        particleTracer.particleTracerId = pid++;
    }
}

std::vector<std::string> MolflowSimulation::SanityCheckModel(bool strictCheck) {
    std::vector<std::string> errLog;

    if (!model->initialized) {
        errLog.push_back("Model not initialized");
    }
    if (model->vertices3.empty()) {
        errLog.push_back("Loaded empty vertex list");
    }
    if (model->facets.empty()) {
        errLog.push_back("Loaded empty facet list");
    }
    if(model->sh.nbFacet != model->facets.size()) {
        char tmp[256];
        snprintf(tmp, 256, "Facet structure not properly initialized, size mismatch: %zu / %zu\n", model->sh.nbFacet, model->facets.size());
        errLog.push_back(tmp);
    }
    for(auto& fac : model->facets){
        bool hasAnyTexture = fac->sh.countDes || fac->sh.countAbs || fac->sh.countRefl || fac->sh.countTrans || fac->sh.countACD || fac->sh.countDirection;
        if (!fac->sh.isTextured && (fac->sh.texHeight * fac->sh.texHeight > 0)) {
            char tmp[256];
            snprintf(tmp, 256, "[Fac #%zu] Untextured facet with texture size\n", fac->globalId+1);
            errLog.push_back(tmp);
        }
        else if (!fac->sh.isTextured && (hasAnyTexture)) {
            fac->sh.countDes = false;
            fac->sh.countAbs = false;
            fac->sh.countRefl = false;
            fac->sh.countTrans = false;
            fac->sh.countACD = false;
            fac->sh.countDirection = false;
            char tmp[256];
            snprintf(tmp, 256, "[Fac #%zu] Untextured facet with texture counters\n", fac->globalId+1);
            errLog.push_back(tmp);
        }
    }

    //Molflow unique
    if (model->wp.enableDecay && model->wp.halfLife <= 0.0) {
        char tmp[255];
        sprintf(tmp, "Particle decay is set, but half life was not set [= %e]\n", model->wp.halfLife);
        errLog.push_back(tmp);
    }

    if(!globalState){
        errLog.push_back("No global simulation state set\n");
    }
    else if(!globalState->initialized){
        errLog.push_back("Global simulation state not initialized\n");
    }

    if(!errLog.empty()){
        Log::console_error("Error log on model sanity check:\n{}", FlattenLines(errLog));
    }
    return errLog;
}



int MolflowSimulation::RebuildAccelStructure() {
    Chronometer timer;
    timer.Start();

    if(model->BuildAccelStructure(globalState, AccelType::BVH, BVHAccel::SplitMethod::SAH, 2))
        return 1;

    for(auto& particleTracer : particleTracers)
        particleTracer.model = (MolflowSimulationModel*) model.get();

    timer.Stop();

    return 0;
}



size_t MolflowSimulation::LoadSimulation(ProcCommData& procInfo, LoadStatus_abstract* loadStatus) {
    Chronometer timer;
    timer.Start();

    procInfo.UpdateControllerStatus(ControllerState::Resetting,"Resetting simulation...", loadStatus);
    ResetSimulation();
    
    procInfo.UpdateControllerStatus(ControllerState::Loading,"Constructing thread hit counters...", loadStatus);
    auto* simModelPtr = (MolflowSimulationModel*) model.get();

    size_t finished = 0;
#pragma omp parallel for
    // New GlobalSimuState structure for threads
    for(int i=0;i<particleTracers.size();i++)
    {
        auto& particleTracer = particleTracers[i];
        auto& tmpResults = particleTracer.tmpState;
        tmpResults.Resize(model);

        // Init tmp vars per thread
        particleTracer.tmpFacetVars.assign(simModelPtr->sh.nbFacet, SimulationFacetTempVar());

        // Update the progress string in a thread-safe manner
#pragma omp critical
        {
            finished++;
            procInfo.UpdateControllerStatus(fmt::format("Constructing thread hit counters... [{}/{} done]",finished, particleTracers.size()), loadStatus);
        }
    }

    //Reserve particle log
    procInfo.UpdateControllerStatus("Setting up particle logs...", loadStatus);
    ReinitializeParticleLog();

    // Build ADS
    procInfo.UpdateControllerStatus("Building ray-tracing accel. structure...", loadStatus);
    RebuildAccelStructure();

    // Initialize simulation
    timer.Stop();

    Log::console_msg_master(3, "  Load {} successful\n", simModelPtr->sh.name);
    Log::console_msg_master(3, "  Geometry: {} vertex {} facets\n", simModelPtr->vertices3.size(), simModelPtr->sh.nbFacet);

    Log::console_msg_master(3, "  Geom size: {} bytes\n", simModelPtr->size());
    Log::console_msg_master(3, "  Number of structure: {}\n", simModelPtr->sh.nbSuper);
    Log::console_msg_master(3, "  Global Hit: {} bytes\n", sizeof(GlobalHitBuffer));
    Log::console_msg_master(3, "  Facet Hit : {} bytes\n", simModelPtr->sh.nbFacet * sizeof(FacetHitBuffer));

    Log::console_msg_master(3, "  Total     : {} bytes\n", GetHitsSize());
    for(auto& particleTracer : particleTracers)
        Log::console_msg_master(5, "  Seed for {}: {}\n", particleTracer.particleTracerId, particleTracer.randomGenerator.GetSeed());
    Log::console_msg_master(3, "  Loading time: {:.2f} ms\n", timer.ElapsedMs());
    
    procInfo.UpdateControllerStatus(ControllerState::Ready,"", loadStatus);
    return 0;
}

size_t MolflowSimulation::GetHitsSize() {
    MolflowSimulationModel* simModelPtr = (MolflowSimulationModel*) model.get();
    return sizeof(GlobalHitBuffer) + model->wp.globalHistogramParams.GetDataSize() +
           + model->sh.nbFacet * sizeof(FacetHitBuffer) * (1+simModelPtr->tdParams.moments.size());
}



void MolflowSimulation::ResetSimulation() {

#pragma omp parallel for
// New GlobalSimuState structure for threads
    for (int i = 0; i < particleTracers.size(); i++)
    {
        auto& particleTracer = particleTracers[i];
        particleTracer.tmpFacetVars.assign(model->sh.nbFacet, SimulationFacetTempVar());
        particleTracer.Reset();
        particleTracer.model = (MolflowSimulationModel*) model.get();
        particleTracer.totalDesorbed = 0;

        particleTracer.tmpParticleLog.clear();
    }
    totalDesorbed = 0;
}

//Commented out: will use ResetSimulation instead, almost identical
/*
void MolflowSimulation::ClearSimulation() {

    //loadOK = false;

#pragma omp parallel for
    // New GlobalSimuState structure for threads
    for (int i = 0; i < particleTracers.size(); i++)
    {
        auto& particleTracer = particleTracers[i];
        particleTracer.tmpFacetVars.assign(model->sh.nbFacet, SimulationFacetTempVar());
        particleTracer.tmpState.Reset();
        particleTracer.model = (MolflowSimulationModel*)model.get();
        particleTracer.totalDesorbed = 0;

        particleTracer.tmpParticleLog.clear();

    }
    totalDesorbed = 0;
}
*/