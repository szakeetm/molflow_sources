#include "MolflowSimulation.h"
#include "IntersectAABB_shared.h"
#include <cstring>
#include <cereal/archives/binary.hpp>
#include <Helper/Chronometer.h>
#include <Helper/ConsoleLogger.h>
#include "Helper/StringHelper.h"
#include "Particle.h"

#include <RayTracing/KDTree.h>
#include <RayTracing/BVH.h>

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
        if (!particleTracer->tmpParticleLog->particleLogMutex.try_lock_for(std::chrono::seconds(10))) {
#pragma omp critical            
            result = -1;
        }
        particleTracer->tmpParticleLog->clear();
        particleTracer->tmpParticleLog->pLog.shrink_to_fit();
        if (model->otfParams.enableLogging) {
           particleTracer->tmpParticleLog->pLog.reserve(model->otfParams.logLimit);
        }
        particleTracer->tmpParticleLog->particleLogMutex.unlock();
    }
    return result;
}

std::shared_ptr<MFSim::ParticleTracer> MolflowSimulation::GetParticleTracerPtr(size_t i) {
    if(i < particleTracers.size())
        return particleTracers[i];
    else
        return nullptr;
}

void MolflowSimulation::ConstructParticleTracers(size_t n, bool fixedSeed) {
    particleTracers.resize(n);
    size_t pid = 0;
    for(auto& particleTracer : particleTracers){
        particleTracer = std::make_shared<MFSim::ParticleTracer>(); //make new
        if(fixedSeed)
         particleTracer->randomGenerator.SetSeed(42424242 + pid);
        else
         particleTracer->randomGenerator.SetSeed(GenerateSeed(pid));
        particleTracer->particleTracerId = pid++;
    }
}

std::vector<std::string> MolflowSimulation::SanityCheckModel(bool strictCheck) {
    std::vector<std::string> errLog = model->SanityCheck();

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
        particleTracer->model = std::static_pointer_cast<MolflowSimulationModel>(model);

    timer.Stop();

    return 0;
}



size_t MolflowSimulation::LoadSimulation(ProcCommData& procInfo, LoadStatus_abstract* loadStatus) {
    Chronometer timer;
    timer.Start();

    procInfo.UpdateControllerStatus({ ControllerState::Resetting }, { "Resetting simulation..." }, loadStatus);
    ResetSimulation();
    
    procInfo.UpdateControllerStatus({ ControllerState::Loading }, { "Constructing thread hit counters..." }, loadStatus);
    //auto simModel = std::static_pointer_cast<MolflowSimulationModel>(model);

    size_t finished = 0;
#pragma omp parallel for
    // New GlobalSimuState structure for threads
    for(int i=0;i<particleTracers.size();i++)
    {
        auto& particleTracer = particleTracers[i];
        auto& tmpResults = particleTracer->tmpState;
        tmpResults->Resize(model);

        // Init tmp vars per thread
        //particleTracer.tmpFacetVars.assign(simModelPtr->sh.nbFacet, SimulationFacetTempVar());

        // Update the progress string in a thread-safe manner
#pragma omp critical
        {
            finished++;
            procInfo.UpdateControllerStatus(std::nullopt, { fmt::format("Constructing thread hit counters... [{}/{} done]",finished, particleTracers.size()) }, loadStatus);
        }
    }

    std::vector<size_t> counterSizes;
    for (const auto& pt : particleTracers) {
        counterSizes.push_back(pt->GetMemSize());
    }
    procInfo.UpdateCounterSizes(counterSizes);

    //Reserve particle log
    procInfo.UpdateControllerStatus(std::nullopt, { "Setting up particle logs..." }, loadStatus);
    ReinitializeParticleLog();

    // Build ADS
    procInfo.UpdateControllerStatus(std::nullopt, { "Building ray-tracing accel. structure..." }, loadStatus);
    RebuildAccelStructure();

    // Initialize simulation
    timer.Stop();

    Log::console_msg_master(3, "  Load {} successful\n", model->sh.name);
    Log::console_msg_master(3, "  Geometry: {} vertex {} facets\n", model->vertices3.size(), model->sh.nbFacet);

    Log::console_msg_master(3, "  Geom size: {} bytes\n", model->memSizeCache);
    Log::console_msg_master(3, "  Number of structure: {}\n", model->sh.nbSuper);
    Log::console_msg_master(3, "  Global Hit: {} bytes\n", sizeof(GlobalHitBuffer));
    Log::console_msg_master(3, "  Facet Hit : {} bytes\n", model->sh.nbFacet * sizeof(FacetHitBuffer));

    Log::console_msg_master(3, "  Total     : {} bytes\n", GetHitsSize());
    for(auto& particleTracer : particleTracers)
        Log::console_msg_master(5, "  Seed for {}: {}\n", particleTracer->particleTracerId, particleTracer->randomGenerator.GetSeed());
    Log::console_msg_master(3, "  Loading time: {:.2f} ms\n", timer.ElapsedMs());
    
    //procInfo.UpdateControllerStatus(ControllerState::Ready,"", loadStatus);
    return 0;
}

size_t MolflowSimulation::GetHitsSize() {
    MolflowSimulationModel* simModelPtr = (MolflowSimulationModel*) model.get();
    return sizeof(GlobalHitBuffer) + model->sp.globalHistogramParams.GetDataSize() +
           + model->sh.nbFacet * sizeof(FacetHitBuffer) * (1+simModelPtr->tdParams.moments.size());
}



void MolflowSimulation::ResetSimulation() {

#pragma omp parallel for
// New GlobalSimuState structure for threads
    for (int i = 0; i < particleTracers.size(); i++)
    {
        auto& particleTracer = particleTracers[i];
        particleTracer->Reset();
        particleTracer->tmpFacetVars.assign(model->sh.nbFacet, SimulationFacetTempVar());
        particleTracer->model = std::static_pointer_cast<MolflowSimulationModel>(model);
        particleTracer->totalDesorbed = 0;

        particleTracer->tmpParticleLog->clear();
    }
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

        particleTracer.tmpParticleLog->clear();

    }
    totalDesorbed = 0;
}
*/