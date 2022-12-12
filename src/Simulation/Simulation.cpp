#include "Simulation.h"
#include "IntersectAABB_shared.h"
#include <cstring>
#include <cereal/archives/binary.hpp>
#include <Helper/Chronometer.h>
#include <Helper/ConsoleLogger.h>
#if defined(USE_OLD_BVH)
// currently always have SuperStructure
#else
#include <RayTracing/KDTree.h>
#include <RayTracing/BVH.h>
#endif

/*SuperStructure::SuperStructure()
{
	aabbTree = NULL;
}

SuperStructure::~SuperStructure()
{
	SAFE_DELETE(aabbTree);
}*/

Simulation::Simulation() : tMutex()
{
	totalDesorbed = 0;

    lastLogUpdateOK = true;

    for(auto& particleTracer : particleTracers) {
        particleTracer.lastHitFacet = nullptr;
        particleTracer.ray.lastIntersected = -1;
    }

    hasVolatile = false;

    globState = nullptr;
    globParticleLog = nullptr;

    //currentParticles.resize(1, CurrentParticleStatus());// = CurrentParticleStatus();
}
Simulation::Simulation(Simulation&& o) noexcept : tMutex() {

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

    globState = o.globState;
    globParticleLog = o.globParticleLog;

}

int Simulation::ReinitializeParticleLog() {
    /*tmpParticleLog.clear();
    tmpParticleLog.shrink_to_fit();
    if (model->otfParams.enableLogging) {
        tmpParticleLog.reserve(model->otfParams.logLimit*//* / model->otfParams.nbProcess*//*);
    }*/

    for(auto& particleTracer : particleTracers) {
        if(!particleTracer.tmpParticleLog.tMutex.try_lock_for(std::chrono::seconds(10)))
           return -1;
        particleTracer.tmpParticleLog.clear();
        particleTracer.tmpParticleLog.pLog.shrink_to_fit();
        if (model->otfParams.enableLogging) {
           particleTracer.tmpParticleLog.pLog.reserve(model->otfParams.logLimit/* / model->otfParams.nbProcess*/);
        }
        particleTracer.tmpParticleLog.tMutex.unlock();
    }
    return 0;
}

MFSim::ParticleTracer * Simulation::GetParticleTracer(size_t i) {
    if(i < particleTracers.size())
        return &particleTracers.at(i);
    else
        return nullptr;
}

void Simulation::SetNParticle(size_t n, bool fixedSeed) {
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

std::pair<int, std::optional<std::string>> Simulation::SanityCheckModel(bool strictCheck) {
    std::string errLog = "[Error Log on Check]\n";
    int errorsOnCheck = 0;

    if (!model->initialized) {
        errLog.append("Model not initialized\n");
        errorsOnCheck++;
    }
    if (model->vertices3.empty()) {
        errLog.append("Loaded empty vertex list\n");
        errorsOnCheck++;
    }
    if (model->facets.empty()) {
        errLog.append("Loaded empty facet list\n");
        errorsOnCheck++;
    }
    if(model->sh.nbFacet != model->facets.size()) {
        char tmp[256];
        snprintf(tmp, 256, "Facet structure not properly initialized, size mismatch: %zu / %zu\n", model->sh.nbFacet, model->facets.size());
        errLog.append(tmp);
        errorsOnCheck++;
    }
    for(auto& fac : model->facets){
        bool hasAnyTexture = fac->sh.countDes || fac->sh.countAbs || fac->sh.countRefl || fac->sh.countTrans || fac->sh.countACD || fac->sh.countDirection;
        if (!fac->sh.isTextured && (fac->sh.texHeight * fac->sh.texHeight > 0)) {
            char tmp[256];
            snprintf(tmp, 256, "[Fac #%zu] Untextured facet with texture size\n", fac->globalId+1);
            errLog.append(tmp);
            if(errLog.size() > 1280) errLog.resize(1280);
            errorsOnCheck++;
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
            errLog.append(tmp);
            if(errLog.size() > 1920) errLog.resize(1920);
            errorsOnCheck++;
        }
    }

    //Molflow unique
    if (model->wp.enableDecay && model->wp.halfLife <= 0.0) {
        char tmp[255];
        sprintf(tmp, "Particle decay is set, but half life was not set [= %e]\n", model->wp.halfLife);
        errLog.append(tmp);
        errorsOnCheck++;
    }

    if(!globState){
        errLog.append("No global simulation state set\n");
        errorsOnCheck++;
    }
    else if(!globState->initialized){
        errLog.append("Global simulation state not initialized\n");
        errorsOnCheck++;
    }

    if(errorsOnCheck){
        Log::console_error("{}", errLog);
    }
    return std::make_pair(errorsOnCheck, (errorsOnCheck > 0 ? std::make_optional(errLog) : std::nullopt)); // 0 = all ok
}

void Simulation::ClearSimulation() {

    //loadOK = false;

    //this->currentParticles.clear();// = CurrentParticleStatus();
    //std::vector<CurrentParticleStatus>(this->nbThreads).swap(this->currentParticles);
    for(auto& particleTracer : particleTracers) {
        particleTracer.tmpFacetVars.assign(model->sh.nbFacet, SimulationFacetTempVar());
        particleTracer.tmpState.Reset();
        particleTracer.model = (MolflowSimulationModel*) model.get();
        particleTracer.totalDesorbed = 0;

        particleTracer.tmpParticleLog.clear();

    }
    totalDesorbed = 0;
    //ResetTmpCounters();
    /*for(auto& tmpResults : tmpGlobalResults)
        tmpResults.Reset();*/
    //tmpParticleLog.clear();

    /*this->model->structures.clear();
    this->model->tdParams.CDFs.clear();
    this->model->tdParams.IDs.clear();
    this->model->tdParams.moments.clear();
    this->model->tdParams.parameters.clear();
    //this->temperatures.clear();
    this->model->vertices3.clear();*/
}

int Simulation::RebuildAccelStructure() {
    Chronometer timer;
    timer.Start();

    if(model->BuildAccelStructure(globState, BVH, BVHAccel::SplitMethod::SAH, 2))
        return 1;

    for(auto& particleTracer : particleTracers)
        particleTracer.model = (MolflowSimulationModel*) model.get();

    timer.Stop();

    return 0;
}



size_t Simulation::LoadSimulation(std::string& loadStatus) {
    Chronometer timer;
    timer.Start();
    loadStatus="Clearing previous simulation";
    ClearSimulation();
    loadStatus="Loading simulation";
    
    auto* simModel = (MolflowSimulationModel*) model.get();

    // New GlobalSimuState structure for threads
    for(auto& particleTracer : particleTracers)
    {
        auto& tmpResults = particleTracer.tmpState;
        tmpResults.Resize(model);

        // Init tmp vars per thread
        particleTracer.tmpFacetVars.assign(simModel->sh.nbFacet, SimulationFacetTempVar());

        //currentParticle.tmpState = *tmpResults;
        //delete tmpResults;
    }

    //Reserve particle log
    ReinitializeParticleLog();

    // Build ADS
    RebuildAccelStructure();

    // Initialise simulation


    //if(!model->sh.name.empty())
    //loadOK = true;
    timer.Stop();

    Log::console_msg_master(3, "  Load {} successful\n", simModel->sh.name);
    Log::console_msg_master(3, "  Geometry: {} vertex {} facets\n", simModel->vertices3.size(), simModel->sh.nbFacet);

    Log::console_msg_master(3, "  Geom size: {} bytes\n", simModel->size());
    Log::console_msg_master(3, "  Number of structure: {}\n", simModel->sh.nbSuper);
    Log::console_msg_master(3, "  Global Hit: {} bytes\n", sizeof(GlobalHitBuffer));
    Log::console_msg_master(3, "  Facet Hit : {} bytes\n", simModel->sh.nbFacet * sizeof(FacetHitBuffer));
/*        printf("  Texture   : %zd bytes\n", textTotalSize);
    printf("  Profile   : %zd bytes\n", profTotalSize);
    printf("  Direction : %zd bytes\n", dirTotalSize);*/

    Log::console_msg_master(3, "  Total     : {} bytes\n", GetHitsSize());
    for(auto& particleTracer : particleTracers)
        Log::console_msg_master(5, "  Seed for {}: {}\n", particleTracer.particleTracerId, particleTracer.randomGenerator.GetSeed());
    Log::console_msg_master(3, "  Loading time: {:.2f} ms\n", timer.ElapsedMs());

    return 0;
}

size_t Simulation::GetHitsSize() {
    MolflowSimulationModel* simModel = (MolflowSimulationModel*) model.get();
    return sizeof(GlobalHitBuffer) + model->wp.globalHistogramParams.GetDataSize() +
           + model->sh.nbFacet * sizeof(FacetHitBuffer) * (1+simModel->tdParams.moments.size());
}



void Simulation::ResetSimulation() {
    //currentParticles.clear();// = CurrentParticleStatus();
    //std::vector<CurrentParticleStatus>(this->nbThreads).swap(this->currentParticles);

    for(auto& particleTracer : particleTracers) {
        particleTracer.Reset();
        particleTracer.tmpFacetVars.assign(model->sh.nbFacet, SimulationFacetTempVar());
        particleTracer.model = (MolflowSimulationModel*) model.get();
        particleTracer.totalDesorbed = 0;

        particleTracer.tmpParticleLog.clear();
    }

    totalDesorbed = 0;
    //tmpParticleLog.clear();
}

/*bool Simulation::StartSimulation() {
    if (!currentParticles.lastHitFacet) StartFromSource();
    return (currentParticles.lastHitFacet != nullptr);
}*/