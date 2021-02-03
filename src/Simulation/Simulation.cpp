#include "Simulation.h"
#include "IntersectAABB_shared.h"
#include "Parameter.h"
#include <cstring>
#include <sstream>
#include <cereal/archives/binary.hpp>
#include <omp.h>
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

    textTotalSize =
    profTotalSize =
    dirTotalSize =
    angleMapTotalSize =
    histogramTotalSize = 0;

    lastLogUpdateOK = true;

    for(auto& particle : particles)
        particle.lastHitFacet = nullptr;

    hasVolatile = false;

	model.sh.nbSuper = 0;
    globState = nullptr;

    //currentParticles.resize(1, CurrentParticleStatus());// = CurrentParticleStatus();
}
Simulation::Simulation(Simulation&& o) noexcept : tMutex() {

    totalDesorbed = o.totalDesorbed;

    textTotalSize = o.textTotalSize;
    profTotalSize = o.profTotalSize;
    dirTotalSize = o.dirTotalSize;
    angleMapTotalSize = o.angleMapTotalSize;
    histogramTotalSize = o.histogramTotalSize;

    lastLogUpdateOK = o.lastLogUpdateOK;

    model = o.model;

    particles = o.particles;
    for(auto& particle : particles) {
        particle.lastHitFacet = nullptr;
        particle.model = &model;
    }

    hasVolatile =  o.hasVolatile;

    globState = o.globState;
}

int Simulation::ReinitializeParticleLog() {
    tmpParticleLog.clear();
    tmpParticleLog.shrink_to_fit();
    if (model.otfParams.enableLogging)
        tmpParticleLog.reserve(model.otfParams.logLimit / model.otfParams.nbProcess);

    return 0;
}

/*bool Simulation::UpdateOntheflySimuParams(Dataport *loader) {
    // Connect the dataport


    if (!AccessDataportTimed(loader, 2000)) {
        //SetErrorSub("Failed to connect to loader DP");
        std::cerr << "Failed to connect to loader DP" << std::endl;
        return false;
    }
    std::string inputString(loader->size,'\0');
    BYTE* buffer = (BYTE*)loader->buff;
    std::copy(buffer, buffer + loader->size, inputString.begin());
    std::stringstream inputStream;
    inputStream << inputString;
    cereal::BinaryInputArchive inputArchive(inputStream);

    inputArchive(model.otfParams);

    ReleaseDataport(loader);

    return true;
}*/

// Global handles
//extern Simulation* sHandle; //Declared at molflowSub.cpp



/*void InitSimulation() {

	// Global handle allocation
	sHandle = new Simulation();
	InitTick();
}*/

int Simulation::SanityCheckGeom() {

    return 0; // all ok
}

void Simulation::ClearSimulation() {

    //loadOK = false;

    textTotalSize =
    profTotalSize =
    dirTotalSize =
    angleMapTotalSize =
    histogramTotalSize = 0;

    //this->currentParticles.clear();// = CurrentParticleStatus();
    //std::vector<CurrentParticleStatus>(this->nbThreads).swap(this->currentParticles);
    for(auto& particle : particles) {
        std::vector<SubProcessFacetTempVar>(model.sh.nbFacet).swap(particle.tmpFacetVars);
        particle.tmpState.Reset();
        particle.model = &model;
        particle.totalDesorbed = 0;
    }
    totalDesorbed = 0;
    //ResetTmpCounters();
    /*for(auto& tmpResults : tmpGlobalResults)
        tmpResults.Reset();*/
    tmpParticleLog.clear();

    /*this->model.structures.clear();
    this->model.tdParams.CDFs.clear();
    this->model.tdParams.IDs.clear();
    this->model.tdParams.moments.clear();
    this->model.tdParams.parameters.clear();
    //this->temperatures.clear();
    this->model.vertices3.clear();*/
}

size_t Simulation::LoadSimulation(SimulationModel *simModel, char *loadStatus) {
    double t0 = GetTick();

    strncpy(loadStatus, "Clearing previous simulation", 127);
    ClearSimulation();
    strncpy(loadStatus, "Loading simulation", 127);
    
    if(!simModel) simModel = &this->model;
    
    // New GlobalSimuState structure for threads
    for(auto& particle : particles)
    {
        auto& tmpResults = particle.tmpState;

        std::vector<FacetState>(simModel->sh.nbFacet).swap(tmpResults.facetStates);
        for(auto& s : simModel->structures){
            for(auto& sFac : s.facets){
                size_t i = sFac.globalId;
                if(!tmpResults.facetStates[i].momentResults.empty())
                    continue; // Skip multiple init when facets exist in all structures
                FacetMomentSnapshot facetMomentTemplate;
                facetMomentTemplate.histogram.Resize(sFac.sh.facetHistogramParams);
                facetMomentTemplate.direction = std::vector<DirectionCell>(sFac.sh.countDirection ? sFac.sh.texWidth*sFac.sh.texHeight : 0);
                facetMomentTemplate.profile = std::vector<ProfileSlice>(sFac.sh.isProfile ? PROFILE_SIZE : 0);
                facetMomentTemplate.texture = std::vector<TextureCell>(sFac.sh.isTextured ? sFac.sh.texWidth*sFac.sh.texHeight : 0);
                //No init for hits
                tmpResults.facetStates[i].momentResults = std::vector<FacetMomentSnapshot>(1 + simModel->tdParams.moments.size(), facetMomentTemplate);
                if (sFac.sh.anglemapParams.record)
                  tmpResults.facetStates[i].recordedAngleMapPdf = std::vector<size_t>(sFac.sh.anglemapParams.GetMapSize());
            }
        }

        //Global histogram
        FacetHistogramBuffer globalHistTemplate;
        globalHistTemplate.Resize(simModel->wp.globalHistogramParams);
        tmpResults.globalHistograms = std::vector<FacetHistogramBuffer>(1 + simModel->tdParams.moments.size(), globalHistTemplate);
        tmpResults.initialized = true;


        // Init tmp vars per thread
        std::vector<SubProcessFacetTempVar>(simModel->sh.nbFacet).swap(particle.tmpFacetVars);
        //currentParticle.tmpState = *tmpResults;
        //delete tmpResults;
    }

    //Reserve particle log
    if (simModel->otfParams.enableLogging)
        tmpParticleLog.reserve(simModel->otfParams.logLimit / simModel->otfParams.nbProcess);

    // Build all AABBTrees
    size_t maxDepth=0;
    for (auto& s : simModel->structures) {
        if(s.aabbTree) s.aabbTree.reset();
        std::vector<SubprocessFacet*> facetPointers; facetPointers.reserve(s.facets.size());
        for (auto& f : s.facets) {
            facetPointers.push_back(&f);
        }
        AABBNODE* tree = BuildAABBTree(facetPointers, 0, maxDepth);
        s.aabbTree = std::make_shared<AABBNODE>(*tree);
        tree = nullptr;
        //delete tree; // pointer unnecessary because of make_shared
    }
    for(auto& particle : particles)
        particle.model = simModel;

    // Initialise simulation


    //if(!model.sh.name.empty())
    //loadOK = true;
    double t1 = GetTick();
    if(omp_get_thread_num() == 0) {
        printf("  Load %s successful\n", simModel->sh.name.c_str());
        printf("  Geometry: %zd vertex %zd facets\n", simModel->vertices3.size(), simModel->sh.nbFacet);

        printf("  Geom size: %zu bytes\n", simModel->size());
        printf("  Number of structure: %zd\n", simModel->sh.nbSuper);
        printf("  Global Hit: %zd bytes\n", sizeof(GlobalHitBuffer));
        printf("  Facet Hit : %zd bytes\n", simModel->sh.nbFacet * sizeof(FacetHitBuffer));
/*        printf("  Texture   : %zd bytes\n", textTotalSize);
        printf("  Profile   : %zd bytes\n", profTotalSize);
        printf("  Direction : %zd bytes\n", dirTotalSize);*/

        printf("  Total     : %zd bytes\n", GetHitsSize());
        //printf("  Seed: %lu\n", randomGenerator.GetSeed());
        printf("  Loading time: %.3f ms\n", (t1 - t0) * 1000.0);
    }
    return 0;
}

bool CurrentParticleStatus::UpdateHits(GlobalSimuState* globState, DWORD timeout) {
    if(!globState) {
        return false;
    }

    //globState = tmpGlobalResults[0];
    bool lastHitUpdateOK = UpdateMCHits(*globState, model->tdParams.moments.size(), timeout);
    // only 1 , so no reduce necessary
    /*ParticleLoggerItem& globParticleLog = tmpParticleLog[0];
    if (dpLog) UpdateLog(dpLog, timeout);*/

    // At last delete tmpCache
    tmpState.Reset();
    //ResetTmpCounters();
    // only reset buffers 1..N-1
    // 0 = global buffer for reduce
    /*for(auto & tmpGlobalResult : tmpGlobalResults)
        tmpGlobalResult.Reset();*/

    return lastHitUpdateOK;
}

size_t Simulation::GetHitsSize() {
    return sizeof(GlobalHitBuffer) + model.wp.globalHistogramParams.GetDataSize() +
           textTotalSize + profTotalSize + dirTotalSize + angleMapTotalSize + histogramTotalSize
           + model.sh.nbFacet * sizeof(FacetHitBuffer) * (1+model.tdParams.moments.size());
}

void CurrentParticleStatus::Reset() {
    position = Vector3d();
    direction = Vector3d();
    oriRatio = 0.0;

    nbBounces = 0;
    lastMomentIndex = 0;
    particleId = 0;
    distanceTraveled = 0;
    generationTime = 0;
    particleTime = 0;
    teleportedFrom = -1;

    velocity = 0.0;
    expectedDecayMoment = 0.0;
    structureId = -1;

    tmpState.Reset();
    lastHitFacet = nullptr;
    randomGenerator.SetSeed(GetSeed());
    model = nullptr;
    transparentHitBuffer.clear();
    tmpFacetVars.clear();
}

void Simulation::ResetSimulation() {
    //currentParticles.clear();// = CurrentParticleStatus();
    //std::vector<CurrentParticleStatus>(this->nbThreads).swap(this->currentParticles);

    for(auto& particle : particles) {
        particle.Reset();
        std::vector<SubProcessFacetTempVar>(model.sh.nbFacet).swap(particle.tmpFacetVars);
        particle.model = &model;
        particle.totalDesorbed = 0;
    }

    totalDesorbed = 0;
    tmpParticleLog.clear();
}

/*bool Simulation::StartSimulation() {
    if (!currentParticles.lastHitFacet) StartFromSource();
    return (currentParticles.lastHitFacet != nullptr);
}*/

void CurrentParticleStatus::RecordHit(const int &type) {
    if (tmpState.globalHits.hitCacheSize < HITCACHESIZE) {
        tmpState.globalHits.hitCache[tmpState.globalHits.hitCacheSize].pos = position;
        tmpState.globalHits.hitCache[tmpState.globalHits.hitCacheSize].type = type;
        ++tmpState.globalHits.hitCacheSize;
    }
}

void CurrentParticleStatus::RecordLeakPos() {
    // Source region check performed when calling this routine
    // Record leak for debugging
    RecordHit(HIT_REF);
    RecordHit(HIT_LAST);
    if (tmpState.globalHits.leakCacheSize < LEAKCACHESIZE) {
        tmpState.globalHits.leakCache[tmpState.globalHits.leakCacheSize].pos = position;
        tmpState.globalHits.leakCache[tmpState.globalHits.leakCacheSize].dir = direction;
        ++tmpState.globalHits.leakCacheSize;
    }
}