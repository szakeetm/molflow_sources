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


    currentParticle.lastHitFacet = nullptr;

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

    currentParticle = o.currentParticle;
    currentParticle.lastHitFacet = nullptr;
    currentParticle.model = &model;

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

    std::vector<SubProcessFacetTempVar>(model.sh.nbFacet).swap(currentParticle.tmpFacetVars);
    currentParticle.tmpState.Reset();
    currentParticle.model = &model;

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

// returns hit size or 0 on error
size_t Simulation::LoadSimulation(char *loadStatus) {
    double t0 = GetTick();

    //SetState(PROCESS_STARTING, "Clearing previous simulation");
    strncpy(loadStatus, "Clearing previous simulation", 127);

    ClearSimulation();

    //SetState(PROCESS_STARTING, "Loading simulation");
    strncpy(loadStatus, "Loading simulation", 127);

    // New GlobalSimuState structure for threads
    {
        auto& tmpResults = currentParticle.tmpState;

        std::vector<FacetState>(model.sh.nbFacet).swap(tmpResults.facetStates);
        for(auto& s : model.structures){
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
                tmpResults.facetStates[i].momentResults = std::vector<FacetMomentSnapshot>(1 + model.tdParams.moments.size(), facetMomentTemplate);
                if (sFac.sh.anglemapParams.record)
                  tmpResults.facetStates[i].recordedAngleMapPdf = std::vector<size_t>(sFac.sh.anglemapParams.GetMapSize());
            }
        }

        //Global histogram
        FacetHistogramBuffer globalHistTemplate;
        globalHistTemplate.Resize(model.wp.globalHistogramParams);
        tmpResults.globalHistograms = std::vector<FacetHistogramBuffer>(1 + model.tdParams.moments.size(), globalHistTemplate);
        tmpResults.initialized = true;


        // Init tmp vars per thread
        std::vector<SubProcessFacetTempVar>(model.sh.nbFacet).swap(currentParticle.tmpFacetVars);
        //currentParticle.tmpState = *tmpResults;
        //delete tmpResults;
    }

    //Reserve particle log
    if (model.otfParams.enableLogging)
        tmpParticleLog.reserve(model.otfParams.logLimit / model.otfParams.nbProcess);

    // Build all AABBTrees
    size_t maxDepth=0;
    for (auto& s : model.structures) {
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
    currentParticle.model = &model;

    // Initialise simulation


    //if(!model.sh.name.empty())
    //loadOK = true;
    double t1 = GetTick();
    printf("  Load %s successful\n", model.sh.name.c_str());
    printf("  Geometry: %zd vertex %zd facets\n", model.vertices3.size(), model.sh.nbFacet);

    printf("  Geom size: %d bytes\n", /*(size_t)(buffer - bufferStart)*/0);
    printf("  Number of structure: %zd\n", model.sh.nbSuper);
    printf("  Global Hit: %zd bytes\n", sizeof(GlobalHitBuffer));
    printf("  Facet Hit : %zd bytes\n", model.sh.nbFacet * sizeof(FacetHitBuffer));
    printf("  Texture   : %zd bytes\n", textTotalSize);
    printf("  Profile   : %zd bytes\n", profTotalSize);
    printf("  Direction : %zd bytes\n", dirTotalSize);

    printf("  Total     : %zd bytes\n", GetHitsSize());
    //printf("  Seed: %lu\n", randomGenerator.GetSeed());
    printf("  Loading time: %.3f ms\n", (t1 - t0)*1000.0);

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
    currentParticle.Reset();
    std::vector<SubProcessFacetTempVar>(model.sh.nbFacet).swap(currentParticle.tmpFacetVars);
    currentParticle.model = &model;

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