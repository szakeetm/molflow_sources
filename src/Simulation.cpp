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

Simulation::Simulation(size_t nbThreads) : tMutex()
{
	totalDesorbed = 0;

    textTotalSize =
    profTotalSize =
    dirTotalSize =
    angleMapTotalSize =
    histogramTotalSize = 0;

    lastLogUpdateOK = true;


    for(auto& particle : currentParticles)
        particle.lastHitFacet = nullptr;

    hasVolatile = false;

	model.sh.nbSuper = 0;

	if(nbThreads == 0)
        this->nbThreads = omp_get_max_threads();
	else
	    this->nbThreads = nbThreads;
    currentParticles.resize(this->nbThreads, CurrentParticleStatus());// = CurrentParticleStatus();
}
Simulation::Simulation(Simulation&& o) noexcept : tMutex() {
    totalDesorbed = o.totalDesorbed;

    textTotalSize = o.textTotalSize;
    profTotalSize = o.profTotalSize;
    dirTotalSize = o.dirTotalSize;
    angleMapTotalSize = o.angleMapTotalSize;
    histogramTotalSize = o.histogramTotalSize;

    lastLogUpdateOK = o.lastLogUpdateOK;

    currentParticles.swap(o.currentParticles);
    for(auto& particle : o.currentParticles) {
        particle.lastHitFacet = nullptr;
    }
    hasVolatile =  o.hasVolatile;

    model =  o.model;

    nbThreads =  o.nbThreads;
}

Simulation::~Simulation()= default;

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
    std::vector<CurrentParticleStatus>(this->nbThreads).swap(this->currentParticles);
    /*this->model.structures.clear();
    this->model.tdParams.CDFs.clear();
    this->model.tdParams.IDs.clear();
    this->model.tdParams.moments.clear();
    this->model.tdParams.parameters.clear();
    //this->temperatures.clear();
    this->model.vertices3.clear();*/
}

// returns hit size or 0 on error
size_t Simulation::LoadSimulation() {
    double t0 = GetTick();

    //SetState(PROCESS_STARTING, "Clearing previous simulation");
    ClearSimulation();

    //SetState(PROCESS_STARTING, "Loading simulation");
/*
    {
        std::string inputString(loader->size,'\0');
        BYTE* buffer = (BYTE*)loader->buff;
        std::copy(buffer, buffer + loader->size, inputString.begin());
        std::stringstream inputStream;
        inputStream << inputString;
        cereal::BinaryInputArchive inputArchive(inputStream);

        //Worker params
        inputArchive(model.wp);
        inputArchive(model.otfParams);
        inputArchive(model.tdParams.CDFs);
        inputArchive(model.tdParams.IDs);
        inputArchive(model.tdParams.parameters);
        //inputArchive(temperatures);
        inputArchive(model.tdParams.moments);
        //inputArchive(desorptionParameterIDs);

        //Geometry
        inputArchive(model.sh);
        inputArchive(model.vertices3);

        model.structures.resize(model.sh.nbSuper); //Create structures
        //model.tdParams.moments = moments;

        //Facets
        for (size_t i = 0; i < model.sh.nbFacet; i++) { //Necessary because facets is not (yet) a vector in the interface
            SubprocessFacet f;
            inputArchive(
                    f.sh,
                    f.indices,
                    f.vertices2,
                    f.outgassingMap,
                    f.angleMap.pdf,
                    f.textureCellIncrements
            );

            //Some initialization
            if (!f.InitializeOnLoad(i, model.tdParams.moments.size(), histogramTotalSize)) return false;
            // Increase size counters
            //histogramTotalSize += 0;
            angleMapTotalSize += f.angleMapSize;
            dirTotalSize += f.directionSize* (1 + model.tdParams.moments.size());
            profTotalSize += f.profileSize* (1 + model.tdParams.moments.size());
            textTotalSize += f.textureSize* (1 + model.tdParams.moments.size());

            hasVolatile |= f.sh.isVolatile;
            if ((f.sh.superDest || f.sh.isVolatile) && ((f.sh.superDest - 1) >= model.sh.nbSuper || f.sh.superDest < 0)) {
                // Geometry error
                //ClearSimulation();
                //ReleaseDataport(loader);
                std::ostringstream err;
                err << "Invalid structure (wrong link on F#" << i + 1 << ")";
                //SetErrorSub(err.str().c_str());
                std::cerr << err.str() << std::endl;
                return 0;
            }

            if (f.sh.superIdx == -1) { //Facet in all structures
                for (auto& s : model.structures) {
                    s.facets.push_back(f);
                }
            }
            else {
                model.structures[f.sh.superIdx].facets.push_back(f); //Assign to structure
            }
        }
    }//inputarchive goes out of scope, file released
*/
    // New GlobalSimuState structure for threads
    {
        GlobalSimuState tmpResults = GlobalSimuState();

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

        tmpGlobalResults.resize(nbThreads);
        for(auto& tmpGlobal : tmpGlobalResults)
            tmpGlobal = tmpResults;
    }

    // Init tmp vars per thread
    for(auto& particle : currentParticles)
        std::vector<SubProcessFacetTempVar>(model.sh.nbFacet).swap(particle.tmpFacetVars);

    //Reserve particle log
    if (model.otfParams.enableLogging)
        tmpParticleLog.reserve(model.otfParams.logLimit / model.otfParams.nbProcess);

    // Build all AABBTrees
    size_t maxDepth=0;
    for (auto& s : model.structures) {
        std::vector<SubprocessFacet*> facetPointers; facetPointers.reserve(s.facets.size());
        for (auto& f : s.facets) {
            facetPointers.push_back(&f);
        }
        s.aabbTree = BuildAABBTree(facetPointers, 0, maxDepth);
    }

    // Initialise simulation
    /*SimulationModel model;
    model.structures = structures;
    model.vertices3 = vertices3;
    model.otfParams = ontheflyParams;
    model.sh = sh;
    model.wp = wp;
    model.tdParams.CDFs = CDFs;
    model.tdParams.IDs = IDs;
    model.tdParams.parameters = parameters;
    model.tdParams.moments = moments;*/


    //if(!model.sh.name.empty())
    //loadOK = true;
    double t1 = GetTick();
    printf("  Load %s successful\n", model.sh.name.c_str());
    printf("  Geometry: %zd vertex %zd facets\n", model.vertices3.size(), model.sh.nbFacet);

    printf("  Geom size: %d bytes\n", /*(size_t)(buffer - bufferStart)*/0);
    printf("  Number of stucture: %zd\n", model.sh.nbSuper);
    printf("  Global Hit: %zd bytes\n", sizeof(GlobalHitBuffer));
    printf("  Facet Hit : %zd bytes\n", model.sh.nbFacet * sizeof(FacetHitBuffer));
    printf("  Texture   : %zd bytes\n", textTotalSize);
    printf("  Profile   : %zd bytes\n", profTotalSize);
    printf("  Direction : %zd bytes\n", dirTotalSize);

    printf("  Total     : %zd bytes\n", GetHitsSize());
    //printf("  Seed: %lu\n", randomGenerator.GetSeed());
    printf("  Loading time: %.3f ms\n", (t1 - t0)*1000.0);

    // calc hit port size
    /*std::ostringstream result;
    {
        cereal::BinaryOutputArchive outputArchive(result);

        outputArchive(
                cereal::make_nvp("GlobalHits", this->tmpGlobalResults[0].globalHits),
                cereal::make_nvp("GlobalHistograms", this->tmpGlobalResults[0].globalHistograms),
                cereal::make_nvp("FacetStates", this->tmpGlobalResults[0].facetStates)
        );
    }

    size_t hSize = result.str().size();//simulation->GetHitsSize();*/
    return 0;

}

void Simulation::UpdateHits(int prIdx, DWORD timeout) {
    if(!globState)
        return;
    if (!tMutex.try_lock_for(std::chrono::seconds (10))) {
        return;
    }
    //globState = tmpGlobalResults[0];
    UpdateMCHits(*globState, prIdx, model.tdParams.moments.size(), timeout);
    // only 1 , so no reduce necessary
    /*ParticleLoggerItem& globParticleLog = tmpParticleLog[0];
    if (dpLog) UpdateLog(dpLog, timeout);*/

    tMutex.unlock();
    //ResetTmpCounters();
    // only reset buffers 1..N-1
    // 0 = global buffer for reduce
    for(auto & tmpGlobalResult : tmpGlobalResults)
        tmpGlobalResult.Reset();
}

size_t Simulation::GetHitsSize() {
    return sizeof(GlobalHitBuffer) + model.wp.globalHistogramParams.GetDataSize() +
           textTotalSize + profTotalSize + dirTotalSize + angleMapTotalSize + histogramTotalSize
           + model.sh.nbFacet * sizeof(FacetHitBuffer) * (1+model.tdParams.moments.size());
}

void Simulation::ResetSimulation() {
    //currentParticles.clear();// = CurrentParticleStatus();
    std::vector<CurrentParticleStatus>(this->nbThreads).swap(this->currentParticles);

    totalDesorbed = 0;
    //ResetTmpCounters();
    for(auto& tmpResults : tmpGlobalResults)
        tmpResults.Reset();
    tmpParticleLog.clear();
}

/*bool Simulation::StartSimulation() {
    if (!currentParticles.lastHitFacet) StartFromSource();
    return (currentParticles.lastHitFacet != nullptr);
}*/

void Simulation::RecordHit(const int &type, const CurrentParticleStatus &currentParticle) {
    const int index = omp_get_thread_num();
    if(index == 0) {
        GlobalSimuState &tmpResults = *currentParticle.tmpState;
        if (tmpResults.globalHits.hitCacheSize < HITCACHESIZE) {
            tmpResults.globalHits.hitCache[tmpResults.globalHits.hitCacheSize].pos = currentParticle.position;
            tmpResults.globalHits.hitCache[tmpResults.globalHits.hitCacheSize].type = type;
            tmpResults.globalHits.hitCacheSize++;
        }
    }
}

void Simulation::RecordLeakPos(const CurrentParticleStatus &currentParticle) {
    // Source region check performed when calling this routine
    // Record leak for debugging

    const int index = omp_get_thread_num();
    if(index == 0) {
        RecordHit(HIT_REF, currentParticle);
        RecordHit(HIT_LAST, currentParticle);
        GlobalSimuState &tmpResults = *currentParticle.tmpState;
        if (tmpResults.globalHits.leakCacheSize < LEAKCACHESIZE) {
            tmpResults.globalHits.leakCache[tmpResults.globalHits.leakCacheSize].pos = currentParticle.position;
            tmpResults.globalHits.leakCache[tmpResults.globalHits.leakCacheSize].dir = currentParticle.direction;
            tmpResults.globalHits.leakCacheSize++;
        }
    }
}