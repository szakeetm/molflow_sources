#include "Simulation.h"
#include "IntersectAABB_shared.h"
#include <cstring>
#include <sstream>
#include <cereal/archives/binary.hpp>
#include <Helper/Chronometer.h>
#include <Helper/ConsoleLogger.h>
#if defined(USE_OLD_BVH)
// currently always have SuperStructure
#include "AABB.h"
#endif

#include <RayTracing/KDTree.h>
#include <RayTracing/BVH.h>
#include <Helper/FormatHelper.h>
#include <omp.h>
#include <cmath>

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

    for(auto& particle : particles) {
        particle.lastHitFacet = nullptr;
        particle.particle.lastIntersected = -1;
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

    particles = o.particles;
    for(auto& particle : particles) {
        particle.lastHitFacet = nullptr;
        particle.particle.lastIntersected = -1;
        particle.model = model.get();
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

    for(auto& particle : particles) {
        if(!particle.tmpParticleLog.tMutex.try_lock_for(std::chrono::seconds(10)))
           return -1;
        particle.tmpParticleLog.clear();
        particle.tmpParticleLog.pLog.shrink_to_fit();
        if (model->otfParams.enableLogging) {
           particle.tmpParticleLog.pLog.reserve(model->otfParams.logLimit/* / model->otfParams.nbProcess*/);
        }
        particle.tmpParticleLog.tMutex.unlock();
    }
    return 0;
}

MFSim::Particle * Simulation::GetParticle(size_t i) {
    if(i < particles.size())
        return &particles.at(i);
    else
        return nullptr;
}

void Simulation::SetNParticle(size_t n, bool fixedSeed) {
    particles.clear();
    particles.resize(n);
    size_t pid = 0;
    for(auto& particle : particles){
        if(fixedSeed)
         particle.randomGenerator.SetSeed(42424242 + pid);
        else
         particle.randomGenerator.SetSeed(GenerateSeed(pid));
        particle.particleId = pid++;
    }
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

    inputArchive(model->otfParams);

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
            snprintf(tmp, 256, "[Fac #%zu] Untextured facet with texture size\n", fac->globalId);
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
            snprintf(tmp, 256, "[Fac #%zu] Untextured facet with texture counters\n", fac->globalId);
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
        printf("%s", errLog.c_str());
    }
    return std::make_pair(errorsOnCheck, (errorsOnCheck > 0 ? std::make_optional(errLog) : std::nullopt)); // 0 = all ok
}

void Simulation::ClearSimulation() {

    //loadOK = false;

    //this->currentParticles.clear();// = CurrentParticleStatus();
    //std::vector<CurrentParticleStatus>(this->nbThreads).swap(this->currentParticles);
    for(auto& particle : particles) {
        particle.tmpFacetVars.assign(model->sh.nbFacet, SubProcessFacetTempVar());
        particle.tmpState.Reset();
        particle.model = model.get();
        particle.totalDesorbed = 0;

        particle.tmpParticleLog.clear();

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
    Chronometer timer(false);
    timer.Start();

    if(model->BuildAccelStructure(globState, model->wp.accel_type, model->wp.splitMethod, model->wp.bvhMaxPrimsInNode,
                                  model->wp.hybridWeight))
        return 1;

    for(auto& particle : particles)
        particle.model = model.get();

    timer.Stop();

    Log::console_msg(4, "Rebuilt Acceleration Structure in %lfs\n", timer.Elapsed());
    return 0;
}



size_t Simulation::LoadSimulation(char *loadStatus) {
    Chronometer timer(false);
    timer.Start();
    strncpy(loadStatus, "Clearing previous simulation", 127);
    ClearSimulation();
    strncpy(loadStatus, "Loading simulation", 127);
    
    auto simModel = this->model;
    // New GlobalSimuState structure for threads
    for(auto& particle : particles)
    {
        auto& tmpResults = particle.tmpState;

        tmpResults.facetStates.assign(model->sh.nbFacet, FacetState());
        for(auto& fac : simModel->facets){
            auto& sFac = *fac;
            size_t i = sFac.globalId;
            if(!tmpResults.facetStates[i].momentResults.empty())
                continue; // Skip multiple init when facets exist in all structures
            FacetMomentSnapshot facetMomentTemplate{};
            facetMomentTemplate.histogram.Resize(sFac.sh.facetHistogramParams);
            facetMomentTemplate.direction.assign((sFac.sh.countDirection ? sFac.sh.texWidth*sFac.sh.texHeight : 0), DirectionCell());
            facetMomentTemplate.profile.assign((sFac.sh.isProfile ? PROFILE_SIZE : 0), ProfileSlice());
            facetMomentTemplate.texture.assign((sFac.sh.isTextured ? sFac.sh.texWidth*sFac.sh.texHeight : 0), TextureCell());

            //No init for hits
            tmpResults.facetStates[i].momentResults.assign(1 + simModel->tdParams.moments.size(), facetMomentTemplate);
            if (sFac.sh.anglemapParams.record)
                tmpResults.facetStates[i].recordedAngleMapPdf.assign(sFac.sh.anglemapParams.GetMapSize(), 0);
        }

        //Global histogram
        FacetHistogramBuffer globalHistTemplate{};
        globalHistTemplate.Resize(simModel->wp.globalHistogramParams);
        tmpResults.globalHistograms.assign(1 + simModel->tdParams.moments.size(), globalHistTemplate);
        tmpResults.initialized = true;


        // Init tmp vars per thread
        particle.tmpFacetVars.assign(simModel->sh.nbFacet, SubProcessFacetTempVar());
        particle.tmpState.hitBattery.resize(simModel->sh.nbFacet);

        //currentParticle.tmpState = *tmpResults;
        //delete tmpResults;
    }

    //Reserve particle log
    ReinitializeParticleLog();

    std::vector<std::vector<std::shared_ptr<Primitive>>> primPointers;
    primPointers.resize(simModel->sh.nbSuper);
    for(auto& sFac : simModel->facets){
        if (sFac->sh.superIdx == -1) { //Facet in all structures
            for (auto& fp_vec : primPointers) {
                fp_vec.push_back(sFac);
            }
        }
        else {
            primPointers[sFac->sh.superIdx].push_back(sFac); //Assign to structure
        }
    }

    for(auto& sFac : simModel->facets){
        if (sFac->sh.opacity_paramId == -1){ //constant sticking
            sFac->sh.opacity = std::clamp(sFac->sh.opacity, 0.0, 1.0);
            sFac->surf = simModel->GetSurface(sFac->sh.opacity);
        }
        else {
            auto* par = &simModel->tdParams.parameters[sFac->sh.opacity_paramId];
            sFac->surf = simModel->GetParameterSurface(sFac->sh.opacity_paramId, par);
        }
    }

#if defined(USE_OLD_BVH)
    std::vector<std::vector<SubprocessFacet*>> facetPointers;
    facetPointers.resize(simModel->sh.nbSuper);
    for(auto& sFac : simModel->facets){
        // TODO: Build structures
        if (sFac->sh.superIdx == -1) { //Facet in all structures
            for (auto& fp_vec : facetPointers) {
                fp_vec.push_back(sFac.get());
            }
        }
        else {
            facetPointers[sFac->sh.superIdx].push_back(sFac.get()); //Assign to structure
        }
    }

    // Build all AABBTrees
    size_t maxDepth=0;
    for (size_t s = 0; s < simModel->sh.nbSuper; ++s) {
        auto& structure = simModel->structures[s];
        if(structure.aabbTree)
            structure.aabbTree.reset();
        AABBNODE* tree = BuildAABBTree(facetPointers[s], 0, maxDepth);
        structure.aabbTree = std::make_shared<AABBNODE>(*tree);
        //delete tree; // pointer unnecessary because of make_shared
    }

#endif

    RebuildAccelStructure();/*
    simModel->accel.clear();
    for (size_t s = 0; s < simModel->sh.nbSuper; ++s) {
        if(model->wp.accel_type == 1)
            simModel->accel.emplace_back(std::make_shared<KdTreeAccel>(KdTreeAccel::SplitMethod::SAH, primPointers[s]));
        else
            simModel->accel.emplace_back(std::make_shared<BVHAccel>(primPointers[s], 2, BVHAccel::SplitMethod::SAH));
    }
// old_bvb
    for(auto& particle : particles)
        particle.model = model.get();*/

    // Initialise simulation


    //if(!model->sh.name.empty())
    //loadOK = true;
    timer.Stop();

    Log::console_msg_master(3, "  Load %s successful\n", simModel->sh.name.c_str());
    Log::console_msg_master(3, "  Geometry: %zd vertex %zd facets\n", simModel->vertices3.size(), simModel->sh.nbFacet);

    Log::console_msg_master(3, "  Geom size: %zu bytes\n", simModel->size());
    Log::console_msg_master(3, "  Number of structure: %zd\n", simModel->sh.nbSuper);
    Log::console_msg_master(3, "  Global Hit: %zd bytes\n", sizeof(GlobalHitBuffer));
    Log::console_msg_master(3, "  Facet Hit : %zd bytes\n", simModel->sh.nbFacet * sizeof(FacetHitBuffer));
/*        printf("  Texture   : %zd bytes\n", textTotalSize);
    printf("  Profile   : %zd bytes\n", profTotalSize);
    printf("  Direction : %zd bytes\n", dirTotalSize);*/

    Log::console_msg_master(3, "  Total     : %zd bytes\n", GetHitsSize());
    for(auto& particle : particles)
        Log::console_msg_master(5, "  Seed for %2zu: %lu\n", particle.particleId, particle.randomGenerator.GetSeed());
    Log::console_msg_master(3, "  Loading time: %.3f ms\n", timer.ElapsedMs());

    return 0;
}

size_t Simulation::GetHitsSize() {
    return sizeof(GlobalHitBuffer) + model->wp.globalHistogramParams.GetDataSize() +
           + model->sh.nbFacet * sizeof(FacetHitBuffer) * (1+model->tdParams.moments.size());
}



void Simulation::ResetSimulation() {
    //currentParticles.clear();// = CurrentParticleStatus();
    //std::vector<CurrentParticleStatus>(this->nbThreads).swap(this->currentParticles);

    // PROFILING
    if(model && model->initialized && !model->accel.empty()) {
        for (auto &accel: model->accel)
            accel->ResetStats();
    }
    // PROFILING -----

    for(auto& particle : particles) {
        particle.Reset();
        particle.tmpFacetVars.assign(model->sh.nbFacet, SubProcessFacetTempVar());
        particle.model = model.get();
        particle.totalDesorbed = 0;

        particle.tmpParticleLog.clear();
    }

    totalDesorbed = 0;
    //tmpParticleLog.clear();
}

bool Simulation::RunParallel(size_t nSteps) {
    bool eos;
    bool lastUpdateOk = false;

    //printf("Lim[%zu] %lu --> %lu\n",threadNum, localDesLimit, simulation->globState->globalHits.globalHits.hit.nbDesorbed);
// Calculate remaining work
    size_t desPerThread = 0;
    size_t remainder = 0;
    if(model->otfParams.desorptionLimit > 0){
        if(model->otfParams.desorptionLimit > (globState->globalHits.globalHits.nbDesorbed)) {
            size_t limitDes_global = model->otfParams.desorptionLimit;
            desPerThread = limitDes_global / particles.size();
            remainder = limitDes_global % particles.size();
        }
    }

    std::vector<size_t> localDesLimits(particles.size());
    for(auto& particle : particles){
        size_t localDes = (desPerThread > particle.totalDesorbed) ? desPerThread - particle.totalDesorbed : 0;
        localDesLimits[particle.particleId] = (particle.particleId < remainder) ? localDes + 1 : localDes;
    }

    bool simEos = false;
#pragma omp parallel for default(none) shared(simEos, localDesLimits, nSteps)
    for(int particleId = 0; particleId < particles.size(); particleId++){
        double timeStart = omp_get_wtime();
        double timeLoopStart = timeStart;
        double timeEnd;
        bool eos = false;
        bool lastUpdateOk = true;
        auto& particle = particles[particleId];
        double stepsPerSec = 1.0;
        do {
            size_t desorptions = localDesLimits[particleId];
            bool end = false;
            {
                size_t nbStep = (nSteps > 0) ? (nSteps / particles.size()) : ((stepsPerSec <= 0.0) ? 250.0 : std::ceil(stepsPerSec + 0.5));
// Check end of simulation
                bool goOn = true;
                size_t remainingDes = 0;

                if (particle.model->otfParams.desorptionLimit > 0) {
                    if (desorptions <= particle.tmpState.globalHits.globalHits.nbDesorbed){
                        //lastHitFacet = nullptr; // reset full particle status or go on from where we left
                        goOn = false;
                    }
                    else {
                        //if(particle->tmpState.globalHits.globalHits.hit.nbDesorbed <= (particle->model->otfParams.desorptionLimit - simulation->globState->globalHits.globalHits.hit.nbDesorbed)/ particle->model->otfParams.nbProcess)
                        remainingDes = desorptions - particle.tmpState.globalHits.globalHits.nbDesorbed;
                    }
                }
                //auto start_time = std::chrono::high_resolution_clock::now();
                if(goOn) {
                    double start_time = omp_get_wtime();
                    goOn = particle.SimulationMCStep(nbStep, particleId, remainingDes);
                    double end_time = omp_get_wtime();

                    //auto end_time = std::chrono::high_resolution_clock::now();

                    if(nSteps > 0)
                        goOn = false;
                    if (goOn) // don't update on end, this will give a false ratio (SimMCStep could return actual steps instead of plain "false"
                    {
                        const double elapsedTimeMs = (end_time - start_time); //std::chrono::duration<float, std::ratio<1, 1>>(end_time - start_time).count();
                        if (elapsedTimeMs != 0.0)
                            stepsPerSec = (1.0 * nbStep) / (elapsedTimeMs); // every 1.0 second
                        else
                            stepsPerSec = (100.0 * nbStep); // in case of fast initial run
                    }
                }
                end = !goOn;      // Run during 1 sec
            }
            if(end) {
#pragma omp critical
                {
                    simEos = true;
                }
            }

            timeEnd = omp_get_wtime();

            bool forceQueue = timeEnd - timeLoopStart > 60 ||
                    particleId == 0; // update after 60s of no update or when thread 0 is called
            if (true || forceQueue) {
                if (model->otfParams.desorptionLimit > 0) {
                    if (localDesLimits[particleId] > particle.tmpState.globalHits.globalHits.nbDesorbed)
                        localDesLimits[particleId] -= particle.tmpState.globalHits.globalHits.nbDesorbed;
                    else localDesLimits[particleId] = 0;
                }

                size_t timeOut = lastUpdateOk ? 0 : 100; //ms
                lastUpdateOk = particle.UpdateHits(globState, globParticleLog, timeOut); // Update hit with 20ms timeout. If fails, probably an other subprocess is updating, so we'll keep calculating and try it later (latest when the simulation is stopped).
                timeLoopStart = omp_get_wtime();
            } else {
                lastUpdateOk = false;
            }
            //printf("[%zu] PUP: %lu , %lu , %lu\n",threadNum, desorptions,localDesLimit, particle->tmpState.globalHits.globalHits.hit.nbDesorbed);
            eos = simEos || (particle.model->otfParams.timeLimit != 0 ? timeEnd - timeStart >=
                                                                               particle.model->otfParams.timeLimit
                                                                             : false);
        } while (!eos);

        if (!lastUpdateOk) {
            particle.UpdateHits(globState, globParticleLog, 20000); // Update hit with 20ms timeout. If fails, probably an other subprocess is updating, so we'll keep calculating and try it later (latest when the simulation is stopped).)
        }
    }

    return simEos;
}

void Simulation::FindBestADS() {

    const std::array<BVHAccel::SplitMethod,7> all_splits = {
            BVHAccel::SplitMethod::SAH, BVHAccel::SplitMethod::HLBVH, BVHAccel::SplitMethod::Middle,
            BVHAccel::SplitMethod::EqualCounts, BVHAccel::SplitMethod::MolflowSplit, BVHAccel::SplitMethod::ProbSplit,
            BVHAccel::SplitMethod::RDH
    };
    const auto runTimePerTest = 2e7;
    // 0. Test runs to gather test battery
    model->otfParams.raySampling = true;
    RunParallel(1e5);
    /*globState->UpdateBatteryFrequencies();

    RunParallel(1e6);*/
    model->otfParams.raySampling = false;
    for (auto &thr: particles) {
        thr.tmpState.StopBatteryChange();
    }
    // 1. BVH
    model->wp.accel_type = 0;
    model->wp.bvhMaxPrimsInNode = 2;

    Chronometer bench(false);

    for(BVHAccel::SplitMethod split : all_splits) {
        if(split == BVHAccel::SplitMethod::RDH)
            continue;
        bench.ReInit();
        model->wp.splitMethod = static_cast<int>(split);

        std::string split_str;
        {
            std::stringstream split_ss;
            split_ss << split;
            split_str = split_ss.str();
        }
        const char* splitName = split_str.c_str();
        Log::console_header(3, "Benchmarking BVH x %s\n", splitName);
        bench.Start();
        if(RebuildAccelStructure()) {
            Log::console_footer(3, "Built invalid BVH, ... Skip!\n", splitName, bench.Elapsed());
            continue;
        }
        Log::console_msg_master(2, "Built BVH x %s: %lf s\n", splitName, bench.Elapsed());
        bench.ReInit();
        bench.Start();
        RunParallel(runTimePerTest);
        Log::console_footer(2, "Benchmarked BVH x %s: %lf s\n", splitName, bench.Elapsed());
    }

    // 2. KD Trees
    const std::array<KdTreeAccel::SplitMethod,6> all_splits_kd = {
            KdTreeAccel::SplitMethod::SAH, KdTreeAccel::SplitMethod::ProbSplit, KdTreeAccel::SplitMethod::ProbHybrid, KdTreeAccel::SplitMethod::TestSplit,
            KdTreeAccel::SplitMethod::HybridSplit, KdTreeAccel::SplitMethod::HybridBin
    };

    model->wp.accel_type = 1;
    for(KdTreeAccel::SplitMethod split : all_splits_kd) {
        if(split == KdTreeAccel::SplitMethod::TestSplit)
            continue;
        bench.ReInit();
        model->wp.splitMethod = static_cast<int>(split);

        std::string split_str;
        {
            std::stringstream split_ss;
            split_ss << split;
            split_str = split_ss.str();
        }
        const char* splitName = split_str.c_str();
        Log::console_header(3, "Benchmarking KD x %s\n", splitName);
        bench.Start();
        if(RebuildAccelStructure()) {
            Log::console_footer(3, "Built invalid KD, ... Skip!\n", splitName, bench.Elapsed());
            continue;
        }
        Log::console_msg_master(2, "Built KD x %s: %lf s\n", splitName, bench.Elapsed());
        bench.ReInit();
        bench.Start();
        RunParallel(runTimePerTest);
        Log::console_footer(2, "Benchmarked KD x %s: %lf s\n", splitName, bench.Elapsed());
    }
}