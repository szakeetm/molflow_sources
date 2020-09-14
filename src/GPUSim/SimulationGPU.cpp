//
// Created by pbahr on 18/05/2020.
//

#include <cereal/archives/binary.hpp>
#include "SimulationGPU.h"
#include "ModelReader.h" // TempFacet

#define LAUNCHSIZE 1920*64*32//1024*64*16//1024*128*64

SimulationGPU::SimulationGPU()
        : SimulationUnit() {
    model = nullptr;
    totalDesorbed = 0;
    memset(&tmpGlobalResult, 0, sizeof(GlobalHitBuffer));
}

SimulationGPU::~SimulationGPU() {
    delete model;
}

int SimulationGPU::SanityCheckGeom() {
    return 0;
}

void SimulationGPU::ClearSimulation() {
    gpuSim.CloseSimulation();
}

bool SimulationGPU::LoadSimulation(Dataport *loader) {
    double t0 = GetTick();
    //SetState(PROCESS_STARTING, "Clearing previous simulation");
    ClearSimulation();

    //SetState(PROCESS_STARTING, "Loading simulation");

    {

        std::string inputString(loader->size,'\0');
        BYTE* buffer = (BYTE*)loader->buff;
        std::copy(buffer, buffer + loader->size, inputString.begin());

        model = flowgpu::loadFromSerialization(inputString);
        if(!model) return false;
        this->ontheflyParams = this->model->ontheflyParams;
    }//inputarchive goes out of scope, file released

    // Initialise simulation

    //TODO: Better sanity check for startup
    if(model->nbFacets_total>0)
        gpuSim.LoadSimulation(model,LAUNCHSIZE);
    //if(!sh.name.empty())

    double t1 = GetTick();
    printf("  Load %s successful\n", model->geomProperties.name.c_str());
    printf("  Loading time: %.3f ms\n", (t1 - t0)*1000.0);
    return true;
}

void SimulationGPU::ResetSimulation() {
    totalDesorbed = 0;
    ResetTmpCounters();
    gpuSim.ResetSimulation();
}

bool SimulationGPU::UpdateOntheflySimuParams(Dataport *loader) {
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

    inputArchive(ontheflyParams);
    model->ontheflyParams = ontheflyParams;

    ReleaseDataport(loader);

    return true;
}

void SimulationGPU::UpdateHits(Dataport *dpHit, Dataport* dpLog, int prIdx, DWORD timeout) {
    //UpdateMCHits(dpHit, prIdx, moments.size(), timeout);
    //if (dpLog) UpdateLog(dpLog, timeout);
    //std::cout << "#SimulationGPU: Updating hits"<<std::endl;
    gpuSim.GetSimulationData(false);
    GlobalCounter* globalCount = gpuSim.GetGlobalCounter();

    BYTE *buffer;
    GlobalHitBuffer *gHits;
    TEXTURE_MIN_MAX texture_limits_old[3];

#if defined(_DEBUG)
    double t0, t1;
    t0 = GetTick();
#endif
    //SetState(PROCESS_STARTING, "Waiting for 'hits' dataport access...", false, true);
    bool lastHitUpdateOK = AccessDataportTimed(dpHit, timeout);
    //SetState(PROCESS_STARTING, "Updating MC hits...", false, true);
    if (!lastHitUpdateOK) return; //Timeout, will try again later

    buffer = (BYTE *) dpHit->buff;
    gHits = (GlobalHitBuffer *) buffer;

    for(unsigned int i = 0; i < globalCount->facetHitCounters.size(); i++) {
        gHits->globalHits.hit.nbMCHit += globalCount->facetHitCounters[i].nbMCHit; // let misses count as 0 (-1+1)
        gHits->globalHits.hit.nbDesorbed += globalCount->facetHitCounters[i].nbDesorbed; // let misses count as 0 (-1+1)
        gHits->globalHits.hit.nbAbsEquiv += globalCount->facetHitCounters[i].nbAbsEquiv; // let misses count as 0 (-1+1)
        gHits->globalHits.hit.nbHitEquiv += globalCount->facetHitCounters[i].nbHitEquiv; // let misses count as 0 (-1+1)
        gHits->globalHits.hit.sum_v_ort += globalCount->facetHitCounters[i].sum_v_ort; // let misses count as 0 (-1+1)
        gHits->globalHits.hit.sum_1_per_velocity += globalCount->facetHitCounters[i].sum_1_per_velocity; // let misses count as 0 (-1+1)
        gHits->globalHits.hit.sum_1_per_ort_velocity += globalCount->facetHitCounters[i].sum_1_per_ort_velocity; // let misses count as 0 (-1+1)

        //if(data.facetHitCounters[i].nbMCHit>0 || data.facetHitCounters[i].nbDesorbed> 0 || data.facetHitCounters[i].nbAbsEquiv>0)
        // facetCounterEveryFile << (i/this->model->nbFacets_total) << " " << (i%this->model->nbFacets_total)+1 << " " << data.facetHitCounters[i].nbMCHit << " " << data.facetHitCounters[i].nbDesorbed << " " << static_cast<unsigned int>(data.facetHitCounters[i].nbAbsEquiv) << std::endl;
    }

    //Memorize current limits, then do a min/max search
    for (int i = 0; i < 3; i++) {
        texture_limits_old[i] = gHits->texture_limits[i];
        gHits->texture_limits[i].min.all = gHits->texture_limits[i].min.moments_only = HITMAX;
        gHits->texture_limits[i].max.all = gHits->texture_limits[i].max.moments_only = 0;
    }

#ifdef DEBUGPOS
    // HHit (Only prIdx 0)
    for (size_t hitIndex = 0; hitIndex < globalCount->positions.size(); hitIndex++) {
        gHits->hitCache[(hitIndex + gHits->lastHitIndex) % HITCACHESIZE].pos.x = globalCount->positions[hitIndex].x;
        gHits->hitCache[(hitIndex + gHits->lastHitIndex) % HITCACHESIZE].pos.y = globalCount->positions[hitIndex].y;
        gHits->hitCache[(hitIndex + gHits->lastHitIndex) % HITCACHESIZE].pos.z = globalCount->positions[hitIndex].z;
        gHits->hitCache[(hitIndex + gHits->lastHitIndex) % HITCACHESIZE].type = HIT_REF;
    }
    gHits->hitCache[gHits->lastHitIndex].type = HIT_REF; //Penup (border between blocks of consecutive hits in the hit cache)
    gHits->lastHitIndex = (gHits->lastHitIndex + globalCount->positions.size()) % HITCACHESIZE;
    gHits->hitCache[gHits->lastHitIndex].type = HIT_LAST; //Penup (border between blocks of consecutive hits in the hit cache)
    gHits->hitCacheSize = std::min(HITCACHESIZE, gHits->hitCacheSize + globalCount->positions.size());
#endif // DEBUGPOS

    // Leak
    gHits->nbLeakTotal += globalCount->leakCounter[0];
#ifdef DEBUGLEAKPOS
    for (size_t leakIndex = 0; leakIndex < globalCount->leakPositions.size(); leakIndex++) {
        gHits->leakCache[(leakIndex + gHits->lastLeakIndex) %
                         LEAKCACHESIZE].pos.x = globalCount->leakPositions[leakIndex].x;
        gHits->leakCache[(leakIndex + gHits->lastLeakIndex) %
                         LEAKCACHESIZE].pos.y = globalCount->leakPositions[leakIndex].y;
        gHits->leakCache[(leakIndex + gHits->lastLeakIndex) %
                         LEAKCACHESIZE].pos.z = globalCount->leakPositions[leakIndex].z;
        gHits->leakCache[(leakIndex + gHits->lastLeakIndex) %
                         LEAKCACHESIZE].dir.x = globalCount->leakDirections[leakIndex].x;
        gHits->leakCache[(leakIndex + gHits->lastLeakIndex) %
                         LEAKCACHESIZE].dir.y = globalCount->leakDirections[leakIndex].y;
        gHits->leakCache[(leakIndex + gHits->lastLeakIndex) %
                         LEAKCACHESIZE].dir.z = globalCount->leakDirections[leakIndex].z;
    }
    gHits->lastLeakIndex = (gHits->lastLeakIndex + globalCount->leakPositions.size()) % LEAKCACHESIZE;
    gHits->leakCacheSize = std::min(LEAKCACHESIZE, gHits->leakCacheSize + globalCount->leakPositions.size());
#endif // DEBUGLEAKPOS

    // Facets

    for(unsigned int i = 0; i < globalCount->facetHitCounters.size(); i++) {
        int realIndex = model->triangle_meshes[0]->poly[i].parentIndex;
        FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *) (buffer + model->tri_facetOffset[realIndex] /*+ sizeof(FacetHitBuffer)*/);
        facetHitBuffer->hit.nbAbsEquiv += globalCount->facetHitCounters[i].nbAbsEquiv;
        facetHitBuffer->hit.nbDesorbed += globalCount->facetHitCounters[i].nbDesorbed;
        facetHitBuffer->hit.nbMCHit += globalCount->facetHitCounters[i].nbMCHit;
        facetHitBuffer->hit.nbHitEquiv += globalCount->facetHitCounters[i].nbHitEquiv;
        facetHitBuffer->hit.sum_v_ort += globalCount->facetHitCounters[i].sum_v_ort;
        facetHitBuffer->hit.sum_1_per_velocity += globalCount->facetHitCounters[i].sum_1_per_velocity;
        facetHitBuffer->hit.sum_1_per_ort_velocity += globalCount->facetHitCounters[i].sum_1_per_ort_velocity;
    } // End nbFacet

    //profiles
    if(!globalCount->profiles.empty()) {
        double timeCorrection = model->wp.finalOutgassingRate;
        for (auto&[id, profiles] : globalCount->profiles) {
            for (auto &mesh : model->triangle_meshes) {
                int previousId = 0;
                for (auto &facet : mesh->poly) {
                    if ((facet.profProps.profileType != flowgpu::PROFILE_FLAGS::noProfile) && (id == facet.parentIndex)) {
                        ProfileSlice *shProfile = (ProfileSlice *) (buffer + (model->tri_facetOffset[id] + sizeof(FacetHitBuffer)));
                        for (unsigned int s = 0; s < PROFILE_SIZE; ++s) {
                            shProfile[s].countEquiv += profiles[s].countEquiv;
                            shProfile[s].sum_v_ort += profiles[s].sum_v_ort_per_area;
                            shProfile[s].sum_1_per_ort_velocity += profiles[s].sum_1_per_ort_velocity;
                        }

                        break; //Only need 1 facet for texture position data
                    }
                }
            }
        }
    }

    //textures
    if(!globalCount->textures.empty()) {
        double timeCorrection = model->wp.finalOutgassingRate;
        for (auto&[id, texels] : globalCount->textures) {
            for (auto &mesh : model->triangle_meshes) {
                int previousId = 0;
                for (auto &facet : mesh->poly) {
                    if ((facet.texProps.textureFlags) && (id == facet.parentIndex)) {
                        int bufferOffset_profSize = facet.profProps.profileType ? PROFILE_SIZE * sizeof(ProfileSlice) : 0;
                        TextureCell *shTexture = (TextureCell *) (buffer + (model->tri_facetOffset[id] + sizeof(FacetHitBuffer) + bufferOffset_profSize));
                        unsigned int width = model->facetTex[facet.texProps.textureOffset].texWidth;
                        unsigned int height = model->facetTex[facet.texProps.textureOffset].texHeight;
                        for (unsigned int h = 0; h < height; ++h) {
                            for (unsigned int w = 0; w < width; ++w) {
                                unsigned int index_glob =
                                        w + h * model->facetTex[facet.texProps.textureOffset].texWidth;

                                shTexture[index_glob].countEquiv += texels[index_glob].countEquiv;
                                shTexture[index_glob].sum_v_ort_per_area += texels[index_glob].sum_v_ort_per_area;
                                shTexture[index_glob].sum_1_per_ort_velocity += texels[index_glob].sum_1_per_ort_velocity;

                                double val[3];  //pre-calculated autoscaling values (Pressure, imp.rate, density)
                                val[0] = shTexture[index_glob].sum_v_ort_per_area *
                                         timeCorrection; //pressure without dCoef_pressure
                                val[1] = shTexture[index_glob].countEquiv * model->texInc[index_glob + facet.texProps.textureOffset] *
                                         timeCorrection; //imp.rate without dCoef
                                val[2] = shTexture[index_glob].sum_1_per_ort_velocity * model->texInc[index_glob + facet.texProps.textureOffset] *
                                         timeCorrection; //particle density without dCoef
                                //Global autoscale
                                for (int v = 0; v < 3; v++) {
                                    if (val[v] > gHits->texture_limits[v].max.all)
                                        gHits->texture_limits[v].max.all = val[v];
                                    if (val[v] > 0.0 && val[v] < gHits->texture_limits[v].min.all)
                                        gHits->texture_limits[v].min.all = val[v];
                                }
                            }
                        }

                        break; //Only need 1 facet for texture position data
                    }
                }
            }
        }
    }

    //if there were no textures:
    for (int v = 0; v < 3; v++) {
        if (gHits->texture_limits[v].min.all == HITMAX)
            gHits->texture_limits[v].min.all = texture_limits_old[v].min.all;
        if (gHits->texture_limits[v].min.moments_only == HITMAX)
            gHits->texture_limits[v].min.moments_only = texture_limits_old[v].min.moments_only;
        if (gHits->texture_limits[v].max.all == 0.0)
            gHits->texture_limits[v].max.all = texture_limits_old[v].max.all;
        if (gHits->texture_limits[v].max.moments_only == 0.0)
            gHits->texture_limits[v].max.moments_only = texture_limits_old[v].max.moments_only;
    }

    ReleaseDataport(dpHit);

    ResetTmpCounters();
    //extern char *GetSimuStatus();
    //SetState(PROCESS_STARTING, GetSimuStatus(), false, true);

#if defined(_DEBUG)
    t1 = GetTick();
    printf("Update hits: %f us\n", (t1 - t0) * 1000000.0);
#endif
}

size_t SimulationGPU::GetHitsSize() {
    return sizeof(GlobalHitBuffer) //+ model->wp.globalHistogramParams.GetDataSize()
           //+ textTotalSize + profTotalSize + dirTotalSize + angleMapTotalSize + histogramTotalSize
           + model->nbFacets_total * sizeof(FacetHitBuffer);
}

void SimulationGPU::ResetTmpCounters() {
    //SetState(0, "Resetting local cache...", false, true);

    memset(&tmpGlobalResult, 0, sizeof(GlobalHitBuffer));
    GlobalCounter* globalCount = gpuSim.GetGlobalCounter();
    std::fill(globalCount->facetHitCounters.begin(),globalCount->facetHitCounters.end(), flowgpu::CuFacetHitCounter64());
    std::fill(globalCount->leakCounter.begin(),globalCount->leakCounter.end(), 0);
    for(auto& tex : globalCount->textures){
        std::fill(tex.second.begin(),tex.second.end(), Texel64());
    }
    for(auto& prof : globalCount->profiles){
        std::fill(prof.second.begin(),prof.second.end(), Texel64());
    }
}

static uint64_t currentDes = 0;
bool SimulationGPU::SimulationMCStep(size_t nbStep){
    for(int i=0;i<nbStep;++i)
        currentDes = gpuSim.RunSimulation();
    //currentDes += nbStep * LAUNCHSIZE;
    bool goOn = this->model->ontheflyParams.desorptionLimit > currentDes;
#ifdef WITHDESORPEXIT
    goOn = true;
#endif
    if(this->model->ontheflyParams.desorptionLimit == 0)
        goOn = true;
    return goOn;
}

int SimulationGPU::ReinitializeParticleLog() {
    // GPU simulations has no support for particle log
    return 0;
}