//
// Created by pbahr on 18/05/2020.
//

#include <cereal/archives/binary.hpp>
#include "SimControllerGPU.h"
#include "ModelReader.h" // TempFacet

bool SimControllerGPU::Load() {
    return false;
}

SimControllerGPU::SimControllerGPU(std::string appName, std::string dpName, size_t parentPID, size_t procIdx)
        : SimulationController(appName, dpName, parentPID, procIdx) {

    totalDesorbed = 0;
    loadOK = false;

    memset(&tmpGlobalResult, 0, sizeof(GlobalHitBuffer));
}

SimControllerGPU::~SimControllerGPU() {
    delete model;
}

int SimControllerGPU::SanityCheckGeom() {
    return 0;
}

void SimControllerGPU::ClearSimulation() {

}

bool SimControllerGPU::LoadSimulation(Dataport *loader) {
    double t0 = GetTick();

    SetState(PROCESS_STARTING, "Clearing previous simulation");
    ClearSimulation();

    SetState(PROCESS_STARTING, "Loading simulation");

    {

        std::string inputString(loader->size,'\0');
        BYTE* buffer = (BYTE*)loader->buff;
        std::copy(buffer, buffer + loader->size, inputString.begin());

        flowgeom::loadFromSerialization(inputString);
    }//inputarchive goes out of scope, file released

    // Initialise simulation

    //if(!sh.name.empty())
    loadOK = true;
    double t1 = GetTick();
    printf("  Load %s successful\n", model->geomProperties.name.c_str());
    printf("  Loading time: %.3f ms\n", (t1 - t0)*1000.0);
    return true;
}

void SimControllerGPU::ResetSimulation() {
    totalDesorbed = 0;
    ResetTmpCounters();
}

bool SimControllerGPU::UpdateOntheflySimuParams(Dataport *loader) {
    // Connect the dataport


    if (!AccessDataportTimed(loader, 2000)) {
        SetErrorSub("Failed to connect to loader DP");
        return false;
    }
    std::string inputString(loader->size,'\0');
    BYTE* buffer = (BYTE*)loader->buff;
    std::copy(buffer, buffer + loader->size, inputString.begin());
    std::stringstream inputStream;
    inputStream << inputString;
    cereal::BinaryInputArchive inputArchive(inputStream);

    inputArchive(ontheflyParams);

    ReleaseDataport(loader);

    return true;
}

void SimControllerGPU::UpdateHits(Dataport *dpHit, Dataport* dpLog,int prIdx, DWORD timeout) {
    //UpdateMCHits(dpHit, prIdx, moments.size(), timeout);
    //if (dpLog) UpdateLog(dpLog, timeout);
}

size_t SimControllerGPU::GetHitsSize() {
    return sizeof(GlobalHitBuffer) //+ model->wp.globalHistogramParams.GetDataSize()
           //+ textTotalSize + profTotalSize + dirTotalSize + angleMapTotalSize + histogramTotalSize
           + model->nbFacets_total * sizeof(FacetHitBuffer);
}

void SimControllerGPU::ResetTmpCounters() {
    SetState(0, "Resetting local cache...", false, true);

    memset(&tmpGlobalResult, 0, sizeof(GlobalHitBuffer));
}