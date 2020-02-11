//
// Created by pbahr on 06/02/2020.
//

#ifndef MOLFLOW_PROJ_HOSTDATA_H
#define MOLFLOW_PROJ_HOSTDATA_H

#include <cstdint>
#include <vector>
#include "LaunchParams.h" //CuFacetHitCounter

struct HostData {
    //std::vector<uint32_t> pixels;
    //std::vector<MolPRD> hitData;
    std::vector<flowgpu::CuFacetHitCounter> facetHitCounters;
    std::vector<uint32_t> leakCounter;

#ifdef DEBUGCOUNT
    std::vector<uint32_t> detCounter;
    std::vector<uint32_t> uCounter;
    std::vector<uint32_t> vCounter;
#endif
#ifdef DEBUGPOS
    std::vector<float3> positions;
    std::vector<uint32_t> posOffset;
#endif
};

#endif //MOLFLOW_PROJ_HOSTDATA_H
