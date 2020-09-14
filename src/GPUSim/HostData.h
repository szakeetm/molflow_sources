//
// Created by pbahr on 06/02/2020.
//

#ifndef MOLFLOW_PROJ_HOSTDATA_H
#define MOLFLOW_PROJ_HOSTDATA_H

#include <cstdint>
#include <vector>
#include <map>
#include "LaunchParams.h" //CuFacetHitCounter



struct Texel64 {
    Texel64() : countEquiv(0), sum_v_ort_per_area(0.0), sum_1_per_ort_velocity(0.0){}

    unsigned long long int countEquiv;
    double sum_v_ort_per_area;
    double sum_1_per_ort_velocity;
};

struct HostData {
    //std::vector<uint32_t> pixels;
#ifdef WITHDESORPEXIT
    std::vector<flowgpu::MolPRD> hitData;
#endif
    std::vector<flowgpu::CuFacetHitCounter> facetHitCounters;
    std::vector<flowgpu::Texel> texels;
    std::vector<flowgpu::Texel> profileSlices;

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
#ifdef DEBUGLEAKPOS
    std::vector<float3> leakPositions;
    std::vector<float3> leakDirections;
    std::vector<uint32_t> leakPosOffset;
#endif
};

struct GlobalCounter {
    //std::vector<uint32_t> pixels;
    //std::vector<MolPRD> hitData;
    std::vector<flowgpu::CuFacetHitCounter64> facetHitCounters;
    std::vector<uint64_t> leakCounter;
    std::map<uint32_t,std::vector<Texel64>> textures;
    std::map<uint32_t,std::vector<Texel64>> profiles;

#ifdef DEBUGCOUNT
    std::vector<uint64_t> detCounter;
    std::vector<uint64_t> uCounter;
    std::vector<uint64_t> vCounter;
#endif
#ifdef DEBUGPOS
    std::vector<float3> positions;
    std::vector<uint32_t> posOffset;
#endif
#ifdef DEBUGLEAKPOS
    std::vector<float3> leakPositions;
    std::vector<float3> leakDirections;
    std::vector<uint32_t> leakPosOffset;
#endif
};

#endif //MOLFLOW_PROJ_HOSTDATA_H
