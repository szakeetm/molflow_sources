//
// Created by pbahr on 06/02/2020.
//

#ifndef MOLFLOW_PROJ_HOSTDATA_H
#define MOLFLOW_PROJ_HOSTDATA_H

#include <cstdint>
#include <vector>
#include <map>
#include "LaunchParams.h" //CuFacetHitCounter

struct CuFacetHitCounter64{
    CuFacetHitCounter64() :
            nbDesorbed(0) , nbMCHit(0) , nbHitEquiv(0.0) ,
            nbAbsEquiv(0.0) , sum_1_per_ort_velocity(0.0) ,
            sum_1_per_velocity(0.0) , sum_v_ort(0.0)
    {};
    // Counts
    uint64_t nbMCHit;                   // Number of hits
    uint64_t nbDesorbed;                // Number of desorbed molec
    double nbAbsEquiv;                   // Equivalent number of absorbed molecules
    double nbHitEquiv;                   //Equivalent number of hits, used for low-flux impingement rate and density calculation
    double sum_1_per_ort_velocity;       // sum of reciprocials of orthogonal velocity components, used to determine the density, regardless of facet orientation
    double sum_1_per_velocity;           //For average molecule speed calculation
    double sum_v_ort;                    // sum of orthogonal speeds of incident velocities, used to determine the pressure
};

struct Texel64 {
    Texel64() : countEquiv(0), sum_v_ort_per_area(0.0), sum_1_per_ort_velocity(0.0){}

    unsigned long long int countEquiv;
    double sum_v_ort_per_area;
    double sum_1_per_ort_velocity;
};

struct HostData {
    //std::vector<uint32_t> pixels;
    //std::vector<MolPRD> hitData;
    std::vector<flowgpu::CuFacetHitCounter> facetHitCounters;
    std::vector<flowgeom::Texel> texels;
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

struct GlobalCounter {
    //std::vector<uint32_t> pixels;
    //std::vector<MolPRD> hitData;
    std::vector<CuFacetHitCounter64> facetHitCounters;
    std::vector<uint64_t> leakCounter;
    std::map<uint32_t,std::vector<Texel64>> textures;

#ifdef DEBUGCOUNT
    std::vector<uint64_t> detCounter;
    std::vector<uint64_t> uCounter;
    std::vector<uint64_t> vCounter;
#endif
#ifdef DEBUGPOS
    std::vector<float3> positions;
    std::vector<uint32_t> posOffset;
#endif
};

#endif //MOLFLOW_PROJ_HOSTDATA_H
