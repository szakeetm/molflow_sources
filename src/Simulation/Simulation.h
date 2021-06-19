/*
Program:     MolFlow+ / Synrad+
Description: Monte Carlo simulator for ultra-high vacuum and synchrotron radiation
Authors:     Jean-Luc PONS / Roberto KERSEVAN / Marton ADY
Copyright:   E.S.R.F / CERN
Website:     https://cern.ch/molflow

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

Full license text: https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html
*/
#pragma once

#include <tuple>
#include <vector>
#include <mutex>

#include "SimulationUnit.h"
#include "Buffer_shared.h" //Facetproperties
#include "SMP.h"
#include "Vector.h"
#include "Random.h"
#include "ProcessControl.h"
#include "SimulationController.h"
#include <../src/GeometrySimu.h>
#include "Particle.h"
#include "RayTracing/RTHelper.h"

class Parameter;

// Local simulation structure
class Simulation : public SimulationUnit {
public:

    Simulation();
    Simulation(Simulation&& o) noexcept ;
    virtual ~Simulation() = default;

    std::pair<int, std::optional<std::string>> SanityCheckModel(bool strictCheck) override;
    void ClearSimulation() override;
    size_t LoadSimulation(char *loadStatus) override;
    int RebuildAccelStructure() override;

    void ResetSimulation() override;

    size_t GetHitsSize() override;

    int ReinitializeParticleLog() override;
    MFSim::Particle * GetParticle(size_t i) override {
        if(i < particles.size())
            return &particles.at(i);
        else
            return nullptr;
    };
    void SetNParticle(size_t n, bool fixedSeed) override {
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
    };

	//size_t totalDesorbed;           // Total desorption number (for this process, not reset on UpdateMCHits)

	// Geometry

	//double stepPerSec=0.0;  // Avg number of step per sec
	//bool loadOK;        // Load OK flag
	//bool lastHitUpdateOK;  // Last hit update timeout
	bool lastLogUpdateOK; // Last log update timeout
	bool hasVolatile;   // Contains volatile facet

    // Particle coordinates (MC)//std::vector<FacetHistogramBuffer> tmpGlobalHistograms; //Recorded histogram since last UpdateMCHits, 1+nbMoment copies
    //ParticleLog tmpParticleLog; //Recorded particle log since last UpdateMCHits
    std::vector<MFSim::Particle> particles;
    mutable std::timed_mutex tMutex;

};
// -- Methods ---------------------------------------------------