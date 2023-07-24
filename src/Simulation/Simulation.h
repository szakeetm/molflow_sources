/*
Program:     MolFlow+ / Synrad+
Description: Monte Carlo simulator for ultra-high vacuum and synchrotron radiation
Authors:     Jean-Luc PONS / Roberto KERSEVAN / Marton ADY / Pascal BAEHR
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

#include "Buffer_shared.h" //Facetproperties
#include "SMP.h"
#include "Vector.h"
#include "Random.h"
#include "ProcessControl.h"
#include "SimulationController.h"
#include "MolflowSimGeom.h"
#include "Particle.h"
#include "RayTracing/RTHelper.h"
#include "Simulation.h"

// Local simulation structure
class Simulation /*: public Simulation_Abstract*/ {
public:

    Simulation();
    //Simulation(Simulation&& o) noexcept ;
    virtual ~Simulation() = default;

    std::pair<int, std::optional<std::string>> SanityCheckModel(bool strictCheck) ;
    void ClearSimulation() ;
    size_t LoadSimulation(std::string& loadStatus) ;
    int RebuildAccelStructure() ;

    void ResetSimulation() ;

    size_t GetHitsSize() ;

    int ReinitializeParticleLog() ;
    MFSim::ParticleTracer * GetParticleTracerPtr(size_t i) ;
    void ConstructParticleTracers(size_t n, bool fixedSeed) ;
	bool lastLogUpdateOK; // Last log update timeout
	bool hasVolatile;   // Contains volatile facet
    std::vector<MFSim::ParticleTracer> particleTracers; //they have their own tmp counters
    mutable std::timed_mutex simuStateMutex;

    std::shared_ptr<SimulationModel> model;
    GlobalSimuState* globStatePtr;
    ParticleLog* globParticleLogPtr; //Recorded particle log since last UpdateMCHits

    size_t totalDesorbed; // todo: should be a "sim counter"

};