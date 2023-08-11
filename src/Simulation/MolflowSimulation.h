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
#include <string>
#include <mutex>

#include "SimulationUnit.h"
#include "Buffer_shared.h" //Facetproperties
#include "SMP.h"
#include "Vector.h"
#include "Random.h"
#include "ProcessControl.h"
#include "SimulationController.h"
#include "MolflowSimGeom.h"
#include "RayTracing/RTHelper.h"

class ParticleTracer;

// Local simulation structure
class MolflowSimulation : public Simulation_Abstract {
public:

    MolflowSimulation()=default;
    //MolflowSimulation(MolflowSimulation&& o) noexcept;
    virtual ~MolflowSimulation() = default;

    std::vector<std::string> SanityCheckModel(bool strictCheck) override;
    //void ClearSimulation() override;
    int LoadSimulation(ProcCommData& procInfo, LoadStatus_abstract* loadStatus=nullptr) override;
    int RebuildAccelStructure() override;

    void ResetSimulation() override;

    int GetHitsSize() override;

    int ReinitializeParticleLog() override;
    std::shared_ptr<MFSim::ParticleTracer> GetParticleTracerPtr(int i) override;
    void ConstructParticleTracers(int n, bool fixedSeed) override;
    bool lastLogUpdateOK=true; // Last log update timeout
    bool hasVolatile=false;   // Contains volatile facet
    std::vector<std::shared_ptr<MFSim::ParticleTracer>> particleTracers; //they have their own tmp counters
    mutable std::timed_mutex simuStateMutex;

};
// -- Methods ---------------------------------------------------