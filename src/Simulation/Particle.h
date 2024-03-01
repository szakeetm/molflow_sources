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

//#include "MolflowSimGeom.h"
//#include "SimulationUnit.h"
#include <Random.h>
#include <mutex>
#include "RayTracing/Ray.h"

struct SimulationFacetTempVar;
class MolflowSimulationModel;
class GlobalSimuState;
struct ParticleLog;
class SimulationFacet;

enum ThreadState:int;

/**
* \brief Namespace containing various simulation only classes and methods
 */
namespace MFSim {

    enum MCStepResult {
        Success,
        DesorptionError,
        MaxReached
    };

    enum ParticleEventType {
        ParticleEvent_FacetHit,
        ParticleEvent_Overtime,
        ParticleEvent_Decay,
        ParticleEvent_Scatter
    };

/**
* \brief Implements particle state and corresponding pre-/post-processing methods (source position, hit recording etc.)
 */
    class ParticleTracer {
    public:

        void IncreaseDistanceCounters(double distanceIncrement);

        MCStepResult SimulationMCStep(size_t nbStep, size_t threadNum, size_t remainingDes);

        bool StartFromSource(Ray& ray);

        bool UpdateMCHits(const std::shared_ptr<GlobalSimuState> globalState, size_t nbMoments,
            std::string& myStatus, std::mutex& statusMutex, size_t timeout_ms);

        void RecordHitOnTexture(const SimulationFacet *f, int m, bool countHit, double velocity_factor,
                                double ortSpeedFactor);

        void ProfileFacet(const SimulationFacet *f, int m, bool countHit, double velocity_factor,
                          double ortSpeedFactor);

        void RecordHit(const int type);

        void RecordLeakPos();

		void IncreaseFacetCounter(const SimulationFacet* f, int m, const size_t hit, const size_t desorb, const size_t absorb,
			const double sum_1_per_v, const double sum_v_ort,
			const Vector3d& impulse, const Vector3d& impulse_square, const Vector3d& impulse_momentum);

        void UpdateVelocity(const SimulationFacet *collidedFacet);

        void LogHit(SimulationFacet *f);

        void RecordDirectionVector(const SimulationFacet *f, int m);

        void RecordAngleMap(const SimulationFacet *collidedFacet);

        void PerformTeleport(SimulationFacet *iFacet);

        void RegisterTransparentPass(SimulationFacet *facet);

        void RecordAbsorb(SimulationFacet *iFacet);

        void PerformBounce(SimulationFacet *iFacet);

        bool PerformScatter();

        void RecordHistograms(SimulationFacet *iFacet, int m);

        bool UpdateHitsAndLog(const std::shared_ptr<GlobalSimuState> globalState, const std::shared_ptr<ParticleLog> particleLog,
            ThreadState& myState, std::string& myStatus, std::mutex& statusMutex, size_t timeout_ms);
        bool UpdateLog(const std::shared_ptr<ParticleLog> globalLog,
            std::string& myStatus, std::mutex& statusMutex, size_t timeout);

        void Reset();

        size_t GetMemSize() const;

        Ray ray; // an object purely for the ray tracing related intersection tests
        double oriRatio; //Represented ratio of desorbed, used for low flux mode

        //Recordings for histogram
        uint64_t totalDesorbed;
        size_t nbBounces; // Number of hit (current particle) since desorption
        size_t lastMomentIndex; // Speedup for binary search
        size_t particleTracerId; //For parallel computing, each core gets a particle id. Lines and leak recording only for id=0
        double distanceTraveled;
        double generationTime; //Time it was created, constant
        //double particleTime; //Actual time, incremented after every hit. (Flight time = actual time - generation time)
        int teleportedFrom;   // We memorize where the particle came from: we can teleport back

        double velocity;
        double initialVelocity; //Used to check for reaching Brownian motion if background collisions are enabled
        double expectedDecayMoment; //for radioactive gases
        //double expectedFreePath; //for background collisions
        //size_t structureId;        // Current structure
        std::unique_ptr<GlobalSimuState> tmpState=std::make_unique<GlobalSimuState>(); //Thread-local "unadded" results, that are reset to 0 when added to global state. Pointer to break circular includes
        std::unique_ptr<ParticleLog> tmpParticleLog=std::make_unique<ParticleLog>(); //Pointer to break circular includes
        SimulationFacet* lastHitFacet=nullptr;     // Last hitted facet, nullptr by default
        MersenneTwister randomGenerator;
        std::shared_ptr<MolflowSimulationModel> model;
        std::vector<SimulationFacet*> transparentHitBuffer; //Storing this buffer simulation-wide is cheaper than recreating it at every Intersect() call
        std::vector <SimulationFacetTempVar> tmpFacetVars; //One per SimulationFacet, for intersect routine

        bool exitRequested{false};

        Vector3d nullVector; //so we don't have to allocate and destroy for dummy uses

        private:
        int LookupMomentIndex(const double time, const size_t startIndex);

    };
}
