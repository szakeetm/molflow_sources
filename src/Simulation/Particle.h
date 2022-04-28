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

#ifndef MOLFLOW_PROJ_PARTICLE_H
#define MOLFLOW_PROJ_PARTICLE_H


#include "../GeometrySimu.h"
#include "SimulationUnit.h"
#include <Random.h>

struct SubProcessFacetTempVar;

namespace MFSim {
    class Particle {
    public:
        //double GenerateRandomVelocity(int CDFId, const double rndVal);

        //double GenerateDesorptionTime(const SubprocessFacet *src, const double rndVal);

        void IncreaseDistanceCounters(double distanceIncrement);

        bool SimulationMCStep(size_t nbStep, size_t threadNum, size_t remainingDes);

        bool StartFromSource();

        bool UpdateMCHits(GlobalSimuState &globSimuState, size_t nbMoments, DWORD timeout);

        void RecordHitOnTexture(const SubprocessFacet *f, int m, bool countHit, double velocity_factor,
                                double ortSpeedFactor);

        void ProfileFacet(const SubprocessFacet *f, int m, bool countHit, double velocity_factor,
                          double ortSpeedFactor);

        void RecordHit(const int &type);

        void RecordLeakPos();

		void IncreaseFacetCounter(const SubprocessFacet* f, int m, const size_t& hit, const size_t& desorb, const size_t& absorb,
			const double& sum_1_per_v, const double& sum_v_ort,
			const Vector3d& impulse = Vector3d(0.0, 0.0, 0.0),
			const Vector3d& impulse_square = Vector3d(0.0, 0.0, 0.0),
			const Vector3d& impulse_momentum = Vector3d(0.0, 0.0, 0.0));

        void UpdateVelocity(const SubprocessFacet *collidedFacet);

        void LogHit(SubprocessFacet *f);

        void RecordDirectionVector(const SubprocessFacet *f, int m);

        void RecordAngleMap(const SubprocessFacet *collidedFacet);

        void PerformTeleport(SubprocessFacet *iFacet);

        void RegisterTransparentPass(SubprocessFacet *facet);

        void RecordAbsorb(SubprocessFacet *iFacet);

        void PerformBounce(SubprocessFacet *iFacet);

        void RecordHistograms(SubprocessFacet *iFacet, int m);

        bool UpdateHits(GlobalSimuState *globState, ParticleLog *particleLog, size_t timeout);
        bool UpdateLog(ParticleLog *globalLog, size_t timeout);

        void Reset();

        Ray particle;
        //Vector3d position;    // Position
        //Vector3d direction;    // Direction
        double oriRatio; //Represented ratio of desorbed, used for low flux mode

        //Recordings for histogram
        uint64_t totalDesorbed;
        size_t nbBounces; // Number of hit (current particle) since desorption
        size_t lastMomentIndex; // Speedup for binary search
        size_t particleId;
        double distanceTraveled;
        double generationTime; //Time it was created, constant
        //double particleTime; //Actual time, incremented after every hit. (Flight time = actual time - generation time)
        int teleportedFrom;   // We memorize where the particle came from: we can teleport back

        double velocity;
        double expectedDecayMoment; //for radioactive gases
        //size_t structureId;        // Current structure
        GlobalSimuState tmpState;
        ParticleLog tmpParticleLog;
        SubprocessFacet *lastHitFacet;     // Last hitted facet
        MersenneTwister randomGenerator;
        SimulationModel *model;
        std::vector<SubprocessFacet*> transparentHitBuffer; //Storing this buffer simulation-wide is cheaper than recreating it at every Intersect() call
        std::vector <SubProcessFacetTempVar> tmpFacetVars; //One per subprocessfacet, for intersect routine

        bool allQuit{false};

        bool StartFromSource(Ray &ray);
    };
}

#endif //MOLFLOW_PROJ_PARTICLE_H
