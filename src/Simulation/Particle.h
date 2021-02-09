//
// Created by pascal on 2/5/21.
//

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

        void RecordHitOnTexture(const SubprocessFacet *f, double time, bool countHit, double velocity_factor,
                                double ortSpeedFactor);

        void ProfileFacet(const SubprocessFacet *f, double time, bool countHit, double velocity_factor,
                          double ortSpeedFactor);

        void RecordHit(const int &type);

        void RecordLeakPos();

        void IncreaseFacetCounter(const SubprocessFacet *f, double time, size_t hit, size_t desorb, size_t absorb,
                                  double sum_1_per_v, double sum_v_ort);

        void UpdateVelocity(const SubprocessFacet *collidedFacet);

        void LogHit(SubprocessFacet *f, std::vector <ParticleLoggerItem> &tmpParticleLog);

        void RecordDirectionVector(const SubprocessFacet *f, double time);

        void RecordAngleMap(const SubprocessFacet *collidedFacet);

        void PerformTeleport(SubprocessFacet *iFacet);

        void RegisterTransparentPass(SubprocessFacet *facet);

        void RecordAbsorb(SubprocessFacet *iFacet);

        void PerformBounce(SubprocessFacet *iFacet);

        void RecordHistograms(SubprocessFacet *iFacet);

        bool UpdateHits(GlobalSimuState *globState, size_t timeout);

        void Reset();

        Vector3d position;    // Position
        Vector3d direction;    // Direction
        double oriRatio; //Represented ratio of desorbed, used for low flux mode

        //Recordings for histogram
        uint64_t totalDesorbed;
        size_t nbBounces; // Number of hit (current particle) since desorption
        size_t lastMomentIndex; // Speedup for binary search
        size_t particleId;
        double distanceTraveled;
        double generationTime; //Time it was created, constant
        double particleTime; //Actual time, incremented after every hit. (Flight time = actual time - generation time)
        int teleportedFrom;   // We memorize where the particle came from: we can teleport back

        double velocity;
        double expectedDecayMoment; //for radioactive gases
        size_t structureId;        // Current structure
        GlobalSimuState tmpState;
        SubprocessFacet *lastHitFacet;     // Last hitted facet
        MersenneTwister randomGenerator;
        SimulationModel *model;
        std::vector<SubprocessFacet *> transparentHitBuffer; //Storing this buffer simulation-wide is cheaper than recreating it at every Intersect() call
        std::vector <SubProcessFacetTempVar> tmpFacetVars; //One per subprocessfacet, for intersect routine
    };
}

#endif //MOLFLOW_PROJ_PARTICLE_H
