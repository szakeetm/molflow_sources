

#pragma once

//#include "MolflowSimGeom.h"
//#include "SimulationUnit.h"
#include <Random.h>
#include <mutex>
#include "RayTracing/Ray.h"

struct FacetHitDetail;
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

        void PerformTeleport(SimulationFacet *teleportSourceFacet);

        void RegisterTransparentPass(SimulationFacet *facet);

        void RecordAbsorb(SimulationFacet *absorbFacet);

        void PerformBounce(SimulationFacet *bounceFacet);

        bool PerformScatter();

        void RecordHistograms(SimulationFacet *histogramFacet, int m);

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
        double expectedScatterPath; //for background collisions
        //size_t structureId;        // Current structure
        std::unique_ptr<GlobalSimuState> tmpState=std::make_unique<GlobalSimuState>(); //Thread-local "unadded" results, that are reset to 0 when added to global state. Pointer to break circular includes
        std::unique_ptr<ParticleLog> tmpParticleLog=std::make_unique<ParticleLog>(); //Pointer to break circular includes
        int lastHitFacetId = -1; //id of last hit facet, that ray tracer should avoid. Remembered to persist across multiple SimulationNBStep(n) calls
        bool insertNewParticleAtNextStep = true; //Remembered to persist across multiple SimulationNBStep(n) calls
        MersenneTwister randomGenerator;
        std::shared_ptr<MolflowSimulationModel> model;
        std::vector<SimulationFacet*> transparentHitBuffer; //Storing this buffer simulation-wide is cheaper than recreating it at every Intersect() call
        std::vector <FacetHitDetail> facetHitDetails; //One per SimulationFacet, hit details for intersect routine

        bool exitRequested{false};

        Vector3d nullVector; //so we don't have to allocate and destroy for dummy uses

        private:
        int LookupMomentIndex(const double time, const size_t startIndex);

    };
}
