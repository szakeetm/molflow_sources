

#pragma once

#include <vector>
#include "../MolflowTypes.h"
#include "../Parameter.h"
#include "Buffer_shared.h"
#include "Helper/ConsoleLogger.h"

#include <mutex>
#include "FacetData.h"
#include "SimulationModel.h"
#include "RayTracing/BVH.h"

#include <cereal/cereal.hpp>
#include <cereal/types/vector.hpp>
#include "RayTracing/KDTree.h"
#include <map>


class SimulationFacet;
class MolflowSimFacet;

class ParameterSurface : public Surface {
    Parameter* timeDepParam; //time-dependent opacity distribution
public:
    ParameterSurface(Parameter *param) : timeDepParam(param) {};

    ~ParameterSurface() = default;

    bool IsHardHit(const Ray &r) override;
};

//Data structure specific to MolflowSimulationModel to store a collection of time-dependent parameters, similar to "Formulas"
//instance: model->tdParams
struct TimeDependentParameters {
    TimeDependentParameters() = default;

    std::vector<Parameter> parameters;
    std::vector<std::vector<IntegratedDesorptionEntry>> IDs; //integrated distribution function for each time-dependent desorption type
    std::vector<Moment> moments;

    static void LoadParameterCatalog(std::vector<Parameter>& parameters);
    static void ClearParameters(std::vector<Parameter>& parameters);
    static std::vector<Parameter> GetCatalogParameters(const std::vector<Parameter>& parameters);
    static size_t InsertParametersBeforeCatalog(std::vector<Parameter>& targetParameters, const std::vector<Parameter>& newParams);
    int GetParamId(const std::string& name) const;
    
    size_t GetMemSize();
};

/*!
 * @brief One instance is the state for one facet for a single moment
 */
class FacetMomentSnapshot {
public:
    FacetMomentSnapshot();
    FacetMomentSnapshot& operator+=(const FacetMomentSnapshot& rhs);

    size_t GetMemSize() const;

    FacetHitBuffer hits;
    std::vector<ProfileSlice> profile;
    std::vector<TextureCell> texture;
    std::vector<DirectionCell> direction;
    FacetHistogramBuffer histogram;

    template<class Archive>
    void serialize(Archive& archive) {
        archive(
            hits,
            profile,
            texture,
            direction,
            histogram
        );
    }
};
FacetMomentSnapshot operator+(const FacetMomentSnapshot& lhs, const FacetMomentSnapshot& rhs);

/*!
 * @brief Object containing all simulation results of an individual facet
 */
class FacetState {
public:
    FacetState& operator+=(const FacetState& rhs);

    size_t GetMemSize() const;

    std::vector<size_t> recordedAngleMapPdf; //Not time-dependent
    std::vector<FacetMomentSnapshot> momentResults; //1+nbMoment
    template<class Archive>
    void serialize(Archive& archive) {
        archive(
            recordedAngleMapPdf,
            momentResults
        );
    }
};
FacetState operator+(const FacetState& lhs, const FacetState& rhs);

/*!
 * @brief Object containing all simulation results, global and per facet
 */
class GlobalSimuState {
public:
    GlobalSimuState& operator=(const GlobalSimuState& src);
    GlobalSimuState& operator+=(const GlobalSimuState& rhs);

    GlobalSimuState() = default; //so a vector can be made of this for particletracer local results
    GlobalSimuState(GlobalSimuState&& rhs) noexcept {
        globalHistograms = std::move(rhs.globalHistograms);
        facetStates = std::move(rhs.facetStates);
        globalStats = rhs.globalStats;
        initialized = rhs.initialized;
    };

    GlobalSimuState(const std::shared_ptr<GlobalSimuState> rhs) {
        *this = rhs; //Uses overloaded operator=
    };

    //size_t memSizeCache=0;
    bool initialized = false;

    void Clear();

    void Resize(std::shared_ptr<SimulationModel> model);

    void Reset();

    size_t GetMemSize() const;

    static std::tuple<size_t, size_t, size_t>
        Compare(const std::shared_ptr<GlobalSimuState> lhsGlobHit, const std::shared_ptr<GlobalSimuState> rhsGlobHit, double globThreshold,
            double locThreshold);

#if defined(MOLFLOW)
    GlobalHitBuffer globalStats;
    std::vector<FacetHistogramBuffer> globalHistograms; //1+nbMoment
    std::vector<FacetState> facetStates; //nbFacet
#endif

    template<class Archive>
    void serialize(Archive& archive) {
        archive(
            CEREAL_NVP(globalStats),
            CEREAL_NVP(globalHistograms),
            CEREAL_NVP(facetStates)
        );
    }

    mutable std::timed_mutex simuStateMutex;
};

// Local simulation structure for Molflow specific simulations
// Appends general SimulationModel by time dependent parameters
class MolflowSimulationModel : public SimulationModel {
  
public:
    /*
    MolflowSimulationModel() : SimulationModel() {};

    ~MolflowSimulationModel();

    MolflowSimulationModel(MolflowSimulationModel &&o) noexcept;

    MolflowSimulationModel(const MolflowSimulationModel &o);
    */
    size_t GetMemSize() override;
    /*
    MolflowSimulationModel &operator=(const MolflowSimulationModel &o) {
        facets = o.facets;
        structures = o.structures;

        rayTracingStructures.insert(rayTracingStructures.begin(), o.rayTracingStructures.begin(), o.rayTracingStructures.end());
        vertices3 = o.vertices3;
        otfParams = o.otfParams;
        tdParams = o.tdParams;
        sp = o.sp;
        sh = o.sh;
        initialized = o.initialized;

        return *this;
    };

    MolflowSimulationModel &operator=(MolflowSimulationModel &&o) noexcept {
        facets = std::move(o.facets);
        structures = std::move(o.structures);

        rayTracingStructures = std::move(o.rayTracingStructures);
        vertices3 = std::move(o.vertices3);
        tdParams = std::move(o.tdParams);
        otfParams = o.otfParams;
        sp = o.sp;
        sh = o.sh;
        initialized = o.initialized;

        return *this;
    };
    */

    //! Do calculations necessary before launching simulation
    void PrepareToRun() override;

    //! Construct acceleration structure with a given splitting method
    int BuildAccelStructure(const std::shared_ptr<GlobalSimuState> globalState, AccelType accel_type, BVHAccel::SplitMethod split,
                            int maxPrimsInNode) override;

    //int InitializeFacets();

    void CalcTotalOutgassing();
    std::vector<std::string> SanityCheck();
    /**
    * \brief Returns an existing or a new surface corresponding to a facet's properties
    * \param facet facet for which a Surface should be found or created
     * \return new or existing Surface corresponding to the chosen parameters of the facet
    */
    std::shared_ptr<Surface> GetSurface(std::shared_ptr<SimulationFacet> facet) override;

    // Sim functions
    double GetOpacityAt(const MolflowSimFacet *f, const double time) const;
    double GetStickingAt(const MolflowSimFacet *f, const double time) const;
    double GetTemperatureAt(const MolflowSimFacet *f, const double time) const;

    TimeDependentParameters tdParams;
    std::vector<Interval> intervalCache; //speedup to store moments as [start_time,end_time], calculated in PrepareToRun();
    std::vector<IntegratedVelocityEntry> maxwell_CDF_1K; //Integrated "surface" maxwell-boltzmann distribution at 1K. TODO: Make global
    void CalcIntervalCache();

    //void BuildPrisma(double L, double R, double angle, double s, int step);
};




[[nodiscard]] std::optional<std::unique_lock<std::timed_mutex>> GetHitLock(GlobalSimuState* simStatePtr, size_t waitMs);

/*!
 * @brief Particle Log structure containing all individual log entries
 */
struct ParticleLog {
public:
    //Methods below are required because it's part of ParticleTracer which is stored in a vector
    ParticleLog &operator=(const ParticleLog &src) {
        pLog = src.pLog;
        return *this;
    }

    ParticleLog(ParticleLog &&rhs) noexcept { //move constructor
        pLog = std::move(pLog);
    };

    ParticleLog(const ParticleLog &rhs) {
        pLog = rhs.pLog;
    };

    ParticleLog() = default;

    void resize(size_t nbLogs) {
        std::vector<ParticleLoggerItem>(nbLogs).swap(pLog);
    };

    void clear() {
        pLog.clear();
    };

    std::vector<ParticleLoggerItem> pLog;
    mutable std::timed_mutex particleLogMutex;
};
