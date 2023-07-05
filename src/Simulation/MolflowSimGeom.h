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

#ifndef MOLFLOW_PROJ_MOLFLOWSIMGEOM_H
#define MOLFLOW_PROJ_MOLFLOWSIMGEOM_H

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
class GlobalSimuState;

class ParameterSurface : public Surface {
    Distribution2D *dist;
public:
    ParameterSurface(Distribution2D *distribution) : dist(distribution) {};

    ~ParameterSurface() = default;

    bool IsHardHit(const Ray &r) override;
};

struct TimeDependentParameters {
    TimeDependentParameters() = default;

    std::vector<Parameter> parameters;
    std::vector<std::vector<IntegratedDesorptionEntry>> IDs; //integrated distribution function for each time-dependent desorption type
    std::vector<Moment> moments;
    
    size_t GetMemSize() {
        size_t sum = 0;
        for (auto &par : parameters) {
            sum += par.GetMemSize();
        }
        auto nbId = this->IDs.size();
        if (nbId > 0) sum += nbId * IDs[0].size() * sizeof(IntegratedDesorptionEntry);
        sum += sizeof(Moment) * moments.capacity();

        return sum;
    }
};

// Local simulation structure for Molflow specific simulations
// Appends general SimulationModel by time dependent parameters
class MolflowSimulationModel : public SimulationModel {
  
public:
    MolflowSimulationModel() : SimulationModel(), /*otfParams(),*/ tdParams()/*, wp(), sh(),*/ {};

    ~MolflowSimulationModel();

    MolflowSimulationModel(MolflowSimulationModel &&o) noexcept;

    MolflowSimulationModel(const MolflowSimulationModel &o);

    size_t size() override;

    MolflowSimulationModel &operator=(const MolflowSimulationModel &o) {
        facets = o.facets;
        structures = o.structures;

        accel.insert(accel.begin(), o.accel.begin(), o.accel.end());
        vertices3 = o.vertices3;
        otfParams = o.otfParams;
        tdParams = o.tdParams;
        wp = o.wp;
        sh = o.sh;
        initialized = o.initialized;

        return *this;
    };

    MolflowSimulationModel &operator=(MolflowSimulationModel &&o) noexcept {
        facets = std::move(o.facets);
        structures = std::move(o.structures);

        accel = std::move(o.accel);
        vertices3 = std::move(o.vertices3);
        tdParams = std::move(o.tdParams);
        otfParams = o.otfParams;
        wp = o.wp;
        sh = o.sh;
        initialized = o.initialized;

        return *this;
    };

    //! Do calculations necessary before launching simulation
    void PrepareToRun() override;

    //! Construct acceleration structure with a given splitting method
    int BuildAccelStructure(GlobalSimuState *globState, AccelType accel_type, BVHAccel::SplitMethod split,
                            int maxPrimsInNode) override;

    //int InitializeFacets();

    void CalcTotalOutgassing();

    /**
    * \brief Returns an existing or a new surface corresponding to a facet's properties
    * \param facet facet for which a Surface should be found or created
     * \return new or existing Surface corresponding to the choosen parameters of the facet
    */
    Surface *GetSurface(SimulationFacet* facet) override {

        if (facet->sh.opacity_paramId == -1){ //constant sticking
            //facet->sh.opacity = std::clamp(facet->sh.opacity, 0.0, 1.0);
            double opacity = std::clamp(facet->sh.opacity, 0.0, 1.0);

            std::shared_ptr<Surface> surface;
            if(opacity == 1.0) {
                surface = std::make_shared<Surface>();
            }
            else if(opacity == 0.0) {
                surface = std::make_shared<TransparentSurface>();
            }
            else {
                surface = std::make_shared<AlphaSurface>(opacity);
            }
            surfaces.insert(std::make_pair(opacity, surface));

            return surface.get();
        }
        else {
            auto opacity_paramId = facet->sh.opacity_paramId;
            auto* par = &this->tdParams.parameters[facet->sh.opacity_paramId];

            double indexed_id = 10.0 + opacity_paramId;
            if (!surfaces.empty()) {
                auto surf = surfaces.find(indexed_id);
                if (surf != surfaces.end())
                    return surf->second.get();
            }

            std::shared_ptr<ParameterSurface> surface;
            surface = std::make_shared<ParameterSurface>(par);
            surfaces.insert(std::make_pair(indexed_id, surface));
            Log::console_msg_master(5, "Insert surface with param id: {}\n", indexed_id);
            return surface.get();
        }
    };

    // Sim functions
    double GetOpacityAt(SimulationFacet *f, double time) const;
    double GetStickingAt(SimulationFacet *f, double time) const;

    TimeDependentParameters tdParams;
    std::vector<Interval> intervalCache; //speedup to store moments as [start_time,end_time], calculated in PrepareToRun();
    std::vector<IntegratedVelocityEntry> maxwell_CDF_1K; //Integrated "surface" maxwell-boltzmann distribution at 1K. TODO: Make global
    void CalcIntervalCache();

    void BuildPrisma(double L, double R, double angle, double s, int step);
};

/*!
 * @brief One instance is the state for one facet for a single moment
 */
class FacetMomentSnapshot {
public:
    FacetMomentSnapshot();
    FacetMomentSnapshot &operator+=(const FacetMomentSnapshot &rhs);
   

    FacetHitBuffer hits;
    std::vector<ProfileSlice> profile;
    std::vector<TextureCell> texture;
    std::vector<DirectionCell> direction;
    FacetHistogramBuffer histogram;

    template<class Archive>
    void serialize(Archive &archive) {
        archive(
                hits,
                profile,
                texture,
                direction,
                histogram
        );
    }
};
FacetMomentSnapshot operator+(const FacetMomentSnapshot &lhs,const FacetMomentSnapshot &rhs);

/*!
 * @brief Object containing all simulation results of an individual facet
 */
class FacetState {
public:
    FacetState &operator+=(const FacetState &rhs);
    

    std::vector<size_t> recordedAngleMapPdf; //Not time-dependent
    std::vector<FacetMomentSnapshot> momentResults; //1+nbMoment
    template<class Archive>
    void serialize(Archive &archive) {
        archive(
                recordedAngleMapPdf,
                momentResults
        );
    }
};
FacetState operator+(const FacetState& lhs,const FacetState& rhs);

/*!
 * @brief Object containing all simulation results, global and per facet
 */
class GlobalSimuState { //replaces old hits dataport
public:
    GlobalSimuState &operator=(const GlobalSimuState &src);
    GlobalSimuState &operator+=(const GlobalSimuState &src);

    GlobalSimuState(GlobalSimuState &&rhs) noexcept: tMutex() {
        globalHistograms = std::move(rhs.globalHistograms);
        facetStates = std::move(rhs.facetStates);
        globalStats = rhs.globalStats;
        initialized = rhs.initialized;
    };

    GlobalSimuState(const GlobalSimuState &rhs) {
        globalStats = rhs.globalStats;
        globalHistograms = rhs.globalHistograms;
        facetStates = rhs.facetStates;
        initialized = rhs.initialized;
    };

    GlobalSimuState() : globalStats(), tMutex() {

    };

    ~GlobalSimuState() {
#if defined(MOLFLOW)
        globalHistograms.clear();
        facetStates.clear();
#endif
    }

    bool initialized = false;

    void Clear();

    void Resize(std::shared_ptr<SimulationModel> model);

    void Reset();

    static std::tuple<int, int, int>
    Compare(const GlobalSimuState &lhsGlobHit, const GlobalSimuState &rhsGlobHit, double globThreshold,
            double locThreshold);

#if defined(MOLFLOW)
    GlobalHitBuffer globalStats;
    std::vector<FacetHistogramBuffer> globalHistograms; //1+nbMoment
    std::vector<FacetState> facetStates; //nbFacet
#endif

    template<class Archive>
    void serialize(Archive &archive) {
        archive(
                CEREAL_NVP(globalStats),
                CEREAL_NVP(globalHistograms),
                CEREAL_NVP(facetStates)
        );
    }

    mutable std::timed_mutex tMutex;
};


[[nodiscard]] std::optional<std::unique_lock<std::timed_mutex>> GetHitLock(GlobalSimuState* simStatePtr, size_t waitMs);

/*!
 * @brief Particle Log structure containing all individual log entries
 */
struct ParticleLog {
public:
    ParticleLog &operator=(const ParticleLog &src) {
        pLog = src.pLog;
        return *this;
    }

    ParticleLog(ParticleLog &&rhs) noexcept: tMutex() {
        pLog = std::move(pLog);
    };

    ParticleLog(const ParticleLog &rhs) {
        pLog = rhs.pLog;
    };

    ParticleLog() : tMutex() {

    };

    void resize(size_t nbLogs) {
        std::vector<ParticleLoggerItem>(nbLogs).swap(pLog);
    };

    void clear() {
        pLog.clear();
    };
    std::vector<ParticleLoggerItem> pLog;
    mutable std::timed_mutex tMutex;
};

#endif //MOLFLOW_PROJ_MOLFLOWSIMGEOM_H
