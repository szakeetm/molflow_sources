//
// Created by Pascal Baehr on 28.04.20.
//

#ifndef MOLFLOW_PROJ_GEOMETRYSIMU_H
#define MOLFLOW_PROJ_GEOMETRYSIMU_H

#include <vector>
#include "MolflowTypes.h"
#include "Buffer_shared.h"
#include "Parameter.h"

#include <mutex>
#include <FacetData.h>
#include <RayTracing/BVH.h>

#include <cereal/cereal.hpp>
#include <cereal/types/vector.hpp>
#include <RayTracing/KDTree.h>
#include <map>

struct SubprocessFacet;

struct ParticleLog;
class SuperStructure;

class ParameterSurface : public Surface {
    Distribution2D *dist;
public:
    ParameterSurface(Distribution2D *distribution) : dist(distribution) {};

    ~ParameterSurface() override = default;

    bool IsHardHit(const Ray &r) override;
};

struct TimeDependentParamters {
    TimeDependentParamters() = default;

    std::vector<Distribution2D> parameters;

    std::vector<std::vector<CDF_p>> CDFs; //cumulative distribution function for each temperature
    std::vector<std::vector<ID_p>> IDs; //integrated distribution function for each time-dependent desorption type
    std::vector<Moment> moments;             //moments when a time-dependent simulation state is recorded
    /*std::vector<UserMoment> userMoments;    //user-defined text values for defining time moments (can be time or time series)
    std::vector<double> temperatures; //keeping track of all temperatures that have a CDF already generated
    std::vector<size_t> desorptionParameterIDs; //time-dependent parameters which are used as desorptions, therefore need to be integrated
*/
    size_t GetMemSize() {
        size_t sum = 0;
        for (auto &par : parameters) {
            sum += par.GetMemSize();
        }
        sum += sizeof(std::vector<std::vector<CDF_p>>);
        for (auto &vec : CDFs) {
            sum += sizeof(std::vector<CDF_p>);
            sum += sizeof(std::pair<double, double>) * vec.capacity();
        }
        sum += sizeof(std::vector<std::vector<ID_p>>);
        for (auto &vec : IDs) {
            sum += sizeof(std::vector<ID_p>);
            sum += sizeof(std::pair<double, double>) * vec.capacity();
        }
        sum += sizeof(std::vector<Moment>);
        sum += sizeof(Moment) * moments.capacity();

        return sum;
    }
};

struct Anglemap {
public:
    std::vector<size_t> pdf;          // Incident angle distribution, phi and theta, not normalized. Used either for recording or for 2nd order interpolation
    std::vector<double> phi_CDFs;    // A table containing phi distributions for each theta, starting from 0 for every line (1 line = 1 theta value). For speed we keep it in one memory block, 1 pointer
    std::vector<size_t> phi_CDFsums; // since CDF runs only to the middle of the last segment, for each theta a line sum is stored here. Also a pdf for theta
    std::vector<double> theta_CDF;      // Theta CDF, not normalized. nth value is the CDF at the end of region n (beginning of first section is always 0)
    size_t theta_CDFsum; // since theta CDF only runs till the middle of the last segment, the map sum is here

    [[nodiscard]] size_t GetMemSize() const {
        size_t sum = 0;
        sum += sizeof(Anglemap);
        sum += sizeof(size_t) * pdf.capacity();
        sum += sizeof(double) * phi_CDFs.capacity();
        sum += sizeof(size_t) * phi_CDFsums.capacity();
        sum += sizeof(double) * theta_CDF.capacity();
        return sum;
    }
};

// Local facet structure
struct SubprocessFacet : public Facet {
    SubprocessFacet();

    explicit SubprocessFacet(size_t nbIndex);

    SubprocessFacet(const SubprocessFacet &o);

    SubprocessFacet(SubprocessFacet &&cpy) noexcept;

    SubprocessFacet &operator=(const SubprocessFacet &o);

    SubprocessFacet &operator=(SubprocessFacet &&o) noexcept;

    std::vector<double> textureCellIncrements;              // Texure increment
    std::vector<bool> largeEnough;      // cells that are NOT too small for autoscaling
    OutgassingMap ogMap;

    Anglemap angleMap; //TODO: -> GeneratingAngleMap from 2.7

    // Temporary var (used in FillHit for hit recording)
    bool isReady{};         // Volatile state

    // Facet hit counters
    //std::vector<FacetHitBuffer> tmpCounter; //1+nbMoment
    //std::vector<FacetHistogramBuffer> tmpHistograms; //1+nbMoment

    //void ResetCounter();
    //void ResizeCounter(size_t nbMoments);
    bool InitializeOnLoad(const size_t &id, const size_t &nbMoments);

    size_t InitializeHistogram(const size_t &nbMoments) const;

    size_t InitializeDirectionTexture(const size_t &nbMoments);

    size_t InitializeProfile(const size_t &nbMoments);

    size_t InitializeTexture(const size_t &nbMoments);

    size_t InitializeAngleMap();

    void InitializeOutgassingMap();

    bool InitializeLinkAndVolatile(const size_t &id);

    [[nodiscard]] size_t GetHitsSize(size_t nbMoments) const;

    [[nodiscard]] size_t GetMemSize() const;

    //void RegisterTransparentPass(SubprocessFacet *facet); //Allows one shared Intersect routine between MolFlow and Synrad

    std::vector<double> InitTextureMesh();
};

// Local simulation structure

class AABBNODE;
class GlobalSimuState;

class SuperStructure {
public:
    SuperStructure() = default;

    ~SuperStructure();

    //std::vector<SubprocessFacet>  facets;   // Facet handles
    std::shared_ptr<AABBNODE> aabbTree; // Structure AABB tree
    std::string strName;
    std::string strFileName;

    size_t GetMemSize() {
        size_t sum = 0;
        /*sum += sizeof (facets);
        for(auto& fac : facets)
            sum += fac.GetMemSize();*/
        sum += sizeof(aabbTree);
        return sum;
    }
};

struct SimulationModel {
public:
    SimulationModel() : otfParams(), tdParams(), wp(), sh(), m(), initialized(false) {};

    ~SimulationModel();

    SimulationModel(SimulationModel &&o) noexcept: m(), initialized(false) {
        *this = std::move(o);
    };

    SimulationModel(const SimulationModel &o) : m(), initialized(false) {
        *this = o;
    };

    size_t size() {
        size_t modelSize = 0;
        modelSize += facets.capacity();
        for (auto &fac : facets)
            modelSize += fac->GetMemSize();
        modelSize += structures.capacity();
        for (auto &struc : structures)
            modelSize += struc.GetMemSize();
        modelSize += sizeof(std::vector<Vector3d>) + sizeof(Vector3d) * vertices3.capacity();
        modelSize += tdParams.GetMemSize();
        modelSize += sizeof(otfParams);
        modelSize += sizeof(wp);
        modelSize += sizeof(sh);
        modelSize += sizeof(m);
        modelSize += sizeof(initialized);

        return modelSize;
    }

    SimulationModel &operator=(const SimulationModel &o) {
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

    SimulationModel &operator=(SimulationModel &&o) noexcept {
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

    int PrepareToRun();

    int BuildAccelStructure(GlobalSimuState *globState, int accel_type, int split,
                            int maxPrimsInNode);

    int InitialiseFacets();

    void CalcTotalOutgassing();

    void CalculateFacetParams(SubprocessFacet *f);

    Surface *GetSurface(double opacity) {

        if (!surfaces.empty()) {
            auto surf = surfaces.find(opacity);
            if (surf != surfaces.end())
                return surf->second.get();
        }
        std::shared_ptr<Surface> surface;
        if (opacity == 1.0) {
            surface = std::make_shared<Surface>();
        } else if (opacity == 0.0) {
            surface = std::make_shared<TransparentSurface>();
        } else {
            surface = std::make_shared<AlphaSurface>(opacity);
        }
        surfaces.insert(std::make_pair(opacity, surface));
        return surface.get();
    };

    Surface *GetParameterSurface(int opacity_paramId, Distribution2D *dist) {

        double indexed_id = 10.0 + opacity_paramId;
        if (!surfaces.empty()) {
            auto surf = surfaces.find(indexed_id);
            if (surf != surfaces.end())
                return surf->second.get();
        }

        std::shared_ptr<ParameterSurface> surface;
        surface = std::make_shared<ParameterSurface>(dist);
        surfaces.insert(std::make_pair(indexed_id, surface));
        printf("Insert param id: %f\n", indexed_id);
        return surface.get();
    };

    // Sim functions
    double GetOpacityAt(SubprocessFacet *f, double time) const;

    double GetStickingAt(SubprocessFacet *f, double time) const;

    // Geometry Description
    std::vector<std::shared_ptr<SubprocessFacet>> facets;    // All facets of this geometry

    std::vector<SuperStructure> structures;
    std::vector<Vector3d> vertices3; // Vertices (3D space)

    std::vector<std::shared_ptr<RTAccel>> accel;
    std::map<double, std::shared_ptr<Surface>> surfaces;

    // Simulation Properties
    OntheflySimulationParams otfParams;
    TimeDependentParamters tdParams;
    WorkerParams wp;

    // Geometry Properties
    GeomProperties sh;

    bool initialized;
    std::mutex m;

    static std::vector<double>
    ComputeHitChances(const std::vector<TestRay> &battery, const std::vector<std::shared_ptr<Facet>> &primitives);
    int ComputeHitStats(const std::vector<TestRay> &battery);

    bool StartFromSource(Ray &ray);

    void PerformBounce(Ray &ray, SubprocessFacet *iFacet);
};

/*!
 * @brief One instance is the state for one facet for a single moment
 */
class FacetMomentSnapshot {
public:
    FacetMomentSnapshot();
    FacetMomentSnapshot &operator+=(const FacetMomentSnapshot &rhs);

    FacetMomentSnapshot &operator+(const FacetMomentSnapshot &rhs);

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

class GlobalSimuState { //replaces old hits dataport
public:
    GlobalSimuState &operator=(const GlobalSimuState &src);

    GlobalSimuState &operator+=(const GlobalSimuState &src);

    GlobalSimuState(GlobalSimuState &&rhs) noexcept: tMutex() {
        globalHistograms = std::move(rhs.globalHistograms);
        facetStates = std::move(rhs.facetStates);
        globalHits = rhs.globalHits;
        initialized = rhs.initialized;
    };

    GlobalSimuState(const GlobalSimuState &rhs) {
        globalHits = rhs.globalHits;
        globalHistograms = rhs.globalHistograms;
        facetStates = rhs.facetStates;
        initialized = rhs.initialized;
    };

    GlobalSimuState() : globalHits(), tMutex() {

    };

    ~GlobalSimuState() {
#if defined(MOLFLOW)
        globalHistograms.clear();
        facetStates.clear();
#endif
    }

    bool initialized = false;
    bool stateChanged = false;

    void clear();

    void Resize(const SimulationModel &model);

    void Reset();

    static std::tuple<int, int, int>
    Compare(const GlobalSimuState &lhsGlobHit, const GlobalSimuState &rhsGlobHit, double globThreshold,
            double locThreshold);

    std::vector<TestRay> PrepareHitBattery();
    int UpdateBatteryFrequencies();
    SampleBattery hitBattery;

#if defined(MOLFLOW)
    GlobalHitBuffer globalHits;
    std::vector<FacetHistogramBuffer> globalHistograms; //1+nbMoment
    std::vector<FacetState> facetStates; //nbFacet
#endif

    template<class Archive>
    void serialize(Archive &archive) {
        archive(
                CEREAL_NVP(globalHits),
                CEREAL_NVP(globalHistograms),
                CEREAL_NVP(facetStates)
        );
    }

    mutable std::timed_mutex tMutex;

    int StopBatteryChange();
};

struct ParticleLog {
public:
    ParticleLog &operator=(const ParticleLog &src) {
        pLog = src.pLog;
        return *this;
    }

    ParticleLog(ParticleLog &&rhs) noexcept: tMutex() {
        pLog = std::move(rhs.pLog);
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

#endif //MOLFLOW_PROJ_GEOMETRYSIMU_H
