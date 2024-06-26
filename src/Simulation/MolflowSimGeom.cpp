#include <set>
#include <sstream>
#include <filesystem>
#include <memory>

#include "Helper/MathTools.h"
#include "CDFGeneration.h"
#include "IDGeneration.h"
#include "Helper/Chronometer.h"
#include "Helper/ConsoleLogger.h"
#include "Polygon.h"
#include "MolflowSimGeom.h"
#include "MolflowSimFacet.h"
#include "IntersectAABB_shared.h" // include needed for recursive delete of AABBNODE
#include "File.h" //FileReader for LoadParamCatalog
#include "GLApp/GLTypes.h"

/**
* \brief Computes with a distribution function and a random number whether a hit on a surface is "hard" or not
* \param r instance of the handled ray for the intersection test
 * \return true on hard hit
*/
bool ParameterSurface::IsHardHit(const Ray &r) {
    const double td_opacity = timeDepParam->InterpY(r.time, false);
    if(td_opacity >= 1.0)
        return true;
    else
        return (r.rng->rnd() < td_opacity);
};

/**
* \brief Testing purpose function, construct an angled PRISMA / parallelepiped
* \param L length
* \param R radius
* \param angle angle between Y-Z
* \param s sticking value
* \param step number of facets used to construct the circular hull
*/
/*
void MolflowSimulationModel::BuildPrisma(double L, double R, double angle, double s, int step) {

    int nbTF = 9 * nbDecade;
    int nbTV = 4 * nbTF;

    sh.nbVertex = 2 * step + nbTV;
    std::vector<Vector3d>(sh.nbVertex).swap(vertices3);

    sh.nbFacet = step + 2 + nbTF;


    sh.nbSuper = 1;
    structures.resize(1);
    structures.begin()->name = "Prisma";

    try{
        facets.resize(sh.nbFacet, nullptr);
    }
    catch(const std::exception &) {
        throw Error("Couldn't allocate memory for facets");
    }

    // Vertices
    for (int i = 0; i < step; i++) {
        double step_angle = (double)i / (double)step * 2 * PI;
        vertices3[2 * i + nbTV].x = R * cos(step_angle);
        vertices3[2 * i + nbTV].y = R * sin(step_angle);
        vertices3[2 * i + nbTV].z = 0;
        vertices3[2 * i + 1 + nbTV].x = R * cos(step_angle);
        vertices3[2 * i + 1 + nbTV].y = R * sin(step_angle) + L * cos(M_PI_2 - angle);
        vertices3[2 * i + 1 + nbTV].z = L * cos(angle);
    }

    try {
        // Cap facet
        facets[0 + nbTF] = std::make_shared<MolflowSimFacet>(step);
        facets[0 + nbTF]->sh.sticking = 1.0;
        facets[0 + nbTF]->sh.desorbType = DES_COSINE;
        facets[0 + nbTF]->sh.outgassing = 1.0;
        for (int i = 0; i < step; i++)
            facets[0 + nbTF]->indices[i] = 2 * i + nbTV;

        facets[1 + nbTF] = std::make_shared<MolflowSimFacet>(step);
        facets[1 + nbTF]->sh.sticking = 1.0;
        facets[1 + nbTF]->sh.desorbType = DES_NONE;
        for (int i = 0; i < step; i++)
            facets[1 + nbTF]->indices[step - i - 1] = 2 * i + 1 + nbTV;

        // Wall facet
        for (int i = 0; i < step; i++) {
            facets[i + 2 + nbTF] = std::make_shared<MolflowSimFacet>(4);
            //facets[i + 2 + nbTF]->sp.reflection.diffusePart = 1.0; //constructor does this already
            //facets[i + 2 + nbTF]->sp.reflection.specularPart = 0.0; //constructor does this already
            facets[i + 2 + nbTF]->sh.sticking = s;
            facets[i + 2 + nbTF]->indices[0] = 2 * i + nbTV;
            facets[i + 2 + nbTF]->indices[1] = 2 * i + 1 + nbTV;
            if (i < step - 1) {
                facets[i + 2 + nbTF]->indices[2] = 2 * (i + 1) + 1 + nbTV;
                facets[i + 2 + nbTF]->indices[3] = 2 * (i + 1) + nbTV;
            }
            else {

                facets[i + 2 + nbTF]->indices[2] = 1 + nbTV;
                facets[i + 2 + nbTF]->indices[3] = 0 + nbTV;
            }
        }
    }
    catch (std::bad_alloc) {
        throw Error("Couldn't reserve memory for the facets");
    }
    catch (...) {
        throw Error("Unspecified Error while building pipe");
    }
}
*/

/**
* \brief Builds ADS given certain parameters
* \param globalState global simulation state for splitting techniques requiring statistical data
* \param accel_type BVH or KD tree
* \param split splitting technique corresponding to the selected AccelType
* \param bvh_width for BVH, the amount of leaves per end node
 * \return 0> for error codes, 0 when no problems
*/
int MolflowSimulationModel::BuildAccelStructure(const std::shared_ptr<GlobalSimuState> globalState, AccelType accel_type, BVHAccel::SplitMethod split,
                                                int bvh_width) {

    initialized = false;
    Chronometer timer;
    timer.Start();

    std::lock_guard<std::mutex> lock(modelMutex);

    std::vector<std::vector<std::shared_ptr<Primitive>>> primPointers;
    primPointers.resize(this->sh.nbSuper);
    for(auto& sFac : this->facets){
        if (sFac->sh.superIdx == -1) { //Facet in all structures
            for (auto& fp_vec : primPointers) {
                fp_vec.push_back(sFac);
            }
        }
        else {
            primPointers[sFac->sh.superIdx].push_back(sFac); //Assign to structure
        }
    }

    for(auto& sFac : this->facets){
        auto mfFac = std::static_pointer_cast<MolflowSimFacet>(sFac);
        if (mfFac->opacity_paramId == -1) { //constant sticking
            sFac->sh.opacity = std::clamp(sFac->sh.opacity, 0.0, 1.0); //sanitize
        }
        sFac->surf = GetSurface(mfFac);
    }

    this->rayTracingStructures.clear();
    if(BVHAccel::SplitMethod::ProbSplit == split && globalState && globalState->initialized && globalState->globalStats.globalHits.nbDesorbed > 0){
        if(globalState->facetStates.size() != this->facets.size())
            return 1;
        std::vector<double> probabilities;
        probabilities.reserve(globalState->facetStates.size());
        for(auto& state : globalState->facetStates) {
            probabilities.emplace_back(state.momentResults[0].hits.nbHitEquiv / globalState->globalStats.globalHits.nbHitEquiv);
        }
        for (size_t s = 0; s < this->sh.nbSuper; ++s) {
            if(accel_type == AccelType::KD)
                this->rayTracingStructures.emplace_back(std::make_unique<KdTreeAccel>(primPointers[s], probabilities));
            else
                this->rayTracingStructures.emplace_back(std::make_unique<BVHAccel>(primPointers[s], bvh_width, BVHAccel::SplitMethod::ProbSplit, probabilities));
        }
    }
    else {
        for (size_t s = 0; s < this->sh.nbSuper; ++s) {
            if(accel_type == AccelType::KD)
                this->rayTracingStructures.emplace_back(std::make_unique<KdTreeAccel>(primPointers[s]));
            else
                this->rayTracingStructures.emplace_back(std::make_unique<BVHAccel>(primPointers[s], bvh_width, split));
        }
    }

    timer.Stop();

    initialized = true;

    return 0;
}

/**
* \brief Do calculations necessary before launching simulation
* determine latest moment
* Generate integrated desorption functions
* match parameters
* Generate speed distribution functions
* Angle map
*/
void MolflowSimulationModel::PrepareToRun() {
    std::lock_guard<std::mutex> lock(modelMutex);

    initialized = false;

    std::string errLog;

    //determine latest moment
    if(!tdParams.moments.empty()) {
        sp.latestMoment = tdParams.moments.back().time + .5 * tdParams.moments.back().window;
    } else {
        sp.latestMoment = sp.timeWindowSize * .5;
    }

    CalcIntervalCache();
    this->maxwell_CDF_1K = CDFGeneration::Generate_CDF(1.0, sp.gasMass);

    std::set<size_t> desorptionParameterIDs;
    std::vector<double> temperatureList;

    //Check and calculate various facet properties for time dependent simulations (CDF, ID )
    for (size_t i = 0; i < sh.nbFacet; i++) {
        const auto facet = std::static_pointer_cast<MolflowSimFacet>(facets[i]);
        // TODO: Find a solution to integrate catalog parameters
        if(facet->outgassing_paramId >= (int) tdParams.parameters.size()){
            throw Error("Facet {} outgassing refers to time-dependent parameter no. {}, but there are only {}.", i + 1, facet->outgassing_paramId + 1, tdParams.parameters.size());
        }
        if(facet->opacity_paramId >= (int) tdParams.parameters.size()){
            throw Error("Facet {} opacity refers to time-dependent parameter no. {}, but there are only {}.", i + 1, facet->opacity_paramId + 1, tdParams.parameters.size());
        }
        if(facet->sticking_paramId >= (int) tdParams.parameters.size()){
            throw Error("Facet {} sticking refers to time-dependent parameter no. {}, but there are only {}.", i + 1, facet->sticking_paramId + 1, tdParams.parameters.size());
        }

        if (facet->outgassing_paramId >= 0) { //if time-dependent desorption
            int id = IDGeneration::GetIDId(desorptionParameterIDs, facet->outgassing_paramId);
            if (id >= 0)
                facet->sh.IDid = id; //we've already generated an ID for this parameter
            else {
                auto[id_new, id_vec] = IDGeneration::GenerateNewID(desorptionParameterIDs, facet->outgassing_paramId, this);
                facet->sh.IDid = id_new;
                tdParams.IDs.emplace_back(std::move(id_vec));
            }
        }
        
        //Angle map
        if (facet->sh.desorbType == DES_ANGLEMAP) {
            auto mfFacet = std::static_pointer_cast<MolflowSimFacet>(facets[i]);
            if (mfFacet->angleMap.pdf.empty()) {
                throw Error("Facet #{} uses angle map desorption but doesn't have a recorded angle map.", i + 1);
            }
            if (mfFacet->sh.anglemapParams.record) {
                throw Error("Facet #{}: Can't RECORD and USE angle map desorption at the same time.", i + 1);
            }
        }
    }

    CalcTotalOutgassing();
    initialized = true;
}

void MolflowSimulationModel::CalcIntervalCache() {
    intervalCache.resize(tdParams.moments.size());
    for (size_t i=0;i<tdParams.moments.size();i++) {
        intervalCache[i].startTime=tdParams.moments[i].time - .5 * tdParams.moments[i].window;
        intervalCache[i].endTime=tdParams.moments[i].time + .5 * tdParams.moments[i].window;
    }
}

std::shared_ptr<Surface> MolflowSimulationModel::GetSurface(std::shared_ptr<SimulationFacet> facet) {
    const auto mfFac = std::static_pointer_cast<MolflowSimFacet>(facet);

    if (mfFac->opacity_paramId == -1) { //constant sticking
        return SimulationModel::GetSurface(facet); //simply based on opacity
    }
    else { //time-dependent opacity
        auto opacity_paramId = mfFac->opacity_paramId;
        Parameter* parPtr = &(this->tdParams.parameters[mfFac->opacity_paramId]);

        double indexed_id = 10.0 + opacity_paramId;
        if (!surfaces.empty()) {
            auto surf = surfaces.find(indexed_id);
            if (surf != surfaces.end())
                return surf->second;
        }
        //not found, make new
        auto surface = std::make_shared<ParameterSurface>(parPtr);
        surfaces.insert(std::make_pair(indexed_id, surface));
        Log::console_msg_master(5, "Insert surface with param id: {}\n", indexed_id);
        return surface;
    }
};

/**
* \brief Compute the outgassing of all source facet depending on the mode (file, regular, time-dependent) and set it to the global settings
*/
void MolflowSimulationModel::CalcTotalOutgassing() {
    // Compute the outgassing of all source facet
    double totalDesorbedMolecules = 0.0;
    double finalOutgassingRate_Pa_m3_sec = 0.0;
    double finalOutgassingRate = 0.0;

    const double latestMoment = sp.latestMoment;


    for (size_t i = 0; i < facets.size(); i++) {
        const auto facet = facets[i];
        const auto mfFacet = std::static_pointer_cast<MolflowSimFacet>(facet);
        if (facet->sh.desorbType != DES_NONE) { //there is a kind of desorption
            if (facet->sh.useOutgassingFile) { //outgassing file
                const auto mfFacet = std::static_pointer_cast<MolflowSimFacet>(facets[i]);
                auto &ogMap = mfFacet->ogMap;
                for (size_t l = 0; l < (ogMap.outgassingMapWidth * ogMap.outgassingMapHeight); l++) {
                    totalDesorbedMolecules +=
                            latestMoment * ogMap.outgassingMap[l] / (1.38E-23 * facet->sh.temperature);
                    finalOutgassingRate += ogMap.outgassingMap[l] / (1.38E-23 * facet->sh.temperature);
                    finalOutgassingRate_Pa_m3_sec += ogMap.outgassingMap[l];
                }
            } else { //regular outgassing
                if (mfFacet->sh.outgassingParam.empty()) { //constant outgassing
                    totalDesorbedMolecules +=
                            latestMoment * facet->sh.outgassing / (1.38E-23 * facet->sh.temperature);
                    finalOutgassingRate +=
                            facet->sh.outgassing / (1.38E-23 * facet->sh.temperature);  //Outgassing molecules/sec
                    finalOutgassingRate_Pa_m3_sec += facet->sh.outgassing;
                } else { //time-dependent outgassing
                    totalDesorbedMolecules +=
                            tdParams.IDs[facet->sh.IDid].back().cumulativeDesValue / (1.38E-23 * facet->sh.temperature);
                    size_t lastIndex = tdParams.parameters[mfFacet->outgassing_paramId].GetSize() - 1;
                    double finalRate_mbar_l_s = tdParams.parameters[mfFacet->outgassing_paramId].GetY(lastIndex);
                    finalOutgassingRate +=
                            finalRate_mbar_l_s * MBARLS_TO_PAM3S / (1.38E-23 * facet->sh.temperature);
                    finalOutgassingRate_Pa_m3_sec += finalRate_mbar_l_s * MBARLS_TO_PAM3S;
                }
            }
        }
    }

    sp.totalDesorbedMolecules = totalDesorbedMolecules;
    sp.finalOutgassingRate_Pa_m3_sec = finalOutgassingRate_Pa_m3_sec;
    sp.finalOutgassingRate = finalOutgassingRate;
}

std::vector<std::string> MolflowSimulationModel::SanityCheck() {
    //Must be called after totaloutgassing is calculated

    std::vector<std::string> errLog;

    if (!initialized) {
        errLog.push_back("Model not initialized");
    }
    if (vertices3.empty()) {
        errLog.push_back("Loaded empty vertex list");
    }
    if (facets.empty()) {
        errLog.push_back("Loaded empty facet list");
    }
    if (sh.nbFacet != facets.size()) {
        char tmp[256];
        snprintf(tmp, 256, "Facet structure not properly initialized, size mismatch: %zu / %zu\n", sh.nbFacet, facets.size());
        errLog.push_back(tmp);
    }
    for (auto& fac : facets) {
        bool hasAnyTexture = fac->sh.countDes || fac->sh.countAbs || fac->sh.countRefl || fac->sh.countTrans || fac->sh.countACD || fac->sh.countDirection;
        if (!fac->sh.isTextured && (fac->sh.texHeight * fac->sh.texHeight > 0)) {
            char tmp[256];
            snprintf(tmp, 256, "[Facet #%zu] Untextured facet with texture size\n", fac->globalId + 1);
            errLog.push_back(tmp);
        }
        else if (!fac->sh.isTextured && (hasAnyTexture)) {
            fac->sh.countDes = false;
            fac->sh.countAbs = false;
            fac->sh.countRefl = false;
            fac->sh.countTrans = false;
            fac->sh.countACD = false;
            fac->sh.countDirection = false;
            char tmp[256];
            snprintf(tmp, 256, "[Fac #%zu] Untextured facet with texture counters\n", fac->globalId + 1);
            errLog.push_back(tmp);
        }
        if (fac->sh.desorbType != DES_NONE && !fac->sh.temperatureParam.empty()) {
            errLog.push_back(fmt::format("[Facet {}]: Time-dependent temperature not allowed on facets with outgassing", fac->globalId + 1));
        }
    }

    //Molflow unique
    if (sp.enableDecay && sp.halfLife <= 0.0) {
        char tmp[255];
        sprintf(tmp, "Particle decay is set, but half life was not set [= %e]\n", sp.halfLife);
        errLog.push_back(tmp);
    }

    // Is there some desorption in the system? (depends on pre calculation)
    if (sp.finalOutgassingRate_Pa_m3_sec <= 0.0) {
        // Do another check for existing desorp facets, needed in case a desorp parameter's final value is 0
        bool found = false;
        size_t nbF = facets.size();
        size_t i = 0;
        while (i < nbF && !found) {
            found = (facets[i]->sh.desorbType != DES_NONE);
            if (!found) i++;
        }

        if (!found)
            throw Error("No desorption facet found");
    }
    if (sp.totalDesorbedMolecules <= 0.0)
        throw Error("Total outgassing is zero.");

    return errLog;
}

/*

MolflowSimulationModel::MolflowSimulationModel(MolflowSimulationModel &&o) noexcept {
    *this = std::move(o);
};

MolflowSimulationModel::MolflowSimulationModel(const MolflowSimulationModel &o) : SimulationModel(o) {
    *this = o;
};

MolflowSimulationModel::~MolflowSimulationModel() = default;
*/

/**
* \brief Calculates the used memory used by the whole simulation model
 * \return memory size used by the whole simulation model
*/
size_t MolflowSimulationModel::GetMemSize() {
    size_t modelSize = 0;
    modelSize += SimulationModel::GetMemSize(); //base class members
    //Molflow-specific members:
    modelSize += tdParams.GetMemSize();
    modelSize += sizeof(Interval) * intervalCache.capacity();
    modelSize += sizeof(IntegratedVelocityEntry) * maxwell_CDF_1K.capacity();
    return modelSize;
}

/**
* \brief Assign operator
* \param src reference to source object
* \return address of this
*/
GlobalSimuState& GlobalSimuState::operator=(const GlobalSimuState&  src) {
    //Copy all but mutex
    facetStates = src.facetStates;
    globalHistograms = src.globalHistograms;
    globalStats = src.globalStats;
    initialized = src.initialized;
    return *this;
}

/**
* \brief Assign operator
* \param src reference to source object
* \return address of this
*/
GlobalSimuState& GlobalSimuState::operator+=(const GlobalSimuState& rhs) {
    //Copy all but mutex
    facetStates += rhs.facetStates;
    globalHistograms += rhs.globalHistograms;
    globalStats += rhs.globalStats;
    return *this;
}

size_t GlobalSimuState::GetMemSize() const
{
    size_t sum = 0;
    sum += sizeof(globalStats); //plain old data
    for (const auto& fs : facetStates) sum += fs.GetMemSize();
    for (const auto& hist : globalHistograms) sum += hist.GetMemSize();
    return sum;
}

/**
* \brief Clears simulation state
*/

void GlobalSimuState::Clear() {
    auto lock = GetHitLock(this, 10000);
    if (!lock) return;
    globalStats = GlobalHitBuffer();
    globalHistograms.clear();
    facetStates.clear();
    initialized = false;
}

/**
* \brief Constructs the 'Global Hit counter structure' structure to hold all results, zero-init
* \param model Contains all related parameters
*/
void GlobalSimuState::Resize(std::shared_ptr<SimulationModel> model) {

    auto lock = GetHitLock(this, 10000);
    if (!lock) return;
    auto mf_model = std::static_pointer_cast<MolflowSimulationModel>(model);
    size_t nbF = model->sh.nbFacet;
    size_t nbMoments = mf_model->tdParams.moments.size();
    facetStates.assign(nbF, FacetState());

    if(!model->facets.empty()) {
        for (size_t i = 0; i < nbF; i++) {
            auto sFac = model->facets[i];
            if (sFac->globalId != i) {
                fmt::print(stderr, "Facet ID mismatch! : {} vs {}\n", sFac->globalId, i);
                exit(0);
            }

            FacetMomentSnapshot facetMomentTemplate{};
            facetMomentTemplate.histogram.Resize(sFac->sh.facetHistogramParams);
            facetMomentTemplate.direction.assign((sFac->sh.countDirection ? sFac->sh.texWidth*sFac->sh.texHeight : 0), DirectionCell());
            facetMomentTemplate.profile.assign((sFac->sh.isProfile ? PROFILE_SIZE : 0), ProfileSlice());
            facetMomentTemplate.texture.assign((sFac->sh.isTextured ? sFac->sh.texWidth*sFac->sh.texHeight : 0), TextureCell());

            //No init for hits
            facetStates[i].momentResults.assign(1 + nbMoments, facetMomentTemplate);
            if (sFac->sh.anglemapParams.record)
                facetStates[i].recordedAngleMapPdf.assign(sFac->sh.anglemapParams.GetMapSize(), 0);
        }
    }

    //Global histogram

    FacetHistogramBuffer globalHistTemplate{};
    globalHistTemplate.Resize(model->sp.globalHistogramParams);
    globalHistograms.assign(1 + nbMoments, globalHistTemplate);
    initialized = true;
}

/**
* \brief zero-init for all structures
*/
void GlobalSimuState::Reset() {
    auto lock = GetHitLock(this, 10000);
    if (!lock) return;
    for (auto& h : globalHistograms) {
        ZEROVECTOR(h.distanceHistogram);
        ZEROVECTOR(h.nbHitsHistogram);
        ZEROVECTOR(h.timeHistogram);
    }
    memset(&globalStats, 0, sizeof(globalStats)); //Plain old data
    for (auto& state : facetStates) {
        ZEROVECTOR(state.recordedAngleMapPdf);
        for (auto& m : state.momentResults) {
            ZEROVECTOR(m.histogram.distanceHistogram);
            ZEROVECTOR(m.histogram.nbHitsHistogram);
            ZEROVECTOR(m.histogram.timeHistogram);
            std::fill(m.direction.begin(), m.direction.end(), DirectionCell());
            std::fill(m.texture.begin(), m.texture.end(), TextureCell());
            std::fill(m.profile.begin(), m.profile.end(), ProfileSlice());
            memset(&(m.hits), 0, sizeof(m.hits));
        }
    }
}

/**
 * @brief Compare function for two simulation states
 * @param lhsGlobHit first simulation state
 * @param rhsGlobHit second simulation state
 * @param globThreshold threshold for relative difference on global counters
 * @param locThreshold threshold for relative difference on facet local counters
 * @return a tuple containing the number of global, local (facet), and fine (facet profile/texture) errors
 */
std::tuple<size_t, size_t, size_t>
GlobalSimuState::Compare(const std::shared_ptr<GlobalSimuState> lhsGlobHit, const std::shared_ptr<GlobalSimuState> rhsGlobHit, double globThreshold,
                         double locThreshold) {

    const double velocityThresholdFactor = 40.0;
    //std::ofstream cmpFile("cmpFile.txt");
    size_t globalErrNb = 0;
    size_t facetErrNb = 0;
    size_t fineErrNb = 0;

    size_t nbFacetSkips = 0;
    size_t nbProfileSkips = 0;
    size_t nbTextureSkips = 0;
    size_t nbDirSkips = 0;
    size_t nbHistSkips_glob = 0;
    size_t nbHistSkips_loc = 0;

    std::string cmpFile;
    std::string cmpFileFine; // extra stream to silence after important outputs

    // Sanity check
    {
        if(lhsGlobHit->globalStats.globalHits.nbDesorbed == 0 && rhsGlobHit->globalStats.globalHits.nbDesorbed == 0){
            cmpFile += fmt::format("[Global][desorp] Neither state has recorded desorptions\n");
            ++globalErrNb;
        }
        else if (lhsGlobHit->globalStats.globalHits.nbDesorbed == 0){
            cmpFile += fmt::format("[Global][desorp] First state has no recorded desorptions\n");
            ++globalErrNb;
        }
        else if (rhsGlobHit->globalStats.globalHits.nbDesorbed == 0){
            cmpFile += fmt::format("[Global][desorp] Second state has no recorded desorptions\n");
            ++globalErrNb;
        }

        if(globalErrNb){ // preemptive return, as global errors make local errors irrelevant
            fmt::print("{}\n", cmpFile);
            return std::make_tuple(static_cast<int>(globalErrNb), -1, -1);
        }
    }

    {
        double absRatio = lhsGlobHit->globalStats.globalHits.nbAbsEquiv / static_cast<double>(lhsGlobHit->globalStats.globalHits.nbDesorbed);
        double absRatio_rhs = rhsGlobHit->globalStats.globalHits.nbAbsEquiv / static_cast<double>(rhsGlobHit->globalStats.globalHits.nbDesorbed);
        if (!IsEqual(absRatio, absRatio_rhs, globThreshold)) {
            cmpFile += fmt::format("[Global][absRatio abs/des] has large difference: {} vs {}\n"
                                   "    {}abs / {}des vs {}abs / {}des\n",
                                   absRatio, absRatio_rhs,
                                lhsGlobHit->globalStats.globalHits.nbAbsEquiv,lhsGlobHit->globalStats.globalHits.nbDesorbed,
                                rhsGlobHit->globalStats.globalHits.nbAbsEquiv,rhsGlobHit->globalStats.globalHits.nbDesorbed);
            ++globalErrNb;
        }
    }

    {
        double hitRatio = static_cast<double>(lhsGlobHit->globalStats.globalHits.nbMCHit) / static_cast<double>(lhsGlobHit->globalStats.globalHits.nbDesorbed);
        double hitRatio_rhs = static_cast<double>(rhsGlobHit->globalStats.globalHits.nbMCHit) / static_cast<double>(rhsGlobHit->globalStats.globalHits.nbDesorbed);
        if (!IsEqual(hitRatio, hitRatio_rhs, globThreshold)) {
            cmpFile += fmt::format("[Global][hits/des] has large difference: {} vs {}\n"
                                   "    {}hit / {}des vs {}hit / {}des\n",
                                   hitRatio,hitRatio_rhs,
                                   lhsGlobHit->globalStats.globalHits.nbMCHit, lhsGlobHit->globalStats.globalHits.nbDesorbed,
                                   rhsGlobHit->globalStats.globalHits.nbMCHit, rhsGlobHit->globalStats.globalHits.nbDesorbed);

            ++globalErrNb;
        }
    }

    if (!IsEqual(lhsGlobHit->globalStats.globalHits.sum_v_ort, rhsGlobHit->globalStats.globalHits.sum_v_ort, globThreshold)) {
        cmpFile += fmt::format("[Global][sum_v_ort] has large difference: {} vs {}\n",
                               lhsGlobHit->globalStats.globalHits.sum_v_ort, rhsGlobHit->globalStats.globalHits.sum_v_ort);
        ++globalErrNb;
    }
    if (!IsEqual(lhsGlobHit->globalStats.globalHits.sum_1_per_velocity, rhsGlobHit->globalStats.globalHits.sum_1_per_velocity, globThreshold)) {
        cmpFile += fmt::format("[Global][sum_1_per_velocity] has large difference: {} vs {}\n",
                               lhsGlobHit->globalStats.globalHits.sum_1_per_velocity, rhsGlobHit->globalStats.globalHits.sum_1_per_velocity);
        ++globalErrNb;
    }
    if (!IsEqual(lhsGlobHit->globalStats.globalHits.sum_1_per_ort_velocity, rhsGlobHit->globalStats.globalHits.sum_1_per_ort_velocity, globThreshold)) {
        cmpFile += fmt::format("[Global][sum_1_per_ort_velocity] has large difference: {} vs {}\n",
                               lhsGlobHit->globalStats.globalHits.sum_1_per_ort_velocity, rhsGlobHit->globalStats.globalHits.sum_1_per_ort_velocity);
        ++globalErrNb;
    }


    // Histogram
    {
        auto& hist_lhs = lhsGlobHit->globalHistograms;
        auto& hist_rhs = rhsGlobHit->globalHistograms;
        for(size_t tHist = 0; tHist < hist_lhs.size(); tHist++) {
            for (size_t hIndex = 0; hIndex < hist_lhs[tHist].nbHitsHistogram.size(); ++hIndex) {
                if(std::sqrt(std::max(1.0,std::min(hist_lhs[tHist].nbHitsHistogram[hIndex], hist_rhs[tHist].nbHitsHistogram[hIndex]))) < 80) {
                    // Sample size not large enough
                    continue;
                }
                double lhRatio = hist_lhs[tHist].nbHitsHistogram[hIndex] / static_cast<double>(lhsGlobHit->globalStats.globalHits.nbMCHit);
                double rhRatio =  hist_rhs[tHist].nbHitsHistogram[hIndex] / static_cast<double>(rhsGlobHit->globalStats.globalHits.nbMCHit);
                if (!IsEqual(lhRatio, rhRatio, locThreshold)) {
                    cmpFile += fmt::format("[Global][Hist][Bounces/global hits][Moment={}] has large difference: {} vs {}\n",
                                           hIndex, lhRatio, rhRatio);
                    ++globalErrNb;
                }
            }
            for (size_t hIndex = 0; hIndex < hist_lhs[tHist].distanceHistogram.size(); ++hIndex) {
                if(std::sqrt(std::max(1.0,std::min(hist_lhs[tHist].distanceHistogram[hIndex], hist_rhs[tHist].distanceHistogram[hIndex]))) < 80) {
                    // Sample size not large enough
                    continue;
                }
                if (!IsEqual(hist_lhs[tHist].distanceHistogram[hIndex] / static_cast<double>(lhsGlobHit->globalStats.globalHits.nbMCHit), hist_rhs[tHist].distanceHistogram[hIndex] / static_cast<double>(rhsGlobHit->globalStats.globalHits.nbMCHit), locThreshold)) {
                    cmpFile += fmt::format("[Global][Hist][Dist][Ind={}] has large difference: {}\n",
                                           hIndex, std::abs(hist_lhs[tHist].distanceHistogram[hIndex] / static_cast<double>(lhsGlobHit->globalStats.globalHits.nbMCHit) - hist_rhs[tHist].distanceHistogram[hIndex] / static_cast<double>(rhsGlobHit->globalStats.globalHits.nbMCHit)));
                    ++globalErrNb;
                }
            }
            for (size_t hIndex = 0; hIndex < hist_lhs[tHist].timeHistogram.size(); ++hIndex) {
                if(std::sqrt(std::max(1.0,std::min(hist_lhs[tHist].timeHistogram[hIndex], hist_rhs[tHist].timeHistogram[hIndex]))) < 80) {
                    // Sample size not large enough
                    continue;
                }
                if (!IsEqual(hist_lhs[tHist].timeHistogram[hIndex] / static_cast<double>(lhsGlobHit->globalStats.globalHits.nbMCHit), hist_rhs[tHist].timeHistogram[hIndex] / static_cast<double>(rhsGlobHit->globalStats.globalHits.nbMCHit), locThreshold)) {
                    cmpFile += fmt::format("[Global][Hist][Time][Ind={}] has large difference: {}\n",
                                           hIndex, std::abs(hist_lhs[tHist].timeHistogram[hIndex] / static_cast<double>(lhsGlobHit->globalStats.globalHits.nbMCHit) - hist_rhs[tHist].timeHistogram[hIndex] / static_cast<double>(rhsGlobHit->globalStats.globalHits.nbMCHit)));

                    ++globalErrNb;
                }
            }
        }
    }
    // facets

    auto locThreshold_bak = locThreshold;
    if (lhsGlobHit->facetStates.size() != rhsGlobHit->facetStates.size()) {
        cmpFile += fmt::format("[FacetStates] Different number of facet states: {} vs {}\n", lhsGlobHit->facetStates.size(), rhsGlobHit->facetStates.size());
        ++facetErrNb;
    }

    size_t nbCmp = std::min(lhsGlobHit->facetStates.size(), rhsGlobHit->facetStates.size());
    for(int facetId = 0; facetId < nbCmp; ++facetId)
    {//cmp
        if(lhsGlobHit->facetStates[facetId].momentResults.size() != rhsGlobHit->facetStates[facetId].momentResults.size()){
            cmpFile += fmt::format("[Facet][{}] Different amount of moments for each state: {} vs {}\n", facetId, lhsGlobHit->facetStates[facetId].momentResults.size(), rhsGlobHit->facetStates[facetId].momentResults.size());
            ++facetErrNb;
            continue;
        }

        for(int m = 0; m < lhsGlobHit->facetStates[facetId].momentResults.size(); ++m) {
            auto &facetCounter_lhs = lhsGlobHit->facetStates[facetId].momentResults[m];
            auto &facetCounter_rhs = rhsGlobHit->facetStates[facetId].momentResults[m];

            // If one facet doesn't have any hits recorded, comparison is pointless, so just skip to next facet
            if ((facetCounter_lhs.hits.nbMCHit == 0 && facetCounter_rhs.hits.nbMCHit == 0)
            || (facetCounter_lhs.hits.nbMCHit < 40000 || facetCounter_rhs.hits.nbMCHit < 40000)) {
                // Skip facet comparison if not enough hits have been recorded for both states
                ++nbFacetSkips;
                continue;
            } else if (facetCounter_lhs.hits.nbMCHit == 0 && facetCounter_rhs.hits.nbMCHit > 0) {
                cmpFile += fmt::format("[Facet][{}][hits][{}] First state has no recorded hits for this facet\n", facetId, m);
                ++facetErrNb;
                ++nbFacetSkips;
                continue;
            } else if (facetCounter_lhs.hits.nbMCHit > 0 && facetCounter_rhs.hits.nbMCHit == 0) {
                cmpFile += fmt::format("[Facet][{}][hits][{}] Second state has no recorded hits for this facet\n", facetId, m);
                ++facetErrNb;
                ++nbFacetSkips;
                continue;
            }

            // Adjust threshold to reasonable limit for moments due to lower amount of hits
            if(m >= 0){
                auto nMC = std::min(lhsGlobHit->facetStates[facetId].momentResults[m].hits.nbMCHit, rhsGlobHit->facetStates[facetId].momentResults[m].hits.nbMCHit);
                locThreshold = std::max(locThreshold_bak, 1.0 / std::sqrt(nMC));
            }
            else{
                locThreshold = locThreshold_bak;
            }

            double scale = 1.0 / static_cast<double>(lhsGlobHit->globalStats.globalHits.nbDesorbed); // getmolpertp
            double scale_rhs = 1.0 / static_cast<double>(rhsGlobHit->globalStats.globalHits.nbDesorbed);
            double fullScale = 1.0;

            // density correction value for scale
            if (facetCounter_lhs.hits.nbMCHit > 0 || facetCounter_lhs.hits.nbDesorbed > 0) {
                if (facetCounter_lhs.hits.nbAbsEquiv > 0.0 ||
                    facetCounter_lhs.hits.nbDesorbed > 0) { //otherwise, save calculation time
                    fullScale = 1.0 - (facetCounter_lhs.hits.nbAbsEquiv + (double) facetCounter_lhs.hits.nbDesorbed) /
                                      (facetCounter_lhs.hits.nbHitEquiv + (double) facetCounter_lhs.hits.nbDesorbed) /
                                      2.0;
                }
            }

            double fullScale_rhs = 1.0;
            if (facetCounter_rhs.hits.nbMCHit > 0 || facetCounter_rhs.hits.nbDesorbed > 0) {
                if (facetCounter_rhs.hits.nbAbsEquiv > 0.0 ||
                    facetCounter_rhs.hits.nbDesorbed > 0) { //otherwise, save calculation time
                    fullScale_rhs = 1.0 -
                                    (facetCounter_rhs.hits.nbAbsEquiv + (double) facetCounter_rhs.hits.nbDesorbed) /
                                    (facetCounter_rhs.hits.nbHitEquiv + (double) facetCounter_rhs.hits.nbDesorbed) /
                                    2.0;
                }
            }

            fullScale *= scale;
            fullScale_rhs *= scale_rhs;

            scale = 1.0 / lhsGlobHit->globalStats.globalHits.nbHitEquiv;
            scale_rhs = 1.0 / rhsGlobHit->globalStats.globalHits.nbHitEquiv;
            fullScale = 1.0 /
                        (lhsGlobHit->globalStats.globalHits.nbHitEquiv + lhsGlobHit->globalStats.globalHits.nbAbsEquiv +
                         static_cast<double>(lhsGlobHit->globalStats.globalHits.nbDesorbed));
            fullScale_rhs = 1.0 /
                            (rhsGlobHit->globalStats.globalHits.nbHitEquiv + rhsGlobHit->globalStats.globalHits.nbAbsEquiv +
                             static_cast<double>(rhsGlobHit->globalStats.globalHits.nbDesorbed));
            double sumHitDes = facetCounter_lhs.hits.nbHitEquiv + static_cast<double>(facetCounter_lhs.hits.nbDesorbed);
            double sumHitDes_rhs =
                    facetCounter_rhs.hits.nbHitEquiv + static_cast<double>(facetCounter_rhs.hits.nbDesorbed);

            if (!(std::sqrt(
                    std::max(1.0, std::min(facetCounter_lhs.hits.nbHitEquiv, facetCounter_rhs.hits.nbHitEquiv))) <
                  80)) {
                double hitRatio = facetCounter_lhs.hits.nbHitEquiv * scale;
                double hitRatio_rhs = facetCounter_rhs.hits.nbHitEquiv * scale_rhs;
                if (!IsEqual(hitRatio, hitRatio_rhs, locThreshold)) {
                    cmpFile += fmt::format("[Facet][{}][facet_hits / global_hits][moment{}] has large difference: {} vs {}\n"
                                           "    {}fh / {}gh vs {}fh / {}gh\n",
                                           facetId, m,
                                           hitRatio,hitRatio_rhs,
                                           facetCounter_lhs.hits.nbHitEquiv, lhsGlobHit->globalStats.globalHits.nbHitEquiv,
                                           facetCounter_rhs.hits.nbHitEquiv, lhsGlobHit->globalStats.globalHits.nbHitEquiv);
                    ++facetErrNb;
                }
                if (!IsEqual(facetCounter_lhs.hits.sum_v_ort * scale, facetCounter_rhs.hits.sum_v_ort * scale_rhs,
                             locThreshold)) {
                    cmpFile += fmt::format("[Facet][{}][sum_v_ort][moment{}] has large difference: {} vs {}\n"
                                           "    {}fv / {}gh vs {}fv / {}gh\n",
                                           facetId, m,
                                           facetCounter_lhs.hits.sum_v_ort * scale, facetCounter_rhs.hits.sum_v_ort * scale_rhs,
                                           facetCounter_lhs.hits.sum_v_ort, lhsGlobHit->globalStats.globalHits.nbHitEquiv,
                                           facetCounter_rhs.hits.sum_v_ort, lhsGlobHit->globalStats.globalHits.nbHitEquiv);
                    ++facetErrNb;
                }
                if (!IsEqual(facetCounter_lhs.hits.sum_1_per_velocity * fullScale,
                             facetCounter_rhs.hits.sum_1_per_velocity * fullScale_rhs,
                             locThreshold * velocityThresholdFactor)) {
                    cmpFile += fmt::format("[Facet][{}][sum_1_per_velocity][moment{}] has large difference: {} (normalized: {})\n",
                                           facetId, m,
                                           std::abs(facetCounter_lhs.hits.sum_1_per_velocity * fullScale - facetCounter_rhs.hits.sum_1_per_velocity * fullScale_rhs),
                                           std::abs(facetCounter_lhs.hits.sum_1_per_velocity * fullScale - facetCounter_rhs.hits.sum_1_per_velocity * fullScale_rhs) / std::max(facetCounter_lhs.hits.sum_1_per_velocity * fullScale, facetCounter_rhs.hits.sum_1_per_velocity * fullScale_rhs));
                    ++facetErrNb;
                }
                if (!IsEqual(facetCounter_lhs.hits.sum_1_per_ort_velocity * fullScale,
                             facetCounter_rhs.hits.sum_1_per_ort_velocity * fullScale_rhs,
                             locThreshold * velocityThresholdFactor)) {
                    cmpFile += fmt::format("[Facet][{}][sum_1_per_ort_velocity][moment{}] has large difference: {} (normalized: {})\n",
                                           facetId, m,
                                           std::abs(facetCounter_lhs.hits.sum_1_per_ort_velocity * fullScale - facetCounter_rhs.hits.sum_1_per_ort_velocity * fullScale_rhs),
                                           std::abs(facetCounter_lhs.hits.sum_1_per_ort_velocity * fullScale - facetCounter_rhs.hits.sum_1_per_ort_velocity * fullScale_rhs) / std::max(facetCounter_lhs.hits.sum_1_per_ort_velocity * fullScale, facetCounter_rhs.hits.sum_1_per_ort_velocity * fullScale_rhs));
                    ++facetErrNb;
                }
            }

            if (!(std::sqrt(
                    std::max(1.0, std::min(facetCounter_lhs.hits.nbAbsEquiv, facetCounter_rhs.hits.nbAbsEquiv))) <
                  80)) {
                double absRatio = facetCounter_lhs.hits.nbAbsEquiv / static_cast<double>(facetCounter_lhs.hits.nbMCHit);
                double absRatio_rhs = facetCounter_rhs.hits.nbAbsEquiv / static_cast<double>(facetCounter_rhs.hits.nbMCHit);
                if (!IsEqual(absRatio, absRatio_rhs, locThreshold)) {
                    cmpFile += fmt::format("[Facet][{}][abs/hits ratio][moment{}] has large difference: {} vs {}\n"
                                           "    {}fa / {}fh vs {}fa / {}fh\n",
                                           facetId, m,
                                           absRatio, absRatio_rhs,
                                            facetCounter_lhs.hits.nbAbsEquiv,facetCounter_lhs.hits.nbMCHit,
                                            facetCounter_rhs.hits.nbAbsEquiv,facetCounter_rhs.hits.nbMCHit);
                    ++facetErrNb;
                }
            }

            if (!(std::sqrt(std::max((size_t) 1,
                                     std::min(facetCounter_lhs.hits.nbDesorbed, facetCounter_rhs.hits.nbDesorbed))) <
                  80)) {
                double desRatio = (double) facetCounter_lhs.hits.nbDesorbed / static_cast<double>(facetCounter_lhs.hits.nbMCHit);
                double desRatio_rhs = (double) facetCounter_rhs.hits.nbDesorbed / static_cast<double>(facetCounter_rhs.hits.nbMCHit);
                if (!IsEqual(desRatio, desRatio_rhs, locThreshold)) {
                    cmpFile += fmt::format("[Facet][{}][desRatio][moment{}] has large difference: {} (normalized: {})\n",
                                           facetId, m,
                                           std::abs(desRatio - desRatio_rhs), std::abs(desRatio - desRatio_rhs) / std::max(desRatio, desRatio_rhs));
                    ++facetErrNb;
                }
            }

            //profile
            {
                auto &prof_lhs = facetCounter_lhs.profile;
                auto &prof_rhs = facetCounter_rhs.profile;

                // use 1,4,6,4,1 as smoothing kernel to get more results
                if(1) { // gaussian 1x5 smoothing kernel
                    for (int id = 2; id < (int)(prof_lhs.size()) - 2; ++id) {
                        auto smooth_countEquiv_lhs =
                                1.0/16.0 * (prof_lhs[id-2].countEquiv + 4.0 * prof_lhs[id-1].countEquiv + 6.0 * prof_lhs[id].countEquiv
                                + 4.0 * prof_lhs[id+1].countEquiv + prof_lhs[id+2].countEquiv);
                        auto smooth_countEquiv_rhs =
                                1.0/16.0 * (prof_rhs[id-2].countEquiv + 4.0 * prof_rhs[id-1].countEquiv + 6.0 * prof_rhs[id].countEquiv
                                        + 4.0 * prof_rhs[id+1].countEquiv + prof_rhs[id+2].countEquiv);
                        if (std::max(1.0, std::min(smooth_countEquiv_lhs, smooth_countEquiv_rhs)) < 100) {
                            // Sample size not large enough
                            ++nbProfileSkips;
                            continue;
                        }
                        if (!IsEqual(smooth_countEquiv_lhs / sumHitDes, smooth_countEquiv_rhs / sumHitDes_rhs,
                                     locThreshold)) {
                            cmpFileFine += fmt::format(
                                    "[Facet][{}][Profile][Moment{}][countEquiv][{}] has large difference: {} : {} - {}\n",
                                    facetId, id, m,
                                    std::abs(smooth_countEquiv_lhs / sumHitDes -
                                                     smooth_countEquiv_rhs / sumHitDes_rhs) /
                                    (smooth_countEquiv_lhs / sumHitDes),
                                    std::abs(smooth_countEquiv_lhs / sumHitDes),
                                    (smooth_countEquiv_rhs / sumHitDes_rhs));

                            ++fineErrNb;
                        }

                        auto smooth_sum_1_per_ort_velocity_lhs =
                                1.0/16.0 * (prof_lhs[id-2].countEquiv + 4.0 * prof_lhs[id-1].countEquiv + 6.0 * prof_lhs[id].countEquiv
                                        + 4.0 * prof_lhs[id+1].countEquiv + prof_lhs[id+2].countEquiv);
                        auto smooth_sum_1_per_ort_velocity_rhs =
                                1.0/16.0 * (prof_rhs[id-2].countEquiv + 4.0 * prof_rhs[id-1].countEquiv + 6.0 * prof_rhs[id].countEquiv
                                        + 4.0 * prof_rhs[id+1].countEquiv + prof_rhs[id+2].countEquiv);

                        if (!IsEqual(smooth_sum_1_per_ort_velocity_lhs * scale,
                                     smooth_sum_1_per_ort_velocity_rhs * scale_rhs,
                                     locThreshold * velocityThresholdFactor)) {
                            cmpFileFine += fmt::format(
                                    "[Facet][{}][Profile][Ind={}][sum_1_per_ort_velocity][{}] has large rel difference: "
                                    "{} : {} - {}\n",
                                    facetId, id, m,
                                    std::abs(smooth_sum_1_per_ort_velocity_lhs * scale -
                                             smooth_sum_1_per_ort_velocity_rhs * scale_rhs) /
                                    (smooth_sum_1_per_ort_velocity_lhs * scale),
                                    std::abs(smooth_sum_1_per_ort_velocity_lhs * scale),
                                    (smooth_sum_1_per_ort_velocity_rhs * scale_rhs));

                            ++fineErrNb;
                        }
                        auto smooth_sum_v_ort_lhs =
                                1.0/16.0 * (prof_lhs[id-2].countEquiv + 4.0 * prof_lhs[id-1].countEquiv + 6.0 * prof_lhs[id].countEquiv
                                        + 4.0 * prof_lhs[id+1].countEquiv + prof_lhs[id+2].countEquiv);
                        auto smooth_sum_v_ort_rhs =
                                1.0/16.0 * (prof_rhs[id-2].countEquiv + 4.0 * prof_rhs[id-1].countEquiv + 6.0 * prof_rhs[id].countEquiv
                                        + 4.0 * prof_rhs[id+1].countEquiv + prof_rhs[id+2].countEquiv);

                        if (!IsEqual(smooth_sum_v_ort_lhs * scale, smooth_sum_v_ort_rhs * scale_rhs,
                                     locThreshold * velocityThresholdFactor)) {
                            cmpFileFine += fmt::format(
                                    "[Facet][{}][Profile][Ind={}][sum_v_ort][{}] has large rel difference: "
                                    "{} : {} - {}\n",
                                    facetId, id, m,
                                    std::abs(smooth_sum_v_ort_lhs * scale - smooth_sum_v_ort_rhs * scale_rhs) /
                                    (smooth_sum_v_ort_lhs * scale),
                                    std::abs(smooth_sum_v_ort_lhs * scale), (smooth_sum_v_ort_rhs * scale_rhs));

                            ++fineErrNb;
                        }
                    }
                }
                else {
                    for (int id = 0; id < prof_lhs.size(); ++id) {
                        if (std::sqrt(std::max(1.0, std::min(prof_lhs[id].countEquiv, prof_rhs[id].countEquiv))) < 10) {
                            // Sample size not large enough
                            ++nbProfileSkips;
                            continue;
                        }
                        if (!IsEqual(prof_lhs[id].countEquiv / sumHitDes, prof_rhs[id].countEquiv / sumHitDes_rhs,
                                     locThreshold)) {
                            cmpFileFine += fmt::format("[Facet][{}][Profile][Ind={}][countEquiv][{}] has large difference: "
                                                   "{} : {} - {}\n",
                                                   facetId, id, m,
                                                   std::abs(prof_lhs[id].countEquiv / sumHitDes - prof_rhs[id].countEquiv / sumHitDes_rhs) / (prof_lhs[id].countEquiv / sumHitDes),
                                                   std::abs(prof_lhs[id].countEquiv / sumHitDes), (prof_rhs[id].countEquiv / sumHitDes_rhs));

                            ++fineErrNb;
                        }
                        if (!IsEqual(prof_lhs[id].sum_1_per_ort_velocity * scale,
                                     prof_rhs[id].sum_1_per_ort_velocity * scale_rhs,
                                     locThreshold * velocityThresholdFactor)) {
                            cmpFileFine += fmt::format("[Facet][{}][Profile][Ind={}][sum_1_per_ort_velocity][{}] has large rel difference: "
                                                   "{} : {} - {}\n",
                                                   facetId, id, m,
                                                   std::abs(prof_lhs[id].sum_1_per_ort_velocity * scale - prof_rhs[id].sum_1_per_ort_velocity * scale_rhs) / (prof_lhs[id].sum_1_per_ort_velocity * scale),
                                                   std::abs(prof_lhs[id].sum_1_per_ort_velocity * scale), (prof_rhs[id].sum_1_per_ort_velocity * scale_rhs));

                            ++fineErrNb;
                        }
                        if (!IsEqual(prof_lhs[id].sum_v_ort * scale, prof_rhs[id].sum_v_ort * scale_rhs,
                                     locThreshold * velocityThresholdFactor)) {
                            cmpFileFine += fmt::format("[Facet][{}][Profile][Ind={}][sum_v_ort][{}] has large rel difference: "
                                                   "{} : {} - {}\n",
                                                   facetId, id, m,
                                                   std::abs(prof_lhs[id].sum_v_ort * scale - prof_rhs[id].sum_v_ort * scale_rhs) / (prof_lhs[id].sum_v_ort * scale),
                                                   std::abs(prof_lhs[id].sum_v_ort * scale), (prof_rhs[id].sum_v_ort * scale_rhs));

                            ++fineErrNb;
                        }
                    }
                }
            }

            //texture
            {
                auto &tex_lhs = facetCounter_lhs.texture;
                auto &tex_rhs = facetCounter_rhs.texture;
                int ix = 0;
                for (int iy = 0; iy < tex_lhs.size(); iy++) {
                    //for (int ix = 0; ix < texWidth_file; ix++) {
                    if (std::max(1.0, std::min(tex_lhs[iy].countEquiv, tex_rhs[iy].countEquiv)) < 640) {
                        // Sample size not large enough
                        ++nbTextureSkips;
                        continue;
                    }
                    if (!IsEqual(tex_lhs[iy].countEquiv / sumHitDes, tex_rhs[iy].countEquiv / sumHitDes_rhs,
                                 locThreshold)) {
                        cmpFileFine += fmt::format("[Facet][{}][Texture][{},{}][countEquiv][{}] has large rel difference: "
                                               "{} : {} - {}\n",
                                               facetId, ix, iy, m,
                                               std::abs( tex_lhs[iy].countEquiv / sumHitDes - tex_rhs[iy].countEquiv / sumHitDes_rhs) / (tex_lhs[iy].countEquiv / sumHitDes),
                                               std::abs(tex_lhs[iy].countEquiv / sumHitDes), (tex_rhs[iy].countEquiv / sumHitDes_rhs));

                        ++fineErrNb;
                    }
                    if (!IsEqual(tex_lhs[iy].sum_1_per_ort_velocity * fullScale,
                                 tex_rhs[iy].sum_1_per_ort_velocity * fullScale_rhs,
                                 locThreshold * velocityThresholdFactor)) {
                        cmpFileFine += fmt::format("[Facet][{}][Texture][{},{}][sum_1_per_ort_velocity][{}] has large rel difference: "
                                                   "{} : {} - {}\n",
                                                   facetId, ix, iy, m,
                                                   std::abs(tex_lhs[iy].sum_1_per_ort_velocity * fullScale - tex_rhs[iy].sum_1_per_ort_velocity * fullScale_rhs) / (tex_lhs[iy].sum_1_per_ort_velocity * fullScale),
                                                   std::abs(tex_lhs[iy].sum_1_per_ort_velocity * fullScale), (tex_rhs[iy].sum_1_per_ort_velocity * fullScale_rhs));

                        ++fineErrNb;
                    }
                    if (!IsEqual(tex_lhs[iy].sum_v_ort_per_area * scale, tex_rhs[iy].sum_v_ort_per_area * scale_rhs,
                                 locThreshold * velocityThresholdFactor)) {
                        cmpFileFine += fmt::format("[Facet][{}][Texture][{},{}][sum_v_ort_per_area][{}] has large rel difference: "
                                                   "{} : {} - {}\n",
                                                   facetId, ix, iy, m,
                                                   std::abs(tex_lhs[iy].sum_v_ort_per_area * scale - tex_rhs[iy].sum_v_ort_per_area * scale_rhs) / (tex_lhs[iy].sum_v_ort_per_area * scale),
                                                   std::abs(tex_lhs[iy].sum_v_ort_per_area * scale), (tex_rhs[iy].sum_v_ort_per_area * scale_rhs));

                        ++fineErrNb;
                    }
                    //}
                } // end for comp texture
            }

            //Directions
            {
                auto &dir_lhs = facetCounter_lhs.direction;
                auto &dir_rhs = facetCounter_rhs.direction;
                int ix = 0;
                for (int iy = 0; iy < dir_lhs.size(); iy++) {
                    //for (int ix = 0; ix < dirWidth_file; ix++) {
                    if (std::sqrt(std::max(1.0, (double) std::min(dir_lhs[iy].count, dir_rhs[iy].count))) < 80) {
                        // Sample size not large enough
                        ++nbDirSkips;
                        continue;
                    }
                    if (!IsEqual((double)dir_lhs[iy].count, (double)dir_rhs[iy].count, locThreshold)) {
                        cmpFileFine += fmt::format("[Facet][{}][dirs][{},{}][count][{}] has large difference: "
                                                   "{} : {} - {}\n",
                                                   facetId, ix, iy, m,
                                                   std::abs((int) dir_lhs[iy].count - (int) dir_rhs[iy].count),
                                                   (int) dir_lhs[iy].count, (int) dir_rhs[iy].count);

                        ++fineErrNb;
                    }
                    if (!IsEqual(dir_lhs[iy].dir.x, dir_rhs[iy].dir.x, locThreshold)) {
                        cmpFileFine += fmt::format("[Facet][{}][dirs][{},{}][dir.x][{}] has large difference: "
                                                   "{} : {} - {}\n",
                                                   facetId, ix, iy, m,
                                                   std::abs((int) dir_lhs[iy].dir.x - (int) dir_rhs[iy].dir.x),
                                                   (int) dir_lhs[iy].dir.x, (int) dir_rhs[iy].dir.x);
                        ++fineErrNb;
                    }
                    if (!IsEqual(dir_lhs[iy].dir.y, dir_rhs[iy].dir.y, locThreshold)) {
                        cmpFileFine += fmt::format("[Facet][{}][dirs][{},{}][dir.y][{}] has large difference: "
                                                   "{} : {} - {}\n",
                                                   facetId, ix, iy, m,
                                                   std::abs((int) dir_lhs[iy].dir.y - (int) dir_rhs[iy].dir.y),
                                                   (int) dir_lhs[iy].dir.y, (int) dir_rhs[iy].dir.y);
                        ++fineErrNb;
                    }
                    if (!IsEqual(dir_lhs[iy].dir.z, dir_rhs[iy].dir.z, locThreshold)) {
                        cmpFileFine += fmt::format("[Facet][{}][dirs][{},{}][dir.z][{}] has large difference: "
                                                   "{} : {} - {}\n",
                                                   facetId, ix, iy, m,
                                                   std::abs((int) dir_lhs[iy].dir.z - (int) dir_rhs[iy].dir.z),
                                                   (int) dir_lhs[iy].dir.z, (int) dir_rhs[iy].dir.z);
                        ++fineErrNb;
                    }
                    //}
                } // end for comp dir
            }

            //facet hist
            {
                auto &hist_lhs = facetCounter_lhs.histogram;
                auto &hist_rhs = facetCounter_rhs.histogram;
                for (size_t hIndex = 0; hIndex < hist_lhs.nbHitsHistogram.size(); ++hIndex) {
                    if (std::sqrt(std::max(1.0, std::min(hist_lhs.nbHitsHistogram[hIndex],
                                                         hist_rhs.nbHitsHistogram[hIndex]))) < 80) {
                        // Sample size not large enough
                        ++nbHistSkips_loc;
                        continue;
                    }
                    if (!IsEqual(hist_lhs.nbHitsHistogram[hIndex] / static_cast<double>(facetCounter_lhs.hits.nbMCHit),
                                 hist_rhs.nbHitsHistogram[hIndex] / static_cast<double>(facetCounter_rhs.hits.nbMCHit), locThreshold)) {
                        cmpFileFine += fmt::format("[Facet][{}][Hist][Bounces][Ind={}][{}] has large difference: "
                                                   "{} : {} - {}\n",
                                                   facetId, hIndex, m,
                                                   std::abs(hist_lhs.nbHitsHistogram[hIndex] / static_cast<double>(facetCounter_lhs.hits.nbMCHit) -
                                                            hist_rhs.nbHitsHistogram[hIndex] / static_cast<double>(facetCounter_rhs.hits.nbMCHit)),
                                                   hist_lhs.nbHitsHistogram[hIndex] / static_cast<double>(facetCounter_lhs.hits.nbMCHit), hist_rhs.nbHitsHistogram[hIndex] / static_cast<double>(facetCounter_rhs.hits.nbMCHit));
                        ++fineErrNb;
                    }
                }

                for (size_t hIndex = 0; hIndex < hist_lhs.distanceHistogram.size(); ++hIndex) {
                    if (std::sqrt(std::max(1.0, std::min(hist_lhs.distanceHistogram[hIndex],
                                                         hist_rhs.distanceHistogram[hIndex]))) < 80) {
                        // Sample size not large enough
                        ++nbHistSkips_loc;
                        continue;
                    }
                    if (!IsEqual(hist_lhs.distanceHistogram[hIndex] / static_cast<double>(facetCounter_lhs.hits.nbMCHit),
                                 hist_rhs.distanceHistogram[hIndex] / static_cast<double>(facetCounter_rhs.hits.nbMCHit), locThreshold)) {
                        cmpFileFine += fmt::format("[Facet][{}][Hist][Dist][Ind={}][{}] has large difference: "
                                                   "{} : {} - {}\n",
                                                   facetId, hIndex, m,
                                                   std::abs(hist_lhs.distanceHistogram[hIndex] / static_cast<double>(facetCounter_lhs.hits.nbMCHit) -
                                                            hist_rhs.distanceHistogram[hIndex] / static_cast<double>(facetCounter_rhs.hits.nbMCHit)),
                                                   hist_lhs.distanceHistogram[hIndex] / static_cast<double>(facetCounter_lhs.hits.nbMCHit), hist_rhs.distanceHistogram[hIndex] / static_cast<double>(facetCounter_rhs.hits.nbMCHit));

                        ++fineErrNb;
                    }
                }

                for (size_t hIndex = 0; hIndex < hist_lhs.timeHistogram.size(); ++hIndex) {
                    if (std::sqrt(
                            std::max(1.0, std::min(hist_lhs.timeHistogram[hIndex], hist_rhs.timeHistogram[hIndex]))) <
                        80) {
                        // Sample size not large enough
                        ++nbHistSkips_loc;
                        continue;
                    }
                    if (!IsEqual(hist_lhs.timeHistogram[hIndex] / static_cast<double>(facetCounter_lhs.hits.nbMCHit),
                                 hist_rhs.timeHistogram[hIndex] / static_cast<double>(facetCounter_rhs.hits.nbMCHit), locThreshold)) {
                        cmpFileFine += fmt::format("[Facet][{}][Hist][Time][Ind={}][{}] has large difference: "
                                                   "{} : {} - {}\n",
                                                   facetId, hIndex, m,
                                                   std::abs(hist_lhs.timeHistogram[hIndex] / static_cast<double>(facetCounter_lhs.hits.nbMCHit) -
                                                            hist_rhs.timeHistogram[hIndex] / static_cast<double>(facetCounter_rhs.hits.nbMCHit)),
                                                   hist_lhs.timeHistogram[hIndex] / static_cast<double>(facetCounter_lhs.hits.nbMCHit), hist_rhs.timeHistogram[hIndex] / static_cast<double>(facetCounter_rhs.hits.nbMCHit));

                        ++fineErrNb;
                    }
                }
            }
        }
    }

    std::string cmp_string;
    int i = 0;
    {
        std::istringstream compStream(cmpFile);
        for (; i < 32 && std::getline(compStream, cmp_string, '\n'); ++i) {
            Log::console_error("{}\n", cmp_string);
        }
    }
    if (AppSettings::verbosity >=4){
        std::istringstream compStreamFine(cmpFileFine);
        for (; i < 32 && std::getline(compStreamFine, cmp_string, '\n'); ++i) {
            Log::console_msg_master(4, "{}\n", cmp_string);
        }
    }

    if(i >= 32) {
        Log::console_error("[Warning] List of differences too long to print: {} glob, {} loc, {} fine\n", globalErrNb, facetErrNb, fineErrNb);
    }

    return std::make_tuple(globalErrNb, facetErrNb, fineErrNb);
}

/**
* \brief += operator, with simple += of underlying structures
* \param rhs reference object on the right hand
* \return address of this (lhs)
*/
FacetHistogramBuffer& FacetHistogramBuffer::operator+=(const FacetHistogramBuffer & rhs) {
    this->nbHitsHistogram += rhs.nbHitsHistogram;
    this->distanceHistogram += rhs.distanceHistogram;
    this->timeHistogram += rhs.timeHistogram;
    return *this;
}

FacetHistogramBuffer operator+(const FacetHistogramBuffer & lhs, const FacetHistogramBuffer & rhs) {
    FacetHistogramBuffer result(lhs);
    result += rhs;
    return result;
}

FacetMomentSnapshot::FacetMomentSnapshot() : hits(), histogram(){

}

/**
* \brief += operator, with simple += of underlying structures
* \param rhs reference object on the right hand
* \return address of this (lhs)
*/
FacetMomentSnapshot& FacetMomentSnapshot::operator+=(const FacetMomentSnapshot & rhs) {
    this->hits += rhs.hits;
    this->profile += rhs.profile;
    this->texture += rhs.texture;
    this->direction += rhs.direction;
    this->histogram += rhs.histogram;
    return *this;
}

size_t FacetMomentSnapshot::GetMemSize() const
{
    size_t sum = 0;
    sum += sizeof(hits);
    sum += profile.capacity() * sizeof(ProfileSlice);
    sum += texture.capacity() * sizeof(TextureCell);
    sum += direction.capacity() * sizeof(DirectionCell);
    sum += histogram.GetMemSize();
    return sum;
}

FacetMomentSnapshot operator+(const FacetMomentSnapshot & lhs, const FacetMomentSnapshot & rhs) {
    FacetMomentSnapshot result(lhs);
    result += rhs;
    return result;
}

/**
* \brief += operator, with simple += of underlying structures
* \param rhs reference object on the right hand
* \return address of this (lhs)
*/
size_t FacetState::GetMemSize() const
{
    size_t sum = 0;
    sum += recordedAngleMapPdf.capacity() * sizeof(size_t);
    for (const auto& facetSnapshot : momentResults) sum += facetSnapshot.GetMemSize();
    return sum;
}

FacetState& FacetState::operator+=(const FacetState & rhs) {
    // Check in case simulation pdf is empty (record==false) but global pdf is not (hasRecorded==true)
    if(this->recordedAngleMapPdf.size() == rhs.recordedAngleMapPdf.size())
        this->recordedAngleMapPdf += rhs.recordedAngleMapPdf;
    this->momentResults += rhs.momentResults;
    return *this;
}

FacetState operator+(const FacetState& lhs, const FacetState& rhs) {
    FacetState result(lhs);
    result += rhs;
    return result;
}

DirectionCell operator+(const DirectionCell& lhs, const DirectionCell& rhs) {
    DirectionCell result(lhs);
    result += rhs;
    return result;
}

//Constructs a lock guard for timed mutex
//Returns false if lock not successful
//Otherwise returns the unique_lock that auto-unlocks the mutex when destroyed
[[nodiscard]] std::optional<std::unique_lock<std::timed_mutex>> GetHitLock(GlobalSimuState* simStatePtr, size_t waitMs) //not shared_ptr as called from within class as well
{
    std::unique_lock<std::timed_mutex> lock(simStatePtr->simuStateMutex, std::chrono::milliseconds(waitMs));
    if (!lock.owns_lock()) {
        return std::nullopt; //couldn't acquire lock, timed out
    }
    else {
        return lock;
    }
}

void TimeDependentParameters::LoadParameterCatalog(std::vector<Parameter>& parameters) {
    char* parse_end;

    std::filesystem::path catalogPath = "parameter_catalog"; //string (POSIX) or wstring (Windows)
    if (!std::filesystem::exists(catalogPath)) return; //No param_catalog directory, don't do anything
    for (const auto& p : std::filesystem::directory_iterator(catalogPath)) {
        if (p.path().extension() == ".csv") {

            std::string csvPath = p.path().u8string();
            std::string csvName = p.path().filename().u8string();

            Parameter newParam{};
            newParam.fromCatalog = true;
            newParam.name = "[catalog] " + csvName;

            std::vector<std::vector<std::string>> table;
            try {

                FileReader file(csvPath);
                table = file.ImportCSV_string();
            }
            catch (const std::exception& e) {
                char errMsg[512];
                sprintf(errMsg, "Failed to load CSV file.\n%s", e.what());
                //GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR); //Can't display dialog window: interface not yet initialized
                std::cerr << errMsg << "\n";
                continue;
            }
            //Parse
            for (size_t i = 0; i < table.size(); i++) {
                std::vector<std::string> row = table[i];
                if (row.size() != 2) {
                    std::stringstream errMsg;
                    errMsg << p.path().filename() << " Row " << i + 1 << "has " << row.size() << " values instead of 2.";
                    //GLMessageBox::Display(errMsg.str().c_str(), "Error", GLDLG_OK, GLDLG_ICONERROR);
                    std::cerr << errMsg.str().c_str() << "\n";
                    break;
                }
                else {
                    double valueX, valueY;
                    try {
                        valueX = std::strtod(row[0].c_str(), &parse_end);
                        if (*parse_end != '\0')
                            throw Error("Malformed input in CSV file.");
                    }
                    catch (const std::exception& e) {
                        char tmp[256];
                        sprintf(tmp, "Can't parse value \"%s\" in row %zd, first column:\n%s", row[0].c_str(), i + 1, e.what());
                        //GLMessageBox::Display(tmp, "Invalid parameter definition", GLDLG_OK, GLDLG_ICONWARNING);
                        std::cerr << tmp << "\n";
                        break;
                    }
                    try {
                        valueY = std::strtod(row[1].c_str(), &parse_end);
                        if (*parse_end != '\0')
                            throw Error("Malformed input in CSV file.");
                    }
                    catch (const std::exception& e) {
                        char tmp[256];
                        sprintf(tmp, "Can't parse value \"%s\" in row %zd, second column:\n%s", row[1].c_str(), i + 1, e.what());
                        //GLMessageBox::Display(tmp, "Invalid parameter definition", GLDLG_OK, GLDLG_ICONWARNING);
                        std::cerr << tmp << "\n";
                        break;
                    }
                    newParam.AddPair(valueX, valueY, true); //insert in correct position
                }
            }
            parameters.push_back(newParam);
        }
    }
}

void TimeDependentParameters::ClearParameters(std::vector<Parameter>& parameters) {

    auto iter = parameters.begin();
    while (iter != parameters.end()) {
        //Delete non-catalog parameters
        if (!iter->fromCatalog)
            iter = parameters.erase(iter);
        else
            ++iter;
    }
}

std::vector<Parameter> TimeDependentParameters::GetCatalogParameters(const std::vector<Parameter>& parameters) {
    std::vector<Parameter> result;
    for (const auto& param : parameters) {
        if (param.fromCatalog) {
            result.push_back(param);
        }
    }
    return result;
}

/**
* \brief Function that inserts a list of new paramters at the beginning of the catalog parameters
* \param newParams vector containing new parameters to be inserted
* \return index to insert position
*/
size_t TimeDependentParameters::InsertParametersBeforeCatalog(std::vector<Parameter>& targetParameters, const std::vector<Parameter>& newParams) {
    size_t index = 0;
    for (; index != targetParameters.size() && targetParameters[index].fromCatalog == false; index++);
    targetParameters.insert(targetParameters.begin() + index, newParams.begin(),
        newParams.end()); //Insert to front (before catalog parameters)
    return index; //returns insert position
}

/**
* \brief Get ID of a parameter (if it exists) for a corresponding name
* \param name name of the parameter that shall be looked up
* \return ID corresponding to the found parameter
*/
int TimeDependentParameters::GetParamId(const std::string& name) const {
    for (int i = 0; i < (int)parameters.size(); i++) {
        if (name == parameters[i].name) return i;
    }
    //Not found
    throw Error("Parameter \"{}\" not found.",name);
}

size_t TimeDependentParameters::GetMemSize() {
    size_t sum = 0;
    for (auto& par : parameters) {
        sum += par.GetMemSize();
    }
    auto nbId = this->IDs.size();
    if (nbId > 0) sum += nbId * IDs[0].size() * sizeof(IntegratedDesorptionEntry);
    sum += sizeof(Moment) * moments.capacity();

    return sum;
}