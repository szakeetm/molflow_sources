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

// M_PI define
#ifdef _WIN32
#define _USE_MATH_DEFINES // activate defines, e.g. M_PI_2
#endif
#include <cmath>
#include <set>
#include <sstream>

#include "Helper/MathTools.h"
#include "CDFGeneration.h"
#include "IDGeneration.h"
#include "Helper/Chronometer.h"
#include "Helper/ConsoleLogger.h"
#include "Polygon.h"
#include "MolflowSimGeom.h"
#include "MolflowSimFacet.h"
#include "IntersectAABB_shared.h" // include needed for recursive delete of AABBNODE

/**
* \brief Computes with a distribution function and a random number whether a hit on a surface is "hard" or not
* \param r instance of the handled ray for the intersection test
 * \return true on hard hit
*/
bool ParameterSurface::IsHardHit(const Ray &r) {
    const double td_opacity = dist->InterpolateY(r.time, false);
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
void  MolflowSimulationModel::BuildPrisma(double L, double R, double angle, double s, int step) {

    int nbDecade = 0;
    int nbTF = 9 * nbDecade;
    int nbTV = 4 * nbTF;

    sh.nbVertex = 2 * step + nbTV;
    std::vector<Vector3d>(sh.nbVertex).swap(vertices3);

    sh.nbFacet = step + 2 + nbTF;


    sh.nbSuper = 1;
    structures.resize(1);
    structures.begin()->strName = strdup("Prisma");

    try{
        facets.resize(sh.nbFacet, nullptr);
    }
    catch(const std::exception &e) {
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
            //facets[i + 2 + nbTF]->wp.reflection.diffusePart = 1.0; //constructor does this already
            //facets[i + 2 + nbTF]->wp.reflection.specularPart = 0.0; //constructor does this already
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

        // Volatile facet
        for (int d = 0; d < nbDecade; d++) {
            for (int i = 0; i < 9; i++) {

                double z = (double)(i + 1) * pow(10, (double)d);
                int idx = d * 36 + i * 4;

                vertices3[idx + 0].x = -R;
                vertices3[idx + 0].y = R;
                vertices3[idx + 0].z = z;
                vertices3[idx + 1].x = R;
                vertices3[idx + 1].y = R;
                vertices3[idx + 1].z = z;
                vertices3[idx + 2].x = R;
                vertices3[idx + 2].y = -R;
                vertices3[idx + 2].z = z;
                vertices3[idx + 3].x = -R;
                vertices3[idx + 3].y = -R;
                vertices3[idx + 3].z = z;

                facets[9 * d + i] = std::make_shared<MolflowSimFacet>(4);
                facets[9 * d + i]->sh.sticking = 0.0;
                facets[9 * d + i]->sh.opacity = 0.0;
                facets[9 * d + i]->sh.isVolatile = true;
                facets[9 * d + i]->indices[0] = idx + 0;
                facets[9 * d + i]->indices[1] = idx + 1;
                facets[9 * d + i]->indices[2] = idx + 2;
                facets[9 * d + i]->indices[3] = idx + 3;

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

/**
* \brief Builds ADS given certain parameters
* \param globState global simulation state for splitting techniques requiring statistical data
* \param accel_type BVH or KD tree
* \param split splitting technique corresponding to the selected AccelType
* \param bvh_width for BVH, the amount of leaves per end node
 * \return 0> for error codes, 0 when no problems
*/
int MolflowSimulationModel::BuildAccelStructure(GlobalSimuState *globState, AccelType accel_type, BVHAccel::SplitMethod split,
                                                int bvh_width) {

    initialized = false;
    Chronometer timer;
    timer.Start();

    if (!m.try_lock()) {
        return 1;
    }

#if defined(USE_OLD_BVH)
    std::vector<std::vector<SimulationFacet*>> facetPointers;
    facetPointers.resize(this->sh.nbSuper);
    for(auto& sFac : this->facets){
        // TODO: Build structures
        if (sFac->sh.superIdx == -1) { //Facet in all structures
            for (auto& fp_vec : facetPointers) {
                fp_vec.push_back(sFac.get());
            }
        }
        else {
            facetPointers[sFac->sh.superIdx].push_back(sFac.get()); //Assign to structure
        }
    }

    // Build all AABBTrees
    size_t maxDepth=0;
    for (size_t s = 0; s < this->sh.nbSuper; ++s) {
        auto& structure = this->structures[s];
        if(structure.aabbTree)
            structure.aabbTree.reset();
        AABBNODE* tree = BuildAABBTree(facetPointers[s], 0, maxDepth);
        structure.aabbTree = std::make_shared<AABBNODE>(*tree);
        //delete tree; // pointer unnecessary because of make_shared
    }

#else
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
        if (sFac->sh.opacity_paramId == -1){ //constant sticking
            sFac->sh.opacity = std::clamp(sFac->sh.opacity, 0.0, 1.0);
            //sFac->surf = simModel->GetSurface(sFac.get());
        }
        /*else {
            auto* par = &simModel->tdParams.parameters[sFac->sh.opacity_paramId];
            sFac->surf = simModel->GetParameterSurface(sFac->sh.opacity_paramId, par);
        }*/
        sFac->surf = GetSurface(sFac.get());
    }

    this->accel.clear();
    if(BVHAccel::SplitMethod::ProbSplit == split && globState && globState->initialized && globState->globalHits.globalHits.nbDesorbed > 0){
        if(globState->facetStates.size() != this->facets.size())
            return 1;
        std::vector<double> probabilities;
        probabilities.reserve(globState->facetStates.size());
        for(auto& state : globState->facetStates) {
            probabilities.emplace_back(state.momentResults[0].hits.nbHitEquiv / globState->globalHits.globalHits.nbHitEquiv);
        }
        for (size_t s = 0; s < this->sh.nbSuper; ++s) {
            if(accel_type == 1)
                this->accel.emplace_back(std::make_shared<KdTreeAccel>(primPointers[s], probabilities));
            else
                this->accel.emplace_back(std::make_shared<BVHAccel>(primPointers[s], bvh_width, BVHAccel::SplitMethod::ProbSplit, probabilities));
        }
    }
    else {
        for (size_t s = 0; s < this->sh.nbSuper; ++s) {
            if(accel_type == 1)
                this->accel.emplace_back(std::make_shared<KdTreeAccel>(primPointers[s]));
            else
                this->accel.emplace_back(std::make_shared<BVHAccel>(primPointers[s], bvh_width, split));
        }
    }
#endif // old_bvb

    timer.Stop();
    m.unlock();

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
 * \return 0> for error codes, 0 when all ok
*/
int MolflowSimulationModel::PrepareToRun() {
    if (!m.try_lock()) {
        return 1;
    }
    initialized = false;

    std::string errLog;

    //determine latest moment
    wp.latestMoment = 1E-10;
    if(!tdParams.moments.empty())
        wp.latestMoment = (tdParams.moments.end()-1)->second;
        //wp.latestMoment = (tdParams.moments.end()-1)->first + (tdParams.moments.end()-1)->second / 2.0;

    std::set<size_t> desorptionParameterIDs;
    std::vector<double> temperatureList;

    //Check and calculate various facet properties for time dependent simulations (CDF, ID )
    for (size_t i = 0; i < sh.nbFacet; i++) {
        const auto facet = facets[i];
        // TODO: Find a solution to integrate catalog parameters
        if(facet->sh.outgassing_paramId >= (int) tdParams.parameters.size()){
            char tmp[256];
            sprintf(tmp, "Facet #%zd: Outgassing parameter \"%d\" isn't defined.", i + 1, facet->sh.outgassing_paramId);
            errLog.append(tmp);
        }
        if(facet->sh.opacity_paramId >= (int) tdParams.parameters.size()){
            char tmp[256];
            sprintf(tmp, "Facet #%zd: Opacity parameter \"%d\" isn't defined.", i + 1, facet->sh.opacity_paramId);
            errLog.append(tmp);
        }
        if(facet->sh.sticking_paramId >= (int) tdParams.parameters.size()){
            char tmp[256];
            sprintf(tmp, "Facet #%zd: Sticking parameter \"%d\" isn't defined.", i + 1, facet->sh.sticking_paramId);
            errLog.append(tmp);
        }

        if (facet->sh.outgassing_paramId >= 0) { //if time-dependent desorption
            int id = IDGeneration::GetIDId(desorptionParameterIDs, facet->sh.outgassing_paramId);
            if (id >= 0)
                facet->sh.IDid = id; //we've already generated an ID for this temperature
            else {
                auto[id_new, id_vec] = IDGeneration::GenerateNewID(desorptionParameterIDs, facet->sh.outgassing_paramId, this);
                facet->sh.IDid = id_new;
                tdParams.IDs.emplace_back(std::move(id_vec));
            }
        }

        // Generate speed distribution functions
        int id = CDFGeneration::GetCDFId(temperatureList, facet->sh.temperature);
        if (id >= 0)
            facet->sh.CDFid = id; //we've already generated a CDF for this temperature
        else {
            auto[cdf_id, cdf_vec] = CDFGeneration::GenerateNewCDF(temperatureList, facet->sh.temperature, wp.gasMass);
            facet->sh.CDFid = cdf_id;
            tdParams.CDFs.emplace_back(cdf_vec);
        }
        //Angle map
        if (facet->sh.desorbType == DES_ANGLEMAP) {
            auto mfFacet = std::dynamic_pointer_cast<MolflowSimFacet>(facets[i]);
            if (mfFacet->angleMap.pdf.empty()) {
                char tmp[256];
                sprintf(tmp, "Facet #%zd: Uses angle map desorption but doesn't have a recorded angle map.", i + 1);
                errLog.append(tmp);
            }
            if (mfFacet->sh.anglemapParams.record) {
                char tmp[256];
                sprintf(tmp, "Facet #%zd: Can't RECORD and USE angle map desorption at the same time.", i + 1);
                errLog.append(tmp);
            }
        }
    }

    if(!errLog.empty()){
        m.unlock();
        return 1;
    }

    CalcTotalOutgassing();

    initialized = true;
    m.unlock();

    return 0;
}

/**
* \brief Compute the outgassing of all source facet depending on the mode (file, regular, time-dependent) and set it to the global settings
*/
void MolflowSimulationModel::CalcTotalOutgassing() {
    // Compute the outgassing of all source facet
    double totalDesorbedMolecules = 0.0;
    double finalOutgassingRate_Pa_m3_sec = 0.0;
    double finalOutgassingRate = 0.0;

    const double latestMoment = wp.latestMoment;


    for (size_t i = 0; i < facets.size(); i++) {
        const auto facet = facets[i];
        if (facet->sh.desorbType != DES_NONE) { //there is a kind of desorption
            if (facet->sh.useOutgassingFile) { //outgassing file
                const auto mfFacet = std::dynamic_pointer_cast<MolflowSimFacet>(facets[i]);
                auto &ogMap = mfFacet->ogMap;
                for (size_t l = 0; l < (ogMap.outgassingMapWidth * ogMap.outgassingMapHeight); l++) {
                    totalDesorbedMolecules +=
                            latestMoment * ogMap.outgassingMap[l] / (1.38E-23 * facet->sh.temperature);
                    finalOutgassingRate += ogMap.outgassingMap[l] / (1.38E-23 * facet->sh.temperature);
                    finalOutgassingRate_Pa_m3_sec += ogMap.outgassingMap[l];
                }
            } else { //regular outgassing
                if (facet->sh.outgassing_paramId == -1) { //constant outgassing
                    totalDesorbedMolecules +=
                            latestMoment * facet->sh.outgassing / (1.38E-23 * facet->sh.temperature);
                    finalOutgassingRate +=
                            facet->sh.outgassing / (1.38E-23 * facet->sh.temperature);  //Outgassing molecules/sec
                    finalOutgassingRate_Pa_m3_sec += facet->sh.outgassing;
                } else { //time-dependent outgassing
                    totalDesorbedMolecules +=
                            tdParams.IDs[facet->sh.IDid].back().second / (1.38E-23 * facet->sh.temperature);
                    size_t lastIndex = tdParams.parameters[facet->sh.outgassing_paramId].GetSize() - 1;
                    double finalRate_mbar_l_s = tdParams.parameters[facet->sh.outgassing_paramId].GetY(lastIndex);
                    finalOutgassingRate +=
                            finalRate_mbar_l_s * 0.100 / (1.38E-23 * facet->sh.temperature); //0.1: mbar*l/s->Pa*m3/s
                    finalOutgassingRate_Pa_m3_sec += finalRate_mbar_l_s * 0.100;
                }
            }
        }
    }

    wp.totalDesorbedMolecules = totalDesorbedMolecules;
    wp.finalOutgassingRate_Pa_m3_sec = finalOutgassingRate_Pa_m3_sec;
    wp.finalOutgassingRate = finalOutgassingRate;
}

MolflowSimulationModel::MolflowSimulationModel(MolflowSimulationModel &&o) noexcept {
    *this = std::move(o);
};

MolflowSimulationModel::MolflowSimulationModel(const MolflowSimulationModel &o) : SimulationModel(o) {
    *this = o;
};

MolflowSimulationModel::~MolflowSimulationModel() = default;

/**
* \brief Calculates the used memory used by the whole simulation model
 * \return memory size used by the whole simulation model
*/
size_t MolflowSimulationModel::size() {
    size_t modelSize = 0;
    modelSize += SimulationModel::size();
    modelSize += tdParams.GetMemSize();
    return modelSize;
}

/**
* \brief Assign operator
* \param src reference to source object
* \return address of this
*/
GlobalSimuState& GlobalSimuState::operator=(const GlobalSimuState & src) {
    //Copy all but mutex
    facetStates = src.facetStates;
    globalHistograms = src.globalHistograms;
    globalHits = src.globalHits;
    initialized = src.initialized;
    return *this;
}

/**
* \brief Assign operator
* \param src reference to source object
* \return address of this
*/
GlobalSimuState& GlobalSimuState::operator+=(const GlobalSimuState & src) {
    //Copy all but mutex
    facetStates += src.facetStates;
    globalHistograms += src.globalHistograms;
    globalHits += src.globalHits;
    return *this;
}

/**
* \brief Clears simulation state
*/
void GlobalSimuState::clear() {
    tMutex.lock();
    globalHits = GlobalHitBuffer();
    globalHistograms.clear();
    facetStates.clear();
    initialized = false;
    tMutex.unlock();
}

/**
* \brief Constructs the 'Global Hit counter structure' structure to hold all results, zero-init
* \param model Contains all related parameters
*/
void GlobalSimuState::Resize(const std::shared_ptr<SimulationModel> &model) {

    tMutex.lock();
    auto mf_model = std::dynamic_pointer_cast<MolflowSimulationModel>(model);
    size_t nbF = model->sh.nbFacet;
    size_t nbMoments = mf_model->tdParams.moments.size();
    facetStates.assign(nbF, FacetState());

    if(!model->facets.empty()) {
        for (size_t i = 0; i < nbF; i++) {
            auto sFac = model->facets[i];
            if (sFac->globalId != i) {
                fmt::print(stderr, "Facet ID mismatch! : {} / {}\n", sFac->globalId, i);
                tMutex.unlock();
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
    globalHistTemplate.Resize(model->wp.globalHistogramParams);
    globalHistograms.assign(1 + nbMoments, globalHistTemplate);
    initialized = true;
    tMutex.unlock();
}

/**
* \brief zero-init for all structures
*/
void GlobalSimuState::Reset() {
    //LockMutex(mutex);
    tMutex.lock();
    for (auto& h : globalHistograms) {
        ZEROVECTOR(h.distanceHistogram);
        ZEROVECTOR(h.nbHitsHistogram);
        ZEROVECTOR(h.timeHistogram);
    }
    memset(&globalHits, 0, sizeof(globalHits)); //Plain old data
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
    tMutex.unlock();
    //ReleaseMutex(mutex);
}

/**
 * @brief Compare function for two simulation states
 * @param lhsGlobHit first simulation state
 * @param rhsGlobHit second simulation state
 * @param globThreshold threshold for relative difference on global counters
 * @param locThreshold threshold for relative difference on facet local counters
 * @return a tuple containing the number of global, local (facet), and fine (facet profile/texture) errors
 */
std::tuple<int, int, int>
GlobalSimuState::Compare(const GlobalSimuState &lhsGlobHit, const GlobalSimuState &rhsGlobHit, double globThreshold,
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
        if(lhsGlobHit.globalHits.globalHits.nbDesorbed == 0 && rhsGlobHit.globalHits.globalHits.nbDesorbed == 0){
            cmpFile += fmt::format("[Global][desorp] Neither state has recorded desorptions\n");
            ++globalErrNb;
        }
        else if (lhsGlobHit.globalHits.globalHits.nbDesorbed == 0){
            cmpFile += fmt::format("[Global][desorp] First state has no recorded desorptions\n");
            ++globalErrNb;
        }
        else if (rhsGlobHit.globalHits.globalHits.nbDesorbed == 0){
            cmpFile += fmt::format("[Global][desorp] Second state has no recorded desorptions\n");
            ++globalErrNb;
        }

        if(globalErrNb){ // preemptive return, as global errors make local errors irrelevant
            fmt::print("{}\n", cmpFile);
            return std::make_tuple(static_cast<int>(globalErrNb), -1, -1);
        }
    }

    {
        double absRatio = lhsGlobHit.globalHits.globalHits.nbAbsEquiv / static_cast<double>(lhsGlobHit.globalHits.globalHits.nbDesorbed);
        double absRatio_rhs = rhsGlobHit.globalHits.globalHits.nbAbsEquiv / static_cast<double>(rhsGlobHit.globalHits.globalHits.nbDesorbed);
        if (!IsEqual(absRatio, absRatio_rhs, globThreshold)) {
            cmpFile += fmt::format("[Global][absRatio] has large difference: {} ({})\n",
                                   std::abs(absRatio - absRatio_rhs), std::abs(absRatio - absRatio_rhs) / std::max(absRatio, absRatio_rhs));
            ++globalErrNb;
        }
    }

    {
        double hitRatio = static_cast<double>(lhsGlobHit.globalHits.globalHits.nbMCHit) / static_cast<double>(lhsGlobHit.globalHits.globalHits.nbDesorbed);
        double hitRatio_rhs = static_cast<double>(rhsGlobHit.globalHits.globalHits.nbMCHit) / static_cast<double>(rhsGlobHit.globalHits.globalHits.nbDesorbed);
        if (!IsEqual(hitRatio, hitRatio_rhs, globThreshold)) {
            cmpFile += fmt::format("[Global][hitRatio] has large difference: {} ({}) --> "
                                   "{} / {} vs {} / {}\n",
                                   std::abs(hitRatio - hitRatio_rhs), std::abs(hitRatio - hitRatio_rhs) / std::max(hitRatio, hitRatio_rhs),
                                   lhsGlobHit.globalHits.globalHits.nbMCHit, lhsGlobHit.globalHits.globalHits.nbDesorbed,
                                   rhsGlobHit.globalHits.globalHits.nbMCHit, rhsGlobHit.globalHits.globalHits.nbDesorbed);

            ++globalErrNb;
        }
    }

    if (!IsEqual(lhsGlobHit.globalHits.globalHits.sum_v_ort, rhsGlobHit.globalHits.globalHits.sum_v_ort, globThreshold)) {
        cmpFile += fmt::format("[Global][sum_v_ort] has large difference: {}\n",
                               std::abs(lhsGlobHit.globalHits.globalHits.sum_v_ort - rhsGlobHit.globalHits.globalHits.sum_v_ort));
        ++globalErrNb;
    }
    if (!IsEqual(lhsGlobHit.globalHits.globalHits.sum_1_per_velocity, rhsGlobHit.globalHits.globalHits.sum_1_per_velocity, globThreshold)) {
        cmpFile += fmt::format("[Global][sum_1_per_velocity] has large difference: {}\n",
                               std::abs(lhsGlobHit.globalHits.globalHits.sum_1_per_velocity - rhsGlobHit.globalHits.globalHits.sum_1_per_velocity));
        ++globalErrNb;
    }
    if (!IsEqual(lhsGlobHit.globalHits.globalHits.sum_1_per_ort_velocity, rhsGlobHit.globalHits.globalHits.sum_1_per_ort_velocity, globThreshold)) {
        cmpFile += fmt::format("[Global][sum_1_per_ort_velocity] has large difference: {}\n",
                               std::abs(lhsGlobHit.globalHits.globalHits.sum_1_per_ort_velocity - rhsGlobHit.globalHits.globalHits.sum_1_per_ort_velocity));
        ++globalErrNb;
    }


    // Histogram
    {
        auto& hist_lhs = lhsGlobHit.globalHistograms;
        auto& hist_rhs = rhsGlobHit.globalHistograms;
        for(size_t tHist = 0; tHist < hist_lhs.size(); tHist++) {
            for (size_t hIndex = 0; hIndex < hist_lhs[tHist].nbHitsHistogram.size(); ++hIndex) {
                if(std::sqrt(std::max(1.0,std::min(hist_lhs[tHist].nbHitsHistogram[hIndex], hist_rhs[tHist].nbHitsHistogram[hIndex]))) < 80) {
                    // Sample size not large enough
                    continue;
                }
                if (!IsEqual(hist_lhs[tHist].nbHitsHistogram[hIndex] / static_cast<double>(lhsGlobHit.globalHits.globalHits.nbMCHit), hist_rhs[tHist].nbHitsHistogram[hIndex] / static_cast<double>(rhsGlobHit.globalHits.globalHits.nbMCHit), locThreshold)) {
                    cmpFile += fmt::format("[Global][Hist][Bounces][Ind={}] has large difference: {}\n",
                                           hIndex, std::abs(hist_lhs[tHist].nbHitsHistogram[hIndex] / static_cast<double>(lhsGlobHit.globalHits.globalHits.nbMCHit) - hist_rhs[tHist].nbHitsHistogram[hIndex] / static_cast<double>(rhsGlobHit.globalHits.globalHits.nbMCHit)));
                    ++globalErrNb;
                }
            }
            for (size_t hIndex = 0; hIndex < hist_lhs[tHist].distanceHistogram.size(); ++hIndex) {
                if(std::sqrt(std::max(1.0,std::min(hist_lhs[tHist].distanceHistogram[hIndex], hist_rhs[tHist].distanceHistogram[hIndex]))) < 80) {
                    // Sample size not large enough
                    continue;
                }
                if (!IsEqual(hist_lhs[tHist].distanceHistogram[hIndex] / static_cast<double>(lhsGlobHit.globalHits.globalHits.nbMCHit), hist_rhs[tHist].distanceHistogram[hIndex] / static_cast<double>(rhsGlobHit.globalHits.globalHits.nbMCHit), locThreshold)) {
                    cmpFile += fmt::format("[Global][Hist][Dist][Ind={}] has large difference: {}\n",
                                           hIndex, std::abs(hist_lhs[tHist].distanceHistogram[hIndex] / static_cast<double>(lhsGlobHit.globalHits.globalHits.nbMCHit) - hist_rhs[tHist].distanceHistogram[hIndex] / static_cast<double>(rhsGlobHit.globalHits.globalHits.nbMCHit)));
                    ++globalErrNb;
                }
            }
            for (size_t hIndex = 0; hIndex < hist_lhs[tHist].timeHistogram.size(); ++hIndex) {
                if(std::sqrt(std::max(1.0,std::min(hist_lhs[tHist].timeHistogram[hIndex], hist_rhs[tHist].timeHistogram[hIndex]))) < 80) {
                    // Sample size not large enough
                    continue;
                }
                if (!IsEqual(hist_lhs[tHist].timeHistogram[hIndex] / static_cast<double>(lhsGlobHit.globalHits.globalHits.nbMCHit), hist_rhs[tHist].timeHistogram[hIndex] / static_cast<double>(rhsGlobHit.globalHits.globalHits.nbMCHit), locThreshold)) {
                    cmpFile += fmt::format("[Global][Hist][Time][Ind={}] has large difference: {}\n",
                                           hIndex, std::abs(hist_lhs[tHist].timeHistogram[hIndex] / static_cast<double>(lhsGlobHit.globalHits.globalHits.nbMCHit) - hist_rhs[tHist].timeHistogram[hIndex] / static_cast<double>(rhsGlobHit.globalHits.globalHits.nbMCHit)));

                    ++globalErrNb;
                }
            }
        }
    }
    // facets

    auto locThreshold_bak = locThreshold;
    for(int facetId = 0; facetId < lhsGlobHit.facetStates.size(); ++facetId)
    {//cmp
        if(lhsGlobHit.facetStates[facetId].momentResults.size() != rhsGlobHit.facetStates[facetId].momentResults.size()){
            cmpFile += fmt::format("[Facet][{}] Different amount of moments for each state: {} vs {}\n", facetId, lhsGlobHit.facetStates[facetId].momentResults.size(), rhsGlobHit.facetStates[facetId].momentResults.size());
            ++facetErrNb;
            continue;
        }

        for(int m = 0; m < lhsGlobHit.facetStates[facetId].momentResults.size(); ++m) {
            auto &facetCounter_lhs = lhsGlobHit.facetStates[facetId].momentResults[m];
            auto &facetCounter_rhs = rhsGlobHit.facetStates[facetId].momentResults[m];

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
                auto nMC = std::min(lhsGlobHit.facetStates[facetId].momentResults[m].hits.nbMCHit, rhsGlobHit.facetStates[facetId].momentResults[m].hits.nbMCHit);
                locThreshold = std::max(locThreshold_bak, 1.0 / std::sqrt(nMC));
            }
            else{
                locThreshold = locThreshold_bak;
            }

            double scale = 1.0 / static_cast<double>(lhsGlobHit.globalHits.globalHits.nbDesorbed); // getmolpertp
            double scale_rhs = 1.0 / static_cast<double>(rhsGlobHit.globalHits.globalHits.nbDesorbed);
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

            scale = 1.0 / lhsGlobHit.globalHits.globalHits.nbHitEquiv;
            scale_rhs = 1.0 / rhsGlobHit.globalHits.globalHits.nbHitEquiv;
            fullScale = 1.0 /
                        (lhsGlobHit.globalHits.globalHits.nbHitEquiv + lhsGlobHit.globalHits.globalHits.nbAbsEquiv +
                         static_cast<double>(lhsGlobHit.globalHits.globalHits.nbDesorbed));
            fullScale_rhs = 1.0 /
                            (rhsGlobHit.globalHits.globalHits.nbHitEquiv + rhsGlobHit.globalHits.globalHits.nbAbsEquiv +
                             static_cast<double>(rhsGlobHit.globalHits.globalHits.nbDesorbed));
            double sumHitDes = facetCounter_lhs.hits.nbHitEquiv + static_cast<double>(facetCounter_lhs.hits.nbDesorbed);
            double sumHitDes_rhs =
                    facetCounter_rhs.hits.nbHitEquiv + static_cast<double>(facetCounter_rhs.hits.nbDesorbed);

            if (!(std::sqrt(
                    std::max(1.0, std::min(facetCounter_lhs.hits.nbHitEquiv, facetCounter_rhs.hits.nbHitEquiv))) <
                  80)) {
                double hitRatio = facetCounter_lhs.hits.nbHitEquiv * scale;
                double hitRatio_rhs = facetCounter_rhs.hits.nbHitEquiv * scale_rhs;
                if (!IsEqual(hitRatio, hitRatio_rhs, locThreshold)) {
                    cmpFile += fmt::format("[Facet][{}][hitRatio][{}] has large difference: "
                                           "{} (normalized: {}) --> "
                                           "{} / {} vs {} / {}\n",
                                           facetId, m,
                                           std::abs(hitRatio - hitRatio_rhs), std::abs(hitRatio - hitRatio_rhs) / std::max(hitRatio, hitRatio_rhs),
                                           facetCounter_lhs.hits.nbHitEquiv, lhsGlobHit.globalHits.globalHits.nbHitEquiv,
                                           facetCounter_rhs.hits.nbHitEquiv, lhsGlobHit.globalHits.globalHits.nbHitEquiv);
                    ++facetErrNb;
                }
                if (!IsEqual(facetCounter_lhs.hits.sum_v_ort * scale, facetCounter_rhs.hits.sum_v_ort * scale_rhs,
                             locThreshold)) {
                    cmpFile += fmt::format("[Facet][{}][sum_v_ort][{}] has large difference: "
                                           "{} (normalized: {}) --> "
                                           "{} / {} vs {} / {}\n",
                                           facetId, m,
                                           std::abs(facetCounter_lhs.hits.sum_v_ort * scale - facetCounter_rhs.hits.sum_v_ort * scale_rhs), std::abs(facetCounter_lhs.hits.sum_v_ort * scale - facetCounter_rhs.hits.sum_v_ort * scale_rhs) / std::max(facetCounter_lhs.hits.sum_v_ort * scale, facetCounter_rhs.hits.sum_v_ort * scale_rhs),
                                           facetCounter_lhs.hits.sum_v_ort, lhsGlobHit.globalHits.globalHits.nbHitEquiv,
                                           facetCounter_rhs.hits.sum_v_ort, lhsGlobHit.globalHits.globalHits.nbHitEquiv);
                    ++facetErrNb;
                }
                if (!IsEqual(facetCounter_lhs.hits.sum_1_per_velocity * fullScale,
                             facetCounter_rhs.hits.sum_1_per_velocity * fullScale_rhs,
                             locThreshold * velocityThresholdFactor)) {
                    cmpFile += fmt::format("[Facet][{}][sum_1_per_velocity][{}] has large difference: "
                                           "{} (normalized: {})\n",
                                           facetId, m,
                                           std::abs(facetCounter_lhs.hits.sum_1_per_velocity * fullScale - facetCounter_rhs.hits.sum_1_per_velocity * fullScale_rhs),
                                           std::abs(facetCounter_lhs.hits.sum_1_per_velocity * fullScale - facetCounter_rhs.hits.sum_1_per_velocity * fullScale_rhs) / std::max(facetCounter_lhs.hits.sum_1_per_velocity * fullScale, facetCounter_rhs.hits.sum_1_per_velocity * fullScale_rhs));
                    ++facetErrNb;
                }
                if (!IsEqual(facetCounter_lhs.hits.sum_1_per_ort_velocity * fullScale,
                             facetCounter_rhs.hits.sum_1_per_ort_velocity * fullScale_rhs,
                             locThreshold * velocityThresholdFactor)) {
                    cmpFile += fmt::format("[Facet][{}][sum_1_per_ort_velocity][{}] has large difference: "
                                           "{} (normalized: {})\n",
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
                    cmpFile += fmt::format("[Facet][{}][absRatio][{}] has large difference: "
                                           "{} (normalized: {})\n",
                                           facetId, m,
                                           std::abs(absRatio - absRatio_rhs), std::abs(absRatio - absRatio_rhs) / std::max(absRatio, absRatio_rhs));
                    ++facetErrNb;
                }
            }

            if (!(std::sqrt(std::max((size_t) 1,
                                     std::min(facetCounter_lhs.hits.nbDesorbed, facetCounter_rhs.hits.nbDesorbed))) <
                  80)) {
                double desRatio = (double) facetCounter_lhs.hits.nbDesorbed / static_cast<double>(facetCounter_lhs.hits.nbMCHit);
                double desRatio_rhs = (double) facetCounter_rhs.hits.nbDesorbed / static_cast<double>(facetCounter_rhs.hits.nbMCHit);
                if (!IsEqual(desRatio, desRatio_rhs, locThreshold)) {
                    cmpFile += fmt::format("[Facet][{}][desRatio][{}] has large difference: "
                                           "{} (normalized: {})\n",
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
                                    "[Facet][{}][Profile][Ind={}][countEquiv][{}] has large difference: "
                                    "{} : {} - {}\n",
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
    {
        std::istringstream compStreamFine(cmpFileFine);
        for (; i < 32 && std::getline(compStreamFine, cmp_string, '\n'); ++i) {
            Log::console_msg_master(4, "{}\n", cmp_string);
        }
    }

    if(i >= 32) {
        Log::console_error("[Warning] List of differences too long: Total = {}\n", globalErrNb + facetErrNb + fineErrNb);
    }

    return std::make_tuple(globalErrNb, facetErrNb, fineErrNb);
}

/**
* \brief += operator, with simple += of underlying structures
* \param rhs reference object on the right hand
* \return address of this (lhs)
*/
FacetHistogramBuffer& FacetHistogramBuffer::operator+=(const FacetHistogramBuffer & rhs) {
    // if (model.wp.globalHistogramParams.recordBounce)
    this->nbHitsHistogram += rhs.nbHitsHistogram;
    this->distanceHistogram += rhs.distanceHistogram;
    this->timeHistogram += rhs.timeHistogram;
    return *this;
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

/**
* \brief + operator, simply calls implemented +=
* \param rhs reference object on the right hand
* \return address of this (lhs)
*/
FacetMomentSnapshot& FacetMomentSnapshot::operator+(const FacetMomentSnapshot & rhs) {
    *this += rhs;
    return *this;
}

/**
* \brief += operator, with simple += of underlying structures
* \param rhs reference object on the right hand
* \return address of this (lhs)
*/
FacetState& FacetState::operator+=(const FacetState & rhs) {
    // Check in case simulation pdf is empty (record==false) but global pdf is not (hasRecorded==true)
    if(this->recordedAngleMapPdf.size() == rhs.recordedAngleMapPdf.size())
        this->recordedAngleMapPdf += rhs.recordedAngleMapPdf;
    this->momentResults += rhs.momentResults;
    return *this;
}