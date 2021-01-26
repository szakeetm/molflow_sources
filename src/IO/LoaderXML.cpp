//
// Created by Pascal Baehr on 20.07.20.
//

#include <sstream>
#include <filesystem>
#include <set>
#include <Helper/MathTools.h>
#include <cereal/archives/binary.hpp>
#include <cmath>
#include "LoaderXML.h"
#include "TimeMoments.h"
#include "File.h"

constexpr size_t cdf_size = 100; //points in a cumulative distribution function

using namespace pugi;
using namespace FlowIO;

static double loadProgress = 0.0;
void setLoadProgress(double newProgress) {
    loadProgress = newProgress;
}

void reportLoadStatus(const std::string& statusString) {
    printf("[Loader at %lf3.2%%] %s", loadProgress , statusString.c_str());
}

/**
* \brief Do calculations necessary before launching simulation
* determine latest moment
* Generate integrated desorption functions
* match parameters
* Generate speed distribution functions
* Angle map
*/
void Loader::PrepareToRun(SimulationModel *model) {
    //determine latest moment
    model->wp.latestMoment = 1E-10;
    if(!model->tdParams.moments.empty())
        model->wp.latestMoment = (model->tdParams.moments.end()-1)->first + (model->tdParams.moments.end()-1)->second / 2.0;

    //Check and calculate various facet properties for time dependent simulations (CDF, ID )
    for (size_t i = 0; i < model->sh.nbFacet; i++) {
        SubprocessFacet& facet = loadFacets[i];
        // TODO: Find a solution to integrate catalog parameters
        if(facet.sh.outgassing_paramId >= (int) model->tdParams.parameters.size()){
            char tmp[256];
            sprintf(tmp, "Facet #%zd: Outgassing parameter \"%d\" isn't defined.", i + 1, facet.sh.outgassing_paramId);
            throw Error(tmp);
        }
        if(facet.sh.opacity_paramId >= (int) model->tdParams.parameters.size()){
            char tmp[256];
            sprintf(tmp, "Facet #%zd: Opacity parameter \"%d\" isn't defined.", i + 1, facet.sh.opacity_paramId);
            throw Error(tmp);
        }
        if(facet.sh.sticking_paramId >= (int) model->tdParams.parameters.size()){
            char tmp[256];
            sprintf(tmp, "Facet #%zd: Sticking parameter \"%d\" isn't defined.", i + 1, facet.sh.sticking_paramId);
            throw Error(tmp);
        }

        if (facet.sh.outgassing_paramId >= 0) { //if time-dependent desorption
            int id = GetIDId(facet.sh.outgassing_paramId);
            if (id >= 0)
                facet.sh.IDid = id; //we've already generated an ID for this temperature
            else
                facet.sh.IDid = GenerateNewID(facet.sh.outgassing_paramId, model);
        }

        // Generate speed distribution functions
        // into Loader-wide std::set<double> temperatureList;
        int id = GetCDFId(facet.sh.temperature);
        if (id >= 0)
            facet.sh.CDFid = id; //we've already generated a CDF for this temperature
        else
            facet.sh.CDFid = GenerateNewCDF(facet.sh.temperature, model->wp.gasMass);

        //Angle map
        if (facet.sh.desorbType == DES_ANGLEMAP) {
            if (!facet.sh.anglemapParams.hasRecorded) {
                char tmp[256];
                sprintf(tmp, "Facet #%zd: Uses angle map desorption but doesn't have a recorded angle map.", i + 1);
                throw Error(tmp);
            }
            if (facet.sh.anglemapParams.record) {
                char tmp[256];
                sprintf(tmp, "Facet #%zd: Can't RECORD and USE angle map desorption at the same time.", i + 1);
                throw Error(tmp);
            }
        }
    }

    CalcTotalOutgassing(model);
}

/**
* \brief Compute the outgassing of all source facet depending on the mode (file, regular, time-dependent) and set it to the global settings
*/
void Loader::CalcTotalOutgassing(SimulationModel* model) {
    // Compute the outgassing of all source facet
    double totalDesorbedMolecules = 0.0;
    double finalOutgassingRate_Pa_m3_sec = 0.0;
    double finalOutgassingRate = 0.0;

    const double latestMoment = model->wp.latestMoment;

    for (int i = 0; i < model->sh.nbFacet; i++) {
        SubprocessFacet& facet = loadFacets[i];
        if (facet.sh.desorbType != DES_NONE) { //there is a kind of desorption
            if (facet.sh.useOutgassingFile) { //outgassing file
                for (int l = 0; l < (facet.sh.outgassingMapWidth * facet.sh.outgassingMapHeight); l++) {
                    totalDesorbedMolecules += latestMoment * facet.outgassingMap[l] / (1.38E-23 * facet.sh.temperature);
                    finalOutgassingRate += facet.outgassingMap[l] / (1.38E-23 * facet.sh.temperature);
                    finalOutgassingRate_Pa_m3_sec += facet.outgassingMap[l];
                }
            } else { //regular outgassing
                if (facet.sh.outgassing_paramId == -1) { //constant outgassing
                    totalDesorbedMolecules += latestMoment * facet.sh.outgassing / (1.38E-23 * facet.sh.temperature);
                    finalOutgassingRate +=
                            facet.sh.outgassing / (1.38E-23 * facet.sh.temperature);  //Outgassing molecules/sec
                    finalOutgassingRate_Pa_m3_sec += facet.sh.outgassing;
                } else { //time-dependent outgassing
                    totalDesorbedMolecules += IDs[facet.sh.IDid].back().second / (1.38E-23 * facet.sh.temperature);
                    size_t lastIndex = model->tdParams.parameters[facet.sh.outgassing_paramId].GetSize() - 1;
                    double finalRate_mbar_l_s = model->tdParams.parameters[facet.sh.outgassing_paramId].GetY(lastIndex);
                    finalOutgassingRate +=
                            finalRate_mbar_l_s * 0.100 / (1.38E-23 * facet.sh.temperature); //0.1: mbar*l/s->Pa*m3/s
                    finalOutgassingRate_Pa_m3_sec += finalRate_mbar_l_s * 0.100;
                }
            }
        }
    }

    model->wp.totalDesorbedMolecules = totalDesorbedMolecules;
    model->wp.finalOutgassingRate_Pa_m3_sec = finalOutgassingRate_Pa_m3_sec;
    model->wp.finalOutgassingRate = finalOutgassingRate;
}

/**
* \brief Get ID (if it exists) of the Commulative Distribution Function (CFD) for a particular temperature (bin)
* \param temperature temperature for the CFD
* \return ID of the CFD
*/
int Loader::GetCDFId(double temperature) {
    if(!temperatureList.empty()) {
        auto lowerBound = std::lower_bound(temperatureList.begin(), temperatureList.end(), temperature);
        if(lowerBound == temperatureList.begin())
            return -1;
        --lowerBound; //even temperatureList.end() can be a bound

        if (std::abs(temperature - *lowerBound) > 1E-5) {
            return std::distance(temperatureList.begin(), lowerBound);
        }
    }
    return -1;
}

/**
* \brief Generate a new Commulative Distribution Function (CFD) for a particular temperature (bin)
* \param temperature for the CFD
* \return Previous size of temperatures vector, which determines new ID
*/
int Loader::GenerateNewCDF(const double temperature, const double gasMass) {
    size_t i = temperatureList.size();
    temperatureList.emplace(temperature);
    CDFs.emplace_back(Generate_CDF(temperature, gasMass, cdf_size));
    return (int) i;
}

/**
* \brief Generate cumulative distribution function (CFD) for the velocity
* \param gasTempKelvins gas temperature in Kelvin
* \param gasMassGramsPerMol molar gas mass in grams per mol
* \param size amount of points/bins of the CFD
* \return CFD as a Vector containing a pair of double values (x value = speed_bin, y value = cumulated value)
*/
std::vector<std::pair<double, double>>
Loader::Generate_CDF(double gasTempKelvins, double gasMassGramsPerMol, size_t size) {
    std::vector<std::pair<double, double>> cdf;
    cdf.reserve(size);
    constexpr double Kb = 1.38E-23;
    constexpr double R = 8.3144621;
    const double a = std::sqrt(Kb * gasTempKelvins /
                    (gasMassGramsPerMol * 1.67E-27)); //distribution a parameter. Converting molar mass to atomic mass

    //Generate cumulative distribution function
    double mostProbableSpeed = std::sqrt(2.0 * R * gasTempKelvins / (gasMassGramsPerMol / 1000.0));
    double binSize = 4.0 * mostProbableSpeed / (double) size; //distribution generated between 0 and 4*V_prob

    for (size_t i = 0; i < size; i++) {
        double x = (double) i * binSize;
        double x_square_per_2_a_square = std::pow(x, 2.0) / (2.0 * std::pow(a, 2.0));
        cdf.emplace_back(std::make_pair(x, 1.0 - std::exp(-x_square_per_2_a_square) * (x_square_per_2_a_square + 1.0)));
    }

    return cdf;
}

/**
* \brief Get ID (if it exists) of the integrated desorption (ID) function for a particular paramId
* \param paramId parameter ID
* \return Id of the integrated desorption function
*/
int Loader::GetIDId(int paramId) {
    if(!desorptionParameterIDs.empty()) {
        auto lowerBound = std::lower_bound(desorptionParameterIDs.begin(), desorptionParameterIDs.end(), paramId);
        if(lowerBound == desorptionParameterIDs.begin())
            return -1;
        --lowerBound; //even temperatureList.end() can be a bound

        if (paramId == *lowerBound) {
            return std::distance(desorptionParameterIDs.begin(), lowerBound);
        }
    }
    return -1;
}

/**
* \brief Generate a new ID (integrated desorption) for desorption parameter for time-dependent simulations
* \param paramId parameter ID
* \return Previous size of IDs vector, which determines new id in the vector
*/
int Loader::GenerateNewID(int paramId, SimulationModel* model) {
    size_t i = desorptionParameterIDs.size();
    desorptionParameterIDs.insert(paramId);
    IDs.push_back(Generate_ID(paramId,model));
    return (int) i;
}

/**
* \brief Generate integrated desorption (ID) function
* \param paramId parameter identifier
* \return ID as a Vector containing a pair of double values (x value = moment, y value = desorption value)
*/
std::vector<std::pair<double, double>> Loader::Generate_ID(int paramId, SimulationModel *model) {
    std::vector<std::pair<double, double>> newID;
    //First, let's check at which index is the latest moment
    size_t indexBeforeLastMoment;
    for (indexBeforeLastMoment = 0; indexBeforeLastMoment < model->tdParams.parameters[paramId].GetSize() &&
                                    (model->tdParams.parameters[paramId].GetX(indexBeforeLastMoment) <
                                     model->wp.latestMoment); indexBeforeLastMoment++);
    if (indexBeforeLastMoment >= model->tdParams.parameters[paramId].GetSize())
        indexBeforeLastMoment = model->tdParams.parameters[paramId].GetSize() - 1; //not found, set as last moment

//Construct integral from 0 to latest moment
//Zero
    newID.emplace_back(0.0, 0.0);

    //First moment
    newID.emplace_back(model->tdParams.parameters[paramId].GetX(0),
                                model->tdParams.parameters[paramId].GetX(0) * model->tdParams.parameters[paramId].GetY(0) *
                                0.100); //for the first moment (0.1: mbar*l/s -> Pa*m3/s)

    //Intermediate moments
    for (size_t pos = 1; pos <= indexBeforeLastMoment; pos++) {
        if (IsEqual(model->tdParams.parameters[paramId].GetY(pos),
                    model->tdParams.parameters[paramId].GetY(pos - 1))) //two equal values follow, simple integration by multiplying
            newID.emplace_back(model->tdParams.parameters[paramId].GetX(pos),
                               newID.back().second +
                                        (model->tdParams.parameters[paramId].GetX(pos) - model->tdParams.parameters[paramId].GetX(pos - 1)) *
                                        model->tdParams.parameters[paramId].GetY(pos) * 0.100);
        else { //difficult case, we'll integrate by dividing to 20 equal sections
            for (double delta = 0.05; delta < 1.0001; delta += 0.05) {
                double delta_t = model->tdParams.parameters[paramId].GetX(pos) - model->tdParams.parameters[paramId].GetX(pos - 1);
                double time = model->tdParams.parameters[paramId].GetX(pos - 1) + delta * delta_t;
                double avg_value = (model->tdParams.parameters[paramId].InterpolateY(time - 0.05 * delta_t, false) +
                                    model->tdParams.parameters[paramId].InterpolateY(time, false)) * 0.100 / 2.0;
                newID.emplace_back(time, newID.back().second + 0.05 * delta_t * avg_value);
            }
        }
    }

    //wp.latestMoment
    double valueAtLatestMoment = model->tdParams.parameters[paramId].InterpolateY(model->wp.latestMoment, false);
    if (IsEqual(valueAtLatestMoment, model->tdParams.parameters[paramId].GetY(
            indexBeforeLastMoment))) //two equal values follow, simple integration by multiplying
        newID.emplace_back(model->wp.latestMoment,
                           newID.back().second +
                                    (model->wp.latestMoment - model->tdParams.parameters[paramId].GetX(indexBeforeLastMoment)) *
                                    model->tdParams.parameters[paramId].GetY(indexBeforeLastMoment) * 0.100);
    else { //difficult case, we'll integrate by dividing two 5equal sections
        for (double delta = 0.0; delta < 1.0001; delta += 0.05) {
            double delta_t = model->wp.latestMoment - model->tdParams.parameters[paramId].GetX(indexBeforeLastMoment);
            double time = model->tdParams.parameters[paramId].GetX(indexBeforeLastMoment) + delta * delta_t;
            double avg_value = (model->tdParams.parameters[paramId].GetY(indexBeforeLastMoment) * 0.100 +
                                model->tdParams.parameters[paramId].InterpolateY(time, false) * 0.100) / 2.0;
            newID.emplace_back(time, newID.back().second + 0.05 * delta_t * avg_value);
        }
    }

    return newID;
}

int LoaderXML::LoadGeometry(const std::string inputFileName, SimulationModel *model) {
    xml_document loadXML;
    auto inputFile = inputFileName.c_str();
    xml_parse_result parseResult = loadXML.load_file(inputFile); //parse xml file directly
    xml_node rootNode = loadXML.child("SimulationEnvironment");
    if(!rootNode){
        std::cerr << "XML file seems to be of older format, please generate a new file with the GUI application!"<<std::endl;
        return 1;
    }

    xml_node geomNode = rootNode.child("Geometry");

    //Vertices
    model->sh.nbVertex = geomNode.child("Vertices").select_nodes("Vertex").size();

    model->vertices3.resize(model->sh.nbVertex); model->vertices3.shrink_to_fit();
    size_t idx = 0;
    for (xml_node vertex : geomNode.child("Vertices").children("Vertex")) {
        model->vertices3[idx].x = vertex.attribute("x").as_double();
        model->vertices3[idx].y = vertex.attribute("y").as_double();
        model->vertices3[idx].z = vertex.attribute("z").as_double();
        idx++;
    }

    //Structures
    model->sh.nbSuper = geomNode.child("Structures").select_nodes("Structure").size();

    //Parameters (needs to precede facets)
    xml_node simuParamNode = rootNode.child("MolflowSimuSettings");
    bool isMolflowFile = (simuParamNode != nullptr); //if no "MolflowSimuSettings" node, it's a Synrad file

    {
        std::vector<Distribution2D> loadedParams;
        if (isMolflowFile) {
            xml_node paramNode = simuParamNode.child("Parameters");
            for (xml_node newParameter : paramNode.children("Parameter")) {
                Distribution2D newPar;
                for (xml_node newMoment : newParameter.children("Moment")) {
                    newPar.AddPair(std::make_pair(newMoment.attribute("t").as_double(),
                                                  newMoment.attribute("value").as_double()));
                }
                loadedParams.push_back(newPar);
            }
        }
        //TODO: Load parameters from catalog explicitly?
        //work->InsertParametersBeforeCatalog(loadedParams);
        model->tdParams.parameters.insert(model->tdParams.parameters.end(),loadedParams.begin(),loadedParams.end());
    }

    //Facets , load for now via temp pointers and convert to vector afterwards
    model->sh.nbFacet = geomNode.child("Facets").select_nodes("Facet").size();
    //SubprocessFacet** loadFacets = (SubprocessFacet **)malloc(model->sh.nbFacet * sizeof(SubprocessFacet *));
    //loadFacets.reserve(model->sh.nbFacet);
    //memset(loadFacets, 0, model->sh.nbFacet * sizeof(SubprocessFacet *));
    idx = 0;
    bool ignoreSumMismatch = false;
    for (xml_node facetNode : geomNode.child("Facets").children("Facet")) {
        size_t nbIndex = facetNode.child("Indices").select_nodes("Indice").size();
        if (nbIndex < 3) {
            char errMsg[128];
            sprintf(errMsg, "Facet %zd has only %zd vertices. ", idx + 1, nbIndex);
            throw Error(errMsg);
        }

        loadFacets.emplace_back(SubprocessFacet(nbIndex));
        LoadFacet(facetNode, &loadFacets[idx], model->sh.nbVertex);
        idx++;
    }

    model->wp.gasMass = simuParamNode.child("Gas").attribute("mass").as_double();
    model->wp.halfLife = simuParamNode.child("Gas").attribute("halfLife").as_double();
    if (simuParamNode.child("Gas").attribute("enableDecay")) {
        model->wp.enableDecay = simuParamNode.child("Gas").attribute("enableDecay").as_bool();
    }
    else {
        model->wp.enableDecay = model->wp.halfLife < 1e100;
    }

    xml_node timeSettingsNode = simuParamNode.child("TimeSettings");
    model->wp.timeWindowSize = timeSettingsNode.attribute("timeWindow").as_double();
    // Default initialization
    /**/
    model->wp.useMaxwellDistribution = timeSettingsNode.attribute("useMaxwellDistr").as_bool();
    model->wp.calcConstantFlow = timeSettingsNode.attribute("calcConstFlow").as_bool();

    std::vector<UserMoment> userMoments;
    xml_node userMomentsNode = timeSettingsNode.child("UserMoments");
    for (xml_node newUserEntry : userMomentsNode.children("UserEntry")) {
        char tmpExpr[512];
        double tmpWindow = 0.0;
        strcpy(tmpExpr, newUserEntry.attribute("content").as_string());
        tmpWindow = newUserEntry.attribute("window").as_double();
        if(tmpWindow==0.0){
            tmpWindow = model->wp.timeWindowSize;
        }

        userMoments.emplace_back(tmpExpr,tmpWindow);
    }
    if(TimeMoments::ParseAndCheckUserMoments(&model->tdParams.moments, userMoments)){
        return 1;
    }

    xml_node motionNode = simuParamNode.child("Motion");
    model->wp.motionType = motionNode.attribute("type").as_int();
    if (model->wp.motionType == 1) { //fixed motion
        xml_node v = motionNode.child("VelocityVector");
        model->wp.motionVector2.x = v.attribute("vx").as_double();
        model->wp.motionVector2.y = v.attribute("vy").as_double();
        model->wp.motionVector2.z = v.attribute("vz").as_double();
    }
    else if (model->wp.motionType == 2) { //rotation
        xml_node v = motionNode.child("AxisBasePoint");
        model->wp.motionVector1.x = v.attribute("x").as_double();
        model->wp.motionVector1.y = v.attribute("y").as_double();
        model->wp.motionVector1.z = v.attribute("z").as_double();
        xml_node v2 = motionNode.child("RotationVector");
        model->wp.motionVector2.x = v2.attribute("x").as_double();
        model->wp.motionVector2.y = v2.attribute("y").as_double();
        model->wp.motionVector2.z = v2.attribute("z").as_double();
    }

    // TODO: InitializeGeometry()
    size_t nbMoments = 1; //Constant flow
#if defined(MOLFLOW)
    nbMoments += model->tdParams.moments.size();
#endif
    size_t fOffset = sizeof(GlobalHitBuffer)+nbMoments * model->wp.globalHistogramParams.GetDataSize();
     for (auto& facet : loadFacets) {
            // Main facet params
            // Current facet
            //SubprocessFacet *f = model->facets[i];
            model->CalculateFacetParams(&facet);
    
            facet.sh.hitOffset = fOffset;
            fOffset += facet.GetHitsSize(model->tdParams.moments.size());
    
    
            // Set some texture parameters
            // bool Facet::SetTexture(double width, double height, bool useMesh)
            if (facet.sh.texWidthD * facet.sh.texHeightD > 0.0000001) {
                const double ceilCutoff = 0.9999999;
                facet.sh.texWidth = (int) std::ceil(facet.sh.texWidthD *
                                            ceilCutoff); //0.9999999: cut the last few digits (convert rounding error 1.00000001 to 1, not 2)
                facet.sh.texHeight = (int) std::ceil(facet.sh.texHeightD * ceilCutoff);
            } else {
                facet.sh.texWidth = 0;
                facet.sh.texHeight = 0;
                facet.sh.texWidthD = 0.0;
                facet.sh.texHeightD = 0.0;
            }
    
        }
    PrepareToRun(model);

    model->tdParams.IDs = this->IDs;
    model->tdParams.CDFs = this->CDFs;

    return 0;
}

int LoaderXML::LoadSimulationState(std::string inputFileName, SimulationModel *model, BYTE* buffer){
    xml_document loadXML;
    xml_parse_result parseResult = loadXML.load_file(inputFileName.c_str()); //parse xml file directly
    xml_node rootNode = loadXML.child("SimulationEnvironment");

    if (!rootNode.child("MolflowResults")) return 1; //simu state not saved with file

    GlobalHitBuffer *gHits = (GlobalHitBuffer *)buffer;
    xml_node resultNode = rootNode.child("MolflowResults");
    xml_node momentsNode = resultNode.child("Moments");
    size_t nbMoments = momentsNode.select_nodes("Moment").size(); //Contains constant flow!
    size_t facetHitsSize = (nbMoments) * sizeof(FacetHitBuffer);
    size_t m = 0;
    for (xml_node newMoment : momentsNode.children("Moment")) {
        setLoadProgress((double) m / (double) nbMoments);
        if (m == 0) { //read global results
            xml_node globalNode = newMoment.child("Global");
            xml_node hitsNode = globalNode.child("Hits");
            gHits->globalHits.hit.nbMCHit = hitsNode.attribute("totalHit").as_llong();
            if (hitsNode.attribute("totalHitEquiv")) {
                gHits->globalHits.hit.nbHitEquiv = hitsNode.attribute("totalHitEquiv").as_double();
            }
            else {
                //Backward compatibility
                gHits->globalHits.hit.nbHitEquiv = static_cast<double>(gHits->globalHits.hit.nbMCHit);
            }
            gHits->globalHits.hit.nbDesorbed = hitsNode.attribute("totalDes").as_llong();
            if (hitsNode.attribute("totalAbsEquiv")) {
                gHits->globalHits.hit.nbAbsEquiv = hitsNode.attribute("totalAbsEquiv").as_double();
            }
            else {
                //Backward compatibility
                gHits->globalHits.hit.nbAbsEquiv = hitsNode.attribute("totalAbs").as_double();
            }
            if (hitsNode.attribute("totalDist_total")) { //if it's in the new format where total/partial are separated
                gHits->distTraveled_total = hitsNode.attribute("totalDist_total").as_double();
                gHits->distTraveledTotal_fullHitsOnly = hitsNode.attribute("totalDist_fullHitsOnly").as_double();
            }
            else
                gHits->distTraveled_total = gHits->distTraveledTotal_fullHitsOnly = hitsNode.attribute("totalDist").as_double();
            gHits->nbLeakTotal = hitsNode.attribute("totalLeak").as_llong();
            //work->desorptionLimit=hitsNode.attribute("maxDesorption").as_llong();

            gHits->hitCacheSize = 0;
            xml_node hitCacheNode = globalNode.child("Hit_Cache");
            for (xml_node newHit : hitCacheNode.children("Hit")) {
                if (gHits->hitCacheSize < HITCACHESIZE) {
                    gHits->hitCache[gHits->hitCacheSize].pos.x = newHit.attribute("posX").as_double();
                    gHits->hitCache[gHits->hitCacheSize].pos.y = newHit.attribute("posY").as_double();
                    gHits->hitCache[gHits->hitCacheSize].pos.z = newHit.attribute("posZ").as_double();
                    gHits->hitCache[gHits->hitCacheSize].type = newHit.attribute("type").as_int();
                    gHits->hitCacheSize++;
                }
            }

            gHits->leakCacheSize = 0;
            xml_node leakCacheNode = globalNode.child("Leak_Cache");
            for (xml_node newLeak : leakCacheNode.children("Leak")) {
                if (gHits->leakCacheSize < LEAKCACHESIZE) {
                    gHits->leakCache[gHits->leakCacheSize].pos.x = newLeak.attribute("posX").as_double();
                    gHits->leakCache[gHits->leakCacheSize].pos.y = newLeak.attribute("posY").as_double();
                    gHits->leakCache[gHits->leakCacheSize].pos.z = newLeak.attribute("posZ").as_double();
                    gHits->leakCache[gHits->leakCacheSize].dir.x = newLeak.attribute("dirX").as_double();
                    gHits->leakCache[gHits->leakCacheSize].dir.y = newLeak.attribute("dirY").as_double();
                    gHits->leakCache[gHits->leakCacheSize].dir.z = newLeak.attribute("dirZ").as_double();
                    gHits->leakCacheSize++;
                }
            }
        } //end global node

        xml_node facetResultsNode = newMoment.child("FacetResults");
        for (xml_node newFacetResult : facetResultsNode.children("Facet")) {
            int facetId = newFacetResult.attribute("id").as_int();
            SubprocessFacet& facet = loadFacets[facetId];
            xml_node facetHitNode = newFacetResult.child("Hits");
            FacetHitBuffer* facetCounter = (FacetHitBuffer *)(buffer + loadFacets[facetId].sh.hitOffset + m * sizeof(FacetHitBuffer));
            if (facetHitNode) { //If there are hit results for the current moment
                facetCounter->hit.nbMCHit = facetHitNode.attribute("nbHit").as_llong();
                if (facetHitNode.attribute("nbHitEquiv")) {
                    facetCounter->hit.nbHitEquiv = facetHitNode.attribute("nbHitEquiv").as_double();
                }
                else {
                    //Backward compatibility
                    facetCounter->hit.nbHitEquiv = static_cast<double>(facetCounter->hit.nbMCHit);
                }
                facetCounter->hit.nbDesorbed = facetHitNode.attribute("nbDes").as_llong();
                if (facetHitNode.attribute("nbAbsEquiv")) {
                    facetCounter->hit.nbAbsEquiv = facetHitNode.attribute("nbAbsEquiv").as_double();
                }
                else {
                    //Backward compatibility
                    facetCounter->hit.nbAbsEquiv = facetHitNode.attribute("nbAbs").as_double();
                }
                facetCounter->hit.sum_v_ort = facetHitNode.attribute("sum_v_ort").as_double();
                facetCounter->hit.sum_1_per_ort_velocity = facetHitNode.attribute("sum_1_per_v").as_double();
                if (facetHitNode.attribute("sum_v")) {
                    facetCounter->hit.sum_1_per_velocity = facetHitNode.attribute("sum_v").as_double();
                }
                else {
                    //Backward compatibility
                    facetCounter->hit.sum_1_per_velocity = 4.0 * Sqr(facetCounter->hit.nbHitEquiv + static_cast<double>(facetCounter->hit.nbDesorbed)) / facetCounter->hit.sum_1_per_ort_velocity;
                }

                /*if (model->displayedMoment == m) { //For immediate display in facet hits list and facet counter
                    facet.facetHitCache.hit = facetCounter->hit;
                }*/
            }
            else { //No hit information, so set to 0
                facetCounter->hit.nbMCHit =
                facetCounter->hit.nbDesorbed =
                        0;
                facetCounter->hit.sum_v_ort =
                facetCounter->hit.nbHitEquiv =
                facetCounter->hit.sum_1_per_ort_velocity =
                facetCounter->hit.sum_1_per_velocity =
                facetCounter->hit.nbAbsEquiv =
                        0.0;
            }

            //Profiles
            if (facet.sh.isProfile) {
                xml_node profileNode = newFacetResult.child("Profile");
                ProfileSlice *profilePtr = (ProfileSlice *)(buffer + facet.sh.hitOffset + facetHitsSize + m * sizeof(ProfileSlice)*PROFILE_SIZE);
                size_t id = 0;
                for (xml_node slice : profileNode.children("Slice")) {
                    if (slice.attribute("countEquiv")) {
                        profilePtr[id].countEquiv = slice.attribute("countEquiv").as_double();
                    }
                    else {
                        //Old format before low-flux
                        profilePtr[id].countEquiv = static_cast<double>(slice.attribute("count").as_llong());
                    }
                    profilePtr[id].sum_1_per_ort_velocity = slice.attribute("sum_1_per_v").as_double();
                    profilePtr[id].sum_v_ort = slice.attribute("sum_v_ort").as_double();
                    id++;
                }
            }

            //Textures
            int ix, iy;
            int profSize = (facet.sh.isProfile) ? ((int)PROFILE_SIZE * (int)sizeof(ProfileSlice)*(1 + (int)model->tdParams.moments.size())) : 0;

            if (facet.sh.texWidth * facet.sh.texHeight > 0) {
                xml_node textureNode = newFacetResult.child("Texture");
                size_t texWidth_file = textureNode.attribute("width").as_llong();
                size_t texHeight_file = textureNode.attribute("height").as_llong();

                /*if (textureNode.attribute("width").as_int() != facet.wp.texWidth ||
                    textureNode.attribute("height").as_int() != facet.wp.texHeight) {
                    std::stringstream msg;
                    msg << "Texture size mismatch on facet " << facetId + 1 << ".\nExpected: " << facet.wp.texWidth << "x" << facet.wp.texHeight << "\n"
                    << "In file: " << textureNode.attribute("width").as_int() << "x" << textureNode.attribute("height").as_int();
                    throw Error(msg.str().c_str());
                    }*/ //We'll treat texture size mismatch, see below

                TextureCell *texture = (TextureCell *)(buffer + facet.sh.hitOffset + facetHitsSize + profSize + m * facet.sh.texWidth*facet.sh.texHeight * sizeof(TextureCell));
                std::stringstream countText, sum1perText, sumvortText;
                if (textureNode.child("countEquiv")) {
                    countText << textureNode.child_value("countEquiv");
                }
                else {
                    countText << textureNode.child_value("count");
                }
                sum1perText << textureNode.child_value("sum_1_per_v");
                sumvortText << textureNode.child_value("sum_v_ort");

                for (iy = 0; iy < (Min(facet.sh.texHeight, texHeight_file)); iy++) { //MIN: If stored texture is larger, don't read extra cells
                    for (ix = 0; ix < (Min(facet.sh.texWidth, texWidth_file)); ix++) { //MIN: If stored texture is larger, don't read extra cells
                        countText >> texture[iy*facet.sh.texWidth + ix].countEquiv;
                        sum1perText >> texture[iy*facet.sh.texWidth + ix].sum_1_per_ort_velocity;
                        sumvortText >> texture[iy*facet.sh.texWidth + ix].sum_v_ort_per_area;

                    }
                    for (int ie = 0; ie < texWidth_file - facet.sh.texWidth; ie++) {//Executed if file texture is bigger than expected texture
                        //Read extra cells from file without doing anything
                        size_t dummy_ll;
                        double dummy_d;
                        countText >> dummy_ll;
                        sum1perText >> dummy_d;
                        sumvortText >> dummy_d;

                    }
                }
                for (int ie = 0; ie < texHeight_file - facet.sh.texHeight; ie++) {//Executed if file texture is bigger than expected texture
                    //Read extra cells ffrom file without doing anything
                    for (int iw = 0; iw < texWidth_file; iw++) {
                        size_t dummy_ll;
                        double dummy_d;
                        countText >> dummy_ll;
                        sum1perText >> dummy_d;
                        sumvortText >> dummy_d;
                    }

                }
            } //end texture

            if (facet.sh.countDirection) {
                xml_node dirNode = newFacetResult.child("Directions");
                if (dirNode.attribute("width").as_int() != facet.sh.texWidth ||
                    dirNode.attribute("height").as_int() != facet.sh.texHeight) {
                    std::stringstream msg;
                    msg << "Direction texture size mismatch on facet " << facetId + 1 << ".\nExpected: " << facet.sh.texWidth << "x" << facet.sh.texHeight << "\n"
                        << "In file: " << dirNode.attribute("width").as_int() << "x" << dirNode.attribute("height").as_int();
                    throw Error(msg.str().c_str());

                }
                DirectionCell *dirs = (DirectionCell *)(buffer + facet.sh.hitOffset + facetHitsSize
                                                        + profSize + (1 + (int)model->tdParams.moments.size())*facet.sh.texWidth*facet.sh.texHeight * sizeof(TextureCell)
                                                        + m * facet.sh.texWidth*facet.sh.texHeight * sizeof(DirectionCell));

                std::stringstream dirText, dirCountText;
                dirText << dirNode.child_value("vel.vectors");
                dirCountText << dirNode.child_value("count");

                for (int iy = 0; iy < facet.sh.texHeight; iy++) {
                    for (int ix = 0; ix < facet.sh.texWidth; ix++) {
                        std::string component;
                        std::getline(dirText, component, ',');
                        dirs[iy*facet.sh.texWidth + ix].dir.x = std::stod(component);
                        std::getline(dirText, component, ',');
                        dirs[iy*facet.sh.texWidth + ix].dir.y = std::stod(component);
                        dirText >> dirs[iy*facet.sh.texWidth + ix].dir.z;
                        dirCountText >> dirs[iy*facet.sh.texWidth + ix].count;
                    }
                }
            } //end directions
        } //end facetResult
        m++;
    } //end moment

    /*
    //Send angle maps //Commented out: CopyGeometryBuffer will send it after LoadXML_geom
    for (size_t i = 0; i < wp.nbFacet; i++) {
        Facet* f = facets[i];
        int profSize = (f->wp.isProfile) ? (PROFILE_SIZE * sizeof(ProfileSlice)*(1 + (int)mApp->worker.moments.size())) : 0;
        size_t *angleMap = (size_t *)((BYTE *)gHits + f->wp.hitOffset + facetHitsSize
            + profSize + (1 + (int)work->moments.size())*f->wp.texWidth*f->wp.texHeight * sizeof(TextureCell)
            + (1 + (int)work->moments.size())*f->wp.texWidth*f->wp.texHeight * sizeof(DirectionCell));
        memcpy(angleMap, f->angleMapCache, f->wp.anglemapParams.phiWidth*(f->wp.anglemapParams.thetaLowerRes + f->wp.anglemapParams.thetaHigherRes) * sizeof(size_t));
    }
    */

    xml_node minMaxNode = resultNode.child("TextureMinMax");
    /* //First write to worker->globState.globalHits, then sync it to ghits(dphit) with SendToHitBuffer()
    gHits->texture_limits[0].min.all = minMaxNode.child("With_constant_flow").child("Pressure").attribute("min").as_double();
    gHits->texture_limits[0].max.all = minMaxNode.child("With_constant_flow").child("Pressure").attribute("max").as_double();
    gHits->texture_limits[1].min.all = minMaxNode.child("With_constant_flow").child("Density").attribute("min").as_double();
    gHits->texture_limits[1].max.all = minMaxNode.child("With_constant_flow").child("Density").attribute("max").as_double();
    gHits->texture_limits[2].min.all = minMaxNode.child("With_constant_flow").child("Imp.rate").attribute("min").as_double();
    gHits->texture_limits[2].max.all = minMaxNode.child("With_constant_flow").child("Imp.rate").attribute("max").as_double();
    gHits->texture_limits[0].min.moments_only = minMaxNode.child("Moments_only").child("Pressure").attribute("min").as_double();
    gHits->texture_limits[0].max.moments_only = minMaxNode.child("Moments_only").child("Pressure").attribute("max").as_double();
    gHits->texture_limits[1].min.moments_only = minMaxNode.child("Moments_only").child("Density").attribute("min").as_double();
    gHits->texture_limits[1].max.moments_only = minMaxNode.child("Moments_only").child("Density").attribute("max").as_double();
    gHits->texture_limits[2].min.moments_only = minMaxNode.child("Moments_only").child("Imp.rate").attribute("min").as_double();
    gHits->texture_limits[2].max.moments_only = minMaxNode.child("Moments_only").child("Imp.rate").attribute("max").as_double();
    */

    gHits->texture_limits[0].min.all = minMaxNode.child("With_constant_flow").child("Pressure").attribute("min").as_double();
    gHits->texture_limits[0].max.all = minMaxNode.child("With_constant_flow").child("Pressure").attribute("max").as_double();
    gHits->texture_limits[1].min.all = minMaxNode.child("With_constant_flow").child("Density").attribute("min").as_double();
    gHits->texture_limits[1].max.all = minMaxNode.child("With_constant_flow").child("Density").attribute("max").as_double();
    gHits->texture_limits[2].min.all = minMaxNode.child("With_constant_flow").child("Imp.rate").attribute("min").as_double();
    gHits->texture_limits[2].max.all = minMaxNode.child("With_constant_flow").child("Imp.rate").attribute("max").as_double();
    gHits->texture_limits[0].min.moments_only = minMaxNode.child("Moments_only").child("Pressure").attribute("min").as_double();
    gHits->texture_limits[0].max.moments_only = minMaxNode.child("Moments_only").child("Pressure").attribute("max").as_double();
    gHits->texture_limits[1].min.moments_only = minMaxNode.child("Moments_only").child("Density").attribute("min").as_double();
    gHits->texture_limits[1].max.moments_only = minMaxNode.child("Moments_only").child("Density").attribute("max").as_double();
    gHits->texture_limits[2].min.moments_only = minMaxNode.child("Moments_only").child("Imp.rate").attribute("min").as_double();
    gHits->texture_limits[2].max.moments_only = minMaxNode.child("Moments_only").child("Imp.rate").attribute("max").as_double();


    //TODO: gHits-> = this->angleMapCache;

    return true;
}

int LoaderXML::LoadSimulationState(const std::string& inputFileName, SimulationModel *model, GlobalSimuState& globState){
    xml_document loadXML;
    xml_parse_result parseResult = loadXML.load_file(inputFileName.c_str()); //parse xml file directly
    xml_node rootNode = loadXML.child("SimulationEnvironment");

    if (!rootNode.child("MolflowResults")) return 1; //simu state not saved with file

    xml_node resultNode = rootNode.child("MolflowResults");
    xml_node momentsNode = resultNode.child("Moments");
    size_t nbMoments = momentsNode.select_nodes("Moment").size(); //Contains constant flow!
    size_t facetHitsSize = (nbMoments) * sizeof(FacetHitBuffer);
    size_t m = 0;
    for (xml_node newMoment : momentsNode.children("Moment")) {
        setLoadProgress((double) m / (double) nbMoments);
        if (m == 0) { //read global results
            xml_node globalNode = newMoment.child("Global");
            xml_node hitsNode = globalNode.child("Hits");
            globState.globalHits.globalHits.hit.nbMCHit = hitsNode.attribute("totalHit").as_llong();
            if (hitsNode.attribute("totalHitEquiv")) {
                globState.globalHits.globalHits.hit.nbHitEquiv = hitsNode.attribute("totalHitEquiv").as_double();
            }
            else {
                //Backward compatibility
                globState.globalHits.globalHits.hit.nbHitEquiv = static_cast<double>(globState.globalHits.globalHits.hit.nbMCHit);
            }
            globState.globalHits.globalHits.hit.nbDesorbed = hitsNode.attribute("totalDes").as_llong();
            if (hitsNode.attribute("totalAbsEquiv")) {
                globState.globalHits.globalHits.hit.nbAbsEquiv = hitsNode.attribute("totalAbsEquiv").as_double();
            }
            else {
                //Backward compatibility
                globState.globalHits.globalHits.hit.nbAbsEquiv = hitsNode.attribute("totalAbs").as_double();
            }
            if (hitsNode.attribute("totalDist_total")) { //if it's in the new format where total/partial are separated
                globState.globalHits.distTraveled_total = hitsNode.attribute("totalDist_total").as_double();
                globState.globalHits.distTraveledTotal_fullHitsOnly = hitsNode.attribute("totalDist_fullHitsOnly").as_double();
            }
            else
                globState.globalHits.distTraveled_total = globState.globalHits.distTraveledTotal_fullHitsOnly = hitsNode.attribute("totalDist").as_double();
            globState.globalHits.nbLeakTotal = hitsNode.attribute("totalLeak").as_llong();
            //work->desorptionLimit=hitsNode.attribute("maxDesorption").as_llong();

            globState.globalHits.hitCacheSize = 0;
            xml_node hitCacheNode = globalNode.child("Hit_Cache");
            for (xml_node newHit : hitCacheNode.children("Hit")) {
                if (globState.globalHits.hitCacheSize < HITCACHESIZE) {
                    globState.globalHits.hitCache[globState.globalHits.hitCacheSize].pos.x = newHit.attribute("posX").as_double();
                    globState.globalHits.hitCache[globState.globalHits.hitCacheSize].pos.y = newHit.attribute("posY").as_double();
                    globState.globalHits.hitCache[globState.globalHits.hitCacheSize].pos.z = newHit.attribute("posZ").as_double();
                    globState.globalHits.hitCache[globState.globalHits.hitCacheSize].type = newHit.attribute("type").as_int();
                    globState.globalHits.hitCacheSize++;
                }
            }

            globState.globalHits.leakCacheSize = 0;
            xml_node leakCacheNode = globalNode.child("Leak_Cache");
            for (xml_node newLeak : leakCacheNode.children("Leak")) {
                if (globState.globalHits.leakCacheSize < LEAKCACHESIZE) {
                    globState.globalHits.leakCache[globState.globalHits.leakCacheSize].pos.x = newLeak.attribute("posX").as_double();
                    globState.globalHits.leakCache[globState.globalHits.leakCacheSize].pos.y = newLeak.attribute("posY").as_double();
                    globState.globalHits.leakCache[globState.globalHits.leakCacheSize].pos.z = newLeak.attribute("posZ").as_double();
                    globState.globalHits.leakCache[globState.globalHits.leakCacheSize].dir.x = newLeak.attribute("dirX").as_double();
                    globState.globalHits.leakCache[globState.globalHits.leakCacheSize].dir.y = newLeak.attribute("dirY").as_double();
                    globState.globalHits.leakCache[globState.globalHits.leakCacheSize].dir.z = newLeak.attribute("dirZ").as_double();
                    globState.globalHits.leakCacheSize++;
                }
            }
        } //end global node

        xml_node facetResultsNode = newMoment.child("FacetResults");
        for (xml_node newFacetResult : facetResultsNode.children("Facet")) {
            int facetId = newFacetResult.attribute("id").as_int();
            SubprocessFacet& facet = model->facets[facetId];
            xml_node facetHitNode = newFacetResult.child("Hits");
            //FacetHitBuffer* facetCounter = (FacetHitBuffer *)(buffer + loadFacets[facetId].sh.hitOffset + m * sizeof(FacetHitBuffer));
            FacetHitBuffer* facetCounter = &globState.facetStates[facetId].momentResults[m].hits;
            if (facetHitNode) { //If there are hit results for the current moment
                facetCounter->hit.nbMCHit = facetHitNode.attribute("nbHit").as_llong();
                if (facetHitNode.attribute("nbHitEquiv")) {
                    facetCounter->hit.nbHitEquiv = facetHitNode.attribute("nbHitEquiv").as_double();
                }
                else {
                    //Backward compatibility
                    facetCounter->hit.nbHitEquiv = static_cast<double>(facetCounter->hit.nbMCHit);
                }
                facetCounter->hit.nbDesorbed = facetHitNode.attribute("nbDes").as_llong();
                if (facetHitNode.attribute("nbAbsEquiv")) {
                    facetCounter->hit.nbAbsEquiv = facetHitNode.attribute("nbAbsEquiv").as_double();
                }
                else {
                    //Backward compatibility
                    facetCounter->hit.nbAbsEquiv = facetHitNode.attribute("nbAbs").as_double();
                }
                facetCounter->hit.sum_v_ort = facetHitNode.attribute("sum_v_ort").as_double();
                facetCounter->hit.sum_1_per_ort_velocity = facetHitNode.attribute("sum_1_per_v").as_double();
                if (facetHitNode.attribute("sum_v")) {
                    facetCounter->hit.sum_1_per_velocity = facetHitNode.attribute("sum_v").as_double();
                }
                else {
                    //Backward compatibility
                    facetCounter->hit.sum_1_per_velocity = 4.0 * Sqr(facetCounter->hit.nbHitEquiv + static_cast<double>(facetCounter->hit.nbDesorbed)) / facetCounter->hit.sum_1_per_ort_velocity;
                }

                /*if (model->displayedMoment == m) { //For immediate display in facet hits list and facet counter
                    facet.facetHitCache.hit = facetCounter->hit;
                }*/
            }
            else { //No hit information, so set to 0
                facetCounter->hit.nbMCHit =
                facetCounter->hit.nbDesorbed =
                        0;
                facetCounter->hit.sum_v_ort =
                facetCounter->hit.nbHitEquiv =
                facetCounter->hit.sum_1_per_ort_velocity =
                facetCounter->hit.sum_1_per_velocity =
                facetCounter->hit.nbAbsEquiv =
                        0.0;
            }

            //Profiles
            if (facet.sh.isProfile) {
                xml_node profileNode = newFacetResult.child("Profile");
                //ProfileSlice *profilePtr = (ProfileSlice *)(buffer + facet.sh.hitOffset + facetHitsSize + m * sizeof(ProfileSlice)*PROFILE_SIZE);
                std::vector<ProfileSlice>& profilePtr = globState.facetStates[facetId].momentResults[m].profile;

                size_t id = 0;
                for (xml_node slice : profileNode.children("Slice")) {
                    if (slice.attribute("countEquiv")) {
                        profilePtr[id].countEquiv = slice.attribute("countEquiv").as_double();
                    }
                    else {
                        //Old format before low-flux
                        profilePtr[id].countEquiv = static_cast<double>(slice.attribute("count").as_llong());
                    }
                    profilePtr[id].sum_1_per_ort_velocity = slice.attribute("sum_1_per_v").as_double();
                    profilePtr[id].sum_v_ort = slice.attribute("sum_v_ort").as_double();
                    id++;
                }
            }

            //Textures
            int ix, iy;
            int profSize = (facet.sh.isProfile) ? ((int)PROFILE_SIZE * (int)sizeof(ProfileSlice)*(1 + (int)model->tdParams.moments.size())) : 0;

            if (facet.sh.texWidth * facet.sh.texHeight > 0) {
                xml_node textureNode = newFacetResult.child("Texture");
                size_t texWidth_file = textureNode.attribute("width").as_llong();
                size_t texHeight_file = textureNode.attribute("height").as_llong();

                /*if (textureNode.attribute("width").as_int() != facet.wp.texWidth ||
                    textureNode.attribute("height").as_int() != facet.wp.texHeight) {
                    std::stringstream msg;
                    msg << "Texture size mismatch on facet " << facetId + 1 << ".\nExpected: " << facet.wp.texWidth << "x" << facet.wp.texHeight << "\n"
                    << "In file: " << textureNode.attribute("width").as_int() << "x" << textureNode.attribute("height").as_int();
                    throw Error(msg.str().c_str());
                    }*/ //We'll treat texture size mismatch, see below

                //TextureCell *texture = (TextureCell *)(buffer + facet.sh.hitOffset + facetHitsSize + profSize + m * facet.sh.texWidth*facet.sh.texHeight * sizeof(TextureCell));
                std::vector<TextureCell>& texture = globState.facetStates[facetId].momentResults[m].texture;

                std::stringstream countText, sum1perText, sumvortText;
                if (textureNode.child("countEquiv")) {
                    countText << textureNode.child_value("countEquiv");
                }
                else {
                    countText << textureNode.child_value("count");
                }
                sum1perText << textureNode.child_value("sum_1_per_v");
                sumvortText << textureNode.child_value("sum_v_ort");

                for (iy = 0; iy < (Min(facet.sh.texHeight, texHeight_file)); iy++) { //MIN: If stored texture is larger, don't read extra cells
                    for (ix = 0; ix < (Min(facet.sh.texWidth, texWidth_file)); ix++) { //MIN: If stored texture is larger, don't read extra cells
                        countText >> texture[iy*facet.sh.texWidth + ix].countEquiv;
                        sum1perText >> texture[iy*facet.sh.texWidth + ix].sum_1_per_ort_velocity;
                        sumvortText >> texture[iy*facet.sh.texWidth + ix].sum_v_ort_per_area;

                    }
                    for (int ie = 0; ie < texWidth_file - facet.sh.texWidth; ie++) {//Executed if file texture is bigger than expected texture
                        //Read extra cells from file without doing anything
                        size_t dummy_ll;
                        double dummy_d;
                        countText >> dummy_ll;
                        sum1perText >> dummy_d;
                        sumvortText >> dummy_d;

                    }
                }
                for (int ie = 0; ie < texHeight_file - facet.sh.texHeight; ie++) {//Executed if file texture is bigger than expected texture
                    //Read extra cells ffrom file without doing anything
                    for (int iw = 0; iw < texWidth_file; iw++) {
                        size_t dummy_ll;
                        double dummy_d;
                        countText >> dummy_ll;
                        sum1perText >> dummy_d;
                        sumvortText >> dummy_d;
                    }

                }
            } //end texture

            if (facet.sh.countDirection) {
                xml_node dirNode = newFacetResult.child("Directions");
                if (dirNode.attribute("width").as_int() != facet.sh.texWidth ||
                    dirNode.attribute("height").as_int() != facet.sh.texHeight) {
                    std::stringstream msg;
                    msg << "Direction texture size mismatch on facet " << facetId + 1 << ".\nExpected: " << facet.sh.texWidth << "x" << facet.sh.texHeight << "\n"
                        << "In file: " << dirNode.attribute("width").as_int() << "x" << dirNode.attribute("height").as_int();
                    throw Error(msg.str().c_str());

                }
                /*DirectionCell *dirs = (DirectionCell *)(buffer + facet.sh.hitOffset + facetHitsSize
                                                        + profSize + (1 + (int)model->tdParams.moments.size())*facet.sh.texWidth*facet.sh.texHeight * sizeof(TextureCell)
                                                        + m * facet.sh.texWidth*facet.sh.texHeight * sizeof(DirectionCell));*/
                std::vector<DirectionCell>& dirs = globState.facetStates[facetId].momentResults[m].direction;

                std::stringstream dirText, dirCountText;
                dirText << dirNode.child_value("vel.vectors");
                dirCountText << dirNode.child_value("count");

                for (int iy = 0; iy < facet.sh.texHeight; iy++) {
                    for (int ix = 0; ix < facet.sh.texWidth; ix++) {
                        std::string component;
                        std::getline(dirText, component, ',');
                        dirs[iy*facet.sh.texWidth + ix].dir.x = std::stod(component);
                        std::getline(dirText, component, ',');
                        dirs[iy*facet.sh.texWidth + ix].dir.y = std::stod(component);
                        dirText >> dirs[iy*facet.sh.texWidth + ix].dir.z;
                        dirCountText >> dirs[iy*facet.sh.texWidth + ix].count;
                    }
                }
            } //end directions
        } //end facetResult
        m++;
    } //end moment

    /*
    //Send angle maps //Commented out: CopyGeometryBuffer will send it after LoadXML_geom
    for (size_t i = 0; i < wp.nbFacet; i++) {
        Facet* f = facets[i];
        int profSize = (f->wp.isProfile) ? (PROFILE_SIZE * sizeof(ProfileSlice)*(1 + (int)mApp->worker.moments.size())) : 0;
        size_t *angleMap = (size_t *)((BYTE *)gHits + f->wp.hitOffset + facetHitsSize
            + profSize + (1 + (int)work->moments.size())*f->wp.texWidth*f->wp.texHeight * sizeof(TextureCell)
            + (1 + (int)work->moments.size())*f->wp.texWidth*f->wp.texHeight * sizeof(DirectionCell));
        memcpy(angleMap, f->angleMapCache, f->wp.anglemapParams.phiWidth*(f->wp.anglemapParams.thetaLowerRes + f->wp.anglemapParams.thetaHigherRes) * sizeof(size_t));
    }
    */

    xml_node minMaxNode = resultNode.child("TextureMinMax");
    /* //First write to worker->globState.globalHits, then sync it to ghits(dphit) with SendToHitBuffer()
    globState.globalHits.texture_limits[0].min.all = minMaxNode.child("With_constant_flow").child("Pressure").attribute("min").as_double();
    globState.globalHits.texture_limits[0].max.all = minMaxNode.child("With_constant_flow").child("Pressure").attribute("max").as_double();
    globState.globalHits.texture_limits[1].min.all = minMaxNode.child("With_constant_flow").child("Density").attribute("min").as_double();
    globState.globalHits.texture_limits[1].max.all = minMaxNode.child("With_constant_flow").child("Density").attribute("max").as_double();
    globState.globalHits.texture_limits[2].min.all = minMaxNode.child("With_constant_flow").child("Imp.rate").attribute("min").as_double();
    globState.globalHits.texture_limits[2].max.all = minMaxNode.child("With_constant_flow").child("Imp.rate").attribute("max").as_double();
    globState.globalHits.texture_limits[0].min.moments_only = minMaxNode.child("Moments_only").child("Pressure").attribute("min").as_double();
    globState.globalHits.texture_limits[0].max.moments_only = minMaxNode.child("Moments_only").child("Pressure").attribute("max").as_double();
    globState.globalHits.texture_limits[1].min.moments_only = minMaxNode.child("Moments_only").child("Density").attribute("min").as_double();
    globState.globalHits.texture_limits[1].max.moments_only = minMaxNode.child("Moments_only").child("Density").attribute("max").as_double();
    globState.globalHits.texture_limits[2].min.moments_only = minMaxNode.child("Moments_only").child("Imp.rate").attribute("min").as_double();
    globState.globalHits.texture_limits[2].max.moments_only = minMaxNode.child("Moments_only").child("Imp.rate").attribute("max").as_double();
    */

    globState.globalHits.texture_limits[0].min.all = minMaxNode.child("With_constant_flow").child("Pressure").attribute("min").as_double();
    globState.globalHits.texture_limits[0].max.all = minMaxNode.child("With_constant_flow").child("Pressure").attribute("max").as_double();
    globState.globalHits.texture_limits[1].min.all = minMaxNode.child("With_constant_flow").child("Density").attribute("min").as_double();
    globState.globalHits.texture_limits[1].max.all = minMaxNode.child("With_constant_flow").child("Density").attribute("max").as_double();
    globState.globalHits.texture_limits[2].min.all = minMaxNode.child("With_constant_flow").child("Imp.rate").attribute("min").as_double();
    globState.globalHits.texture_limits[2].max.all = minMaxNode.child("With_constant_flow").child("Imp.rate").attribute("max").as_double();
    globState.globalHits.texture_limits[0].min.moments_only = minMaxNode.child("Moments_only").child("Pressure").attribute("min").as_double();
    globState.globalHits.texture_limits[0].max.moments_only = minMaxNode.child("Moments_only").child("Pressure").attribute("max").as_double();
    globState.globalHits.texture_limits[1].min.moments_only = minMaxNode.child("Moments_only").child("Density").attribute("min").as_double();
    globState.globalHits.texture_limits[1].max.moments_only = minMaxNode.child("Moments_only").child("Density").attribute("max").as_double();
    globState.globalHits.texture_limits[2].min.moments_only = minMaxNode.child("Moments_only").child("Imp.rate").attribute("min").as_double();
    globState.globalHits.texture_limits[2].max.moments_only = minMaxNode.child("Moments_only").child("Imp.rate").attribute("max").as_double();


    //TODO: globState.globalHits. = this->angleMapCache;

    return true;
}

void LoaderXML::LoadFacet(pugi::xml_node facetNode, SubprocessFacet *facet, size_t nbTotalVertices) {
    int idx = 0;
    int facetId = facetNode.attribute("id").as_int();
    for (xml_node indice : facetNode.child("Indices").children("Indice")) {
        facet->indices[idx] = indice.attribute("vertex").as_int() + 0; //+ vertexOffset;
        if (facet->indices[idx] >= nbTotalVertices) {
            char err[128];
            sprintf(err, "Facet %d refers to vertex %d which doesn't exist", facetId + 1, idx + 1);
            throw Error(err);
        }
        idx++;
    }
    facet->sh.opacity = facetNode.child("Opacity").attribute("constValue").as_double();
    facet->sh.is2sided = facetNode.child("Opacity").attribute("is2sided").as_int();
    facet->sh.superIdx = facetNode.child("Structure").attribute("inStructure").as_int();
    facet->sh.superDest = facetNode.child("Structure").attribute("linksTo").as_int();
    facet->sh.teleportDest = facetNode.child("Teleport").attribute("target").as_int();

    // Only parse Molflow files
    facet->sh.sticking = facetNode.child("Sticking").attribute("constValue").as_double();
    facet->sh.sticking_paramId = facetNode.child("Sticking").attribute("parameterId").as_int();
    facet->sh.opacity_paramId = facetNode.child("Opacity").attribute("parameterId").as_int();
    facet->sh.outgassing = facetNode.child("Outgassing").attribute("constValue").as_double();
    facet->sh.desorbType = facetNode.child("Outgassing").attribute("desType").as_int();
    facet->sh.desorbTypeN = facetNode.child("Outgassing").attribute("desExponent").as_double();
    facet->sh.outgassing_paramId = facetNode.child("Outgassing").attribute("parameterId").as_int();
    bool hasOutgassingFile = facetNode.child("Outgassing").attribute("hasOutgassingFile").as_bool();
    facet->sh.useOutgassingFile = facetNode.child("Outgassing").attribute("useOutgassingFile").as_bool();
    facet->sh.temperature = facetNode.child("Temperature").attribute("value").as_double();
    facet->sh.accomodationFactor = facetNode.child("Temperature").attribute("accFactor").as_double();
    xml_node reflNode = facetNode.child("Reflection");
    if (reflNode.attribute("diffusePart") && reflNode.attribute("specularPart")) { //New format
        facet->sh.reflection.diffusePart = reflNode.attribute("diffusePart").as_double();
        facet->sh.reflection.specularPart = reflNode.attribute("specularPart").as_double();
        if (reflNode.attribute("cosineExponent")) {
            facet->sh.reflection.cosineExponent = reflNode.attribute("cosineExponent").as_double();
        }
        else {
            facet->sh.reflection.cosineExponent = 0.0; //uniform
        }
    }
    else { //old XML format: fully diffuse / specular / uniform reflections
        int oldReflType = reflNode.attribute("type").as_int();
        if (oldReflType == REFLECTION_DIFFUSE) {
            facet->sh.reflection.diffusePart = 1.0;
            facet->sh.reflection.specularPart = 0.0;
        }
        else if (oldReflType == REFLECTION_SPECULAR) {
            facet->sh.reflection.diffusePart = 0.0;
            facet->sh.reflection.specularPart = 1.0;
        }
        else { //Uniform
            facet->sh.reflection.diffusePart = 0.0;
            facet->sh.reflection.specularPart = 0.0;
            facet->sh.reflection.cosineExponent = 0.0;
        }
    }

    if (reflNode.attribute("enableSojournTime")) {
        facet->sh.enableSojournTime = reflNode.attribute("enableSojournTime").as_bool();
        if (!reflNode.attribute("sojournFreq")) {//Backward compatibility with ver. before 2.6.25
            facet->sh.sojournFreq = 1.0 / reflNode.attribute("sojournTheta0").as_double();
            facet->sh.sojournE = 8.31 * reflNode.attribute("sojournE").as_double();
        }
        else {
            facet->sh.sojournFreq = reflNode.attribute("sojournFreq").as_double();
            facet->sh.sojournE = reflNode.attribute("sojournE").as_double();
        }
    }
    else {
        //Already set to default when calling Molflow::LoadFile()
    }
    facet->sh.isMoving = facetNode.child("Motion").attribute("isMoving").as_bool();
    xml_node recNode = facetNode.child("Recordings");
    facet->sh.profileType = recNode.child("Profile").attribute("type").as_int();
    xml_node incidentAngleNode = recNode.child("IncidentAngleMap");
    if (incidentAngleNode) {
        facet->sh.anglemapParams.record = recNode.child("IncidentAngleMap").attribute("record").as_bool();
        facet->sh.anglemapParams.phiWidth = recNode.child("IncidentAngleMap").attribute("phiWidth").as_ullong();
        facet->sh.anglemapParams.thetaLimit = recNode.child("IncidentAngleMap").attribute("thetaLimit").as_double();
        facet->sh.anglemapParams.thetaLowerRes = recNode.child("IncidentAngleMap").attribute("thetaLowerRes").as_ullong();
        facet->sh.anglemapParams.thetaHigherRes = recNode.child("IncidentAngleMap").attribute("thetaHigherRes").as_ullong();
    }
    xml_node texNode = recNode.child("Texture");
    facet->sh.texWidthD = texNode.attribute("texDimX").as_double();
    facet->sh.texHeightD = texNode.attribute("texDimY").as_double();
    facet->sh.countDes = texNode.attribute("countDes").as_bool();
    facet->sh.countAbs = texNode.attribute("countAbs").as_bool();
    facet->sh.countRefl = texNode.attribute("countRefl").as_bool();
    facet->sh.countTrans = texNode.attribute("countTrans").as_bool();
    facet->sh.countDirection = texNode.attribute("countDir").as_bool();
    facet->sh.countACD = texNode.attribute("countAC").as_bool();

    xml_node outgNode = facetNode.child("DynamicOutgassing");
    if ((hasOutgassingFile) && outgNode && outgNode.child("map")) {
        facet->sh.outgassingMapWidth = outgNode.attribute("width").as_int();
        facet->sh.outgassingMapHeight = outgNode.attribute("height").as_int();
        facet->sh.outgassingFileRatio = outgNode.attribute("ratio").as_double();
        double totalDose = outgNode.attribute("totalDose").as_double();
        facet->sh.totalOutgassing = outgNode.attribute("totalOutgassing").as_double();
        double totalFlux = outgNode.attribute("totalFlux").as_double();

        double sum = 0.0;

        std::stringstream outgText;
        outgText << outgNode.child_value("map");
        std::vector<double>(facet->sh.outgassingMapWidth*facet->sh.outgassingMapHeight).swap(facet->outgassingMap);

        for (int iy = 0; iy < facet->sh.outgassingMapHeight; iy++) {
            for (int ix = 0; ix < facet->sh.outgassingMapWidth; ix++) {
                outgText >> facet->outgassingMap[iy*facet->sh.outgassingMapWidth + ix];
                sum += facet->outgassingMap[iy*facet->sh.outgassingMapWidth + ix];
            }
        }
    }
    else hasOutgassingFile = facet->sh.useOutgassingFile = false; //if outgassing map was incorrect, don't use it

    xml_node angleMapNode = facetNode.child("IncidentAngleMap");
    if (angleMapNode && angleMapNode.child("map") && angleMapNode.attribute("angleMapThetaLimit")) {

        facet->sh.anglemapParams.phiWidth = angleMapNode.attribute("angleMapPhiWidth").as_ullong();
        facet->sh.anglemapParams.thetaLimit = angleMapNode.attribute("angleMapThetaLimit").as_double();
        facet->sh.anglemapParams.thetaLowerRes = angleMapNode.attribute("angleMapThetaLowerRes").as_ullong();
        facet->sh.anglemapParams.thetaHigherRes = angleMapNode.attribute("angleMapThetaHigherRes").as_ullong();

        std::stringstream angleText;
        angleText << angleMapNode.child_value("map");

        //size_t* angleMapCache = (size_t*)malloc(facet->sh.anglemapParams.GetDataSize());
        angleMapCache.emplace(std::make_pair(facet->globalId,std::vector<size_t>()));
        auto& angleMap = angleMapCache.at(facet->globalId);
        for (int iy = 0; iy < (facet->sh.anglemapParams.thetaLowerRes + facet->sh.anglemapParams.thetaHigherRes); iy++) {
            for (int ix = 0; ix < facet->sh.anglemapParams.phiWidth; ix++) {
                angleText >> angleMap[iy*facet->sh.anglemapParams.phiWidth + ix];
            }
        }
        facet->sh.anglemapParams.hasRecorded = true;

        //size_t mapSize = facet->sh.anglemapParams.GetMapSize();
        //facet->angleMap.pdf.resize(mapSize);
        //memcpy(facet->angleMap.pdf.data(), angleMapCache, facet->sh.anglemapParams.GetRecordedDataSize());
        //free(angleMapCache);
    }
    else {
        facet->sh.anglemapParams.hasRecorded = false; //if angle map was incorrect, don't use it
        if (facet->sh.desorbType == DES_ANGLEMAP) facet->sh.desorbType = DES_NONE;
    }

    //Update flags
    facet->sh.isProfile = (facet->sh.profileType != PROFILE_NONE);
    //wp.isOpaque = (wp.opacity != 0.0);
    facet->sh.isTextured = ((facet->sh.texWidthD * facet->sh.texHeightD) > 0);
}

void Loader::MoveFacetsToStructures(SimulationModel* model) {
    model->structures.resize(model->sh.nbSuper);
    for (size_t i = 0; i < model->sh.nbFacet; i++) { //Necessary because facets is not (yet) a vector in the interface
        if (loadFacets[i].sh.superIdx == -1) { //Facet in all structures
            for (auto &s : model->structures) {
                s.facets.emplace_back(loadFacets[i]);
                s.facets.back().globalId = i;
            }
        } else {
            model->structures[loadFacets[i].sh.superIdx].facets.emplace_back(loadFacets[i]); //Assign to structure
            model->structures[loadFacets[i].sh.superIdx].facets.back().globalId = i;
        }
    }
}

/**
* \brief Serialization function for a binary cereal archive for the worker attributes
* \return output string stream containing the result of the archiving
*/
std::ostringstream Loader::SerializeForLoader(SimulationModel* model) {
    std::ostringstream result;
    cereal::BinaryOutputArchive outputArchive(result);

    std::vector<Moment> momentIntervals;
    momentIntervals.reserve(model->tdParams.moments.size());
    for(auto& moment : model->tdParams.moments){
        momentIntervals.emplace_back(std::make_pair(moment.first - (0.5 * moment.second), moment.first + (0.5 * moment.second)));
    }
    outputArchive(
            cereal::make_nvp("wp",model->wp),
            cereal::make_nvp("ontheflyParams",model->otfParams),
            cereal::make_nvp("CDFs",model->tdParams.CDFs),
            cereal::make_nvp("IDs",model->tdParams.IDs),
            cereal::make_nvp("parameters",model->tdParams.parameters),
            //CEREAL_NVP(temperatures),
            //CEREAL_NVP(moments)
            cereal::make_nvp("moments",momentIntervals)
            //CEREAL_NVP(desorptionParameterIDs)
    ); //Worker

    //geom->SerializeForLoader(outputArchive);
    outputArchive(
            CEREAL_NVP(model->sh),
            CEREAL_NVP(model->vertices3)
    ); //Geom

    size_t fOffset = sizeof(GlobalHitBuffer) + (1 + model->tdParams.moments.size())*model->wp.globalHistogramParams.GetDataSize(); //calculating offsets for all facets for the hits dataport during the simulation

    for (size_t facIdx = 0; facIdx < model->sh.nbFacet; facIdx++) {
        SubprocessFacet& sFac = loadFacets[facIdx];
        sFac.sh.hitOffset = fOffset; //Marking the offsets for the hits, but here we don't actually send any hits.
        fOffset += sFac.GetHitsSize(model->tdParams.moments.size());

        //model->facets[i]->SerializeForLoader(outputArchive);
        // Anglemap already part of XML load
        //size_t mapSize = sFac->sh.anglemapParams.GetMapSize();
        //std::vector<size_t> angleMapVector(mapSize);
        //memcpy(angleMapVector.data(), sFac->angleMapCache, sFac->sh.anglemapParams.GetRecordedDataSize());

        std::vector<double> textIncVector;

        // Add surface elements area (reciprocal)
        if (sFac.sh.isTextured) {
            textIncVector.resize(sFac.sh.texHeight*sFac.sh.texWidth);

            double rw = sFac.sh.U.Norme() / (double)(sFac.sh.texWidthD);
            double rh = sFac.sh.V.Norme() / (double)(sFac.sh.texHeightD);
            double area = rw * rh;
            size_t add = 0;
            for (int j = 0; j < sFac.sh.texHeight; j++) {
                for (int i = 0; i < sFac.sh.texWidth; i++) {
                    if (area > 0.0) {
                        textIncVector[add] = 1.0 / area;
                    }
                    else {
                        textIncVector[add] = 0.0;
                    }
                    add++;
                }
            }
        }

        outputArchive(
                CEREAL_NVP(sFac.sh), //Contains anglemapParams
                CEREAL_NVP(sFac.indices),
                CEREAL_NVP(sFac.vertices2)
#if defined(MOLFLOW)
                , CEREAL_NVP(sFac.outgassingMap)
                , CEREAL_NVP(sFac.angleMap.pdf) //angleMapVector
                , CEREAL_NVP(textIncVector)
#endif
        );
    }

    return result;
}
/**
* \brief Serialization function for a binary cereal archive for the worker attributes
* \return output string stream containing the result of the archiving
*/
std::ostringstream Loader::SerializeResultsForLoader(GlobalSimuState* globState) {
    std::ostringstream result;
    cereal::BinaryOutputArchive outputArchive(result);

    outputArchive(
            cereal::make_nvp("GlobalHits",globState->globalHits),
            cereal::make_nvp("GlobalHistograms",globState->globalHistograms),
            cereal::make_nvp("FacetStates",globState->facetStates)
    );

    return result;
}