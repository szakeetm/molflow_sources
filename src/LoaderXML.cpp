//
// Created by Pascal Baehr on 20.07.20.
//

#include <sstream>
#include <filesystem>
#include <set>
#include <GLApp/MathTools.h>
#include "LoaderXML.h"
#include "PugiXML/pugixml.hpp"
#include "TimeMoments.h"
#include "File.h"

constexpr size_t cdf_size = 100; //points in a cumulative distribution function

using namespace pugi;
using namespace MFLoad;

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
        SubprocessFacet *f = model->facets[i];
        // TODO: Find a solution to integrate catalog parameters
        if(f->sh.outgassing_paramId >= model->tdParams.parameters.size()){
            char tmp[256];
            sprintf(tmp, "Facet #%zd: Outgassing parameter \"%d\" isn't defined.", i + 1, f->sh.outgassing_paramId);
            throw Error(tmp);
        }
        if(f->sh.opacity_paramId >= model->tdParams.parameters.size()){
            char tmp[256];
            sprintf(tmp, "Facet #%zd: Opacity parameter \"%d\" isn't defined.", i + 1, f->sh.opacity_paramId);
            throw Error(tmp);
        }
        if(f->sh.sticking_paramId >= model->tdParams.parameters.size()){
            char tmp[256];
            sprintf(tmp, "Facet #%zd: Sticking parameter \"%d\" isn't defined.", i + 1, f->sh.sticking_paramId);
            throw Error(tmp);
        }

        if (f->sh.outgassing_paramId >= 0) { //if time-dependent desorption
            int id = GetIDId(f->sh.outgassing_paramId);
            if (id >= 0)
                f->sh.IDid = id; //we've already generated an ID for this temperature
            else
                f->sh.IDid = GenerateNewID(f->sh.outgassing_paramId);
        }

        // Generate speed distribution functions
        // into Loader-wide std::set<double> temperatureList;
        int id = GetCDFId(f->sh.temperature);
        if (id >= 0)
            f->sh.CDFid = id; //we've already generated a CDF for this temperature
        else
            f->sh.CDFid = GenerateNewCDF(f->sh.temperature, model->wp.gasMass);

        //Angle map
        if (f->sh.desorbType == DES_ANGLEMAP) {
            if (!f->sh.anglemapParams.hasRecorded) {
                char tmp[256];
                sprintf(tmp, "Facet #%zd: Uses angle map desorption but doesn't have a recorded angle map.", i + 1);
                throw Error(tmp);
            }
            if (f->sh.anglemapParams.record) {
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
        SubprocessFacet *f = model->facets[i];
        if (f->sh.desorbType != DES_NONE) { //there is a kind of desorption
            if (f->sh.useOutgassingFile) { //outgassing file
                for (int l = 0; l < (f->sh.outgassingMapWidth * f->sh.outgassingMapHeight); l++) {
                    totalDesorbedMolecules += latestMoment * f->outgassingMap[l] / (1.38E-23 * f->sh.temperature);
                    finalOutgassingRate += f->outgassingMap[l] / (1.38E-23 * f->sh.temperature);
                    finalOutgassingRate_Pa_m3_sec += f->outgassingMap[l];
                }
            } else { //regular outgassing
                if (f->sh.outgassing_paramId == -1) { //constant outgassing
                    totalDesorbedMolecules += latestMoment * f->sh.outgassing / (1.38E-23 * f->sh.temperature);
                    finalOutgassingRate +=
                            f->sh.outgassing / (1.38E-23 * f->sh.temperature);  //Outgassing molecules/sec
                    finalOutgassingRate_Pa_m3_sec += f->sh.outgassing;
                } else { //time-dependent outgassing
                    totalDesorbedMolecules += IDs[f->sh.IDid].back().second / (1.38E-23 * f->sh.temperature);
                    size_t lastIndex = model->tdParams.parameters[f->sh.outgassing_paramId].GetSize() - 1;
                    double finalRate_mbar_l_s = model->tdParams.parameters[f->sh.outgassing_paramId].GetY(lastIndex);
                    finalOutgassingRate +=
                            finalRate_mbar_l_s * 0.100 / (1.38E-23 * f->sh.temperature); //0.1: mbar*l/s->Pa*m3/s
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
    int i;
    for (i = 0; i < (int) desorptionParameterIDs.size() &&
                (paramId != desorptionParameterIDs[i]); i++); //check if we already had this parameter Id
    if (i >= (int) desorptionParameterIDs.size()) i = -1; //not found
    return i;

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
int Loader::GenerateNewID(int paramId) {
    size_t i = desorptionParameterIDs.size();
    desorptionParameterIDs.insert(paramId);
    IDs.push_back(Generate_ID(paramId));
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

void LoaderXML::LoadGeometry(const std::string inputFileName, SimulationModel *model) {
    xml_document loadXML;
    xml_parse_result parseResult = loadXML.load_file(inputFileName.c_str()); //parse xml file directly
    xml_node rootNode = loadXML.child("SimulationEnvironment");
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
    xml_node simuParamNode = loadXML.child("MolflowSimuSettings");
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

    //Facets
    model->sh.nbFacet = geomNode.child("Facets").select_nodes("Facet").size();
    model->facets = (SubprocessFacet **)malloc(model->sh.nbFacet * sizeof(SubprocessFacet *));
    memset(model->facets, 0, model->sh.nbFacet * sizeof(SubprocessFacet *));
    idx = 0;
    bool ignoreSumMismatch = false;
    for (xml_node facetNode : geomNode.child("Facets").children("Facet")) {
        size_t nbIndex = facetNode.child("Indices").select_nodes("Indice").size();
        if (nbIndex < 3) {
            char errMsg[128];
            sprintf(errMsg, "Facet %zd has only %zd vertices. ", idx + 1, nbIndex);
            throw Error(errMsg);
        }

        model->facets[idx] = new SubprocessFacet(nbIndex);
        LoadFacet(facetNode, model->facets[idx], model->sh.nbVertex);
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
        return;
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

    for (int i = 0; i < model->sh.nbFacet; i++) {
        // Main facet params
        // Current facet
        SubprocessFacet *f = model->facets[i];
        model->CalculateFacetParams(f);

        f->sh.hitOffset = fOffset;
        fOffset += f->GetHitsSize(model->tdParams.moments.size());
    }

    PrepareToRun(model);

}

void LoaderXML::LoadSimulationState(){

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
        size_t* angleMapCache = (size_t*)malloc(facet->sh.anglemapParams.GetDataSize());

        for (int iy = 0; iy < (facet->sh.anglemapParams.thetaLowerRes + facet->sh.anglemapParams.thetaHigherRes); iy++) {
            for (int ix = 0; ix < facet->sh.anglemapParams.phiWidth; ix++) {
                angleText >> angleMapCache[iy*facet->sh.anglemapParams.phiWidth + ix];
            }
        }
        facet->sh.anglemapParams.hasRecorded = true;

        size_t mapSize = facet->sh.anglemapParams.GetMapSize();
        facet->angleMap.pdf.resize(mapSize);
        memcpy(facet->angleMap.pdf.data(), angleMapCache, facet->sh.anglemapParams.GetRecordedDataSize());

        free(angleMapCache);
    }
    else {
        facet->sh.anglemapParams.hasRecorded = false; //if angle map was incorrect, don't use it
        if (facet->sh.desorbType == DES_ANGLEMAP) facet->sh.desorbType = DES_NONE;
    }

    //Update flags
    facet->sh.isProfile = (facet->sh.profileType != PROFILE_NONE);
    //wp.isOpaque = (wp.opacity != 0.0);
    facet->sh.isTextured = ((facet->sh.texWidth * facet->sh.texHeight) > 0);
}
