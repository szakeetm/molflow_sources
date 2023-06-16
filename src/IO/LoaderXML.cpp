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

#include <sstream>
#include <set>
#include <Helper/MathTools.h>
#include <cmath>
#include <iomanip> // setprecision
#include "LoaderXML.h"
#include "TimeMoments.h"
#include "File.h"
#include "Helper/ConsoleLogger.h"
#include "Helper/StringHelper.h"
#include "Simulation/MolflowSimFacet.h"
#include <fmt/core.h>
#include <Formulas.h>

using namespace pugi;
using namespace FlowIO;

// Use work->InsertParametersBeforeCatalog(loadedParams);
// if loaded from GUI side
int LoaderXML::LoadGeometry(const std::string &inputFileName, std::shared_ptr<MolflowSimulationModel> model, GLProgress_Abstract& prg) {
    if (!model->m.try_lock()) {
        return 1;
    }

    xml_document loadXML;
    auto inputFile = inputFileName.c_str();
    xml_parse_result parseResult = loadXML.load_file(inputFile); //parse xml file directly
    xml_node rootNode = loadXML.child("SimulationEnvironment");
    if(!rootNode){
        std::cerr << "Info: XML file seems to be of older format, you can upgrade by saving with Molflow 2.8+" <<std::endl;
        rootNode = loadXML.root();
        //return 1;
    }
    model->sh.name = FileUtils::GetFilename(inputFile);

    xml_node geomNode = rootNode.child("Geometry");

    prg.SetMessage("Loading vertices...");
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
    idx = 0;
    model->structures.resize(model->sh.nbSuper);
    for (xml_node structure : geomNode.child("Structures").children("Structure")) {
        model->structures[idx].strName = structure.attribute("name").value();
        model->structures[idx].strFileName = fmt::format("{}.txt",model->structures[idx].strName); // For backward compatibilty with STR
        idx++;
    }

    //Parameters (needs to precede facets)
    xml_node simuParamNode = rootNode.child("MolflowSimuSettings");
    bool isMolflowFile = (simuParamNode != nullptr); //if no "MolflowSimuSettings" node, it's a Synrad file

    {
        if (isMolflowFile) {
            xml_node paramNode = simuParamNode.child("Parameters");
            for (xml_node newParameter : paramNode.children("Parameter")) {
                Parameter& newPar = geometrySettings.parameters.emplace_back();
                newPar.name = newParameter.attribute("name").as_string();
                if (newParameter.attribute("logXinterp")) {
                    newPar.logXinterp = newParameter.attribute("logXinterp").as_bool();
                } //else set to false by constructor
                if (newParameter.attribute("logYinterp")) {
                    newPar.logYinterp = newParameter.attribute("logYinterp").as_bool();
                } //else set to false by constructor
                for (xml_node newMoment : newParameter.children("Moment")) {
                    newPar.AddPair(std::make_pair(newMoment.attribute("t").as_double(),
                                                  newMoment.attribute("value").as_double()));
                }
            }
        }
        //TODO: Load parameters from catalog explicitly?
        //work->InsertParametersBeforeCatalog(loadedParams);
        model->tdParams.parameters.insert(model->tdParams.parameters.end(),geometrySettings.parameters.begin(),geometrySettings.parameters.end());
    }

    prg.SetMessage("Loading facets...",false);
    //Facets , load for now via temp pointers and convert to vector afterwards
    model->sh.nbFacet = geomNode.child("Facets").select_nodes("Facet").size();
    //SubprocessFac** loadFacets = (SubprocessFac **)malloc(model->sh.nbFacet * sizeof(SubprocessFac *));
    //loadFacets.reserve(model->sh.nbFacet);
    //memset(loadFacets, 0, model->sh.nbFacet * sizeof(SubprocessFac *));
    idx = 0;
    bool ignoreSumMismatch = false;
    std::vector<std::shared_ptr<SimulationFacet>> loadFacets; // tmp facet holder
    for (xml_node facetNode : geomNode.child("Facets").children("Facet")) {
        size_t nbIndex = facetNode.child("Indices").select_nodes("Indice").size();
        if (nbIndex < 3) {
            char errMsg[128];
            sprintf(errMsg, "Facet %zd has only %zd vertices. ", idx + 1, nbIndex);
            throw Error(errMsg);
        }

        loadFacets.emplace_back(std::make_shared<MolflowSimFacet>(nbIndex));
        LoadFacet(facetNode, (MolflowSimFacet*)loadFacets[idx].get(), model->sh.nbVertex);

        //Set param names for interface
        /*if (facets[idx]->sh.sticking_paramId > -1) facets[idx]->userSticking = work->parameters[facets[idx]->sh.sticking_paramId].name;
        if (facets[idx]->sh.opacity_paramId > -1) facets[idx]->userOpacity = work->parameters[facets[idx]->sh.opacity_paramId].name;
        if (facets[idx]->sh.outgassing_paramId > -1) facets[idx]->userOutgassing = work->parameters[facets[idx]->sh.outgassing_paramId].name;
*/
        idx++;
        prg.SetProgress((double)idx/(double)model->sh.nbFacet);
    }

    model->wp.gasMass = simuParamNode.child("Gas").attribute("mass").as_double();
    model->wp.halfLife = simuParamNode.child("Gas").attribute("halfLife").as_double();
    if (simuParamNode.child("Gas").attribute("enableDecay")) {
        model->wp.enableDecay = simuParamNode.child("Gas").attribute("enableDecay").as_bool();
    }
    else {
        model->wp.enableDecay = model->wp.halfLife < 1e100;
    }

    prg.SetMessage("Loading physics parameters...",false);
    xml_node timeSettingsNode = simuParamNode.child("TimeSettings");
    model->wp.timeWindowSize = timeSettingsNode.attribute("timeWindow").as_double();
    // Default initialization
    /**/
    model->wp.useMaxwellDistribution = timeSettingsNode.attribute("useMaxwellDistr").as_bool();
    model->wp.calcConstantFlow = timeSettingsNode.attribute("calcConstFlow").as_bool();

    //geometrySettings.userMoments.clear();
    xml_node userMomentsNode = timeSettingsNode.child("UserMoments");
    for (xml_node newUserEntry : userMomentsNode.children("UserEntry")) {
        UserMoment um;
        um.content=newUserEntry.attribute("content").as_string();
        um.timeWindow=newUserEntry.attribute("window").as_double();
        if(um.timeWindow==0.0){
            um.timeWindow = model->wp.timeWindowSize;
        }
        
        geometrySettings.userMoments.emplace_back(um);
    }
    try {
        TimeMoments::ParseAndCheckUserMoments(model->tdParams.moments, geometrySettings.userMoments, prg);
    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        model->m.unlock();
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

    auto forcesNode = simuParamNode.child("MeasureForces");
    if (!forcesNode) {
        model->wp.enableForceMeasurement = false;
        model->wp.torqueRefPoint = Vector3d(0.0, 0.0, 0.0);
    }
    else {
        model->wp.enableForceMeasurement = forcesNode.attribute("enabled").as_bool();
        auto torqueNode = forcesNode.child("Torque");
        if (torqueNode) {
            auto v = torqueNode.child("refPoint");
            model->wp.torqueRefPoint.x = v.attribute("x").as_double();
            model->wp.torqueRefPoint.y = v.attribute("y").as_double();
            model->wp.torqueRefPoint.z = v.attribute("z").as_double();
        }
    }

    prg.SetMessage("Loading histograms...",false);
    xml_node globalHistNode = simuParamNode.child("Global_histograms");
    if (globalHistNode) { // Molflow version before 2.8 didn't save histograms
        xml_node nbBounceNode = globalHistNode.child("Bounces");
        if (nbBounceNode) {
            model->wp.globalHistogramParams.recordBounce=true;
            model->wp.globalHistogramParams.nbBounceBinsize=nbBounceNode.attribute("binSize").as_ullong();
            model->wp.globalHistogramParams.nbBounceMax=nbBounceNode.attribute("max").as_ullong();
        }
        xml_node distanceNode = globalHistNode.child("Distance");
        if (distanceNode) {
            model->wp.globalHistogramParams.recordDistance=true;
            model->wp.globalHistogramParams.distanceBinsize=distanceNode.attribute("binSize").as_double();
            model->wp.globalHistogramParams.distanceMax=distanceNode.attribute("max").as_double();
        }
#ifdef MOLFLOW
        xml_node timeNode = globalHistNode.child("Time");
        if (timeNode) {
            model->wp.globalHistogramParams.recordTime=true;
            model->wp.globalHistogramParams.timeBinsize=timeNode.attribute("binSize").as_double();
            model->wp.globalHistogramParams.timeMax=timeNode.attribute("max").as_double();
        }
#endif
    }

    model->tdParams.IDs = this->IDs;
    model->facets = std::move(loadFacets);

    model->m.unlock();

    return 0;
}

std::vector<SelectionGroup> LoaderXML::LoadSelections(const std::string& inputFileName) {
    std::vector<SelectionGroup> selGroup;

    xml_document loadXML;
    auto inputFile = inputFileName.c_str();
    xml_parse_result parseResult = loadXML.load_file(inputFile); //parse xml file directly
    xml_node rootNode = loadXML.child("SimulationEnvironment");
    if(!rootNode){
        std::cerr << "Info: XML file seems to be of older format, you can upgrade by saving with Molflow 2.8+" <<std::endl;
        return selGroup;
    }


    xml_node interfNode = loadXML.child("Interface");
    xml_node selNode = interfNode.child("Selections");
    //int nbS = selNode.select_nodes("Selection").size();

    selGroup.reserve(std::distance(selNode.children("Selection").begin(),selNode.children("Selection").end()));
    for (xml_node sNode : selNode.children("Selection")) {
        SelectionGroup s;
        s.name = sNode.attribute("name").as_string();
        s.selection.reserve(sNode.select_nodes("selItem").size());
        for (xml_node iNode : sNode.children("selItem"))
            s.selection.push_back(iNode.attribute("facet").as_llong());
        selGroup.push_back(s);
    }

    return selGroup;
}

int LoaderXML::LoadSimulationState(const std::string &inputFileName, std::shared_ptr<MolflowSimulationModel> model,
                                   GlobalSimuState *globState, GLProgress_Abstract& prg) {

    try {
        xml_document loadXML;
        xml_parse_result parseResult = loadXML.load_file(inputFileName.c_str()); //parse xml file directly
        xml_node rootNode = loadXML.child("SimulationEnvironment");

        if (!rootNode) {
            std::cerr << "Info: XML file seems to be of older format, you can upgrade by saving with Molflow 2.8+"
                      << std::endl;
            rootNode = loadXML.root();
        }

        if (!rootNode.child("MolflowResults"))
            return 1; //simu state not saved with file

        auto lock = GetHitLock(globState, 10000);
        if (!lock) return 1;

        xml_node resultNode = rootNode.child("MolflowResults");
        xml_node momentsNode = resultNode.child("Moments");
        size_t nbMoments = momentsNode.select_nodes("Moment").size(); //Contains constant flow!
        size_t facetHitsSize = (nbMoments) * sizeof(FacetHitBuffer);
        size_t m = 0;
        for (xml_node newMoment: momentsNode.children("Moment")) {
            
            if (m == 0) { //read global results
                prg.SetMessage(fmt::format("Loading global results...",m));
                xml_node globalNode = newMoment.child("Global");
                xml_node hitsNode = globalNode.child("Hits");
                globState->globalStats.globalHits.nbMCHit = hitsNode.attribute("totalHit").as_llong();
                if (hitsNode.attribute("totalHitEquiv")) {
                    globState->globalStats.globalHits.nbHitEquiv = hitsNode.attribute("totalHitEquiv").as_double();
                } else {
                    //Backward compatibility
                    globState->globalStats.globalHits.nbHitEquiv = static_cast<double>(globState->globalStats.globalHits.nbMCHit);
                }
                globState->globalStats.globalHits.nbDesorbed = hitsNode.attribute("totalDes").as_llong();
                if (hitsNode.attribute("totalAbsEquiv")) {
                    globState->globalStats.globalHits.nbAbsEquiv = hitsNode.attribute("totalAbsEquiv").as_double();
                } else {
                    //Backward compatibility
                    globState->globalStats.globalHits.nbAbsEquiv = hitsNode.attribute("totalAbs").as_double();
                }
                if (hitsNode.attribute(
                        "totalDist_total")) { //if it's in the new format where total/partial are separated
                    globState->globalStats.distTraveled_total = hitsNode.attribute("totalDist_total").as_double();
                    globState->globalStats.distTraveledTotal_fullHitsOnly = hitsNode.attribute(
                            "totalDist_fullHitsOnly").as_double();
                } else
                    globState->globalStats.distTraveled_total = globState->globalStats.distTraveledTotal_fullHitsOnly = hitsNode.attribute(
                            "totalDist").as_double();
                globState->globalStats.nbLeakTotal = hitsNode.attribute("totalLeak").as_llong();
                //work->desorptionLimit=hitsNode.attribute("maxDesorption").as_llong();

                globState->globalStats.hitCacheSize = 0;
                xml_node hitCacheNode = globalNode.child("Hit_Cache");
                for (xml_node newHit: hitCacheNode.children("Hit")) {
                    if (globState->globalStats.hitCacheSize < HITCACHESIZE) {
                        globState->globalStats.hitCache[globState->globalStats.hitCacheSize].pos.x = newHit.attribute(
                                "posX").as_double();
                        globState->globalStats.hitCache[globState->globalStats.hitCacheSize].pos.y = newHit.attribute(
                                "posY").as_double();
                        globState->globalStats.hitCache[globState->globalStats.hitCacheSize].pos.z = newHit.attribute(
                                "posZ").as_double();
                        globState->globalStats.hitCache[globState->globalStats.hitCacheSize].type = newHit.attribute(
                                "type").as_int();
                        globState->globalStats.hitCacheSize++;
                    }
                }

                globState->globalStats.leakCacheSize = 0;
                xml_node leakCacheNode = globalNode.child("Leak_Cache");
                for (xml_node newLeak: leakCacheNode.children("Leak")) {
                    if (globState->globalStats.leakCacheSize < LEAKCACHESIZE) {
                        globState->globalStats.leakCache[globState->globalStats.leakCacheSize].pos.x = newLeak.attribute(
                                "posX").as_double();
                        globState->globalStats.leakCache[globState->globalStats.leakCacheSize].pos.y = newLeak.attribute(
                                "posY").as_double();
                        globState->globalStats.leakCache[globState->globalStats.leakCacheSize].pos.z = newLeak.attribute(
                                "posZ").as_double();
                        globState->globalStats.leakCache[globState->globalStats.leakCacheSize].dir.x = newLeak.attribute(
                                "dirX").as_double();
                        globState->globalStats.leakCache[globState->globalStats.leakCacheSize].dir.y = newLeak.attribute(
                                "dirY").as_double();
                        globState->globalStats.leakCache[globState->globalStats.leakCacheSize].dir.z = newLeak.attribute(
                                "dirZ").as_double();
                        globState->globalStats.leakCacheSize++;
                    }
                }
            } //end global node

            prg.SetMessage(fmt::format("Loading histograms [moment {}]...", m),false);
            bool hasHistogram =
                    model->wp.globalHistogramParams.recordBounce || model->wp.globalHistogramParams.recordDistance;
#ifdef MOLFLOW
            hasHistogram = hasHistogram || model->wp.globalHistogramParams.recordTime;
#endif
            if (hasHistogram) {
                xml_node histNode = newMoment.child("Histograms");
                if (histNode) { //Versions before 2.8 didn't save histograms
                    //Retrieve histogram map from hits dp
                    auto &globalHistogram = globState->globalHistograms[m];
                    if (model->wp.globalHistogramParams.recordBounce) {
                        auto &nbHitsHistogram = globalHistogram.nbHitsHistogram;
                        xml_node hist = histNode.child("Bounces");
                        if (hist) {
                            size_t histSize = model->wp.globalHistogramParams.GetBounceHistogramSize();
                            size_t saveHistSize = hist.attribute("size").as_ullong();
                            if (histSize == saveHistSize) {
                                //Can do: compare saved with expected size
                                size_t h = 0;
                                for (auto bin: hist.children("Bin")) {
                                    if (h < histSize) {
                                        nbHitsHistogram[h++] = bin.attribute("count").as_double();
                                    } else {
                                        //Treat errors
                                    }
                                }
                            } else {
                                //Treat errors
                            }
                        }
                    }
                    if (model->wp.globalHistogramParams.recordDistance) {
                        auto &distanceHistogram = globalHistogram.distanceHistogram;
                        xml_node hist = histNode.child("Distance");
                        if (hist) {
                            size_t histSize = model->wp.globalHistogramParams.GetDistanceHistogramSize();
                            size_t saveHistSize = hist.attribute("size").as_ullong();
                            if (histSize == saveHistSize) {
                                //Can do: compare saved with expected size
                                size_t h = 0;
                                for (auto bin: hist.children("Bin")) {
                                    if (h < histSize) {
                                        distanceHistogram[h++] = bin.attribute("count").as_double();
                                    } else {
                                        //Treat errors
                                    }
                                }
                            } else {
                                //Treat errors
                            }
                        }
                    }
                    if (model->wp.globalHistogramParams.recordTime) {
                        auto &timeHistogram = globalHistogram.timeHistogram;
                        xml_node hist = histNode.child("Time");
                        if (hist) {
                            size_t histSize = model->wp.globalHistogramParams.GetTimeHistogramSize();
                            size_t saveHistSize = hist.attribute("size").as_ullong();
                            if (histSize == saveHistSize) {
                                //Can do: compare saved with expected size
                                size_t h = 0;
                                for (auto bin: hist.children("Bin")) {
                                    if (h < histSize) {
                                        timeHistogram[h++] = bin.attribute("count").as_double();
                                    } else {
                                        //Treat errors
                                    }
                                }
                            } else {
                                //Treat errors
                            }
                        }
                    }
                }
            }

            prg.SetMessage(fmt::format("Loading facet results [moment {}]...", m),false);
            xml_node facetResultsNode = newMoment.child("FacetResults");
            for (xml_node newFacetResult: facetResultsNode.children("Facet")) {
                int facetId = newFacetResult.attribute("id").as_int();
                prg.SetProgress((double)((m * model->sh.nbFacet) + facetId) /
                    ((double)nbMoments * model->sh.nbFacet));
                if(facetId >= model->facets.size()){
                    throw std::runtime_error(fmt::format("Accessing simulation state for facet #{}, but only {} facets have been loaded!\nMaybe the input file is corrupted?",facetId+1, model->facets.size()));
                }
                auto sFac = model->facets[facetId];
                xml_node facetHitNode = newFacetResult.child("Hits");
                //FacetHitBuffer* facetCounter = (FacetHitBuffer *)(buffer + loadFacets[facetId].sh.hitOffset + m * sizeof(FacetHitBuffer));
                FacetHitBuffer *facetCounter = &globState->facetStates[facetId].momentResults[m].hits;
                if (facetHitNode) { //If there are hit results for the current moment
                    facetCounter->nbMCHit = facetHitNode.attribute("nbHit").as_llong();
                    if (facetHitNode.attribute("nbHitEquiv")) {
                        facetCounter->nbHitEquiv = facetHitNode.attribute("nbHitEquiv").as_double();
                    } else {
                        //Backward compatibility
                        facetCounter->nbHitEquiv = static_cast<double>(facetCounter->nbMCHit);
                    }
                    facetCounter->nbDesorbed = facetHitNode.attribute("nbDes").as_llong();
                    if (facetHitNode.attribute("nbAbsEquiv")) {
                        facetCounter->nbAbsEquiv = facetHitNode.attribute("nbAbsEquiv").as_double();
                    } else {
                        //Backward compatibility
                        facetCounter->nbAbsEquiv = facetHitNode.attribute("nbAbs").as_double();
                    }
                    facetCounter->sum_v_ort = facetHitNode.attribute("sum_v_ort").as_double();
                    facetCounter->sum_1_per_ort_velocity = facetHitNode.attribute("sum_1_per_v").as_double();
                    if (facetHitNode.attribute("sum_v")) {
                        facetCounter->sum_1_per_velocity = facetHitNode.attribute("sum_v").as_double();
                    } else {
                        //Backward compatibility
                        facetCounter->sum_1_per_velocity =
                                4.0 * Square(facetCounter->nbHitEquiv + static_cast<double>(facetCounter->nbDesorbed)) /
                                facetCounter->sum_1_per_ort_velocity;
                    }
                    auto forcesNode = newFacetResult.child("Forces");

                    if (forcesNode) { //Load if there's recorded information

                        auto impulseNode = forcesNode.child("Impulse");
                        if (impulseNode) {
                            facetCounter->impulse = Vector3d(
                                impulseNode.attribute("x").as_double(),
                                impulseNode.attribute("y").as_double(),
                                impulseNode.attribute("z").as_double()
                            );
                        }
                        else {
                            facetCounter->impulse = Vector3d(0.0, 0.0, 0.0);
                        }
                        auto impulse_sqr_Node = forcesNode.child("Impulse_square");
                        if (impulse_sqr_Node) {
                            facetCounter->impulse_square = Vector3d(
                                impulse_sqr_Node.attribute("x").as_double(),
                                impulse_sqr_Node.attribute("y").as_double(),
                                impulse_sqr_Node.attribute("z").as_double()
                            );
                        }
                        else {
                            facetCounter->impulse_square = Vector3d(0.0, 0.0, 0.0);
                        }
                        auto impulse_momentum_Node = forcesNode.child("Impulse_momentum");
                        if (impulse_momentum_Node) {
                            facetCounter->impulse_momentum = Vector3d(
                                impulse_momentum_Node.attribute("x").as_double(),
                                impulse_momentum_Node.attribute("y").as_double(),
                                impulse_momentum_Node.attribute("z").as_double()
                            );
                        }
                        else {
                            facetCounter->impulse_momentum = Vector3d(0.0, 0.0, 0.0);
                        }
                    }

                    // Do this after XML load
                    /*if (model->displayedMoment == m) { //For immediate display in facet hits list and facet counter
                        facet.facetHitCache.hit = facetCounter->hit;
                    }*/
                } else { //No hit information, so set to 0
					facetCounter->nbMCHit =
						facetCounter->nbDesorbed =
						0;
					facetCounter->sum_v_ort =
						facetCounter->nbHitEquiv =
						facetCounter->sum_1_per_ort_velocity =
						facetCounter->sum_1_per_velocity =
						facetCounter->nbAbsEquiv =
						0.0;
					facetCounter->impulse =
                        facetCounter->impulse_square =
                        facetCounter->impulse_momentum =
                        Vector3d(0.0, 0.0, 0.0);
                }

                //Profiles
                if (sFac->sh.isProfile) {
                    xml_node profileNode = newFacetResult.child("Profile");
                    //ProfileSlice *profilePtr = (ProfileSlice *)(buffer + facet.sh.hitOffset + facetHitsSize + m * sizeof(ProfileSlice)*PROFILE_SIZE);
                    std::vector<ProfileSlice> &profilePtr = globState->facetStates[facetId].momentResults[m].profile;

                    size_t id = 0;
                    for (xml_node slice: profileNode.children("Slice")) {
                        if (slice.attribute("countEquiv")) {
                            profilePtr[id].countEquiv = slice.attribute("countEquiv").as_double();
                        } else {
                            //Old format before low-flux
                            profilePtr[id].countEquiv = static_cast<double>(slice.attribute("count").as_llong());
                        }
                        profilePtr[id].sum_1_per_ort_velocity = slice.attribute("sum_1_per_v").as_double();
                        profilePtr[id].sum_v_ort = slice.attribute("sum_v_ort").as_double();
                        id++;
                    }
                }

                //Textures
                if (sFac->sh.texWidth * sFac->sh.texHeight > 0) {
                    xml_node textureNode = newFacetResult.child("Texture");
                    size_t texWidth_file = textureNode.attribute("width").as_llong();
                    size_t texHeight_file = textureNode.attribute("height").as_llong();

                    std::vector<TextureCell> &texture = globState->facetStates[facetId].momentResults[m].texture;

                    std::stringstream countText, sum1perText, sumvortText;
                    if (textureNode.child("countEquiv")) {
                        countText << textureNode.child_value("countEquiv");
                    } else {
                        countText << textureNode.child_value("count");
                    }
                    sum1perText << textureNode.child_value("sum_1_per_v");
                    sumvortText << textureNode.child_value("sum_v_ort");

                    for (size_t iy = 0; iy < (std::min(sFac->sh.texHeight,
                                                  texHeight_file)); iy++) { //MIN: If stored texture is larger, don't read extra cells
                        for (size_t ix = 0; ix < (std::min(sFac->sh.texWidth,
                                                      texWidth_file)); ix++) { //MIN: If stored texture is larger, don't read extra cells
                            countText >> texture[iy * sFac->sh.texWidth + ix].countEquiv;
                            sum1perText >> texture[iy * sFac->sh.texWidth + ix].sum_1_per_ort_velocity;
                            sumvortText >> texture[iy * sFac->sh.texWidth + ix].sum_v_ort_per_area;

                        }
                        for (int ie = 0; ie < texWidth_file -
                                              sFac->sh.texWidth; ie++) {//Executed if file texture is bigger than expected texture
                            //Read extra cells from file without doing anything
                            size_t dummy_ll;
                            double dummy_d;
                            countText >> dummy_ll;
                            sum1perText >> dummy_d;
                            sumvortText >> dummy_d;

                        }
                    }
                    for (int ie = 0; ie < texHeight_file -
                                          sFac->sh.texHeight; ie++) {//Executed if file texture is bigger than expected texture
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

                if (sFac->sh.countDirection) {
                    xml_node dirNode = newFacetResult.child("Directions");
                    if (dirNode.attribute("width").as_int() != sFac->sh.texWidth ||
                        dirNode.attribute("height").as_int() != sFac->sh.texHeight) {
                        std::stringstream msg;
                        msg << "Direction texture size mismatch on facet " << facetId + 1 << ".\nExpected: "
                            << sFac->sh.texWidth << "x" << sFac->sh.texHeight << "\n"
                            << "In file: " << dirNode.attribute("width").as_int() << "x"
                            << dirNode.attribute("height").as_int();
                        throw Error(msg.str().c_str());

                    }
                    
                    std::vector<DirectionCell> &dirs = globState->facetStates[facetId].momentResults[m].direction;

                    std::stringstream dirText, dirCountText;
                    dirText << dirNode.child_value("vel.vectors");
                    dirCountText << dirNode.child_value("count");

                    for (size_t iy = 0; iy < sFac->sh.texHeight; iy++) {
                        for (size_t ix = 0; ix < sFac->sh.texWidth; ix++) {
                            std::string component;
                            std::getline(dirText, component, ',');
                            dirs[iy * sFac->sh.texWidth + ix].dir.x = std::stod(component);
                            std::getline(dirText, component, ',');
                            dirs[iy * sFac->sh.texWidth + ix].dir.y = std::stod(component);
                            dirText >> dirs[iy * sFac->sh.texWidth + ix].dir.z;
                            dirCountText >> dirs[iy * sFac->sh.texWidth + ix].count;
                        }
                    }
                } //end directions

                // Facet histogram
                hasHistogram =
                        sFac->sh.facetHistogramParams.recordBounce || sFac->sh.facetHistogramParams.recordDistance;
#ifdef MOLFLOW
                hasHistogram = hasHistogram || sFac->sh.facetHistogramParams.recordTime;
#endif
                if (hasHistogram) {
                    xml_node histNode = newFacetResult.child("Histograms");
                    if (histNode) { //Versions before 2.8 didn't save histograms
                        //Retrieve histogram map from hits dp
                        auto &facetHistogram = globState->facetStates[facetId].momentResults[m].histogram;
                        if (sFac->sh.facetHistogramParams.recordBounce) {
                            auto &nbHitsHistogram = facetHistogram.nbHitsHistogram;
                            xml_node hist = histNode.child("Bounces");
                            if (hist) {
                                size_t histSize = sFac->sh.facetHistogramParams.GetBounceHistogramSize();
                                size_t saveHistSize = hist.attribute("size").as_ullong();
                                if (histSize == saveHistSize) {
                                    //Can do: compare saved with expected size
                                    size_t h = 0;
                                    for (auto bin: hist.children("Bin")) {
                                        if (h < histSize) {
                                            nbHitsHistogram[h++] = bin.attribute("count").as_double();
                                        } else {
                                            //Treat errors
                                        }
                                    }
                                } else {
                                    //Treat errors
                                }
                            }
                        }
                        if (sFac->sh.facetHistogramParams.recordDistance) {
                            auto &distanceHistogram = facetHistogram.distanceHistogram;
                            xml_node hist = histNode.child("Distance");
                            if (hist) {
                                size_t histSize = sFac->sh.facetHistogramParams.GetDistanceHistogramSize();
                                size_t saveHistSize = hist.attribute("size").as_ullong();
                                if (histSize == saveHistSize) {
                                    //Can do: compare saved with expected size
                                    size_t h = 0;
                                    for (auto bin: hist.children("Bin")) {
                                        if (h < histSize) {
                                            distanceHistogram[h++] = bin.attribute("count").as_double();
                                        } else {
                                            //Treat errors
                                        }
                                    }
                                } else {
                                    //Treat errors
                                }
                            }
                        }
                        if (sFac->sh.facetHistogramParams.recordTime) {
                            auto &timeHistogram = facetHistogram.timeHistogram;
                            xml_node hist = histNode.child("Time");
                            if (hist) {
                                size_t histSize = sFac->sh.facetHistogramParams.GetTimeHistogramSize();
                                size_t saveHistSize = hist.attribute("size").as_ullong();
                                if (histSize == saveHistSize) {
                                    //Can do: compare saved with expected size
                                    size_t h = 0;
                                    for (auto bin: hist.children("Bin")) {
                                        if (h < histSize) {
                                            timeHistogram[h++] = bin.attribute("count").as_double();
                                        } else {
                                            //Treat errors
                                        }
                                    }
                                } else {
                                    //Treat errors
                                }
                            }
                        }
                    }
                }
            } //end facetResult
            m++;
        } //end moment

        prg.SetMessage("File load complete.\n",false);
    }
    catch (const std::exception &e) {
        Log::console_error("[LoaderXML] {}", e.what());
        throw;
    }
    return 0;
}

void LoaderXML::LoadFacet(pugi::xml_node facetNode, MolflowSimFacet *facet, size_t nbTotalVertices) {
    int idx = 0;
    bool ignoreSumMismatch = true;
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

    // TODO: Only parse Molflow files
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
    facet->sh.texWidth_precise = texNode.attribute("texDimX").as_double();
    facet->sh.texHeight_precise = texNode.attribute("texDimY").as_double();
    facet->sh.countDes = texNode.attribute("countDes").as_bool();
    facet->sh.countAbs = texNode.attribute("countAbs").as_bool();
    facet->sh.countRefl = texNode.attribute("countRefl").as_bool();
    facet->sh.countTrans = texNode.attribute("countTrans").as_bool();
    facet->sh.countDirection = texNode.attribute("countDir").as_bool();
    facet->sh.countACD = texNode.attribute("countAC").as_bool();
    
    xml_node outgNode = facetNode.child("DynamicOutgassing");
    if ((hasOutgassingFile) && outgNode && outgNode.child("map")) {
        facet->ogMap.outgassingMapWidth = outgNode.attribute("width").as_int();
        facet->ogMap.outgassingMapHeight = outgNode.attribute("height").as_int();
        if (outgNode.attribute("ratioU")) { //New format supporting non-square textures
            facet->ogMap.outgassingFileRatioU = outgNode.attribute("ratioU").as_double();
            facet->ogMap.outgassingFileRatioV = outgNode.attribute("ratioV").as_double();
        }
        else { //Old format for square textures
            facet->ogMap.outgassingFileRatioU = facet->ogMap.outgassingFileRatioV = outgNode.attribute("ratio").as_double();
        }
        facet->ogMap.totalDose = outgNode.attribute("totalDose").as_double();
        facet->sh.totalOutgassing = outgNode.attribute("totalOutgassing").as_double();
        facet->ogMap.totalFlux = outgNode.attribute("totalFlux").as_double();

        double sum = 0.0;

        std::stringstream outgText;
        outgText << outgNode.child_value("map");

        auto& ogMap = facet->ogMap;
        std::vector<double>(ogMap.outgassingMapWidth*ogMap.outgassingMapHeight).swap(ogMap.outgassingMap);

        for (size_t iy = 0; iy < ogMap.outgassingMapHeight; iy++) {
            for (size_t ix = 0; ix < ogMap.outgassingMapWidth; ix++) {
                outgText >> ogMap.outgassingMap[iy*ogMap.outgassingMapWidth + ix];
                sum += ogMap.outgassingMap[iy*ogMap.outgassingMapWidth + ix];
            }
        }
        if (!ignoreSumMismatch && !IsEqual(sum, facet->sh.totalOutgassing)) {
            std::stringstream msg; msg << std::setprecision(8);
            msg << "Facet " << facetId + 1 << ":\n";
            msg << "The total dynamic outgassing (" << 10.0 * facet->sh.totalOutgassing << " mbar.l/s)\n";
            msg << "doesn't match the sum of the dynamic outgassing cells (" << 10.0 * sum << " mbar.l/s).";
            /*if (1 == GLMessageBox::Display(msg.str(), "Dynamic outgassing mismatch", { "OK","Ignore rest" },GLDLG_ICONINFO))
                ignoreSumMismatch = true;*/
            std::cerr << msg.str() << std::endl;
            ignoreSumMismatch = true; // TODO: For now silence after one
        }
    }
    else {
        hasOutgassingFile = facet->sh.useOutgassingFile = false; //if outgassing map was incorrect, don't use it
        facet->ogMap.outgassingMapWidth = 0;
        facet->ogMap.outgassingMapHeight = 0;
        facet->ogMap.outgassingFileRatioU = 0.0;
        facet->ogMap.outgassingFileRatioV = 0.0;
        facet->ogMap.totalDose = 0.0;
        facet->sh.totalOutgassing = 0.0;
        facet->ogMap.totalFlux = 0.0;
    }
    xml_node angleMapNode = facetNode.child("IncidentAngleMap");
    if (angleMapNode && angleMapNode.child("map") && angleMapNode.attribute("angleMapThetaLimit")) {

        facet->sh.anglemapParams.phiWidth = angleMapNode.attribute("angleMapPhiWidth").as_ullong();
        facet->sh.anglemapParams.thetaLimit = angleMapNode.attribute("angleMapThetaLimit").as_double();
        facet->sh.anglemapParams.thetaLowerRes = angleMapNode.attribute("angleMapThetaLowerRes").as_ullong();
        facet->sh.anglemapParams.thetaHigherRes = angleMapNode.attribute("angleMapThetaHigherRes").as_ullong();

        std::stringstream angleText;
        angleText << angleMapNode.child_value("map");


        //std::map<size_t,std::vector<size_t>> angleMapCache;
        //size_t* angleMapCache = (size_t*)malloc(facet->sh.anglemapParams.GetDataSize());
        //angleMapCache.emplace(std::make_pair(facet->globalId,std::vector<size_t>()));
        auto& angleMap = facet->angleMap.pdf;
        try {
            angleMap.clear();
            angleMap.resize( facet->sh.anglemapParams.phiWidth * (facet->sh.anglemapParams.thetaLowerRes + facet->sh.anglemapParams.thetaHigherRes));
        }
        catch(...) {
            std::stringstream err;
            err << "Not enough memory for incident angle map on facet ";
            throw Error(err.str().c_str());
        }

        size_t angleMapSum = 0;
        for (size_t iy = 0; iy < (facet->sh.anglemapParams.thetaLowerRes + facet->sh.anglemapParams.thetaHigherRes); iy++) {
            for (size_t ix = 0; ix < facet->sh.anglemapParams.phiWidth; ix++) {
                angleText >> angleMap[iy*facet->sh.anglemapParams.phiWidth + ix];
                angleMapSum += angleMap[iy*facet->sh.anglemapParams.phiWidth + ix];
            }
        }
        /*if(angleMapSum > 0) // only has recorded if at least one value is set
            facet->sh.anglemapParams.hasRecorded = true;*/
    }
    else {
        //facet->sh.anglemapParams.hasRecorded = false; //if angle map was incorrect, don't use it
        if (facet->sh.desorbType == DES_ANGLEMAP) facet->sh.desorbType = DES_NONE;
    }

    /*
    // Init by default as true
    facet->viewSettings.textureVisible = facetNode.child("ViewSettings").attribute("textureVisible").as_bool(true);
    facet->viewSettings.volumeVisible = facetNode.child("ViewSettings").attribute("volumeVisible").as_bool(true);
    */

    xml_node facetHistNode = facetNode.child("Histograms");
    if (facetHistNode) { // Molflow version before 2.8 didn't save histograms
        xml_node nbBounceNode = facetHistNode.child("Bounces");
        if (nbBounceNode) {
            facet->sh.facetHistogramParams.recordBounce=true;
            facet->sh.facetHistogramParams.nbBounceBinsize=nbBounceNode.attribute("binSize").as_ullong();
            facet->sh.facetHistogramParams.nbBounceMax=nbBounceNode.attribute("max").as_ullong();
        }
        xml_node distanceNode = facetHistNode.child("Distance");
        if (distanceNode) {
            facet->sh.facetHistogramParams.recordDistance=true;
            facet->sh.facetHistogramParams.distanceBinsize=distanceNode.attribute("binSize").as_double();
            facet->sh.facetHistogramParams.distanceMax=distanceNode.attribute("max").as_double();
        }
#ifdef MOLFLOW
        xml_node timeNode = facetHistNode.child("Time");
        if (timeNode) {
            facet->sh.facetHistogramParams.recordTime=true;
            facet->sh.facetHistogramParams.timeBinsize=timeNode.attribute("binSize").as_double();
            facet->sh.facetHistogramParams.timeMax=timeNode.attribute("max").as_double();
        }
#endif
    }

    //Update flags
    facet->sh.isProfile = (facet->sh.profileType != PROFILE_NONE);
    //wp.isOpaque = (wp.opacity != 0.0);
    facet->sh.isTextured = ((facet->sh.texWidth_precise * facet->sh.texHeight_precise) > 0);

    // Do some fixes
    bool hasAnyTexture = facet->sh.countDes || facet->sh.countAbs || facet->sh.countRefl || facet->sh.countTrans || facet->sh.countACD;
    if (!facet->sh.isTextured && (hasAnyTexture)) {
        facet->sh.countDes = false;
        facet->sh.countAbs = false;
        facet->sh.countRefl = false;
        facet->sh.countTrans = false;
        facet->sh.countACD = false;
        
        std::stringstream msg; msg << std::setprecision(8);
        msg << "Facet (#"<< facetId << ") has no valid mesh, but active texture counters: removing...\n";
        std::cerr << msg.str();
    }
}

int LoaderXML::LoadConvergenceValues(const std::string &inputFileName, std::vector<ConvergenceData> *convergenceValues,
                                     GLProgress_Abstract& prg) {

    xml_document loadXML;
    xml_parse_result parseResult = loadXML.load_file(inputFileName.c_str()); //parse xml file directly
    xml_node rootNode = loadXML.child("SimulationEnvironment");

    if (!rootNode) {
        std::cerr << "Info: XML file seems to be of older format, you can upgrade by saving with Molflow 2.8+"
                  << std::endl;
        rootNode = loadXML.root();
    }

    if (!rootNode.child("MolflowResults"))
        return 1; //simu state not saved with file

    xml_node resultNode = rootNode.child("MolflowResults");
    xml_node convNode = resultNode.child("Convergence");

    convergenceValues->resize(0);
    for(auto& convVec : convNode.children()){
        std::stringstream convText;
        ConvergenceData convData;
        std::vector<std::pair<size_t, double>>& vec = convData.conv_vec;
        convText << convVec.child_value();
        // get length of file:
        convText.seekg (0, std::stringstream::end);
        //int length = convText.tellg();
        convText.seekg (0, std::stringstream::beg);
        if(convText.peek() == '\n') {
            char nl;
            convText.get(nl);
        }
        std::string line;
        while(!convText.eof()){
            std::getline(convText, line);
            size_t posOfTab = line.find ('\t');
            //std::string second = pieces.substr(pos + 1);

            if (posOfTab==std::string::npos)
                continue;
            size_t nbDes = 0;
            double convVal = 0.0;
            try{
                nbDes = stringToNumber<size_t>(line.substr(0, posOfTab));
                convVal = stringToNumber<double>(line.substr(posOfTab+1));
            }
            catch (const std::exception &e){
                // Just write an error and move to next line e.g. when fail on inf/nan
                std::cerr << "[XML][Convergence] Parsing error: "<<e.what()<< std::endl;
                continue;
            }

            //if(nbDes < vec[vec.size()-1].first) break; // skip if data is malformed (desorptions should increase)
            vec.emplace_back(std::make_pair(nbDes, convVal));
        }
        convergenceValues->push_back(convData);
    }

    return 0;
}
/*
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
}*/
