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

#include <iomanip>      // std::setprecision
#include <sstream>
#include <Helper/StringHelper.h>
#include <Helper/ConsoleLogger.h>
#include <Helper/GLProgress_abstract.hpp>
#include <pugixml.hpp>

#include "WriterXML.h"
#include "versionId.h"
#include "Simulation/MolflowSimFacet.h"
#include "Simulation/MolflowSimGeom.h"
#include "MolflowTypes.h"
#include "GeometryTypes.h"

using namespace FlowIO;
using namespace pugi;

xml_node XmlWriter::GetRootNode(xml_document &saveDoc) {
    // Check whether we do a old to new format update
    // (CLI doesn't create save files from scratch, to keep interface settings)
    // If yes, move old scheme to new scheme by adding the root node "SimulationEnvironment"
    if(updateRootNode){
        //Determine if original file is in old format
        bool oldFormatUsed = false;
        auto rootNode = saveDoc.document_element();
        if(!saveDoc.child("SimulationEnvironment")){
            rootNode = saveDoc.root();
            if(rootNode.child("Geometry")){
                oldFormatUsed = true;
            }
        }
        // If in old format and an update is needed, create SimulationEnvironment root and copy old children
        // (that were originally at root level)
        if(oldFormatUsed && !useOldXMLFormat) {
            xml_document newDoc;
            auto newHead = newDoc.append_child("SimulationEnvironment");
            newDoc.append_attribute("type") = "molflow";
            newDoc.append_attribute("version") = appVersionId;
            for(auto& node : rootNode.children()){
                newHead.append_copy(node);
            }

            saveDoc.reset(newDoc);
        }
    }

    xml_node rootNode;
    if (useOldXMLFormat) {
        rootNode = saveDoc.root(); //If we keep old format, just return the original root (downgrade not supported)
    } else { //If new format should be used
        if(updateRootNode) {
            rootNode = saveDoc.document_element();
            rootNode = saveDoc.child("SimulationEnvironment");
        }
        if(!rootNode) {
            if(!saveDoc.child("SimulationEnvironment")) {
                rootNode = saveDoc.append_child("SimulationEnvironment");
                rootNode.append_attribute("type") = "molflow";
                rootNode.append_attribute("version") = appVersionId;
            }
            else
                rootNode = saveDoc.document_element();
        }
    }

    return rootNode;
}

void XmlWriter::SaveGeometry(pugi::xml_document &saveDoc, const std::shared_ptr<MolflowSimulationModel> model,
            GLProgress_Abstract& prg, const std::vector<size_t> &selectionToSave) {
    xml_node rootNode = GetRootNode(saveDoc);


    bool saveAllFacets = selectionToSave.empty();

    if(updateRootNode)
        rootNode.remove_child("Geometry");
    prg.SetMessage("Saving vertices...");
    xml_node geomNode = rootNode.prepend_child("Geometry");
    geomNode.append_child("Vertices").append_attribute(
            "nb") = model->vertices3.size(); //creates Vertices node, adds nb attribute and sets its value to sp.nbVertex
    for (size_t i = 0; i < model->vertices3.size(); i++) {
        //prg.SetProgress(0.166*((double)i / (double)model->vertices3.size()));
        xml_node v = geomNode.child("Vertices").append_child("Vertex");
        v.append_attribute("id") = i;
        v.append_attribute("x") = model->vertices3[i].x;
        v.append_attribute("y") = model->vertices3[i].y;
        v.append_attribute("z") = model->vertices3[i].z;
    }
    prg.SetMessage("Writing facets...", false);
    
    geomNode.append_child("Facets");
    
    //Save every facet or only those of the selection
    size_t nbF = saveAllFacets ? model->facets.size() : selectionToSave.size();
    geomNode.child("Facets").append_attribute("nb") = nbF;
    for (size_t i = 0; i < nbF; i++) {
        prg.SetProgress((double)i / (double)nbF);
        xml_node currentFacetNode = geomNode.child("Facets").append_child("Facet");
        currentFacetNode.append_attribute("id") = i;
        size_t facetId = saveAllFacets ? i : selectionToSave[i];
        auto mfFac = std::static_pointer_cast<MolflowSimFacet>(model->facets[facetId]);
        SaveFacet(currentFacetNode, mfFac, model->vertices3.size(),model->tdParams);
    }
    

    prg.SetMessage("Writing structures...");
    geomNode.append_child("Structures").append_attribute("nb") = model->sh.nbSuper;

    if (model->structures.size() == model->sh.nbSuper) {
        for (size_t i = 0; i < model->sh.nbSuper; i++) {
            xml_node s = geomNode.child("Structures").append_child("Structure");
            s.append_attribute("id") = i;
            s.append_attribute("name") = model->structures[i].name.c_str();
        }
    }

    // Simulation Settings
    prg.SetMessage("Writing simulation parameters...");
    if(updateRootNode)
        rootNode.remove_child("MolflowSimuSettings");
    xml_node simuParamNode = rootNode.insert_child_after("MolflowSimuSettings",geomNode);

    simuParamNode.append_child("Gas").append_attribute("mass") = model->sp.gasMass;
    simuParamNode.child("Gas").append_attribute(
            "enableDecay") = (int) model->sp.enableDecay; //backward compatibility: 0 or 1
    simuParamNode.child("Gas").append_attribute("halfLife") = model->sp.halfLife;
    xml_node scatteringNode = simuParamNode.child("Gas").append_child("Scattering");
    scatteringNode.append_attribute("enabled") = model->sp.scattering.enabled;
    scatteringNode.append_attribute("massRatio") = model->sp.scattering.massRatio;
    scatteringNode.append_attribute("meanFreePath_cm") = model->sp.scattering.meanFreePath_cm;

    auto lowFluxNode = simuParamNode.append_child("LowFluxMode");
    lowFluxNode.append_attribute("enabled") = model->otfParams.lowFluxMode;
    lowFluxNode.append_attribute("cutoff") = model->otfParams.lowFluxCutoff;

    xml_node timeSettingsNode = simuParamNode.append_child("TimeSettings");

    xml_node userMomentsNode = timeSettingsNode.append_child("UserMoments");
    userMomentsNode.append_attribute("nb") = interfaceSettings->userMoments.size();
    for (size_t i = 0; i < interfaceSettings->userMoments.size(); i++) {
        xml_node newUserEntry = userMomentsNode.append_child("UserEntry");
        newUserEntry.append_attribute("id") = i;
        newUserEntry.append_attribute("content") = interfaceSettings->userMoments[i].content.c_str();
        newUserEntry.append_attribute("window") = interfaceSettings->userMoments[i].timeWindow;
    }

    timeSettingsNode.append_attribute("timeWindow") = model->sp.timeWindowSize;
    timeSettingsNode.append_attribute(
            "useMaxwellDistr") = (int) model->sp.useMaxwellDistribution; //backward compatibility: 0 or 1
    timeSettingsNode.append_attribute(
            "calcConstFlow") = (int) model->sp.calcConstantFlow; //backward compatibility: 0 or 1

    xml_node motionNode = simuParamNode.append_child("Motion");
    motionNode.append_attribute("type") = model->sp.motionType;
    if (model->sp.motionType == 1) { //fixed motion
        xml_node v = motionNode.append_child("VelocityVector");
        v.append_attribute("vx") = model->sp.motionVector2.x;
        v.append_attribute("vy") = model->sp.motionVector2.y;
        v.append_attribute("vz") = model->sp.motionVector2.z;
    } else if (model->sp.motionType == 2) { //rotation
        xml_node v = motionNode.append_child("AxisBasePoint");
        v.append_attribute("x") = model->sp.motionVector1.x;
        v.append_attribute("y") = model->sp.motionVector1.y;
        v.append_attribute("z") = model->sp.motionVector1.z;
        xml_node v2 = motionNode.append_child("RotationVector");
        v2.append_attribute("x") = model->sp.motionVector2.x;
        v2.append_attribute("y") = model->sp.motionVector2.y;
        v2.append_attribute("z") = model->sp.motionVector2.z;
    }

    auto forcesNode = simuParamNode.append_child("MeasureForces");
    forcesNode.append_attribute("enabled") = model->sp.enableForceMeasurement;
    auto torqueNode = forcesNode.append_child("Torque");
    auto v = torqueNode.append_child("refPoint");
    v.append_attribute("x") = model->sp.torqueRefPoint.x;
    v.append_attribute("y") = model->sp.torqueRefPoint.y;
    v.append_attribute("z") = model->sp.torqueRefPoint.z;

    xml_node paramNode = simuParamNode.append_child("Parameters");
    size_t nbNonCatalogParameters = 0;

    for (auto &parameter : model->tdParams.parameters) {
        if (!parameter.fromCatalog) { //Don't save catalog parameters
            xml_node newParameter = paramNode.append_child("Parameter");
            newParameter.append_attribute("id") = nbNonCatalogParameters;
            newParameter.append_attribute("name") = parameter.name.c_str();
            newParameter.append_attribute("nbMoments") = (int) parameter.GetSize();
            newParameter.append_attribute("logXinterp") = parameter.logXinterp;
            newParameter.append_attribute("logYinterp") = parameter.logYinterp;
            for (size_t m = 0; m < parameter.GetSize(); m++) {
                xml_node newMoment = newParameter.append_child("Moment");
                newMoment.append_attribute("id") = m;
                newMoment.append_attribute("t") = parameter.GetX(m);
                newMoment.append_attribute("value") = parameter.GetY(m);
            }
            nbNonCatalogParameters++;
        }
    }
    paramNode.append_attribute("nb") = nbNonCatalogParameters;

    xml_node globalHistNode = simuParamNode.append_child("Global_histograms");
    if (model->sp.globalHistogramParams.recordBounce) {
        xml_node nbBounceNode = globalHistNode.append_child("Bounces");
        nbBounceNode.append_attribute("binSize") = model->sp.globalHistogramParams.nbBounceBinsize;
        nbBounceNode.append_attribute("max") = model->sp.globalHistogramParams.nbBounceMax;
    }
    if (model->sp.globalHistogramParams.recordDistance) {
        xml_node distanceNode = globalHistNode.append_child("Distance");
        distanceNode.append_attribute("binSize") = model->sp.globalHistogramParams.distanceBinsize;
        distanceNode.append_attribute("max") = model->sp.globalHistogramParams.distanceMax;
    }
#ifdef MOLFLOW
    if (model->sp.globalHistogramParams.recordTime) {
        xml_node timeNode = globalHistNode.append_child("Time");
        timeNode.append_attribute("binSize") = model->sp.globalHistogramParams.timeBinsize;
        timeNode.append_attribute("max") = model->sp.globalHistogramParams.timeMax;
    }
#endif

    prg.SetMessage("Writing interface settings...");
    if (updateRootNode)
        rootNode.remove_child("Interface");
    xml_node interfNode = rootNode.insert_child_after("Interface", simuParamNode);

    xml_node selNode = interfNode.append_child("Selections");
    selNode.append_attribute("nb") = saveAllFacets ? interfaceSettings->selections.size() : 0; //Don't save sel.groups if save file restricted to a subset of facets
    if (saveAllFacets) {
        for (size_t i = 0; i < interfaceSettings->selections.size(); i++) { //don't save selections when exporting part of the geometry (saveSelected)
            xml_node newSel = selNode.append_child("Selection");
            newSel.append_attribute("id") = i;
            newSel.append_attribute("name") = interfaceSettings->selections[i].name.c_str();
            newSel.append_attribute("nb") = interfaceSettings->selections[i].facetIds.size();
            for (size_t j = 0; j < interfaceSettings->selections[i].facetIds.size(); j++) {
                xml_node newItem = newSel.append_child("selItem");
                newItem.append_attribute("id") = j;
                newItem.append_attribute("facet") = interfaceSettings->selections[i].facetIds[j];
            }
        }
    }

	if (saveAllFacets) { //don't save views when exporting part of the geometry (saveSelected)
		xml_node viewNode = interfNode.append_child(useOldXMLFormat ? "Views" : "UserSavedViews"); // New name since 2.9.17
		viewNode.append_attribute("nb") = saveAllFacets ? interfaceSettings->userViews.size() : 0;

		for (int i = 0; i < interfaceSettings->userViews.size(); i++) { 
			const auto& v = interfaceSettings->userViews[i];
			xml_node newView = viewNode.append_child("View");
			newView.append_attribute("id") = i;
			CameraViewToXml(v, newView);
		}
	}

	if (saveAllFacets && !interfaceSettings->viewerCurrentViews.empty()) {
		xml_node viewNode = interfNode.append_child("ViewerCurrentViews");

		for (const auto& currentView : interfaceSettings->viewerCurrentViews) {
			xml_node newView = viewNode.append_child("View");
			CameraViewToXml(currentView, newView);
		}
	}

    xml_node formulaNode = interfNode.append_child("Formulas");
    formulaNode.append_attribute("nb") = saveAllFacets ? interfaceSettings->userFormulas.size() : 0;
    if (saveAllFacets) { //don't save formulas when exporting part of the geometry (saveSelected)
        for (size_t i = 0; i < interfaceSettings->userFormulas.size(); i++) {
            xml_node newFormula = formulaNode.append_child("Formula");
            newFormula.append_attribute("id") = i;
            newFormula.append_attribute("name") = interfaceSettings->userFormulas[i].name.c_str();
            newFormula.append_attribute("expression") = interfaceSettings->userFormulas[i].expression.c_str();
        }
    }

    if (interfaceSettings->profilePlotterSettings.hasData) {
        xml_node profilePlotterNode = interfNode.append_child("ProfilePlotter");
        profilePlotterNode.append_child("Parameters").append_attribute("logScale") = interfaceSettings->profilePlotterSettings.logYscale;
        xml_node viewsNode = profilePlotterNode.append_child("Views");
        for (int v : interfaceSettings->profilePlotterSettings.viewIds) {
            xml_node view = viewsNode.append_child("View");
            view.append_attribute("facetId") = v;
        }
    }

    if (interfaceSettings->convergencePlotterSettings.hasData) {
        xml_node convergencePlotterNode = interfNode.append_child("ConvergencePlotter");
        convergencePlotterNode.append_child("Parameters").append_attribute("logScale") = interfaceSettings->convergencePlotterSettings.logYscale;
        xml_node viewsNode = convergencePlotterNode.append_child("Views");
        for (int v : interfaceSettings->convergencePlotterSettings.viewIds) {
            xml_node view = viewsNode.append_child("View");
            view.append_attribute("formulaHash") = v;
        }
    }

	//Texture Min/Max
	xml_node textureSettingsNode = interfNode.append_child("TextureSettings");

	xml_node textureDisplayNode = textureSettingsNode.append_child("DisplaySettings");
	textureDisplayNode.append_attribute("autoscale") = interfaceSettings->texAutoScale;
	textureDisplayNode.append_attribute("autoscaleMode") = static_cast<int>(interfaceSettings->texAutoscaleMode);
	textureDisplayNode.append_attribute("color") = interfaceSettings->texColormap;
	textureDisplayNode.append_attribute("logScale") = interfaceSettings->texLogScale;
	textureDisplayNode.append_attribute("physicsMode") = interfaceSettings->textureMode;
	{
		xml_node textureLimitsNode = textureSettingsNode.append_child("Limits");
		xml_node autoNode = textureLimitsNode.append_child("Autoscale");
		autoNode.append_child("Constant_flow").append_child("Pressure").append_attribute("min") = interfaceSettings->textureLimits[0].autoscale.min.steady_state;
		autoNode.child("Constant_flow").child("Pressure").append_attribute("max") = interfaceSettings->textureLimits[0].autoscale.max.steady_state;
		autoNode.child("Constant_flow").append_child("Density").append_attribute("min") = interfaceSettings->textureLimits[1].autoscale.min.steady_state;
		autoNode.child("Constant_flow").child("Density").append_attribute("max") = interfaceSettings->textureLimits[1].autoscale.max.steady_state;
		autoNode.child("Constant_flow").append_child("Imp.rate").append_attribute("min") = interfaceSettings->textureLimits[2].autoscale.min.steady_state;
		autoNode.child("Constant_flow").child("Imp.rate").append_attribute("max") = interfaceSettings->textureLimits[2].autoscale.max.steady_state;

		autoNode.append_child("Moments_only").append_child("Pressure").append_attribute("min") = interfaceSettings->textureLimits[0].autoscale.min.moments_only;
		autoNode.child("Moments_only").child("Pressure").append_attribute("max") = interfaceSettings->textureLimits[0].autoscale.max.moments_only;
		autoNode.child("Moments_only").append_child("Density").append_attribute("min") = interfaceSettings->textureLimits[1].autoscale.min.moments_only;
		autoNode.child("Moments_only").child("Density").append_attribute("max") = interfaceSettings->textureLimits[1].autoscale.max.moments_only;
		autoNode.child("Moments_only").append_child("Imp.rate").append_attribute("min") = interfaceSettings->textureLimits[2].autoscale.min.moments_only;
		autoNode.child("Moments_only").child("Imp.rate").append_attribute("max") = interfaceSettings->textureLimits[2].autoscale.max.moments_only;

		xml_node manualNode = textureLimitsNode.append_child("Manual"); //Manual: only "steady state" variable used
		manualNode.append_child("Constant_flow").append_child("Pressure").append_attribute("min") = interfaceSettings->textureLimits[0].manual.min.steady_state;
		manualNode.child("Constant_flow").child("Pressure").append_attribute("max") = interfaceSettings->textureLimits[0].manual.max.steady_state;
		manualNode.child("Constant_flow").append_child("Density").append_attribute("min") = interfaceSettings->textureLimits[1].manual.min.steady_state;
		manualNode.child("Constant_flow").child("Density").append_attribute("max") = interfaceSettings->textureLimits[1].manual.max.steady_state;
		manualNode.child("Constant_flow").append_child("Imp.rate").append_attribute("min") = interfaceSettings->textureLimits[2].manual.min.steady_state;
		manualNode.child("Constant_flow").child("Imp.rate").append_attribute("max") = interfaceSettings->textureLimits[2].manual.max.steady_state;
	}
}

// Save XML document to file
bool XmlWriter::WriteXMLToFile(xml_document &saveDoc, const std::string &outputFileName) {
    if (!saveDoc.save_file(outputFileName.c_str())) {
        std::cerr << "Error writing XML file." << std::endl;
        return false;
    }
    return true;
}

// Directly append to file (load + save)
bool XmlWriter::AppendSimulationStateToFile(const std::string& outputFileName, const std::shared_ptr<MolflowSimulationModel> model, GLProgress_Abstract& prg, const std::shared_ptr<GlobalSimuState> globalState) {                                
    xml_document saveDoc;
    xml_parse_result parseResult = saveDoc.load_file(outputFileName.c_str()); //parse xml file directly

    SaveSimulationState(saveDoc, model, prg, globalState);

    if (!saveDoc.save_file(outputFileName.c_str())) {
        std::cerr << "Error writing XML file." << std::endl; //successful save
        return false;
    }
    return true;
}

// Append to open XML node
bool XmlWriter::SaveSimulationState(xml_document &saveDoc, const std::shared_ptr<MolflowSimulationModel> model, GLProgress_Abstract& prg, const std::shared_ptr<GlobalSimuState> globalState) {

    auto lock = GetHitLock(globalState.get(), 10000);
    if (!lock) return false;

    xml_node rootNode = GetRootNode(saveDoc);
    rootNode.remove_child("MolflowResults"); // clear previous results to replace with new status

    xml_node resultNode = rootNode.append_child("MolflowResults");
    xml_node momentsNode = resultNode.append_child("Moments");
    momentsNode.append_attribute("nb") = model->tdParams.moments.size() + 1;
    //size_t facetHitsSize = (1 + model->tdParams.moments.size()) * sizeof(FacetHitBuffer);

    prg.SetMessage("Writing simulation results...");
    for (size_t m = 0; m <= model->tdParams.moments.size(); m++) {
        prg.SetProgress(0.5 + 0.5 * (double) m / (1.0 + (double) model->tdParams.moments.size()));
        xml_node newMoment = momentsNode.append_child("Moment");
        newMoment.append_attribute("id") = m;
        if (m == 0) {
            newMoment.append_attribute("time") = "Constant flow";
            newMoment.append_attribute("timeWindow") = 0;
        } else {
            newMoment.append_attribute("time") = model->tdParams.moments[m - 1].time;
            newMoment.append_attribute("timeWindow") = model->tdParams.moments[m - 1].window;
        }

        if (m == 0) { //Write global results. Later these results will probably be time-dependent as well.
            xml_node globalNode = newMoment.append_child("Global");

            xml_node hitsNode = globalNode.append_child("Hits");
            hitsNode.append_attribute("totalHit") = globalState->globalStats.globalHits.nbMCHit;
            hitsNode.append_attribute("totalHitEquiv") = globalState->globalStats.globalHits.nbHitEquiv;
            hitsNode.append_attribute("totalDes") = globalState->globalStats.globalHits.nbDesorbed;
            hitsNode.append_attribute("totalAbsEquiv") = globalState->globalStats.globalHits.nbAbsEquiv;
            hitsNode.append_attribute("totalDist_total") = globalState->globalStats.distTraveled_total;
            hitsNode.append_attribute("totalDist_fullHitsOnly") = globalState->globalStats.distTraveledTotal_fullHitsOnly;
            hitsNode.append_attribute("totalLeak") = globalState->globalStats.nbLeakTotal;
            hitsNode.append_attribute("maxDesorption") = model->otfParams.desorptionLimit;

            xml_node hitCacheNode = globalNode.append_child("Hit_Cache");
            hitCacheNode.append_attribute("nb") = globalState->globalStats.hitCacheSize;

            for (size_t i = 0; i < globalState->globalStats.hitCacheSize; i++) {
                xml_node newHit = hitCacheNode.append_child("Hit");
                newHit.append_attribute("id") = i;
                newHit.append_attribute("posX") = globalState->globalStats.hitCache[i].pos.x;
                newHit.append_attribute("posY") = globalState->globalStats.hitCache[i].pos.y;
                newHit.append_attribute("posZ") = globalState->globalStats.hitCache[i].pos.z;
                newHit.append_attribute("type") = globalState->globalStats.hitCache[i].type;
            }

            xml_node leakCacheNode = globalNode.append_child("Leak_Cache");
            leakCacheNode.append_attribute("nb") = globalState->globalStats.leakCacheSize;
            for (size_t i = 0; i < globalState->globalStats.leakCacheSize; i++) {
                xml_node newLeak = leakCacheNode.append_child("Leak");
                newLeak.append_attribute("id") = i;
                newLeak.append_attribute("posX") = globalState->globalStats.leakCache[i].pos.x;
                newLeak.append_attribute("posY") = globalState->globalStats.leakCache[i].pos.y;
                newLeak.append_attribute("posZ") = globalState->globalStats.leakCache[i].pos.z;
                newLeak.append_attribute("dirX") = globalState->globalStats.leakCache[i].dir.x;
                newLeak.append_attribute("dirY") = globalState->globalStats.leakCache[i].dir.y;
                newLeak.append_attribute("dirZ") = globalState->globalStats.leakCache[i].dir.z;
            }
        } //end global node

        bool hasHistogram =
                model->sp.globalHistogramParams.recordBounce || model->sp.globalHistogramParams.recordDistance;
#ifdef MOLFLOW
        hasHistogram = hasHistogram || model->sp.globalHistogramParams.recordTime;
#endif
        if (hasHistogram) {
            xml_node histNode = newMoment.append_child("Histograms");
            //Retrieve histogram map from hits dp
            auto &globalHist = globalState->globalHistograms[m];
            if (model->sp.globalHistogramParams.recordBounce) {
                auto &nbHitsHistogram = globalHist.nbHitsHistogram;
                xml_node hist = histNode.append_child("Bounces");
                size_t histSize = model->sp.globalHistogramParams.GetBounceHistogramSize();
                hist.append_attribute("size") = histSize;
                hist.append_attribute(
                        "binSize") = model->sp.globalHistogramParams.nbBounceBinsize; //redundancy for human-reading or export
                hist.append_attribute(
                        "max") = model->sp.globalHistogramParams.nbBounceMax; //redundancy for human-reading or export
                for (size_t h = 0; h < histSize; h++) {
                    xml_node bin = hist.append_child("Bin");
                    auto value = bin.append_attribute("start");
                    if (h == histSize - 1) value = "overRange";
                    else value = h * model->sp.globalHistogramParams.nbBounceBinsize;
                    bin.append_attribute("count") = nbHitsHistogram[h];
                }
            }
            if (model->sp.globalHistogramParams.recordDistance) {
                auto &distanceHistogram = globalHist.distanceHistogram;
                xml_node hist = histNode.append_child("Distance");
                size_t histSize = model->sp.globalHistogramParams.GetDistanceHistogramSize();
                hist.append_attribute("size") = histSize;
                hist.append_attribute(
                        "binSize") = model->sp.globalHistogramParams.distanceBinsize; //redundancy for human-reading or export
                hist.append_attribute(
                        "max") = model->sp.globalHistogramParams.distanceMax; //redundancy for human-reading or export
                for (size_t h = 0; h < histSize; h++) {
                    xml_node bin = hist.append_child("Bin");
                    auto value = bin.append_attribute("start");
                    if (h == histSize - 1) value = "overRange";
                    else value = h * model->sp.globalHistogramParams.distanceBinsize;
                    bin.append_attribute("count") = distanceHistogram[h];
                }
            }
            if (model->sp.globalHistogramParams.recordTime) {
                auto &timeHistogram = globalHist.timeHistogram;
                xml_node hist = histNode.append_child("Time");
                size_t histSize = model->sp.globalHistogramParams.GetTimeHistogramSize();
                hist.append_attribute("size") = histSize;
                hist.append_attribute(
                        "binSize") = model->sp.globalHistogramParams.timeBinsize; //redundancy for human-reading or export
                hist.append_attribute(
                        "max") = model->sp.globalHistogramParams.timeMax; //redundancy for human-reading or export
                for (size_t h = 0; h < histSize; h++) {
                    xml_node bin = hist.append_child("Bin");
                    auto value = bin.append_attribute("start");
                    if (h == histSize - 1) value = "overRange";
                    else value = (double)h * model->sp.globalHistogramParams.timeBinsize;
                    bin.append_attribute("count") = timeHistogram[h];
                }
            }
        }

        xml_node facetResultsNode = newMoment.append_child("FacetResults");

        for (auto &fac : model->facets) {
            auto &sFac = *fac;
            //SimulationFacet& f = model->structures[0].facets[0].;
            xml_node newFacetResult = facetResultsNode.append_child("Facet");
            newFacetResult.append_attribute("id") = sFac.globalId;

            xml_node facetHitNode = newFacetResult.append_child("Hits");
            //FacetHitBuffer* facetCounter = (FacetHitBuffer *)(buffer + sFac.sh.hitOffset + m * sizeof(FacetHitBuffer));
            const auto &facetCounter = globalState->facetStates[sFac.globalId].momentResults[m].hits;

            facetHitNode.append_attribute("nbHit") = facetCounter.nbMCHit;
            facetHitNode.append_attribute("nbHitEquiv") = facetCounter.nbHitEquiv;
            facetHitNode.append_attribute("nbDes") = facetCounter.nbDesorbed;
            facetHitNode.append_attribute("nbAbsEquiv") = facetCounter.nbAbsEquiv;
            facetHitNode.append_attribute("sum_v_ort") = facetCounter.sum_v_ort;
            facetHitNode.append_attribute("sum_1_per_v") = facetCounter.sum_1_per_ort_velocity;
            facetHitNode.append_attribute("sum_v") = facetCounter.sum_1_per_velocity;

            if (model->sp.enableForceMeasurement) { //don't save all-zero quantities if not measured

                auto forcesNode = newFacetResult.append_child("Forces");

                auto impulseNode = forcesNode.append_child("Impulse");
                impulseNode.append_attribute("x") = facetCounter.impulse.x;
                impulseNode.append_attribute("y") = facetCounter.impulse.y;
                impulseNode.append_attribute("z") = facetCounter.impulse.z;

                auto impulse_square_Node = forcesNode.append_child("Impulse_square");
                impulse_square_Node.append_attribute("x") = facetCounter.impulse_square.x;
                impulse_square_Node.append_attribute("y") = facetCounter.impulse_square.y;
                impulse_square_Node.append_attribute("z") = facetCounter.impulse_square.z;

                auto impulse_momentum_Node = forcesNode.append_child("Impulse_momentum");
                impulse_momentum_Node.append_attribute("x") = facetCounter.impulse_momentum.x;
                impulse_momentum_Node.append_attribute("y") = facetCounter.impulse_momentum.y;
                impulse_momentum_Node.append_attribute("z") = facetCounter.impulse_momentum.z;

            }

            if (sFac.sh.isProfile) {
                xml_node profileNode = newFacetResult.append_child("Profile");
                profileNode.append_attribute("size") = PROFILE_SIZE;
                //ProfileSlice *pr = (ProfileSlice *)(buffer + sFac.sh.hitOffset + facetHitsSize + m * sizeof(ProfileSlice)*PROFILE_SIZE);
                const auto &pr = globalState->facetStates[sFac.globalId].momentResults[m].profile;

                for (int p = 0; p < PROFILE_SIZE; p++) {
                    xml_node slice = profileNode.append_child("Slice");
                    slice.append_attribute("id") = p;
                    slice.append_attribute("countEquiv") = pr[p].countEquiv;
                    slice.append_attribute("sum_1_per_v") = pr[p].sum_1_per_ort_velocity;
                    slice.append_attribute("sum_v_ort") = pr[p].sum_v_ort;
                }
            }

            //size_t profSize = (sFac.sh.isProfile) ? (PROFILE_SIZE * sizeof(ProfileSlice) * (1 + model->tdParams.moments.size())) : 0;
            size_t height = sFac.sh.texHeight;
            size_t width = sFac.sh.texWidth;

            if (sFac.sh.texWidth * sFac.sh.texHeight > 0) {
                xml_node textureNode = newFacetResult.append_child("Texture");
                textureNode.append_attribute("width") = sFac.sh.texWidth;
                textureNode.append_attribute("height") = sFac.sh.texHeight;

                //TextureCell *texture = (TextureCell *)(buffer + sFac.sh.hitOffset + facetHitsSize + profSize + m * w*h * sizeof(TextureCell));
                const auto &texture = globalState->facetStates[sFac.globalId].momentResults[m].texture;

                std::stringstream countText, sum1perText, sumvortText;
                countText << '\n'; //better readability in file
                sum1perText << std::setprecision(8) << '\n';
                sumvortText << std::setprecision(8) << '\n';

                for (size_t iy = 0; iy < height; iy++) {
                    for (size_t ix = 0; ix < width; ix++) {
                        countText << texture[iy * sFac.sh.texWidth + ix].countEquiv << '\t';
                        sum1perText << texture[iy * sFac.sh.texWidth + ix].sum_1_per_ort_velocity << '\t';
                        sumvortText << texture[iy * sFac.sh.texWidth + ix].sum_v_ort_per_area << '\t';
                    }
                    countText << '\n';
                    sum1perText << '\n';
                    sumvortText << '\n';
                }
                textureNode.append_child("count").append_child(node_cdata).set_value(countText.str().c_str());
                textureNode.append_child("sum_1_per_v").append_child(node_cdata).set_value(sum1perText.str().c_str());
                textureNode.append_child("sum_v_ort").append_child(node_cdata).set_value(sumvortText.str().c_str());

            } //end texture

            if (sFac.sh.countDirection/* && sFac.dirCache*/) {
                xml_node dirNode = newFacetResult.append_child("Directions");
                dirNode.append_attribute("width") = sFac.sh.texWidth;
                dirNode.append_attribute("height") = sFac.sh.texHeight;

                const auto &dirs = globalState->facetStates[sFac.globalId].momentResults[m].direction;

                std::stringstream dirText, dirCountText;
                dirText << std::setprecision(8) << '\n'; //better readability in file
                dirCountText << '\n';

                for (size_t iy = 0; iy < height; iy++) {
                    for (size_t ix = 0; ix < width; ix++) {
                        dirText << dirs[iy * sFac.sh.texWidth + ix].dir.x << ",";
                        dirText << dirs[iy * sFac.sh.texWidth + ix].dir.y << ",";
                        dirText << dirs[iy * sFac.sh.texWidth + ix].dir.z << "\t";
                        dirCountText << dirs[iy * sFac.sh.texWidth + ix].count << "\t";

                    }
                    dirText << "\n";
                    dirCountText << "\n";
                }
                dirNode.append_child("vel.vectors").append_child(node_cdata).set_value(dirText.str().c_str());
                dirNode.append_child("count").append_child(node_cdata).set_value(dirCountText.str().c_str());
            } //end directions

            //Facet histograms (1 per moment) comes here
            bool hasFHistogram =
                    sFac.sh.facetHistogramParams.recordBounce || sFac.sh.facetHistogramParams.recordDistance;
#ifdef MOLFLOW
            hasFHistogram = hasFHistogram || sFac.sh.facetHistogramParams.recordTime;
#endif
            if (hasFHistogram) {
                xml_node histNode = newFacetResult.append_child("Histograms");
                //Retrieve histogram map from hits dp
                auto &histogram = globalState->facetStates[sFac.globalId].momentResults[m].histogram;
                if (sFac.sh.facetHistogramParams.recordBounce) {
                    auto &nbHitsHistogram = histogram.nbHitsHistogram;
                    xml_node hist = histNode.append_child("Bounces");
                    size_t histSize = sFac.sh.facetHistogramParams.GetBounceHistogramSize();
                    hist.append_attribute("size") = histSize;
                    hist.append_attribute(
                            "binSize") = sFac.sh.facetHistogramParams.nbBounceBinsize; //redundancy for human-reading or export
                    hist.append_attribute(
                            "max") = sFac.sh.facetHistogramParams.nbBounceMax; //redundancy for human-reading or export
                    for (size_t h = 0; h < histSize; h++) {
                        xml_node bin = hist.append_child("Bin");
                        auto value = bin.append_attribute("start");
                        if (h == histSize - 1) value = "overRange";
                        else value = h * sFac.sh.facetHistogramParams.nbBounceBinsize;
                        bin.append_attribute("count") = nbHitsHistogram[h];
                    }
                }
                if (sFac.sh.facetHistogramParams.recordDistance) {
                    auto &distanceHistogram = histogram.distanceHistogram;
                    xml_node hist = histNode.append_child("Distance");
                    size_t histSize = sFac.sh.facetHistogramParams.GetDistanceHistogramSize();
                    hist.append_attribute("size") = histSize;
                    hist.append_attribute(
                            "binSize") = sFac.sh.facetHistogramParams.distanceBinsize; //redundancy for human-reading or export
                    hist.append_attribute(
                            "max") = sFac.sh.facetHistogramParams.distanceMax; //redundancy for human-reading or export
                    for (size_t h = 0; h < histSize; h++) {
                        xml_node bin = hist.append_child("Bin");
                        auto value = bin.append_attribute("start");
                        if (h == histSize - 1) value = "overRange";
                        else value = h * sFac.sh.facetHistogramParams.distanceBinsize;
                        bin.append_attribute("count") = distanceHistogram[h];
                    }
                }
                if (sFac.sh.facetHistogramParams.recordTime) {
                    auto &timeHistogram = histogram.timeHistogram;
                    xml_node hist = histNode.append_child("Time");
                    size_t histSize = sFac.sh.facetHistogramParams.GetTimeHistogramSize();
                    hist.append_attribute("size") = histSize;
                    hist.append_attribute(
                            "binSize") = sFac.sh.facetHistogramParams.timeBinsize; //redundancy for human-reading or export
                    hist.append_attribute(
                            "max") = sFac.sh.facetHistogramParams.timeMax; //redundancy for human-reading or export
                    for (size_t h = 0; h < histSize; h++) {
                        xml_node bin = hist.append_child("Bin");
                        auto value = bin.append_attribute("start");
                        if (h == histSize - 1) value = "overRange";
                        else value = (double)h * sFac.sh.facetHistogramParams.timeBinsize;
                        bin.append_attribute("count") = timeHistogram[h];
                    }
                }
            }
        }
    }
    return true;
}

/**
* \brief To save facet data for the geometry in XML
* \param facetNode XML node representing a facet
*/
void XmlWriter::SaveFacet(pugi::xml_node facetNode, std::shared_ptr<MolflowSimFacet> facet, size_t nbTotalVertices, const TimeDependentParameters& tdParams) {
    
    {   
        xml_node e = facetNode.append_child("Sticking");
        std::string paramName = facet->sh.stickingParam;
        if (!paramName.empty()) {
            e.append_attribute("stickingParam") = paramName.c_str();
            e.append_attribute("parameterId") = tdParams.GetParamId(paramName);
            e.append_attribute("constValue") = 0.0; //for backward compatibility, won't be read
        }
        else {
            e.append_attribute("constValue") = facet->sh.sticking;
            e.append_attribute("parameterId") = -1; //for backward compatibility, versions<2.9.15 expect it
        }
    }
    {
        xml_node e = facetNode.append_child("Opacity");
        std::string paramName = facet->sh.opacityParam;
        if (!paramName.empty()) {
            e.append_attribute("opacityParam") = paramName.c_str();
            e.append_attribute("parameterId") = tdParams.GetParamId(paramName);
            e.append_attribute("constValue") = 1.0; //for backward compatibility, won't be read
        }
        else {
            e.append_attribute("constValue") = facet->sh.opacity;
            e.append_attribute("parameterId") = -1; //for backward compatibility, versions<2.9.15 expect it
        }
        e.append_attribute("is2sided") = (int) facet->sh.is2sided; //backward compatibility: 0 or 1
    }

    xml_node e = facetNode.append_child("Outgassing");
    {
        
        std::string paramName = facet->sh.outgassingParam;
        if (!paramName.empty()) {
            e.append_attribute("outgassingParam") = paramName.c_str();
            e.append_attribute("parameterId") = tdParams.GetParamId(paramName);
            e.append_attribute("constValue") = 0.0; //for backward compatibility, won't be read
        }
        else {
            e.append_attribute("constValue") = facet->sh.outgassing;
            e.append_attribute("parameterId") = -1; //for backward compatibility, versions<2.9.15 expect it
        }
        e.append_attribute("is2sided") = (int)facet->sh.is2sided; //backward compatibility: 0 or 1
    }

    e.append_attribute("desType") = facet->sh.desorbType;
    e.append_attribute("desExponent") = facet->sh.desorbTypeN;
    e.append_attribute(
            "hasOutgassingFile") = (int) facet->sh.useOutgassingFile;//TODO:(int)hasOutgassingFile; //backward compatibility: 0 or 1
    e.append_attribute("useOutgassingFile") = (int) facet->sh.useOutgassingFile; //backward compatibility: 0 or 1

    e = facetNode.append_child("Temperature");
    if (!facet->sh.temperatureParam.empty()) { //Introduced in 2.9.15, no backward compatibility
        e.append_attribute("tempParam") = facet->sh.temperatureParam.c_str();
    }
    else {
        e.append_attribute("value") = facet->sh.temperature;
    }
    e.append_attribute("accFactor") = facet->sh.accomodationFactor;

    e = facetNode.append_child("Reflection");

    e.append_attribute("diffusePart") = facet->sh.reflection.diffusePart;
    e.append_attribute("specularPart") = facet->sh.reflection.specularPart;
    e.append_attribute("cosineExponent") = facet->sh.reflection.cosineExponent;

    //For backward compatibility
    if (facet->sh.reflection.diffusePart > 0.99) {
        e.append_attribute("type") = REFLECTION_DIFFUSE;
    } else if (facet->sh.reflection.specularPart > 0.99) {
        e.append_attribute("type") = REFLECTION_SPECULAR;
    } else
        e.append_attribute("type") = REFLECTION_UNIFORM;

    e.append_attribute("enableSojournTime") = (int) facet->sh.enableSojournTime; //backward compatibility: 0 or 1
    e.append_attribute("sojournFreq") = facet->sh.sojournFreq;
    e.append_attribute("sojournE") = facet->sh.sojournE;

    e = facetNode.append_child("Structure");
    e.append_attribute("inStructure") = facet->sh.superIdx;
    e.append_attribute("linksTo") = facet->sh.superDest;

    e = facetNode.append_child("Teleport");
    e.append_attribute("target") = facet->sh.teleportDest;

    e = facetNode.append_child("Motion");
    e.append_attribute("isMoving") = (int) facet->sh.isMoving; //backward compatibility: 0 or 1

    e = facetNode.append_child("Recordings");
    xml_node t = e.append_child("Profile");
    t.append_attribute("type") = facet->sh.profileType;
    switch (facet->sh.profileType) {
        case 0:
            t.append_attribute("name") = "none";
            break;
        case 1:
            t.append_attribute("name") = "pressure u";
            break;
        case 2:
            t.append_attribute("name") = "pressure v";
            break;
        case 3:
            t.append_attribute("name") = "angular";
            break;
        case 4:
            t.append_attribute("name") = "speed";
            break;
        case 5:
            t.append_attribute("name") = "ortho.v";
            break;
        case 6:
            t.append_attribute("name") = "tang.v.";
            break;
    }
    t = e.append_child("Texture");
    //assert(!(cellPropertiesIds == NULL && (facet->sh.countAbs || facet->sh.countDes || facet->sh.countRefl || facet->sh.countTrans))); //Count texture on non-existent texture

    t.append_attribute("hasMesh") = (facet->sh.countAbs || facet->sh.countDes || facet->sh.countRefl ||
                                     facet->sh.countTrans);
    t.append_attribute("texDimX") = facet->sh.texWidth_precise;
    t.append_attribute("texDimY") = facet->sh.texHeight_precise;
    t.append_attribute("countDes") = (int) facet->sh.countDes; //backward compatibility: 0 or 1
    t.append_attribute("countAbs") = (int) facet->sh.countAbs; //backward compatibility: 0 or 1
    t.append_attribute("countRefl") = (int) facet->sh.countRefl; //backward compatibility: 0 or 1
    t.append_attribute("countTrans") = (int) facet->sh.countTrans; //backward compatibility: 0 or 1
    t.append_attribute("countDir") = (int) facet->sh.countDirection; //backward compatibility: 0 or 1
    t.append_attribute("countAC") = (int) facet->sh.countACD; //backward compatibility: 0 or 1

    if (facet->sh.anglemapParams.record) {
        t = e.append_child("IncidentAngleMap");
        t.append_attribute("record") = facet->sh.anglemapParams.record;
        t.append_attribute("phiWidth") = facet->sh.anglemapParams.phiWidth;
        t.append_attribute("thetaLimit") = facet->sh.anglemapParams.thetaLimit;
        t.append_attribute("thetaLowerRes") = facet->sh.anglemapParams.thetaLowerRes;
        t.append_attribute("thetaHigherRes") = facet->sh.anglemapParams.thetaHigherRes;
    }

    if (!interfaceSettings->facetSettings.empty()) {
        e = facetNode.append_child("ViewSettings");
        e.append_attribute("textureVisible") = interfaceSettings->facetSettings[facet->globalId].textureVisible; //backward compatibility: 0 or 1
        e.append_attribute("volumeVisible") = interfaceSettings->facetSettings[facet->globalId].volumeVisible; //backward compatibility: 0 or 1
    }

    facetNode.append_child("Indices").append_attribute("nb") = facet->sh.nbIndex;
    for (size_t i = 0; i < facet->sh.nbIndex; i++) {
        xml_node indice = facetNode.child("Indices").append_child("Indice");
        indice.append_attribute("id") = i;
        indice.append_attribute("vertex") = facet->indices[i];
    }

    if (facet->sh.useOutgassingFile) { // hasOutgassingFile
        xml_node textureNode = facetNode.append_child("DynamicOutgassing");
        textureNode.append_attribute("width") = facet->ogMap.outgassingMapWidth;
        textureNode.append_attribute("height") = facet->ogMap.outgassingMapHeight;
        textureNode.append_attribute("ratioU") = facet->ogMap.outgassingFileRatioU;
        textureNode.append_attribute("ratioV") = facet->ogMap.outgassingFileRatioV;
        textureNode.append_attribute("totalDose") = facet->ogMap.totalDose;
        textureNode.append_attribute("totalOutgassing") = facet->sh.totalOutgassing;
        textureNode.append_attribute("totalFlux") = facet->ogMap.totalFlux;

        std::stringstream outgText;
        outgText << std::setprecision(8);
        outgText << '\n'; //better readability in file
        for (int iy = 0; iy < facet->ogMap.outgassingMapHeight; iy++) {
            for (int ix = 0; ix < facet->ogMap.outgassingMapWidth; ix++) {
                outgText << facet->ogMap.outgassingMap[iy * facet->ogMap.outgassingMapWidth + ix] << '\t';
            }
            outgText << '\n';
        }
        textureNode.append_child("map").append_child(node_cdata).set_value(outgText.str().c_str());

    } //end texture

    if (!facet->angleMap.pdf.empty()) {
        xml_node textureNode = facetNode.append_child("IncidentAngleMap");
        textureNode.append_attribute("angleMapPhiWidth") = facet->sh.anglemapParams.phiWidth;
        textureNode.append_attribute("angleMapThetaLimit") = facet->sh.anglemapParams.thetaLimit;
        textureNode.append_attribute("angleMapThetaLowerRes") = facet->sh.anglemapParams.thetaLowerRes;
        textureNode.append_attribute("angleMapThetaHigherRes") = facet->sh.anglemapParams.thetaHigherRes;

        std::stringstream angleText;
        angleText << '\n'; //better readability in file
        for (int iy = 0;
             iy < (facet->sh.anglemapParams.thetaLowerRes + facet->sh.anglemapParams.thetaHigherRes); iy++) {
            for (int ix = 0; ix < facet->sh.anglemapParams.phiWidth; ix++) {
                angleText << facet->angleMap.pdf[iy * facet->sh.anglemapParams.phiWidth + ix] << '\t';
            }
            angleText << '\n';
        }
        textureNode.append_child("map").append_child(node_cdata).set_value(angleText.str().c_str());

    } //end angle map

    xml_node histNode = facetNode.append_child("Histograms");
    if (facet->sh.facetHistogramParams.recordBounce) {
        xml_node nbBounceNode = histNode.append_child("Bounces");
        nbBounceNode.append_attribute("binSize") = facet->sh.facetHistogramParams.nbBounceBinsize;
        nbBounceNode.append_attribute("max") = facet->sh.facetHistogramParams.nbBounceMax;
    }
    if (facet->sh.facetHistogramParams.recordDistance) {
        xml_node distanceNode = histNode.append_child("Distance");
        distanceNode.append_attribute("binSize") = facet->sh.facetHistogramParams.distanceBinsize;
        distanceNode.append_attribute("max") = facet->sh.facetHistogramParams.distanceMax;
    }
#ifdef MOLFLOW
    if (facet->sh.facetHistogramParams.recordTime) {
        xml_node timeNode = histNode.append_child("Time");
        timeNode.append_attribute("binSize") = facet->sh.facetHistogramParams.timeBinsize;
        timeNode.append_attribute("max") = facet->sh.facetHistogramParams.timeMax;
    }
#endif
}

XmlWriter::XmlWriter(bool useOldXMLFormat, bool updateRootNode) : useOldXMLFormat(useOldXMLFormat), updateRootNode(updateRootNode){
    interfaceSettings = std::make_unique<MolflowInterfaceSettings>();
}

void XmlWriter::WriteConvergenceValues(pugi::xml_document& saveDoc, const std::vector<std::vector<FormulaHistoryDatapoint>>& convergenceData) {
    //make sure interfaceSettings->userFormuals is set (for names)
    auto rootNode = GetRootNode(saveDoc);
    xml_node resultNode = rootNode.child("MolflowResults");
    if (!resultNode) {
        resultNode = rootNode.append_child("MolflowResults");
    }
    //Convergence results
    xml_node convNode = resultNode.append_child("Convergence");

    int formulaId = 0;
    for (const auto& formulaVec : convergenceData) {
        std::stringstream convText;
        convText << std::setprecision(8) << '\n';
        for (const auto& convVal : formulaVec) {
            convText << convVal.nbDes << "\t" << convVal.value << "\n";
        }
        xml_node newFormulaNode = convNode.append_child("ConvData");
        if (interfaceSettings->userFormulas.size() > formulaId) {
            newFormulaNode.append_attribute("Formula") = interfaceSettings->userFormulas[formulaId].expression.c_str();
            newFormulaNode.append_attribute("nbEntries") = formulaVec.size();
        }
        xml_node newConv = newFormulaNode.append_child(node_cdata);
        newConv.set_value(convText.str().c_str());
        formulaId++;
    }
}

void XmlWriter::CameraViewToXml(const CameraView& v, xml_node& targetViewNode) {
    targetViewNode.append_attribute("name") = v.name.c_str();
    targetViewNode.append_attribute("projMode") = v.projMode;
    targetViewNode.append_attribute("camAngleOx") = v.camAngleOx;
    targetViewNode.append_attribute("camAngleOy") = v.camAngleOy;
    targetViewNode.append_attribute("camAngleOz") = v.camAngleOz;
    targetViewNode.append_attribute("camDist") = v.camDist;
    targetViewNode.append_attribute("lightAngleOx") = v.lightAngleOx;
    targetViewNode.append_attribute("lightAngleOy") = v.lightAngleOy;
    targetViewNode.append_attribute("camOffset.x") = v.camOffset.x;
    targetViewNode.append_attribute("camOffset.y") = v.camOffset.y;
    targetViewNode.append_attribute("camOffset.z") = v.camOffset.z;
    targetViewNode.append_attribute("performXY") = v.performXY;
    targetViewNode.append_attribute("vLeft") = v.vLeft;
    targetViewNode.append_attribute("vRight") = v.vRight;
    targetViewNode.append_attribute("vTop") = v.vTop;
    targetViewNode.append_attribute("vBottom") = v.vBottom;
    auto clippingNode = targetViewNode.append_child("Clipping");
    clippingNode.append_attribute("enabled") = v.enableClipping;
    clippingNode.append_attribute("a") = v.clipPlane.a;
    clippingNode.append_attribute("b") = v.clipPlane.b;
    clippingNode.append_attribute("c") = v.clipPlane.c;
    clippingNode.append_attribute("d") = v.clipPlane.d;
}