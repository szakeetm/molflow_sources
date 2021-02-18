//
// Created by Pascal Baehr on 30.07.20.
//

#include "WriterXML.h"
#include "PugiXML/pugixml.hpp"
#include "versionId.h"
#include <iomanip>      // std::setprecision
#include <sstream>

using namespace FlowIO;
using namespace pugi;

constexpr short useOldXMLFormat = 0;

static double writeProgress = 0.0;
void setWriteProgress(double newProgress) {
    writeProgress = newProgress;
}

void reportWriteStatus(const std::string& statusString) {
    auto time_point = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(time_point);
    char s[256];
    struct tm * p = localtime(&now_c);
    strftime(s, 256, "%F_%T", p);
    printf("[%s] %s\n", s , statusString.c_str());
}

void WriterXML::SaveGeometry(const std::string& outputFileName, SimulationModel *model) {
    xml_document saveDoc;
    xml_node rootNode;
    if(useOldXMLFormat){
        rootNode = saveDoc;
    }
    else {
        rootNode = saveDoc.append_child("SimulationEnvironment");
        rootNode.append_attribute("type") = "molflow";
        rootNode.append_attribute("version") = appVersionId;
    }
    xml_node geomNode = rootNode.append_child("Geometry");

}

int WriterXML::SaveSimulationState(const std::string& outputFileName, SimulationModel *model, GlobalSimuState& globState) {
    xml_document saveDoc;
    xml_parse_result parseResult = saveDoc.load_file(outputFileName.c_str()); //parse xml file directly

    xml_node rootNode = saveDoc.child("SimulationEnvironment");
    rootNode.remove_child("MolflowResults"); // clear previous results to replace with new status

    xml_node resultNode = rootNode.append_child("MolflowResults");
    reportWriteStatus("Writing simulation results...");
    xml_node momentsNode = resultNode.append_child("Moments");
    momentsNode.append_attribute("nb") = model->tdParams.moments.size() + 1;
    size_t facetHitsSize = (1 + model->tdParams.moments.size()) * sizeof(FacetHitBuffer);
    
    for (size_t m = 0; m <= model->tdParams.moments.size(); m++) {
        setWriteProgress(0.5 + 0.5 * (double) m / (1.0 + (double) model->tdParams.moments.size()));
        xml_node newMoment = momentsNode.append_child("Moment");
        newMoment.append_attribute("id") = m;
        if (m == 0) {
            newMoment.append_attribute("time") = "Constant flow";
            newMoment.append_attribute("timeWindow") = 0;
        }
        else {
            newMoment.append_attribute("time") = model->tdParams.moments[m - 1].first;
            newMoment.append_attribute("timeWindow") = model->tdParams.moments[m - 1].second;
        }

        if (m == 0) { //Write global results. Later these results will probably be time-dependent as well.
            xml_node globalNode = newMoment.append_child("Global");

            xml_node hitsNode = globalNode.append_child("Hits");
            hitsNode.append_attribute("totalHit") = globState.globalHits.globalHits.hit.nbMCHit;
            hitsNode.append_attribute("totalHitEquiv") = globState.globalHits.globalHits.hit.nbHitEquiv;
            hitsNode.append_attribute("totalDes") = globState.globalHits.globalHits.hit.nbDesorbed;
            hitsNode.append_attribute("totalAbsEquiv") = globState.globalHits.globalHits.hit.nbAbsEquiv;
            hitsNode.append_attribute("totalDist_total") = globState.globalHits.distTraveled_total;
            hitsNode.append_attribute("totalDist_fullHitsOnly") = globState.globalHits.distTraveledTotal_fullHitsOnly;
            hitsNode.append_attribute("totalLeak") = globState.globalHits.nbLeakTotal;
            hitsNode.append_attribute("maxDesorption") = model->otfParams.desorptionLimit;

            xml_node hitCacheNode = globalNode.append_child("Hit_Cache");
            hitCacheNode.append_attribute("nb") = globState.globalHits.hitCacheSize;

            for (size_t i = 0; i < globState.globalHits.hitCacheSize; i++) {
                xml_node newHit = hitCacheNode.append_child("Hit");
                newHit.append_attribute("id") = i;
                newHit.append_attribute("posX") = globState.globalHits.hitCache[i].pos.x;
                newHit.append_attribute("posY") = globState.globalHits.hitCache[i].pos.y;
                newHit.append_attribute("posZ") = globState.globalHits.hitCache[i].pos.z;
                newHit.append_attribute("type") = globState.globalHits.hitCache[i].type;
            }

            xml_node leakCacheNode = globalNode.append_child("Leak_Cache");
            leakCacheNode.append_attribute("nb") = globState.globalHits.leakCacheSize;
            for (size_t i = 0; i < globState.globalHits.leakCacheSize; i++) {
                xml_node newLeak = leakCacheNode.append_child("Leak");
                newLeak.append_attribute("id") = i;
                newLeak.append_attribute("posX") = globState.globalHits.leakCache[i].pos.x;
                newLeak.append_attribute("posY") = globState.globalHits.leakCache[i].pos.y;
                newLeak.append_attribute("posZ") = globState.globalHits.leakCache[i].pos.z;
                newLeak.append_attribute("dirX") = globState.globalHits.leakCache[i].dir.x;
                newLeak.append_attribute("dirY") = globState.globalHits.leakCache[i].dir.y;
                newLeak.append_attribute("dirZ") = globState.globalHits.leakCache[i].dir.z;
            }
        } //end global node

        xml_node facetResultsNode = newMoment.append_child("FacetResults");

        for (size_t i = 0; i < model->sh.nbFacet; i++) {
            SubprocessFacet* subFac = nullptr;

            for(auto& fac : model->facets){
                if(i == fac.globalId){
                    subFac = &fac;
                    break;
                }
            }

            if(!subFac){
                std::cerr << "Facet "<<i<< " not found in model!" << std::endl;
                continue;
            }

            //SubprocessFacet& f = model->structures[0].facets[0].;
            xml_node newFacetResult = facetResultsNode.append_child("Facet");
            newFacetResult.append_attribute("id") = i;

            xml_node facetHitNode = newFacetResult.append_child("Hits");
            //FacetHitBuffer* facetCounter = (FacetHitBuffer *)(buffer + subFac->sh.hitOffset + m * sizeof(FacetHitBuffer));
            const auto& facetCounter = globState.facetStates[subFac->globalId].momentResults[m].hits;
            
            facetHitNode.append_attribute("nbHit") = facetCounter.hit.nbMCHit;
            facetHitNode.append_attribute("nbHitEquiv") = facetCounter.hit.nbHitEquiv;
            facetHitNode.append_attribute("nbDes") = facetCounter.hit.nbDesorbed;
            facetHitNode.append_attribute("nbAbsEquiv") = facetCounter.hit.nbAbsEquiv;
            facetHitNode.append_attribute("sum_v_ort") = facetCounter.hit.sum_v_ort;
            facetHitNode.append_attribute("sum_1_per_v") = facetCounter.hit.sum_1_per_ort_velocity;
            facetHitNode.append_attribute("sum_v") = facetCounter.hit.sum_1_per_velocity;

            if (subFac->sh.isProfile) {
                xml_node profileNode = newFacetResult.append_child("Profile");
                profileNode.append_attribute("size") = PROFILE_SIZE;
                //ProfileSlice *pr = (ProfileSlice *)(buffer + subFac->sh.hitOffset + facetHitsSize + m * sizeof(ProfileSlice)*PROFILE_SIZE);
                const auto& pr = globState.facetStates[subFac->globalId].momentResults[m].profile;

                for (int p = 0; p < PROFILE_SIZE; p++) {
                    xml_node slice = profileNode.append_child("Slice");
                    slice.append_attribute("id") = p;
                    slice.append_attribute("countEquiv") = pr[p].countEquiv;
                    slice.append_attribute("sum_1_per_v") = pr[p].sum_1_per_ort_velocity;
                    slice.append_attribute("sum_v_ort") = pr[p].sum_v_ort;
                }
            }

            size_t profSize = (subFac->sh.isProfile) ? (PROFILE_SIZE * sizeof(ProfileSlice)*(1 + model->tdParams.moments.size())) : 0;
            size_t h = subFac->sh.texHeight;
            size_t w = subFac->sh.texWidth;

            if (subFac->sh.texWidth * subFac->sh.texHeight > 0) {
                xml_node textureNode = newFacetResult.append_child("Texture");
                textureNode.append_attribute("width") = subFac->sh.texWidth;
                textureNode.append_attribute("height") = subFac->sh.texHeight;

                //TextureCell *texture = (TextureCell *)(buffer + subFac->sh.hitOffset + facetHitsSize + profSize + m * w*h * sizeof(TextureCell));
                const auto& texture = globState.facetStates[subFac->globalId].momentResults[m].texture;

                std::stringstream countText, sum1perText, sumvortText;
                countText << '\n'; //better readability in file
                sum1perText << std::setprecision(8) << '\n';
                sumvortText << std::setprecision(8) << '\n';

                for (size_t iy = 0; iy < h; iy++) {
                    for (size_t ix = 0; ix < w; ix++) {
                        countText << texture[iy*subFac->sh.texWidth + ix].countEquiv << '\t';
                        sum1perText << texture[iy*subFac->sh.texWidth + ix].sum_1_per_ort_velocity << '\t';
                        sumvortText << texture[iy*subFac->sh.texWidth + ix].sum_v_ort_per_area << '\t';
                    }
                    countText << '\n';
                    sum1perText << '\n';
                    sumvortText << '\n';
                }
                textureNode.append_child("count").append_child(node_cdata).set_value(countText.str().c_str());
                textureNode.append_child("sum_1_per_v").append_child(node_cdata).set_value(sum1perText.str().c_str());
                textureNode.append_child("sum_v_ort").append_child(node_cdata).set_value(sumvortText.str().c_str());

            } //end texture

            if (subFac->sh.countDirection/* && subFac->dirCache*/) {
                xml_node dirNode = newFacetResult.append_child("Directions");
                dirNode.append_attribute("width") = subFac->sh.texWidth;
                dirNode.append_attribute("height") = subFac->sh.texHeight;

                //DirectionCell *dirs = (DirectionCell *)(buffer + subFac->sh.hitOffset + facetHitsSize + profSize + (1 + (int)model->tdParams.moments.size())*w*h * sizeof(TextureCell) + m * w*h * sizeof(DirectionCell));
                const auto& dirs = globState.facetStates[subFac->globalId].momentResults[m].direction;

                std::stringstream dirText, dirCountText;
                dirText << std::setprecision(8) << '\n'; //better readability in file
                dirCountText << '\n';

                for (size_t iy = 0; iy < h; iy++) {
                    for (size_t ix = 0; ix < w; ix++) {
                        dirText << dirs[iy*subFac->sh.texWidth + ix].dir.x << ",";
                        dirText << dirs[iy*subFac->sh.texWidth + ix].dir.y << ",";
                        dirText << dirs[iy*subFac->sh.texWidth + ix].dir.z << "\t";
                        dirCountText << dirs[iy*subFac->sh.texWidth + ix].count << "\t";

                    }
                    dirText << "\n";
                    dirCountText << "\n";
                }
                dirNode.append_child("vel.vectors").append_child(node_cdata).set_value(dirText.str().c_str());
                dirNode.append_child("count").append_child(node_cdata).set_value(dirCountText.str().c_str());
            } //end directions
        }
    }

    //Texture Min/Max
    xml_node minMaxNode = resultNode.append_child("TextureMinMax");
    minMaxNode.append_child("With_constant_flow").append_child("Pressure").append_attribute("min") = 0.0;
    minMaxNode.child("With_constant_flow").child("Pressure").append_attribute("max") = 0.0;
    minMaxNode.child("With_constant_flow").append_child("Density").append_attribute("min") = 0.0;
    minMaxNode.child("With_constant_flow").child("Density").append_attribute("max") = 0.0;
    minMaxNode.child("With_constant_flow").append_child("Imp.rate").append_attribute("min") = 0.0;
    minMaxNode.child("With_constant_flow").child("Imp.rate").append_attribute("max") = 0.0;

    minMaxNode.append_child("Moments_only").append_child("Pressure").append_attribute("min") = 0.0;
    minMaxNode.child("Moments_only").child("Pressure").append_attribute("max") = 0.0;
    minMaxNode.child("Moments_only").append_child("Density").append_attribute("min") = 0.0;
    minMaxNode.child("Moments_only").child("Density").append_attribute("max") = 0.0;
    minMaxNode.child("Moments_only").append_child("Imp.rate").append_attribute("min") = 0.0;
    minMaxNode.child("Moments_only").child("Imp.rate").append_attribute("max") = 0.0;

    if (!saveDoc.save_file(outputFileName.c_str()))
        throw Error("Error writing XML file."); //successful save
    return 0;
}