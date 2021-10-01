//
// Created by pascal on 1/25/21.
//

#include "ParameterParser.h"
#include "Helper/StringHelper.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <tuple>
#include <Helper/ConsoleLogger.h>

namespace Parameters {

    enum FacetParam : size_t {
        opacity,
        temperature,
        sticking,
        outgassing
    };

    enum SimuParam : size_t {
        mass,
        enableDecay,
        halfLife
    };

    static std::unordered_map<std::string,FacetParam> const tableFac = {
            {"opacity",FacetParam::opacity},
            {"temperature",FacetParam::temperature},
            {"sticking",FacetParam::sticking},
            {"outgassing",FacetParam::outgassing}
    };
    static std::unordered_map<std::string,SimuParam> const tableSim = {
            {"mass",SimuParam::mass},
            {"enableDecay",SimuParam::enableDecay},
            {"halfLife",SimuParam::halfLife}
    };
    std::vector<std::tuple<size_t, FacetParam, double>> facetParams;
    std::vector<std::tuple<SimuParam, double>> simuParams;
}

void parseFacet(std::istringstream &facetString, const std::vector<SelectionGroup> &selections) {
    std::string id_str;
    std::string param_str;
    std::string paramVal_str;
    std::getline(facetString, id_str, '.');
    std::getline(facetString, param_str, '=');
    std::getline(facetString, paramVal_str);
    std::vector<size_t> id_range;

    bool fromSelection = false; // for corresponding output message
    if(id_str.find('\"') != std::string::npos){ // selection
        id_str.erase (std::remove(id_str.begin(), id_str.end(), '\"'), id_str.end());
        for(const auto& sel : selections){
            if(sel.name == id_str) {
                id_range = sel.selection;
                fromSelection = true;
                break;
            }
        }
    }
    else { // facet id or range
        // facet list returns indices [0,inf], input is given by [1,inf]

        try {
            // For now get facet list for all combinations (with 3 parameter), check for valid ids later
            splitFacetList(id_range, id_str, 1e7);
        } catch (std::exception& e) {
            Log::console_error("[%s] Could not parse facet id or range:\n", __FUNCTION__);
            Log::console_error("\t%s\n", id_str.c_str());
        }
    }
    auto tablePair = Parameters::tableFac.find(param_str);
    if(tablePair == Parameters::tableFac.end()) {
        Log::console_error("[%s] Invalid option was given:\n", __FUNCTION__);
        Log::console_error("\t%s\n", param_str.c_str());
        return;
    }
    auto param = tablePair->second;
    double paramVal = std::strtod(paramVal_str.c_str(),nullptr);
    for(auto& id : id_range)
        Parameters::facetParams.emplace_back(std::make_tuple(id, param, paramVal));
    if(fromSelection)
        Log::console_msg_master(3, "[ParameterChange][Facet][Group: %s] Changing parameter %s to %s\n", id_str.c_str(), param_str.c_str(), paramVal_str.c_str());
    else
        Log::console_msg_master(3, "[ParameterChange][Facet][ID: %s] Changing parameter %s to %s\n", id_str.c_str(), param_str.c_str(), paramVal_str.c_str());
}

void parseSimu(std::istringstream& facetString){
    std::string param_str;
    std::string paramVal_str;
    std::getline(facetString, param_str, '=');
    std::getline(facetString, paramVal_str);
    auto tablePair = Parameters::tableSim.find(param_str);
    if(tablePair == Parameters::tableSim.end()) {
        Log::console_error("[%s] Invalid option was given:\n", __FUNCTION__);
        Log::console_error("\t%s\n", param_str.c_str());
        return;
    }
    auto param = tablePair->second;
    double paramVal = std::strtod(paramVal_str.c_str(),nullptr);
    Parameters::simuParams.emplace_back(std::make_tuple(param, paramVal));

    Log::console_msg_master(3, "[ParameterChange][Simulation] Changing parameter %s to %s\n", param_str.c_str(), paramVal_str.c_str());

}

void parseFacet(const std::string& facetString){
    auto token = facetString.find('.');
    std::string id = facetString.substr(0, token); // token is "scott"
    auto tokenEq = facetString.find('=');
    std::string param = facetString.substr(token+1, facetString.size()-token); // token is "scott"
    std::string paramVal = facetString.substr(tokenEq+1); // token is "scott"
    Log::console_msg_master(3,"[Facet #%s] %s = %s\n", id.c_str(), param.c_str(), paramVal.c_str());
}

void parseInputStream(std::stringstream& inputLineStream, const std::vector<SelectionGroup> &selections){
    size_t i = 0;
    for (std::string line; inputLineStream >> line; ) {
        std::istringstream lineStream(line);
        std::string optionType;
        std::getline(lineStream, optionType, '.');

        if (optionType == "facet") {
            parseFacet(lineStream, selections);
            //printf("[%zu] Parsing %s\n", i, lineStream.str().c_str());
        } else if (optionType == "simulation") {
            parseSimu(lineStream);
        } else {
            Log::console_error("[Line #%zu] Unknown input %s\n", i, line.c_str());
        }
        ++i;
    }
}


void ParameterParser::ParseFile(const std::string &paramFile, const std::vector<SelectionGroup> &selections) {
    //convDistr=std::vector<std::pair<double,double>>();

    std::ifstream inputFileStream(paramFile);
    std::stringstream inputStream;

    copy(std::istreambuf_iterator<char>(inputFileStream),
         std::istreambuf_iterator<char>(),
         std::ostreambuf_iterator<char>(inputStream));

    parseInputStream(inputStream, selections);

}

void ParameterParser::ParseInput(const std::vector<std::string> &paramSweep, const std::vector<SelectionGroup> &selections) {
    std::stringstream inputStream(std::ios_base::app | std::ios_base::out | std::ios_base::in);
    const char delimiter = ';';
    for(auto sweep : paramSweep){
        size_t pos = 0;
        std::string token;
        while ((pos = sweep.find(delimiter)) != std::string::npos) {
            token = sweep.substr(0, pos);
            inputStream << token << '\n';
            sweep.erase(0, pos + 1);
        }
        inputStream << sweep << '\n';
    }

    parseInputStream(inputStream, selections);
}

void ParameterParser::ChangeSimuParams(WorkerParams& params){
    for(auto& par : Parameters::simuParams){
        auto type = std::get<0>(par);
        switch (type) {
            case (Parameters::SimuParam::enableDecay):
                params.enableDecay = (std::get<1>(par) > 0.0);
                break;
            case (Parameters::SimuParam::mass):
                params.gasMass = std::get<1>(par);
                break;
            case (Parameters::SimuParam::halfLife):
                params.halfLife = std::get<1>(par);
                break;
            default:
                Log::console_error("Unknown SimuParam %s\n", std::get<0>(par));
        }
    }
}

void ParameterParser::ChangeFacetParams(std::vector<std::shared_ptr<SubprocessFacet>> &facets) {
    for(auto& par : Parameters::facetParams){
        size_t id = std::get<0>(par);
        if(id < facets.size()) {
            auto& facet = *facets.at(id);
            auto type = std::get<1>(par);
            switch (type) {
                case (Parameters::FacetParam::opacity):
                    facet.sh.opacity = std::get<2>(par);
                    break;
                case (Parameters::FacetParam::outgassing):
                    facet.sh.outgassing = std::get<2>(par);
                    break;
                case (Parameters::FacetParam::sticking):
                    facet.sh.sticking = std::get<2>(par);
                    break;
                case (Parameters::FacetParam::temperature):
                    facet.sh.temperature = std::get<2>(par);
                    break;
                default:
                    Log::console_error("Unknown FacetParam %s\n", std::get<1>(par));
            }
        }
    }
}
