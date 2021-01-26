//
// Created by pascal on 1/25/21.
//

#include "ParameterParser.h"

#include <fstream>
#include <sstream>
#include <vector>
#include <tuple>

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

    std::vector<std::tuple<size_t, FacetParam, double>> facetParams;
    std::vector<std::tuple<SimuParam, double>> simuParams;
}

void parseFacet(std::istringstream& facetString){
    std::string id_str;
    std::string param_str;
    std::string paramVal_str;
    std::getline(facetString, id_str, '.');
    std::getline(facetString, param_str, '=');
    std::getline(facetString, paramVal_str);
    size_t id= std::strtoul(id_str.c_str(),nullptr,10);
    auto param = (Parameters::FacetParam) std::strtoul(param_str.c_str(),nullptr,10);
    double paramVal = std::strtod(paramVal_str.c_str(),nullptr);
    Parameters::facetParams.emplace_back(std::make_tuple(id, param, paramVal));
    //printf("[Facet #%s] %s = %s\n", id.c_str(), param.c_str(), paramVal.c_str());
}

void parseSimu(std::istringstream& facetString){
    std::string param_str;
    std::string paramVal_str;
    std::getline(facetString, param_str, '=');
    std::getline(facetString, paramVal_str);
    auto param = (Parameters::SimuParam) std::strtoul(param_str.c_str(),nullptr,10);
    double paramVal = std::strtod(paramVal_str.c_str(),nullptr);
    Parameters::simuParams.emplace_back(std::make_tuple(param, paramVal));
    //printf("[Facet #%s] %s = %s\n", id.c_str(), param.c_str(), paramVal.c_str());
}

void parseFacet(const std::string& facetString){
    auto token = facetString.find('.');
    std::string id = facetString.substr(0, token); // token is "scott"
    auto tokenEq = facetString.find('=');
    std::string param = facetString.substr(token+1, facetString.size()-token); // token is "scott"
    std::string paramVal = facetString.substr(tokenEq+1); // token is "scott"
    printf("[Facet #%s] %s = %s\n", id.c_str(), param.c_str(), paramVal.c_str());
}

void ParameterParser::Parse(const std::string& paramFile){
    //convDistr=std::vector<std::pair<double,double>>();

    std::ifstream inputFileStream(paramFile);

    size_t i = 0;
    for (std::string line; std::getline(inputFileStream, line); ) {
        std::istringstream lineStream(line);
        std::string optionType;
        std::getline(lineStream, optionType, '.');

        if( optionType == "facet" ) {
            parseFacet(lineStream);
            printf("[%zu] Parsing %s\n", i, lineStream.str().c_str());
        }
        else {
            printf("[%zu] Unknown input\n", i);
        }
        ++i;
    }
}

void ParameterParser::ChangeSimuParams(WorkerParams& params){
    for(auto& par : Parameters::simuParams){
        if(std::get<0>(par) == Parameters::SimuParam::enableDecay)
            params.enableDecay = (std::get<1>(par) > 0.0);
        else if(std::get<0>(par) == Parameters::SimuParam::mass)
            params.gasMass = std::get<1>(par);
        else if(std::get<0>(par) == Parameters::SimuParam::halfLife)
            params.halfLife = std::get<1>(par);
    }
}

void ParameterParser::ChangeFacetParams(std::vector<SubprocessFacet> &facets) {
    for(auto& par : Parameters::facetParams){
        size_t id = std::get<0>(par);
        if(id < facets.size()) {
            auto& facet = facets[id];
            if (std::get<1>(par) == Parameters::FacetParam::opacity)
                facet.sh.opacity = std::get<2>(par);
            else if (std::get<1>(par) == Parameters::FacetParam::outgassing)
                facet.sh.outgassing = std::get<2>(par);
            else if (std::get<1>(par) == Parameters::FacetParam::sticking)
                facet.sh.sticking = std::get<2>(par);
            else if (std::get<1>(par) == Parameters::FacetParam::temperature)
                facet.sh.temperature = std::get<2>(par);
        }
    }
}
