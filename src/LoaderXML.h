//
// Created by Pascal Baehr on 20.07.20.
//

#ifndef MOLFLOW_PROJ_LOADERXML_H
#define MOLFLOW_PROJ_LOADERXML_H

#include <set>
#include <map>
#include "GeometrySimu.h"
#include "PugiXML/pugixml.hpp"

namespace FlowIO {
    class Loader {
        std::vector<std::string> parameterList;

        std::set<double> temperatureList;
        std::set<size_t> desorptionParameterIDs; //time-dependent parameters which are used as desorptions, therefore need to be integrated

    protected:
        void PrepareToRun(SimulationModel *model);
        void CalcTotalOutgassing(SimulationModel* model);

        // CDF
        int GetCDFId(double temperature);
        int GenerateNewCDF(double temperature, double gasMass);
        static std::vector<std::pair<double, double>> Generate_CDF(double gasTempKelvins, double gasMassGramsPerMol, size_t size);

        // ID
        int GetIDId(int paramId);
        int GenerateNewID(int paramId, SimulationModel* model);
        static std::vector<std::pair<double, double>> Generate_ID(int paramId, SimulationModel *model);

        std::vector<std::vector<std::pair<double, double>>> IDs;         //integrated distribution function for each time-dependent desorption type
        std::vector<std::vector<std::pair<double, double>>> CDFs;        //cumulative distribution function for each temperature
        std::map<size_t,std::vector<size_t>> angleMapCache;
    public:
        std::vector<SubprocessFacet> loadFacets;

        virtual int LoadGeometry(std::string inputFileName, SimulationModel *model) = 0;
        virtual int LoadSimulationState(std::string inputFileName, SimulationModel *model, BYTE* buffer) = 0;

        std::ostringstream SerializeForLoader(SimulationModel* model);
        std::ostringstream SerializeResultsForLoader(GlobalSimuState* globState);
        void MoveFacetsToStructures(SimulationModel* model);
    };

    class LoaderXML : public Loader {

    protected:
        void LoadFacet(pugi::xml_node facetNode, SubprocessFacet *facet, size_t nbTotalVertices);
    public:
        int LoadGeometry(std::string inputFileName, SimulationModel *model) override;
        int LoadSimulationState(std::string inputFileName, SimulationModel *model, BYTE* buffer) override;
        int LoadSimulationState(std::string inputFileName, SimulationModel *model, GlobalSimuState& globState);

    };
}

#endif //MOLFLOW_PROJ_LOADERXML_H
