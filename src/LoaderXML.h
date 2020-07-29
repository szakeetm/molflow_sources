//
// Created by Pascal Baehr on 20.07.20.
//

#ifndef MOLFLOW_PROJ_LOADERXML_H
#define MOLFLOW_PROJ_LOADERXML_H

#include "GeometrySimu.h"
#include "PugiXML/pugixml.hpp"

namespace MFLoad {
    class Loader {
        std::vector<std::string> parameterList;

        std::set<double> temperatureList;
        std::set<size_t> desorptionParameterIDs; //time-dependent parameters which are used as desorptions, therefore need to be integrated

        std::vector<std::vector<std::pair<double, double>>> CDFs; //cumulative distribution function for each temperature
        std::vector<std::vector<std::pair<double, double>>> IDs; //integrated distribution function for each time-dependent desorption type

    protected:
        void PrepareToRun(SimulationModel *model);
        void CalcTotalOutgassing(SimulationModel* model);

        // CDF
        int GetCDFId(double temperature);
        int GenerateNewCDF(double temperature, double gasMass);
        static std::vector<std::pair<double, double>> Generate_CDF(double gasTempKelvins, double gasMassGramsPerMol, size_t size);

        // ID
        int GetIDId(int paramId);
        int GenerateNewID(int paramId);
        static std::vector<std::pair<double, double>> Generate_ID(int paramId, SimulationModel *model);

    public:
        virtual void LoadGeometry(std::string inputFileName, SimulationModel *model) = 0;
        virtual void LoadSimulationState() = 0;
    };

    class LoaderXML : public Loader {

    protected:
        void LoadFacet(pugi::xml_node facetNode, SubprocessFacet *facet, size_t nbTotalVertices);
    public:
        void LoadGeometry(std::string inputFileName, SimulationModel *model);
        void LoadSimulationState();
    };
}

#endif //MOLFLOW_PROJ_LOADERXML_H
