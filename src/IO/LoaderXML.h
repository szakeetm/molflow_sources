//
// Created by Pascal Baehr on 20.07.20.
//

#ifndef MOLFLOW_PROJ_LOADERXML_H
#define MOLFLOW_PROJ_LOADERXML_H

#include <set>
#include <map>
#include <GeometryTypes.h>
#include "GeometrySimu.h"
#include "PugiXML/pugixml.hpp"

namespace FlowIO {
    class Loader {
    protected:
        std::vector<std::vector<std::pair<double, double>>> IDs;         //integrated distribution function for each time-dependent desorption type
        std::vector<std::vector<std::pair<double, double>>> CDFs;        //cumulative distribution function for each temperature
    public:
        virtual int LoadGeometry(std::string inputFileName, SimulationModel *model) = 0;
    };

    class LoaderXML : public Loader {

    protected:
        void LoadFacet(pugi::xml_node facetNode, SubprocessFacet *facet, size_t nbTotalVertices);
    public:
        int LoadGeometry(std::string inputFileName, SimulationModel *model) override;
        static std::vector<SelectionGroup> LoadSelections(const std::string& inputFileName);
        static int LoadSimulationState(const std::string& inputFileName, SimulationModel *model, GlobalSimuState& globState);

        UserInput uInput;
    };
}

#endif //MOLFLOW_PROJ_LOADERXML_H
