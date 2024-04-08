

#ifndef MOLFLOW_PROJ_PARAMETERPARSER_H
#define MOLFLOW_PROJ_PARAMETERPARSER_H


#include <string>
#include <Buffer_shared.h>
#include "Simulation/MolflowSimGeom.h"
#include "GeometryTypes.h"

//! Tools for parsing CLI arguments for the parameter sweep
class ParameterParser {
public:
    static void ParseFile(const std::string &paramFile, const std::vector<SelectionGroup> &selections);

    static void ParseInput(const std::vector<std::string> &paramChanges, const std::vector<SelectionGroup> &selections);

    static void ChangeSimuParams(SimuParams& params);

    static int ChangeFacetParams(std::vector<std::shared_ptr<SimulationFacet>> facets);
};


#endif //MOLFLOW_PROJ_PARAMETERPARSER_H
