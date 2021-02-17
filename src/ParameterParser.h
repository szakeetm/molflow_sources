//
// Created by pascal on 1/25/21.
//

#ifndef MOLFLOW_PROJ_PARAMETERPARSER_H
#define MOLFLOW_PROJ_PARAMETERPARSER_H


#include <string>
#include <Buffer_shared.h>
#include "GeometrySimu.h"
#include "GeometryTypes.h"

class ParameterParser {
public:
    static void Parse(const std::string &paramFile, const std::vector<SelectionGroup> &selections);

    static void ChangeSimuParams(WorkerParams& params);

    static void ChangeFacetParams(std::vector<SubprocessFacet> &facets);
};


#endif //MOLFLOW_PROJ_PARAMETERPARSER_H
