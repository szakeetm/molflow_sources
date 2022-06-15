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

#ifndef MOLFLOW_PROJ_PARAMETERPARSER_H
#define MOLFLOW_PROJ_PARAMETERPARSER_H


#include <string>
#include <Buffer_shared.h>
#include "Simulation/MolflowSimGeom.h"
#include "GeometryTypes.h"

class ParameterParser {
public:
    static void ParseFile(const std::string &paramFile, const std::vector<SelectionGroup> &selections);

    static void ParseInput(const std::vector<std::string> &paramSweep, const std::vector<SelectionGroup> &selections);

    static void ChangeSimuParams(WorkerParams& params);

    static int ChangeFacetParams(std::vector<std::shared_ptr<SimulationFacet>> &facets);
};


#endif //MOLFLOW_PROJ_PARAMETERPARSER_H
