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

#ifndef MOLFLOW_PROJ_IDGENERATION_H
#define MOLFLOW_PROJ_IDGENERATION_H

#include <vector>
#include "MolflowSimGeom.h"
#include <set>

namespace IDGeneration {
    // ID
    int GetIDId(const std::set<size_t>& desorptionParameterIDs, int paramId);
    std::pair<int, std::vector<DesorptionEntry>> GenerateNewID(std::set<size_t>& desorptionParameterIDs, int paramId, MolflowSimulationModel* model);
    std::vector<IntegratedDesorptionEntry> Generate_ID(int paramId, MolflowSimulationModel* model);
};


#endif //MOLFLOW_PROJ_IDGENERATION_H
