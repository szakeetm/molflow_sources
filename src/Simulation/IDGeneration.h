

#ifndef MOLFLOW_PROJ_IDGENERATION_H
#define MOLFLOW_PROJ_IDGENERATION_H

#include <vector>
#include "MolflowSimGeom.h"
#include <set>

namespace IDGeneration {
    // ID
    int GetIDId(const std::set<size_t>& desorptionParameterIDs, int paramId);
    std::pair<int, std::vector<IntegratedDesorptionEntry>> GenerateNewID(std::set<size_t>& desorptionParameterIDs, int paramId, MolflowSimulationModel* model);
    std::vector<IntegratedDesorptionEntry> Generate_ID(int paramId, MolflowSimulationModel* model);
};


#endif //MOLFLOW_PROJ_IDGENERATION_H
