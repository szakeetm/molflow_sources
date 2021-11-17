//
// Created by pascal on 2/9/21.
//

#ifndef MOLFLOW_PROJ_IDGENERATION_H
#define MOLFLOW_PROJ_IDGENERATION_H

#include <vector>
#include <GeometrySimu.h>
#include <set>

namespace IDGeneration {
    // ID
    int GetIDId(const std::set<size_t>& desorptionParameterIDs, int paramId);
    std::pair<int, std::vector<ID_p>> GenerateNewID( std::set<size_t>& desorptionParameterIDs, int paramId, SimulationModel* model);
    std::vector<std::pair<double, double>> Generate_ID(int paramId, SimulationModel *model);
};


#endif //MOLFLOW_PROJ_IDGENERATION_H
