//
// Created by pascal on 5/17/21.
//

#ifndef MOLFLOW_PROJ_SETTINGSIO_MF_H
#define MOLFLOW_PROJ_SETTINGSIO_MF_H

#include "SettingsIO.h"

class GlobalSimuState;
class SimulationModel;

namespace SettingsIO {
    // export utility functions
    void export_facet_details(GlobalSimuState* glob, SimulationModel* model);
    void export_facet_quantities(GlobalSimuState* glob, SimulationModel* model);
}


#endif //MOLFLOW_PROJ_SETTINGSIO_MF_H
