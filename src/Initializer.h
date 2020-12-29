//
// Created by Pascal Baehr on 01.08.20.
//

#ifndef MOLFLOW_PROJ_INITIALIZER_H
#define MOLFLOW_PROJ_INITIALIZER_H

#include <list>
#include <string>
#include <SimulationManager.h>
#include "GeometrySimu.h"

namespace Settings {
    extern SimulationManager simManager;
    extern double nbCPUCores;
    extern size_t nbThreadsPerCore;
    extern uint64_t simDuration;
    extern uint64_t autoSaveDuration;
    extern std::list<uint64_t> desLimit;
    extern bool resetOnStart;
    extern std::string req_real_file;
}

class Initializer {
private:
    static int parseCommands(int argc, char** argv);
    static int loadFromXML(SimulationManager *simManager, SimulationModel *model, GlobalSimuState *globState,
                           bool loadState);
    static int setSharedStorage();

public:
    static int init(int argc, char **argv, SimulationManager *simManager, SimulationModel *model,
                    GlobalSimuState *globState);
};


#endif //MOLFLOW_PROJ_INITIALIZER_H
