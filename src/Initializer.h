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
    extern size_t nbThreads;
    extern uint64_t simDuration;
    extern uint64_t autoSaveDuration;
    extern bool loadAutosave;
    extern std::list<uint64_t> desLimit;
    extern bool resetOnStart;
    extern std::string inputFile;
    extern std::string outputFile;
    extern std::string paramFile;
}

class Initializer {
private:
    static int parseCommands(int argc, char** argv);
    static int loadFromXML(SimulationManager *simManager, SimulationModel *model, GlobalSimuState *globState,
                           bool loadState);
    static int setSharedStorage();
    static int initDesLimit(SimulationModel& model, GlobalSimuState& globState);
    static int initSimModel(SimulationModel* model);
public:
    static std::string getAutosaveFile();
    static int init(int argc, char **argv, SimulationManager *simManager, SimulationModel *model,
                    GlobalSimuState *globState);

    static int initSimUnit(SimulationManager *simManager, SimulationModel *model, GlobalSimuState *globState);
};


#endif //MOLFLOW_PROJ_INITIALIZER_H
