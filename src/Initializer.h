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
    extern std::string paramFile;
    extern std::vector<std::string> paramSweep;
}

class Initializer {
private:
    static int parseCommands(int argc, char** argv);
    static int loadFromXML(const std::string &fileName, bool loadState, std::shared_ptr<SimulationModel> model,
                           GlobalSimuState *globState);
    static int setSharedStorage();
    static int initSimModel(std::shared_ptr<SimulationModel> model);
public:
    static std::string getAutosaveFile();
    static int initFromFile(SimulationManager *simManager, std::shared_ptr<SimulationModel> model, GlobalSimuState *globState);
    static int initDesLimit(std::shared_ptr<SimulationModel> model, GlobalSimuState& globState);

    static int initFromArgv(int argc, char **argv, SimulationManager *simManager, std::shared_ptr<SimulationModel> model);
};


#endif //MOLFLOW_PROJ_INITIALIZER_H
