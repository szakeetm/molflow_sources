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
    extern std::vector<std::string> paramSweep;
    extern std::string outputPath;
}

class Initializer {
private:
    static int parseCommands(int argc, char** argv);
    static int loadFromXML(const std::string &fileName, bool loadState, SimulationModel *model,
                           GlobalSimuState *globState);
    static int setSharedStorage();
    static int initSimModel(SimulationModel* model);
public:
    static std::string getAutosaveFile();
    static int initFromFile(SimulationManager *simManager, SimulationModel *model, GlobalSimuState *globState);
    static int initDesLimit(SimulationModel& model, GlobalSimuState& globState);

    static int initFromArgv(int argc, char **argv, SimulationManager *simManager, SimulationModel *model);
};


#endif //MOLFLOW_PROJ_INITIALIZER_H
