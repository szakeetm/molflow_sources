//
// Created by Pascal Baehr on 01.08.20.
//

#ifndef MOLFLOW_PROJ_INITIALIZERGPU_H
#define MOLFLOW_PROJ_INITIALIZERGPU_H

#include <list>
#include <string>
#include <SimulationManager.h>
class MolflowSimulationModel;

namespace Settings {
    extern size_t nbThreads;
    extern uint64_t simDuration;
    extern uint64_t outputDuration;
    extern uint64_t autoSaveDuration;
    extern bool loadAutosave;
    extern std::list<uint64_t> desLimit;
    extern bool resetOnStart;
    extern std::string paramFile;
    extern std::vector<std::string> paramSweep;
}

class InitializerGPU {
private:
    static int parseCommands(int argc, char** argv);
    static int loadFromXML(const std::string &fileName, bool loadState, const std::shared_ptr<MolflowSimulationModel>& model,
                           GlobalSimuState *globState);
    static int initSimModel(const std::shared_ptr<MolflowSimulationModel>& model);
public:
    static std::string getAutosaveFile();
    static int initFromFile(SimulationManager *simManager, const std::shared_ptr<MolflowSimulationModel>& model, GlobalSimuState *globState);
    static int initAutoGenerated(SimulationManager *simManager, const std::shared_ptr<MolflowSimulationModel> &model,
                                 GlobalSimuState *globState, double ratio, int steps, double angle);
    static int initDesLimit(const std::shared_ptr<MolflowSimulationModel>& model, GlobalSimuState &globState);

    static int initFromArgv(int argc, char **argv, SimulationManager *simManager, const std::shared_ptr<MolflowSimulationModel>& model);

    static int initTimeLimit(const std::shared_ptr<MolflowSimulationModel> &model, double time);

    static int
    loadFromGeneration(const std::shared_ptr<MolflowSimulationModel> &model, GlobalSimuState *globState, double ratio,
                       int step, double angle);
};


#endif //MOLFLOW_PROJ_INITIALIZERGPU_H
