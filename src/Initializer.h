

#pragma once

#include <list>
#include <string>
#include <SimulationManager.h>
#include "Simulation/MolflowSimGeom.h"
#include "SettingsIO.h"

class Initializer {
private:
    static SettingsIO::CLIArguments parseArguments(int argc, char** argv);
    static std::shared_ptr<MolflowSimulationModel> CLILoadFromXML(const std::string &fileName, bool loadSimulationState, const std::shared_ptr<MolflowSimulationModel> model,
        const std::shared_ptr<GlobalSimuState> globalState, MolflowInterfaceSettings& interfaceSettings, SettingsIO::CLIArguments& parsedArgs);
    static void initSimModel(const std::shared_ptr<MolflowSimulationModel> model);
public:
    static std::string getAutosaveFile(SettingsIO::CLIArguments& parsedArgs);
    static std::shared_ptr<MolflowSimulationModel> initFromFile(SimulationManager& simManager, const std::shared_ptr<MolflowSimulationModel> model, const std::shared_ptr<GlobalSimuState> globalState, MolflowInterfaceSettings& interfaceSettings, SettingsIO::CLIArguments& parsedArgs);
    
    static int initDesLimit(const std::shared_ptr<MolflowSimulationModel> model, const std::shared_ptr<GlobalSimuState> globalState, SettingsIO::CLIArguments& parsedArgs);

    [[nodiscard]] static SettingsIO::CLIArguments initFromArgv(int argc, char **argv, SimulationManager& simManager, const std::shared_ptr<MolflowSimulationModel> model);

    static int initTimeLimit(const std::shared_ptr<MolflowSimulationModel> model, SettingsIO::CLIArguments& parsedArgs);
    /*
    static int
    loadFromGeneration(const std::shared_ptr<MolflowSimulationModel> model, const std::shared_ptr<GlobalSimuState> globalState, double ratio,
                       int step, double angle);
                       */
};
