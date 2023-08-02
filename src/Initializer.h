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

#ifndef MOLFLOW_PROJ_INITIALIZER_H
#define MOLFLOW_PROJ_INITIALIZER_H

#include <list>
#include <string>
#include <SimulationManager.h>
#include "Simulation/MolflowSimGeom.h"

namespace CLIArguments {
    extern size_t nbThreads;
    extern uint64_t simDuration;
    extern uint64_t outputInterval;
    extern uint64_t autoSaveInterval;
    extern bool loadAutosave;
    extern std::list<uint64_t> desLimit;
    extern bool resetOnStart;
    extern std::string paramFile;
    extern std::vector<std::string> paramChanges;
}

class Initializer {
private:
    static int parseCommands(int argc, char** argv);
    static std::shared_ptr<MolflowSimulationModel> CLILoadFromXML(const std::string &fileName, bool loadSimulationState, std::shared_ptr<MolflowSimulationModel> model,
                           GlobalSimuState *globState, MolflowUserSettings& userSettings, bool noProgress);
    static int initSimModel(std::shared_ptr<MolflowSimulationModel> model);
public:
    static std::string getAutosaveFile();
    static std::shared_ptr<MolflowSimulationModel> initFromFile(SimulationManager *simManager, std::shared_ptr<MolflowSimulationModel> model, GlobalSimuState *globState, MolflowUserSettings& userSettings);
    
    static int initDesLimit(std::shared_ptr<MolflowSimulationModel> model, GlobalSimuState& globState);

    static void initFromArgv(int argc, char **argv, SimulationManager *simManager, std::shared_ptr<MolflowSimulationModel> model);

    static int initTimeLimit(std::shared_ptr<MolflowSimulationModel> model, double time);

    static int
    loadFromGeneration(std::shared_ptr<MolflowSimulationModel> model, GlobalSimuState *globState, double ratio,
                       int step, double angle);
};


#endif //MOLFLOW_PROJ_INITIALIZER_H
