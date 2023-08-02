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

#pragma once

#include <list>
#include <string>
#include <SimulationManager.h>
#include "Simulation/MolflowSimGeom.h"
#include "SettingsIO.h"

class Initializer {
private:
    static SettingsIO::CLIArguments parseArguments(int argc, char** argv);
    static std::shared_ptr<MolflowSimulationModel> CLILoadFromXML(const std::string &fileName, bool loadSimulationState, std::shared_ptr<MolflowSimulationModel> model,
                           GlobalSimuState *globState, MolflowUserSettings& userSettings, SettingsIO::CLIArguments& parsedArgs);
    static int initSimModel(std::shared_ptr<MolflowSimulationModel> model);
public:
    static std::string getAutosaveFile(SettingsIO::CLIArguments& parsedArgs);
    static std::shared_ptr<MolflowSimulationModel> initFromFile(SimulationManager& simManager, std::shared_ptr<MolflowSimulationModel> model, GlobalSimuState *globState, MolflowUserSettings& userSettings, SettingsIO::CLIArguments& parsedArgs);
    
    static int initDesLimit(std::shared_ptr<MolflowSimulationModel> model, GlobalSimuState& globState, SettingsIO::CLIArguments& parsedArgs);

    [[nodiscard]] static SettingsIO::CLIArguments initFromArgv(int argc, char **argv, SimulationManager& simManager, std::shared_ptr<MolflowSimulationModel> model);

    static int initTimeLimit(std::shared_ptr<MolflowSimulationModel> model, SettingsIO::CLIArguments& parsedArgs);

    static int
    loadFromGeneration(std::shared_ptr<MolflowSimulationModel> model, GlobalSimuState *globState, double ratio,
                       int step, double angle);
};
