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

#ifndef MOLFLOW_PROJ_CONFIG_H
#define MOLFLOW_PROJ_CONFIG_H

#include <string>
#include <list>

namespace FlowIO {
    class Config {
        public:
            Config() {
                nbCPUCores = 0;
                nbThreadsPerCore = 0;
                simDuration = 10;
                autoSaveDuration = 600; // default: autosave every 600s=10min
                loadAutosave = false;
                resetOnStart = false;
            }
            ~Config() = default;

            int ReadFromFile(const std::string& configFileName);
            double nbCPUCores;
            size_t nbThreadsPerCore;
            uint64_t simDuration;
            uint64_t autoSaveDuration; // default: autosave every 600s=10min
            bool loadAutosave;
            std::list<uint64_t> desLimit;
            bool resetOnStart;
            std::string req_real_file;
    };
}

#endif //MOLFLOW_PROJ_CONFIG_H
