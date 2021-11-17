//
// Created by Pascal Baehr on 30.07.20.
//

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
