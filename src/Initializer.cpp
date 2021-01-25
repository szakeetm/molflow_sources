//
// Created by Pascal Baehr on 01.08.20.
//

#include "Initializer.h"
#include "CLI11/CLI11.hpp"
#include "IO/LoaderXML.h"

namespace Settings {
    double nbCPUCores = 0;
    size_t nbThreadsPerCore = 0;
    uint64_t simDuration = 10;
    uint64_t autoSaveDuration = 600; // default: autosave every 600s=10min
    bool loadAutosave = false;
    std::list<uint64_t> desLimit;
    bool resetOnStart = false;
    std::string req_real_file;
}

class FlowFormatter : public CLI::Formatter {
public:
    std::string make_usage(const CLI::App *app, std::string name) const override {
        return "Usage: ./"
               +std::filesystem::path(name).filename().string()
               +" [options]";
    }
};

int Initializer::init(int argc, char **argv, SimulationManager *simManager, SimulationModel *model,
                      GlobalSimuState *globState) {

#if defined(WIN32) || defined(__APPLE__)
    setlocale(LC_ALL, "C");
#else
    std::setlocale(LC_ALL, "C");
#endif
    parseCommands(argc, argv);

    simManager->nbThreads = Settings::nbThreadsPerCore;
    simManager->nbCores = Settings::nbCPUCores;
    simManager->useCPU = true;

    if(simManager->InitSimUnits()) {
        std::cout << "Error: Initialising subprocesses: " << simManager->nbCores << std::endl;
        return 1;
    }
    std::cout << "Active cores: " << simManager->nbCores << std::endl;
    model->otfParams.nbProcess = simManager->nbThreads;
    //model->otfParams.desorptionLimit = Settings::desLimit.front();

    // Set desorption limit if used
    if(initDesLimit(*model,*globState)) {
        exit(0);
    }
    else{
        model->otfParams.desorptionLimit = 0;
    }

    loadFromXML(simManager, model, globState, !Settings::resetOnStart);

    simManager->ForwardGlobalCounter(globState);

    return 0;
}

int Initializer::parseCommands(int argc, char** argv) {
    CLI::App app{"Molflow+/Synrad+ Simulation Management"};
    app.formatter(std::make_shared<FlowFormatter>());

    std::vector<double> limits;
    // Define options
    app.add_option("-j,--threads", Settings::nbThreadsPerCore, "# threads per core");
    app.add_option("-p,--procs", Settings::nbCPUCores, "# CPU cores");
    app.add_option("-t,--time", Settings::simDuration, "Simulation duration in seconds");
    app.add_option("-d,--ndes", limits, "Desorption limit for simulation end");
    app.add_option("-f,--file", Settings::req_real_file, "Require an existing file")
            ->required()
            ->check(CLI::ExistingFile);
    app.add_option("-a,--autosaveDuration", Settings::autoSaveDuration, "Seconds for autosave if not zero");
    app.add_flag("--loadAutosave", Settings::loadAutosave, "Whether autosave_ file should be used if exists");

    app.add_flag("-r,--reset", Settings::resetOnStart, "Resets simulation status loaded from while");
    app.set_config("--config");
    CLI11_PARSE(app, argc, argv);

    //std::cout<<app.config_to_str(true,true);

    std::cout << "Number used CPU cores: " << Settings::nbCPUCores << std::endl;
    for(auto& lim : limits)
        Settings::desLimit.emplace_back(lim);
    return 0;
}

int Initializer::loadFromXML(SimulationManager *simManager, SimulationModel *model, GlobalSimuState *globState,
                             bool loadState) {

    //1. Load Input File (regular XML)
    FlowIO::LoaderXML loader;
    // Easy if XML is split into multiple parts
    // Geometry
    // Settings
    // Previous results
    model->m.lock();
    if(loader.LoadGeometry(Settings::req_real_file, model)){
        std::cerr << "[Error (LoadGeom)] Please check the input file!" << std::endl;
        model->m.unlock();
        exit(0);
    }

    // 2. Create simulation dataports


    try {

        simManager->ResetSimulations();
        //progressDlg->SetMessage("Creating Logger...");
        /*size_t logDpSize = 0;
        if (model->otfParams.enableLogging) {
            logDpSize = sizeof(size_t) + model->otfParams.logLimit * sizeof(ParticleLoggerItem);
        }
        simManager->ReloadLogBuffer(logDpSize, true);*/
        //progressDlg->SetMessage("Creating hit buffer...");
        size_t nbMoments = model->tdParams.moments.size();

        /*// Calc hitsize to init hit buffer
        {
            size_t hitSize = 0;
            hitSize += sizeof(GlobalHitBuffer) + (1 + nbMoments) * model->wp.globalHistogramParams.GetDataSize();
            for (int i = 0; i < model->sh.nbFacet; i++) {
                hitSize += loader.loadFacets[i].GetHitsSize(nbMoments);
            }
            hitSize =
            simManager->ReloadHitBuffer(hitSize);
        }
        BYTE* buffer = simManager->GetLockedHitBuffer();
*/
        loader.InitSimModel(model);
        // temp facets from loader to model 2d (structure, facet)
        //
        globState->Resize(*model);

        // 3. init counters with previous results
        if(loadState) {

            if(Settings::loadAutosave){
                std::string fileName = std::filesystem::path(Settings::req_real_file).filename().string();
                std::string autoSavePrefix = "autosave_";
                fileName = autoSavePrefix + fileName;
                if(std::filesystem::exists(fileName)) {
                    std::cout << "Found autosave file! Loading simulation state..." << std::endl;
                    loader.LoadSimulationState(fileName, model, *globState);
                }
            }
            else {
                loader.LoadSimulationState(Settings::req_real_file, model, *globState);
            }
        }


    }
    catch (std::exception& e) {
        std::cerr << "[Warning (LoadGeom)] " << e.what() << std::endl;
    }

    // Some postprocessing
    //loader.MoveFacetsToStructures(model);
    simManager->ForwardSimModel(model);

    if (simManager->ExecuteAndWait(COMMAND_LOAD, PROCESS_READY, 0, 0)) {
        //CloseLoaderDP();
        model->m.unlock();
        std::string errString = "Failed to send geometry to sub process:\n";
        errString.append(simManager->GetErrorDetails());
        throw std::runtime_error(errString);
    }
    model->m.unlock();
    return 0;
}

int Initializer::initDesLimit(SimulationModel& model, GlobalSimuState& globState){
    model.otfParams.desorptionLimit = 0;

    // Skip desorptions if limit was already reached
    if(!Settings::desLimit.empty())
    {
        size_t oldDesNb = globState.globalHits.globalHits.hit.nbDesorbed;
        size_t listSize = Settings::desLimit.size();
        for(size_t l = 0; l < listSize; ++l) {
            model.otfParams.desorptionLimit = Settings::desLimit.front();
            Settings::desLimit.pop_front();

            if (oldDesNb > model.otfParams.desorptionLimit){
                printf("Skipping desorption limit: %zu\n", model.otfParams.desorptionLimit);
            }
            else{
                printf("Starting with desorption limit: %zu from %zu\n", model.otfParams.desorptionLimit , oldDesNb);
                return 0;
            }
        }
        if(Settings::desLimit.empty()){
            printf("All given desorption limits have been reached. Consider resetting the simulation results from the input file (--reset): Starting desorption %zu\n", oldDesNb);
            return 1;
        }
    }

    return 0;
}

// TODO: Combine with loadXML function
std::string Initializer::getAutosaveFile(){
    // Create copy of input file for autosave
    std::string autoSave;
    if(Settings::autoSaveDuration > 0)
    {
        autoSave = std::filesystem::path(Settings::req_real_file).filename().string();

        std::string autoSavePrefix = "autosave_";
        if(autoSave.size() > autoSavePrefix.size() && std::search(autoSave.begin(), autoSave.begin()+autoSavePrefix.size(), autoSavePrefix.begin(), autoSavePrefix.end()) == autoSave.begin())
        {
            autoSave = std::filesystem::path(Settings::req_real_file).filename().string();
            Settings::req_real_file = autoSave.substr( autoSavePrefix.size(), autoSave.size() - autoSavePrefix.size());
            std::cout << "Using autosave file " << autoSave << " for "<<Settings::req_real_file<<'\n';
        }
        else {
            // create autosavefile from copy of original
            std::stringstream autosaveFile;
            autosaveFile << autoSavePrefix<< autoSave;
            autoSave = autosaveFile.str();

            try {
                std::filesystem::copy_file(Settings::req_real_file, autoSave,
                                           std::filesystem::copy_options::overwrite_existing);
            } catch (std::filesystem::filesystem_error &e) {
                std::cout << "Could not copy file: " << e.what() << '\n';
            }
        }
    }

    return autoSave;
}
