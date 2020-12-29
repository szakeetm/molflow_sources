//
// Created by Pascal Baehr on 01.08.20.
//

#include "Initializer.h"
#include "CLI11/CLI11.hpp"
#include "LoaderXML.h"

namespace Settings {
    SimulationManager simManager("molflow","MFLW");
    double nbCPUCores = 0;
    size_t nbThreadsPerCore = 0;
    uint64_t simDuration = 10;
    uint64_t autoSaveDuration = 600; // default: autosave every 600s=10min
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
    if(!Settings::desLimit.empty()) {
        model->otfParams.desorptionLimit = Settings::desLimit.front();
        Settings::desLimit.pop_front();
    }
    else{
        model->otfParams.desorptionLimit = 0;
    }

    loadFromXML(simManager, model, globState, !Settings::resetOnStart);
    std::cout << "Start with "<< globState->globalHits.globalHits.hit.nbMCHit
              << " : " << globState->globalHits.globalHits.hit.nbDesorbed << " // "<< Settings::resetOnStart << std::endl;

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
    app.add_option("-a", Settings::autoSaveDuration, "Seconds for autosave if not zero.");

    app.add_flag("-r,--reset", Settings::resetOnStart, "Resets simulation status loaded from while");
    CLI11_PARSE(app, argc, argv);

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
        if(loadState)
            loader.LoadSimulationState(Settings::req_real_file, model, *globState);
        //simManager->UnlockHitBuffer();



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
