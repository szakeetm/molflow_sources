//
// Created by Pascal Baehr on 01.08.20.
//

#include "Initializer.h"
#include "CLI11/CLI11.hpp"
#include "LoaderXML.h"

namespace Settings {
    SimulationManager simManager("molflow","MFLW");
    double nbCPUCores = 0;
    uint64_t simDuration = 10;
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

int Initializer::init(int argc, char **argv, SimulationManager *simManager, SimulationModel *model) {
    parseCommands(argc, argv);

    simManager->nbCores = Settings::nbCPUCores;
    simManager->useCPU = true;

    if(simManager->InitSimUnits())
        std::cout << "Error: Initialising subprocesses: " << simManager->simHandles.size() << std::endl;

    std::cout << "Active cores: " << simManager->simHandles.size() << std::endl;
    model->otfParams.nbProcess = simManager->nbCores;

    loadFromXML(simManager, model);





    return 0;
}

int Initializer::parseCommands(int argc, char** argv) {
    CLI::App app{"Molflow+/Synrad+ Simulation Management"};
    app.formatter(std::make_shared<FlowFormatter>());

    // Define options
    app.add_option("-p,--procs", Settings::nbCPUCores, "# CPU cores");
    app.add_option("-t,--time", Settings::simDuration, "Simulation duration in seconds");
    app.add_option("-f,--file", Settings::req_real_file, "Require an existing file")
            ->required()
            ->check(CLI::ExistingFile);
    CLI11_PARSE(app, argc, argv);

    std::cout << "Number used CPU cores: " << Settings::nbCPUCores << std::endl;
    return 0;
}

int Initializer::loadFromXML(SimulationManager *simManager, SimulationModel *model) {

    //1. Load Input File (regular XML)
    FlowIO::LoaderXML loader;
    // Easy if XML is split into multiple parts
    // Geometry
    // Settings
    // Previous results
    if(loader.LoadGeometry(Settings::req_real_file, model)){
        std::cerr << "[Error (LoadGeom)] Please check the input file!" << std::endl;
        exit(0);
    }

    // 2. Create simulation dataports
    try {

        simManager->ResetSimulations();
        //progressDlg->SetMessage("Creating Logger...");
        size_t logDpSize = 0;
        if (model->otfParams.enableLogging) {
            logDpSize = sizeof(size_t) + model->otfParams.logLimit * sizeof(ParticleLoggerItem);
        }
        simManager->ReloadLogBuffer(logDpSize, true);
        //progressDlg->SetMessage("Creating hit buffer...");
        size_t nbMoments = model->tdParams.moments.size();

        // Calc hitsize to init hit buffer
        {
            size_t hitSize = 0;
            hitSize += sizeof(GlobalHitBuffer) + (1 + nbMoments) * model->wp.globalHistogramParams.GetDataSize();
            for (int i = 0; i < model->sh.nbFacet; i++) {
                hitSize += loader.loadFacets[i].GetHitsSize(nbMoments);
            }
            simManager->ReloadHitBuffer(hitSize);
        }
        BYTE* buffer = simManager->GetLockedHitBuffer();

        // 3. init counters with previous results
        loader.LoadSimulationState(Settings::req_real_file, model, buffer);
        simManager->UnlockHitBuffer();

        // temp facets from loader to model 2d (structure, facet)

        std::string loaderString = loader.SerializeForLoader(model).str();

        if (simManager->ShareWithSimUnits((BYTE *) loaderString.c_str(), loaderString.size(), LoadType::LOADGEOM)) {
            std::string errString = "Failed to send params to sub process!\n";
            std::cerr << "[Warning (LoadGeom)] " << errString.c_str() << std::endl;
            exit(0);
        }
    }
    catch (std::exception& e) {
        std::cerr << "[Warning (LoadGeom)] " << e.what() << std::endl;
    }

    // Some postprocessing
    loader.MoveFacetsToStructures(model);
    
    return 0;
}
