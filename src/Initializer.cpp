//
// Created by Pascal Baehr on 01.08.20.
//

#include "Initializer.h"
#include "CLI11/CLI11.hpp"
#include "IO/LoaderXML.h"
#include "ParameterParser.h"

namespace Settings {
    double nbCPUCores = 0;
    size_t nbThreads = 0;
    uint64_t simDuration = 10;
    uint64_t autoSaveDuration = 600; // default: autosave every 600s=10min
    bool loadAutosave = false;
    std::list<uint64_t> desLimit;
    bool resetOnStart = false;
    std::string inputFile;
    std::string outputFile;
    std::string paramFile;
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

    simManager->nbThreads = Settings::nbThreads;
    simManager->useCPU = true;

    if(simManager->InitSimUnits()) {
        std::cout << "Error: Initialising simulation unit: " << simManager->nbThreads << std::endl;
        return 1;
    }
    std::cout << "Active cores: " << simManager->nbThreads << std::endl;
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
    if(!Settings::paramFile.empty()){
        // Sweep parameters from file
        ParameterParser::Parse(Settings::paramFile);
        ParameterParser::ChangeSimuParams(model->wp);
        ParameterParser::ChangeFacetParams(model->facets);
    }

    simManager->ForwardGlobalCounter(globState);

    return 0;
}

int Initializer::parseCommands(int argc, char** argv) {
    CLI::App app{"Molflow+/Synrad+ Simulation Management"};
    app.formatter(std::make_shared<FlowFormatter>());

    std::vector<double> limits;
    // Define options
    app.add_option("-j,--threads", Settings::nbThreads, "# threads per core");
    app.add_option("-p,--procs", Settings::nbCPUCores, "# CPU cores");
    app.add_option("-t,--time", Settings::simDuration, "Simulation duration in seconds");
    app.add_option("-d,--ndes", limits, "Desorption limit for simulation end");
    app.add_option("-f,--file", Settings::inputFile, "Required input file (XML only)")
            ->required()
            ->check(CLI::ExistingFile);
    app.add_option("-o,--output", Settings::outputFile, "Output file if different from input file");
    app.add_option("-a,--autosaveDuration", Settings::autoSaveDuration, "Seconds for autosave if not zero");
    app.add_flag("--loadAutosave", Settings::loadAutosave, "Whether autosave_ file should be used if exists");
    app.add_option("--setParams", Settings::paramFile, "Parameter file for ad hoc change of the given geometry parameters")
            ->check(CLI::ExistingFile);

    app.add_flag("-r,--reset", Settings::resetOnStart, "Resets simulation status loaded from while");
    app.set_config("--config");
    CLI11_PARSE(app, argc, argv);

    //std::cout<<app.config_to_str(true,true);

    std::cout << "Number used CPU cores: " << Settings::nbCPUCores << std::endl;

    // Save to inputfile
    if(Settings::outputFile.empty())
        Settings::outputFile = Settings::inputFile;

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
    if(loader.LoadGeometry(Settings::inputFile, model)){
        std::cerr << "[Error (LoadGeom)] Please check the input file!" << std::endl;
        model->m.unlock();
        exit(0);
    }

    std::cout << "[LoadGeom] Loaded geometry of " << sizeof(*model) << " / "  << model->size() << " bytes!" << std::endl;

    model->facets = std::move(loader.loadFacets);

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
        std::cout << "[LoadGeom] Initializing geometry!" << std::endl;
        initSimModel(model);
        // temp facets from loader to model 2d (structure, facet)
        //
        std::cout << "[LoadGeom] Resizing state!" << std::endl;
        globState->Resize(*model);

        // 3. init counters with previous results
        if(loadState) {
            std::cout << "[LoadGeom] Initializing previous simulation state!" << std::endl;

            if(Settings::loadAutosave){
                std::string fileName = std::filesystem::path(Settings::inputFile).filename().string();
                std::string autoSavePrefix = "autosave_";
                fileName = autoSavePrefix + fileName;
                if(std::filesystem::exists(fileName)) {
                    std::cout << "Found autosave file! Loading simulation state..." << std::endl;
                    loader.LoadSimulationState(fileName, model, *globState);
                }
            }
            else {
                loader.LoadSimulationState(Settings::inputFile, model, *globState);
            }
        }


    }
    catch (std::exception& e) {
        std::cerr << "[Warning (LoadGeom)] " << e.what() << std::endl;
    }

    // Some postprocessing
    //loader.MoveFacetsToStructures(model);
    std::cout << "[LoadGeom] Forwarding model to simulation units!" << std::endl;

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
        autoSave = std::filesystem::path(Settings::inputFile).filename().string();

        std::string autoSavePrefix = "autosave_";
        if(autoSave.size() > autoSavePrefix.size() && std::search(autoSave.begin(), autoSave.begin()+autoSavePrefix.size(), autoSavePrefix.begin(), autoSavePrefix.end()) == autoSave.begin())
        {
            autoSave = std::filesystem::path(Settings::inputFile).filename().string();
            Settings::inputFile = autoSave.substr(autoSavePrefix.size(), autoSave.size() - autoSavePrefix.size());
            std::cout << "Using autosave file " << autoSave << " for " << Settings::inputFile << '\n';
        }
        else {
            // create autosavefile from copy of original
            std::stringstream autosaveFile;
            autosaveFile << autoSavePrefix<< autoSave;
            autoSave = autosaveFile.str();

            try {
                std::filesystem::copy_file(Settings::inputFile, autoSave,
                                           std::filesystem::copy_options::overwrite_existing);
            } catch (std::filesystem::filesystem_error &e) {
                std::cout << "Could not copy file: " << e.what() << '\n';
            }
        }
    }

    return autoSave;
}

/**
* \brief Serialization function for a binary cereal archive for the worker attributes
* \return output string stream containing the result of the archiving
*/
int Initializer::initSimModel(SimulationModel* model) {


    std::vector<Moment> momentIntervals;
    momentIntervals.reserve(model->tdParams.moments.size());
    for(auto& moment : model->tdParams.moments){
        momentIntervals.emplace_back(std::make_pair(moment.first - (0.5 * moment.second), moment.first + (0.5 * moment.second)));
    }

    model->tdParams.moments = momentIntervals;


    model->structures.resize(model->sh.nbSuper); //Create structures


    //TODO: Globalize Size values
    size_t angleMapTotalSize = 0;
    size_t dirTotalSize = 0;
    size_t profTotalSize = 0;
    size_t textTotalSize = 0;
    size_t histogramTotalSize = 0;

    bool hasVolatile = false;

    auto& loadFacets = model->facets;
    for (size_t facIdx = 0; facIdx < model->sh.nbFacet; facIdx++) {
        SubprocessFacet& sFac = loadFacets[facIdx];

        std::vector<double> textIncVector;
        // Add surface elements area (reciprocal)
        if (sFac.sh.isTextured) {
            textIncVector.resize(sFac.sh.texHeight*sFac.sh.texWidth);

            double rw = sFac.sh.U.Norme() / (double)(sFac.sh.texWidthD);
            double rh = sFac.sh.V.Norme() / (double)(sFac.sh.texHeightD);
            double area = rw * rh;
            area *= (sFac.sh.is2sided) ? 2.0 : 1.0;
            size_t add = 0;
            for (size_t j = 0; j < sFac.sh.texHeight; j++) {
                for (size_t i = 0; i < sFac.sh.texWidth; i++) {
                    if (area > 0.0) {
                        textIncVector[add] = 1.0 / area;
                    }
                    else {
                        textIncVector[add] = 0.0;
                    }
                    add++;
                }
            }
        }
        sFac.textureCellIncrements = textIncVector;

        //Some initialization
        if (!sFac.InitializeOnLoad(facIdx, model->tdParams.moments.size(), histogramTotalSize)) return false;
        // Increase size counters
        //histogramTotalSize += 0;
        angleMapTotalSize += sFac.angleMapSize;
        dirTotalSize += sFac.directionSize* (1 + model->tdParams.moments.size());
        profTotalSize += sFac.profileSize* (1 + model->tdParams.moments.size());
        textTotalSize += sFac.textureSize* (1 + model->tdParams.moments.size());

        hasVolatile |= sFac.sh.isVolatile;

        if ((sFac.sh.superDest || sFac.sh.isVolatile) && ((sFac.sh.superDest - 1) >= model->sh.nbSuper || sFac.sh.superDest < 0)) {
            // Geometry error
            //ClearSimulation();
            //ReleaseDataport(loader);
            std::ostringstream err;
            err << "Invalid structure (wrong link on F#" << facIdx + 1 << ")";
            //SetErrorSub(err.str().c_str());
            std::cerr << err.str() << std::endl;
            return 1;
        }

        if (sFac.sh.superIdx == -1) { //Facet in all structures
            for (auto& s : model->structures) {
                s.facets.push_back(sFac);
            }
        }
        else {
            model->structures[sFac.sh.superIdx].facets.push_back(sFac); //Assign to structure
        }
    }

    return 0;
}