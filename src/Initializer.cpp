//
// Created by Pascal Baehr on 01.08.20.
//

#include "Initializer.h"
#include "IO/LoaderXML.h"
#include "ParameterParser.h"

#include <CLI11/CLI11.hpp>
#include <ziplib/ZipArchive.h>
#include <ziplib/ZipFile.h>
#include <File.h>
#include <Helper/StringHelper.h>
#include <Helper/ConsoleLogger.h>
#include <SettingsIO.h>

namespace Settings {
    size_t nbThreads = 0;
    uint64_t simDuration = 10;
    uint64_t autoSaveDuration = 600; // default: autosave every 600s=10min
    bool loadAutosave = false;
    std::list<uint64_t> desLimit;
    bool resetOnStart = false;
    std::string paramFile;
    std::vector<std::string> paramSweep;
}

void initDefaultSettings() {
    Settings::nbThreads = 0;
    Settings::simDuration = 0;
    Settings::autoSaveDuration = 600;
    Settings::loadAutosave = false;
    Settings::desLimit.clear();
    Settings::resetOnStart = false;
    Settings::paramFile.clear();
    Settings::paramSweep.clear();

    SettingsIO::overwrite = false;
    SettingsIO::workFile.clear();
    SettingsIO::inputFile.clear();
    SettingsIO::outputFile.clear();
    SettingsIO::workPath.clear();
    SettingsIO::inputPath.clear();
    SettingsIO::outputPath.clear();
}

class FlowFormatter : public CLI::Formatter {
public:
    std::string make_usage(const CLI::App *app, std::string name) const override {
        return "Usage: ./"
               + std::filesystem::path(name).filename().string()
               + " [options]";
    }
};

int Initializer::parseCommands(int argc, char **argv) {
    CLI::App app{"Molflow+/Synrad+ Simulation Management"};
    app.formatter(std::make_shared<FlowFormatter>());

    // Local variables for parsing and immediate processing
    bool verbose = false;
    std::vector<double> limits;

    // Define options
    app.add_option("-j,--threads", Settings::nbThreads, "# Threads to be deployed");
    app.add_option("-t,--time", Settings::simDuration, "Simulation duration in seconds");
    app.add_option("-d,--ndes", limits, "Desorption limit for simulation end");
    app.add_option("-f,--file", SettingsIO::inputFile, "Required input file (XML only)")
            ->required()
            ->check(CLI::ExistingFile);
    CLI::Option *optOfile = app.add_option("-o,--output", SettingsIO::outputFile,
                                           R"(Output file name (e.g. 'outfile.xml', defaults to 'out_{inputFileName}')");
    CLI::Option *optOpath = app.add_option("--outputPath", SettingsIO::outputPath,
                                           "Output path, defaults to \'Results_{date}\'");
    app.add_option("-a,--autosaveDuration", Settings::autoSaveDuration, "Seconds for autoSave if not zero");
    app.add_option("--setParamsByFile", Settings::paramFile,
                   "Parameter file for ad hoc change of the given geometry parameters")
            ->check(CLI::ExistingFile);
    app.add_option("--setParams", Settings::paramSweep,
                   "Direct parameter input for ad hoc change of the given geometry parameters");
    app.add_option("--verbosity", Settings::verbosity, "Restrict console output to different levels");

    app.add_flag("--loadAutosave", Settings::loadAutosave, "Whether autoSave_ file should be used if exists");
    app.add_flag("-r,--reset", Settings::resetOnStart, "Resets simulation status loaded from file");
    app.add_flag("--verbose", verbose, "Verbose console output (all levels)");
    CLI::Option *optOverwrite = app.add_flag("--overwrite", SettingsIO::overwrite,
                                             "Overwrite input file with new results")->excludes(optOfile, optOpath);
    optOfile->excludes(optOverwrite);
    optOpath->excludes(optOverwrite);
    app.set_config("--config");

    CLI11_PARSE(app, argc, argv);

    if (verbose)
        Settings::verbosity = 42;

    //std::cout<<app.config_to_str(true,true);
    for (auto &lim : limits)
        Settings::desLimit.emplace_back(lim);

    if (Settings::simDuration == 0 && Settings::desLimit.empty()) {
        Log::console_error("No end criterion has been set!\n");
        Log::console_error(" Either use: -t or -d\n");
        return 1;
    }

    return 0;
}

int Initializer::initFromArgv(int argc, char **argv, SimulationManager *simManager,
                              std::shared_ptr<SimulationModel> model) {
    Log::console_header(1, "Commence: Initialising!\n");

#if defined(WIN32) || defined(__APPLE__)
    setlocale(LC_ALL, "C");
#else
    std::setlocale(LC_ALL, "C");
#endif

    initDefaultSettings();

    int err = 0;
    if ((err = parseCommands(argc, argv))) {
        return err;
    }

    simManager->nbThreads = Settings::nbThreads;
    simManager->useCPU = true;

    if (simManager->InitSimUnits()) {
        Log::console_error("Error: Initialising simulation units: %zu\n", simManager->nbThreads);
        return 1;
    }

    model->otfParams.nbProcess = simManager->nbThreads;
    model->otfParams.timeLimit = (double) Settings::simDuration;
    //model->otfParams.desorptionLimit = Settings::desLimit.front();
    Log::console_msg_master(4, "Active cores: %zu\n", simManager->nbThreads);
    Log::console_msg_master(4, "Running simulation for: %zu sec\n", Settings::simDuration);

    return 0;
}

int Initializer::initFromFile(SimulationManager *simManager, std::shared_ptr<SimulationModel> model,
                              GlobalSimuState *globState) {
    if (SettingsIO::prepareIO()) {
        Log::console_error("Error preparing I/O folders\n");
        return 1;
    }

    if (std::filesystem::path(SettingsIO::workFile).extension() == ".xml")
        loadFromXML(SettingsIO::workFile, !Settings::resetOnStart, model, globState);
    else {
        Log::console_error("Invalid file extension for input file detected: %s\n",
                           std::filesystem::path(SettingsIO::workFile).extension().c_str());
        return 1;
    }
    if (!Settings::paramFile.empty() || !Settings::paramSweep.empty()) {
        // 1. Load selection groups in case we need them for parsing
        std::vector<SelectionGroup> selGroups = FlowIO::LoaderXML::LoadSelections(SettingsIO::workFile);
        // 2. Sweep parameters from file
        if (!Settings::paramFile.empty())
            ParameterParser::ParseFile(Settings::paramFile, selGroups);
        if (!Settings::paramSweep.empty())
            ParameterParser::ParseInput(Settings::paramSweep, selGroups);
        ParameterParser::ChangeSimuParams(model->wp);
        ParameterParser::ChangeFacetParams(model->facets);
    }

    // Set desorption limit if used
    if (initDesLimit(model, *globState)) {
        return 1;
    }

    Log::console_msg_master(2, "Forwarding model to simulation units!\n");
    try {
        simManager->InitSimulation(model, globState);
    }
    catch (std::exception &ex) {
        Log::console_error("Failed initialising simulation units:\n%s\n", ex.what());
        return 1;
    }
    Log::console_footer(1, "Finalize: Initialising!\n");

    return 0;
}

int Initializer::loadFromXML(const std::string &fileName, bool loadState, std::shared_ptr<SimulationModel> model,
                             GlobalSimuState *globState) {

    Log::console_header(1, "[ ] Loading geometry from file %s\n", fileName.c_str());

    //1. Load Input File (regular XML)
    FlowIO::LoaderXML loader;
    // Easy if XML is split into multiple parts
    // Geometry
    // Settings
    // Previous results
    model->m.lock();
    double progress = 0.0;
    if (loader.LoadGeometry(fileName, model, &progress)) {
        Log::console_error("Please check the input file!\n");
            model->m.unlock();
        return 1;
    }

    // InsertParamsBeforeCatalog
    std::vector<Parameter> paramCatalog;
    Parameter::LoadParameterCatalog(paramCatalog);
    model->tdParams.parameters.insert(model->tdParams.parameters.end(), paramCatalog.begin(), paramCatalog.end());

    //TODO: Load parameters from catalog explicitly?
    // For GUI
    // work->InsertParametersBeforeCatalog(loadedParams);
    // Load viewsettings for each facet

    Log::console_msg_master(3, " Loaded geometry of %zu bytes!\n", model->size());

    //InitializeGeometry();
    model->InitialiseFacets();

    Log::console_msg_master(3, " Initializing geometry!\n");
    initSimModel(model);
    model->PrepareToRun();

    // 2. Create simulation dataports
    try {
        //progressDlg->SetMessage("Creating Logger...");
        /*size_t logDpSize = 0;
        if (model->otfParams.enableLogging) {
            logDpSize = sizeof(size_t) + model->otfParams.logLimit * sizeof(ParticleLoggerItem);
        }
        simManager->ReloadLogBuffer(logDpSize, true);*/

        Log::console_msg_master(3, " Resizing state!\n");
        globState->Resize(*model);

        // 3. init counters with previous results
        if (loadState) {
            Log::console_msg_master(3, " Initializing previous simulation state!\n");

            if (Settings::loadAutosave) {
                std::string fileName = std::filesystem::path(SettingsIO::workFile).filename().string();
                std::string autoSavePrefix = "autosave_";
                fileName = autoSavePrefix + fileName;
                if (std::filesystem::exists(fileName)) {
                    Log::console_msg_master(2, " Found autosave file! Loading simulation state...\n");
                    FlowIO::LoaderXML::LoadSimulationState(fileName, model, globState, nullptr);
                }
            } else {
                FlowIO::LoaderXML::LoadSimulationState(SettingsIO::workFile, model, globState, nullptr);
            }
        }
    }
    catch (std::exception &e) {
        Log::console_error("[Warning] %s\n", e.what());
    }

    model->m.unlock();
    Log::console_footer(1, "[x] Loaded geometry\n");

    return 0;
}

int Initializer::initDesLimit(std::shared_ptr<SimulationModel> model, GlobalSimuState &globState) {
    model->otfParams.desorptionLimit = 0;

    // Skip desorptions if limit was already reached
    if (!Settings::desLimit.empty()) {
        size_t oldDesNb = globState.globalHits.globalHits.nbDesorbed;
        size_t listSize = Settings::desLimit.size();
        for (size_t l = 0; l < listSize; ++l) {
            model->otfParams.desorptionLimit = Settings::desLimit.front();
            Settings::desLimit.pop_front();

            if (oldDesNb > model->otfParams.desorptionLimit) {
                Log::console_msg_master(1, "Skipping desorption limit: %zu\n", model->otfParams.desorptionLimit);
            } else {
                Log::console_msg_master(1, "Starting with desorption limit: %zu from %zu\n",
                                        model->otfParams.desorptionLimit, oldDesNb);
                return 0;
            }
        }
        if (Settings::desLimit.empty()) {
            Log::console_msg_master(1,
                                    "All given desorption limits have been reached. Consider resetting the simulation results from the input file (--reset): Starting desorption %zu\n",
                                    oldDesNb);
            return 1;
        }
    }

    return 0;
}

// TODO: Combine with loadXML function
std::string Initializer::getAutosaveFile() {
    // Create copy of input file for autosave
    std::string autoSave;
    if (Settings::autoSaveDuration > 0) {
        autoSave = std::filesystem::path(SettingsIO::workFile).filename().string();

        std::string autoSavePrefix = "autosave_";
        // Check if autosave_ is part of the input filename, if yes, generate a new input file without the prefix
        if (autoSave.size() > autoSavePrefix.size() &&
            std::search(autoSave.begin(), autoSave.begin() + autoSavePrefix.size(), autoSavePrefix.begin(),
                        autoSavePrefix.end()) == autoSave.begin()) {
            // TODO: Revisit wether input/output is acceptable here
            autoSave = std::filesystem::path(SettingsIO::workFile).filename().string();
            SettingsIO::inputFile = autoSave.substr(autoSavePrefix.size(), autoSave.size() - autoSavePrefix.size());
            Log::console_msg_master(2, "Using autosave file %s for %s\n", autoSave.c_str(),
                                    SettingsIO::inputFile.c_str());
        } else {
            // create autosavefile from copy of original
            autoSave = std::filesystem::path(SettingsIO::outputPath).append(autoSavePrefix).concat(autoSave).string();
            try {
                std::filesystem::copy_file(SettingsIO::workFile, autoSave,
                                           std::filesystem::copy_options::overwrite_existing);
            } catch (std::filesystem::filesystem_error &e) {
                Log::console_error("Could not copy file: %s\n", e.what());
            }
        }
    }

    return autoSave;
}

/**
* \brief Prepares data structures for use in simulation
* \return error code: 0=no error, 1=error
*/
int Initializer::initSimModel(std::shared_ptr<SimulationModel> model) {

    std::vector<Moment> momentIntervals;
    momentIntervals.reserve(model->tdParams.moments.size());
    for (auto &moment : model->tdParams.moments) {
        momentIntervals.emplace_back(
                std::make_pair(moment.first - (0.5 * moment.second), moment.first + (0.5 * moment.second)));
    }

    model->tdParams.moments = momentIntervals;


    model->structures.resize(model->sh.nbSuper); //Create structures

    bool hasVolatile = false;

    for (size_t facIdx = 0; facIdx < model->sh.nbFacet; facIdx++) {
        auto sFac = model->facets[facIdx];

        std::vector<double> textIncVector;
        // Add surface elements area (reciprocal)
        if (sFac->sh.isTextured) {
            textIncVector.resize(sFac->sh.texHeight * sFac->sh.texWidth);

            double rw = sFac->sh.U.Norme() / (double) (sFac->sh.texWidth_precise);
            double rh = sFac->sh.V.Norme() / (double) (sFac->sh.texHeight_precise);
            double area = rw * rh;
            area *= (sFac->sh.is2sided) ? 2.0 : 1.0;
            size_t add = 0;
            for (size_t j = 0; j < sFac->sh.texHeight; j++) {
                for (size_t i = 0; i < sFac->sh.texWidth; i++) {
                    if (area > 0.0) {
                        textIncVector[add] = 1.0 / area;
                    } else {
                        textIncVector[add] = 0.0;
                    }
                    add++;
                }
            }
        }
        sFac->textureCellIncrements = textIncVector;

        //Some initialization
        if (!sFac->InitializeOnLoad(facIdx, model->tdParams.moments.size())) return false;

        hasVolatile |= sFac->sh.isVolatile;

        if ((sFac->sh.superDest || sFac->sh.isVolatile) &&
            ((sFac->sh.superDest - 1) >= model->sh.nbSuper || sFac->sh.superDest < 0)) {
            // Geometry error
            //ClearSimulation();
            //ReleaseDataport(loader);
            //SetErrorSub(err.str().c_str());
            Log::console_error("Invalid structure (wrong link on F#%d)\n", facIdx + 1);

            return 1;
        }
    }

    return 0;
}