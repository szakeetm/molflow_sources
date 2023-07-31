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

#include "Initializer.h"
#include "IO/LoaderXML.h"
#include "ParameterParser.h"
#include <Helper/GLProgress_CLI.hpp>

#include <CLI11/CLI11.hpp>
#include <Helper/StringHelper.h>
#include <Helper/ConsoleLogger.h>
#include <SettingsIO.h>

namespace CLIArguments {
    size_t nbThreads = 0;
    uint64_t simDuration = 10;
    uint64_t outputDuration = 60;
    uint64_t autoSaveDuration = 600; // default: autosave every 600s=10min
    bool loadAutosave = false;
    std::list<uint64_t> desLimit;
    bool resetOnStart = false;
    std::string paramFile;
    std::vector<std::string> paramChanges;
    bool noProgress = false;  //If true, doesn't print percentage updates for CLI progressbars
}

void initDefaultSettings() {
    CLIArguments::nbThreads = 0;
    CLIArguments::simDuration = 0;
    CLIArguments::outputDuration = 60;
    CLIArguments::autoSaveDuration = 600;
    CLIArguments::loadAutosave = false;
    CLIArguments::desLimit.clear();
    CLIArguments::resetOnStart = false;
    CLIArguments::paramFile.clear();
    CLIArguments::paramChanges.clear();
    CLIArguments::noProgress = false;

    SettingsIO::outputFacetDetails = false;
    SettingsIO::outputFacetQuantities = false;
    SettingsIO::overwrite = false;
    SettingsIO::autogenerateTest = false;

    SettingsIO::workFile.clear();
    SettingsIO::inputFile.clear();
    SettingsIO::outputFile.clear();
    SettingsIO::workPath.clear();
    //SettingsIO::inputPath.clear();
    SettingsIO::outputPath.clear();
}

/**
* \brief Modifies the default --help command output
 */
class FlowFormatter : public CLI::Formatter {
public:
    std::string make_usage(const CLI::App *app, std::string name) const override {
        return "Usage: ./"
               + std::filesystem::path(name).filename().string()
               + " [options]";
    }
};

/**
* \brief Parse all CLI commands to use with the corresponding variables (@Settings / @SettingsIO)
 * \return 0> error code, 0 when ok
 */
int Initializer::parseCommands(int argc, char **argv) {
    CLI::App app{"Molflow+/Synrad+ Simulation Management"};
    app.formatter(std::make_shared<FlowFormatter>());

    // Local variables for parsing and immediate processing
    bool verbose = false;
    std::vector<double> limits;

    // Define options
    app.add_option("-j,--threads", CLIArguments::nbThreads, "# Threads to be deployed");
    app.add_option("-t,--time", CLIArguments::simDuration, "Simulation duration in seconds");
    app.add_option("-d,--ndes", limits, "Desorption limit for simulation end");

    auto group = app.add_option_group("subgroup");
    group->add_option("-f,--file", SettingsIO::inputFile, "Required input file (XML/ZIP only)")
            ->check(CLI::ExistingFile);
    group->add_flag("--auto", SettingsIO::autogenerateTest, "Use auto generated test case");
    group->require_option(1);

    CLI::Option *optOfile = app.add_option("-o,--output", SettingsIO::outputFile,
                                           R"(Output file name (e.g. 'outfile.xml', defaults to 'out_{inputFileName}')");
    CLI::Option *optOpath = app.add_option("--outputPath", SettingsIO::outputPath,
                                           "Output path, defaults to \'Results_{date}\'");
    app.add_option("-s,--outputDuration", CLIArguments::outputDuration, "Seconds between each stat output if not zero");
    app.add_option("-a,--autosaveDuration", CLIArguments::autoSaveDuration, "Seconds for autoSave if not zero");
    app.add_flag("--writeFacetDetails", SettingsIO::outputFacetDetails,
                   "Will write a CSV file containing all facet details including physical quantities");
    app.add_flag("--writeFacetQuantities", SettingsIO::outputFacetQuantities,
                   "Will write a CSV file containing all physical quantities for each facet");

    app.add_option("--setParamsByFile", CLIArguments::paramFile,
                   "Parameter file for ad hoc change of the given geometry parameters")
            ->check(CLI::ExistingFile);
    app.add_option("--setParams", CLIArguments::paramChanges,
                   "Direct parameter input for ad hoc change of the given geometry parameters");
    app.add_option("--verbosity", AppSettings::verbosity, "Restrict console output to different levels");
    app.add_flag("--noProgress", CLIArguments::noProgress, "Log file mode: No percentage updates printed of progress");
    app.add_flag("--loadAutosave", CLIArguments::loadAutosave, "Whether autosave_ file should be used if exists");
    app.add_flag("-r,--reset", CLIArguments::resetOnStart, "Resets simulation status loaded from file");
    app.add_flag("--verbose", verbose, "Verbose console output (all levels)");
    CLI::Option *optOverwrite = app.add_flag("--overwrite", SettingsIO::overwrite,
                                             "Overwrite input file with new results")->excludes(optOfile, optOpath);
    optOfile->excludes(optOverwrite);
    optOpath->excludes(optOverwrite);
    app.set_config("--config");

    CLI11_PARSE(app, argc, argv);

    if (verbose)
        AppSettings::verbosity = 42;

    //std::cout<<app.config_to_str(true,true);
    for (auto& lim : limits)
        CLIArguments::desLimit.push_back(static_cast<size_t>(lim));

    if (CLIArguments::simDuration == 0 && CLIArguments::desLimit.empty()) {
        Log::console_error("No end criterion has been set. Use either -t or -d\n");
        return 0;
    }

    return -1;
}

/**
* \brief Comfort function wrapping around default initialization from command line arguments, taking care of default initializations
 * \return 0> error code, -1 when ok
 */
int Initializer::initFromArgv(int argc, char **argv, SimulationManager *simManager,
                              std::shared_ptr<MolflowSimulationModel> model) {

// Set local to parse input files the same on all systems
//duplicate, in case we called this function from the test suite and not from main()
#if defined(__APPLE__)
    setlocale(LC_ALL, "en_US.UTF-8");
#else
    std::setlocale(LC_ALL, "en_US.UTF-8");
#endif

    initDefaultSettings();

    int err = 1;
    if (-1 < (err = parseCommands(argc, argv))) {
        Log::console_error("Error: Initializing parsing arguments\n");
        return err;
    }

    Log::console_header(1, "Initializing simulation...\n");

    simManager->nbThreads = CLIArguments::nbThreads;
    simManager->noProgress = CLIArguments::noProgress;

    if (simManager->SetUpSimulation()) { //currently only calls CreateCPUHandle()
        Log::console_error("Error: Setting up simulation units [{} threads]...\n", simManager->nbThreads);
        return 1;
    }

    model->otfParams.nbProcess = simManager->nbThreads;
    model->otfParams.timeLimit = (double) CLIArguments::simDuration;
    //model->otfParams.desorptionLimit = CLIArguments::desLimit.front();
    Log::console_msg_master(1, "Active cores: {}\n", simManager->nbThreads);
    if (CLIArguments::simDuration != 0) {
        Log::console_msg_master(1, "Running simulation for: {} sec\n", CLIArguments::simDuration);
    }

    return -1;
}

/**
* \brief Initializes the simulation model from a valid input file and handles parameter changes
 */
void Initializer::initFromFile(SimulationManager *simManager, std::shared_ptr<MolflowSimulationModel> model,
                              GlobalSimuState *globStatePtr, UserSettings& userSettings) {
    if (SettingsIO::prepareIO()) {
        throw Error("Error preparing I/O folders");
    }

    if (std::filesystem::path(SettingsIO::workFile).extension() == ".xml") {
        CLILoadFromXML(SettingsIO::workFile, !CLIArguments::resetOnStart, model, globStatePtr, userSettings, simManager->noProgress);
    }
    else {
        throw Error(fmt::format("Invalid file extension for input file detected: {}\n",
                           std::filesystem::path(SettingsIO::workFile).extension().string()));
    }
    if (!CLIArguments::paramFile.empty() || !CLIArguments::paramChanges.empty()) {
        // Apply parameter changes from file
        if (!CLIArguments::paramFile.empty())
            ParameterParser::ParseFile(CLIArguments::paramFile, userSettings.selections);
        if (!CLIArguments::paramChanges.empty())
            ParameterParser::ParseInput(CLIArguments::paramChanges, userSettings.selections);
        ParameterParser::ChangeSimuParams(model->wp);
        if(ParameterParser::ChangeFacetParams(model->facets)){
            throw Error("Error in ParameterParser::ChangeFacetParams()");
        }
    }

    // Set desorption limit if used
    if (initDesLimit(model, *globStatePtr)) {
        throw Error("Error in Initializer::initDesLimit");
    }

    simManager->simulationChanged = true;
    Log::console_msg_master(2, "Forwarding model to simulation units...\n");
    try {
        simManager->InitSimulation(model, globStatePtr);
    }
    catch (std::exception &ex) {
        Log::console_error("Failed Initializing simulation units:\n{}\n", ex.what());
        throw ex;
    }
    Log::console_footer(1, "Forwarded successfully.\n");
}

/**
* \brief Initialize simulation from automatically generated test case (prism)
 * \return 0> error code, 0 when ok
 */
void Initializer::initAutoGenerated(SimulationManager *simManager, std::shared_ptr<MolflowSimulationModel> model,
                                   GlobalSimuState *globStatePtr, double ratio, int steps, double angle) {
    if (SettingsIO::prepareIO()) {
        throw Error("Error preparing I/O folders\n");
    }

    loadFromGeneration(model, globStatePtr, ratio, steps, angle);

    // Set desorption limit if used
    if (initDesLimit(model, *globStatePtr)) {
        throw Error("Failed to set desorption limit(s).");
    }

    simManager->simulationChanged = true;
    Log::console_msg_master(2, "Forwarding model to simulation units...\n");
    simManager->InitSimulation(model, globStatePtr);
    Log::console_footer(1, "Forwarded successfully.\n");
}

/**
* \brief Generate a test case with an oblique prism given an angle
 * \return 0> error code, 0 when ok
 */
int
Initializer::loadFromGeneration(std::shared_ptr<MolflowSimulationModel> model, GlobalSimuState *globStatePtr, double ratio,
                                int step, double angle) {

    Log::console_header(1, "[ ] Loading geometry : PRISM\n");

    //double ratio = 10.0;
    double R = 1.0;
    double L = ratio * R;
    //int    step = 10;
    //1. Load Input File (regular XML)
    // Geometry
    model->BuildPrisma(L, R, angle, 0.0, step);
    // Settings
    // Previous results

    TimeDependentParameters::LoadParameterCatalog(model->tdParams.parameters);

    Log::console_msg_master(3, "Loaded geometry ({} bytes)\n", model->size());

    //InitializeGeometry();
    model->InitializeFacets();

    initSimModel(model);
    try {
        model->PrepareToRun();
    }
    catch (std::exception& err) {
        Log::console_error("Error in model->PrepareToRun()\n{}",err.what());
        return 1; //to convert to throwing error
    }
    globStatePtr->Resize(model);
    return 0;
}

/**
* \brief Wrapper for CLI's initFromFile for XML loading with LoaderXML
 */
void Initializer::CLILoadFromXML(const std::string &fileName, bool loadSimulationState, std::shared_ptr<MolflowSimulationModel> model,
                             GlobalSimuState *globStatePtr, UserSettings& userSettings, bool noProgress) {

    //Log::console_header(1, "Loading geometry from file {}\n", fileName);
    GLProgress_CLI prg(fmt::format("Loading geometry from file {}", fileName));
    prg.noProgress = noProgress; //Suppress percentages if ran in script

    //1. Load Input File (regular XML)
    {
        FlowIO::LoaderXML loader;
        // Easy if XML is split into multiple parts
        // Geometry
        // Settings
        // Previous results
        loader.LoadGeometry(fileName, model, prg);
        userSettings = loader.userSettings; //persistent for saving
    }

    prg.SetMessage(fmt::format("Loaded geometry ({} bytes)", model->size()));

    //InitializeGeometry();
    model->InitializeFacets();

    prg.SetMessage("Initializing geometry...",true);
    initSimModel(model);
    model->PrepareToRun();

    // 2. Create simulation dataports
    try {

        prg.SetMessage("Resizing global state counters...");
        globStatePtr->Resize(model);

        // 3. init counters with previous results
        if (loadSimulationState) {
            prg.SetMessage("Loading simulation state...");
            if (CLIArguments::loadAutosave) {
                std::string autosaveFileName = std::filesystem::path(SettingsIO::workFile).filename().string();
                std::string autoSavePrefix = "autosave_";
                autosaveFileName = autoSavePrefix + autosaveFileName;
                
                if (std::filesystem::exists(autosaveFileName)) {
                    prg.SetMessage(fmt::format("Found autosave file {}, loading simulation state...\n",autosaveFileName),true);
                    FlowIO::LoaderXML::LoadSimulationState(autosaveFileName, model, globStatePtr, prg);
                }
            } else {
                FlowIO::LoaderXML::LoadSimulationState(SettingsIO::workFile, model, globStatePtr, prg);
            }

            // Update Angle map status
            for(int i = 0; i < model->facets.size(); i++ ) {
#if defined(MOLFLOW)
                auto f = std::dynamic_pointer_cast<MolflowSimFacet>(model->facets[i]);
                if (f->sh.anglemapParams.record) { //Recording, so needs to be updated
                    //Retrieve angle map from hits dp
                    globStatePtr->facetStates[i].recordedAngleMapPdf = f->angleMap.pdf;
                }
#endif
            }
        }
    }
    catch (const std::exception &e) {
        Log::console_error("[Warning] {}\n", e.what());
    }
}

/**
* \brief Initialize a desorption limit by checking against the given input (vector)
 * \return 0> error code, 0 when ok
 */
int Initializer::initDesLimit(std::shared_ptr<MolflowSimulationModel> model, GlobalSimuState &globState) {
    try { //unti initDesLimit will throw error
        std::lock_guard<std::mutex> lock(model->modelMutex);
    }
    catch (...) {
        return 1;
    }

    model->otfParams.desorptionLimit = 0;

    // Skip desorptions if limit was already reached
    if (!CLIArguments::desLimit.empty()) {
        size_t oldDesNb = globState.globalStats.globalHits.nbDesorbed;
        size_t listSize = CLIArguments::desLimit.size();
        for (size_t l = 0; l < listSize; ++l) {
            model->otfParams.desorptionLimit = CLIArguments::desLimit.front();
            CLIArguments::desLimit.pop_front();

            if (oldDesNb > model->otfParams.desorptionLimit) {
                Log::console_msg_master(1, "Skipping desorption limit: {}\n", model->otfParams.desorptionLimit);
            } else {
                Log::console_msg_master(1, "Starting with desorption limit: {} from {}\n",
                                        model->otfParams.desorptionLimit, oldDesNb);

                
                return 0;
            }
        }
        if (CLIArguments::desLimit.empty()) {
            Log::console_msg_master(1,
                                    "All given desorption limits have been reached. Consider resetting the simulation results from the input file (--reset): Starting desorption {}\n",
                                    oldDesNb);
            
            return 1;
        }
    }
    return 0;
}

/**
* \brief Initializes time limit for the simulation
 * \return 0> error code, 0 when ok
 */
int Initializer::initTimeLimit(std::shared_ptr<MolflowSimulationModel> model, double time) {
    try { //unti initTimeLimit will throw error
        std::lock_guard<std::mutex> lock(model->modelMutex);
    }
    catch (...) {
        return 1;
    }

    model->otfParams.timeLimit = time;
    CLIArguments::simDuration = static_cast<size_t>(time);

    return 0;
}

// TODO: Combine with loadXML function
/**
* \brief Prepares autosave file that is created in user specified intervals
 * \return output file name
 */
std::string Initializer::getAutosaveFile() {
    // Create copy of input file for autosave
    std::string autoSave;
    if (CLIArguments::autoSaveDuration > 0) {
        autoSave = std::filesystem::path(SettingsIO::workFile).filename().string();

        std::string autoSavePrefix = "autosave_";
        // Check if autosave_ is part of the input filename, if yes, generate a new input file without the prefix
        if (autoSave.size() > autoSavePrefix.size() &&
            std::search(autoSave.begin(), autoSave.begin() + autoSavePrefix.size(), autoSavePrefix.begin(),
                        autoSavePrefix.end()) == autoSave.begin()) {
            // TODO: Revisit wether input/output is acceptable here
            autoSave = std::filesystem::path(SettingsIO::workPath).append(SettingsIO::workFile).filename().string();
            SettingsIO::inputFile = autoSave.substr(autoSavePrefix.size(), autoSave.size() - autoSavePrefix.size());
            Log::console_msg_master(2, "Using autosave file {} for {}\n", autoSave, SettingsIO::inputFile);
        } else {
            // create autosavefile from copy of original
            autoSave = std::filesystem::path(SettingsIO::workPath).append(autoSavePrefix).concat(autoSave).string();
            try {
                if(!SettingsIO::workFile.empty() && std::filesystem::exists(SettingsIO::workFile))
                    std::filesystem::copy_file(SettingsIO::workFile, autoSave,
                                           std::filesystem::copy_options::overwrite_existing);
            } catch (std::filesystem::filesystem_error &e) {
                Log::console_error("Could not copy file to create autosave file: {}\n", e.what());
            }
        }
    }

    return autoSave;
}

/**
* \brief Prepares data structures for use in simulation
* \return error code: 0=no error, 1=error
*/
int Initializer::initSimModel(std::shared_ptr<MolflowSimulationModel> model) {

    try { //unti initSimModel will throw error
        std::lock_guard<std::mutex> lock(model->modelMutex);
    }
    catch (...) {
        return 1;
    }

    model->structures.resize(model->sh.nbSuper); //Create structures

    bool hasVolatile = false;

    for (size_t facIdx = 0; facIdx < model->sh.nbFacet; facIdx++) {
        auto mfFac = std::dynamic_pointer_cast<MolflowSimFacet>(model->facets[facIdx]);

        std::vector<double> textIncVector;
        // Add surface elements area (reciprocal)
        if (mfFac->sh.isTextured) {
            auto meshAreas = mfFac->InitTextureMesh();
            textIncVector.resize(mfFac->sh.texHeight * mfFac->sh.texWidth);

            double rw = mfFac->sh.U.Norme() / (double) (mfFac->sh.texWidth_precise);
            double rh = mfFac->sh.V.Norme() / (double) (mfFac->sh.texHeight_precise);
            double area = rw * rh;
            area *= (mfFac->sh.is2sided) ? 2.0 : 1.0;
            size_t add = 0;

            for (size_t j = 0; j < mfFac->sh.texHeight; j++) {
                for (size_t i = 0; i < mfFac->sh.texWidth; i++) {
                    if (meshAreas[add] < 0.0) {
                        textIncVector[add] = 1.0 / area;
                    } else {
                        textIncVector[add] = 1.0 / (meshAreas[add] * ((mfFac->sh.is2sided) ? 2.0 : 1.0));
                    }
                    add++;
                }
            }
        }
        mfFac->textureCellIncrements = textIncVector;

        //Some initialization
        try {
            if (!mfFac->InitializeOnLoad(facIdx)) return false;
        }
        catch (const std::exception& err){
            Log::console_error("Failed to initialize facet (F#{})\n{}\n", facIdx + 1, err.what());
            return 1;
        }
        hasVolatile |= mfFac->sh.isVolatile;

        if ((mfFac->sh.superDest || mfFac->sh.isVolatile) &&
            ((mfFac->sh.superDest - 1) >= model->sh.nbSuper || mfFac->sh.superDest < 0)) {
            // Geometry error
            Log::console_error("Invalid structure (wrong link on F#{})\n", facIdx + 1);
            return 1;
        }
    }

    return 0;
}