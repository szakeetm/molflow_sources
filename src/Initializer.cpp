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
SettingsIO::CLIArguments Initializer::parseArguments(int argc, char **argv) {
    CLI::App app{"Molflow+/Synrad+ Simulation Management"};
    app.formatter(std::make_shared<FlowFormatter>());

    // Local variables for parsing and immediate processing
    bool verbose = false;
    std::vector<double> limits;

    SettingsIO::CLIArguments parsedArgs;

    // Define options
    app.add_option("-j,--threads", parsedArgs.nbThreads, "# Threads to be deployed");
    app.add_option("-t,--time", parsedArgs.simDuration, "Simulation duration in seconds");
    app.add_option("-d,--ndes", limits, "Desorption limit for simulation end");

    auto group = app.add_option_group("subgroup");
    group->add_option("-f,--file", parsedArgs.inputFile, "Required input file (XML/ZIP only)")
            ->check(CLI::ExistingFile);
    group->require_option(1);

    CLI::Option *optOfile = app.add_option("-o,--output", parsedArgs.outputFile,
                                           R"(Output file name (e.g. 'outfile.xml', defaults to 'out_{inputFileName}')");
    CLI::Option *optOpath = app.add_option("--outputPath", parsedArgs.outputPath,
                                           "Output path, defaults to \'Results_{date}\'");
    app.add_option("-s,--statprintInterval", parsedArgs.statprintInterval, "Seconds between each stat output, zero=disabled");
    app.add_option("-a,--autosaveInterval", parsedArgs.autoSaveInterval, "Autosave interval in seconds, zero=disabled)");
    app.add_flag("--writeFacetDetails", parsedArgs.outputFacetDetails,
                   "Will write a CSV file containing all facet details including physical quantities");
    app.add_flag("--writeFacetQuantities", parsedArgs.outputFacetQuantities,
                   "Will write a CSV file containing all physical quantities for each facet");

    app.add_option("--setParamsByFile", parsedArgs.paramFile,
                   "Parameter file for ad hoc change of the given geometry parameters")
            ->check(CLI::ExistingFile);
    app.add_option("--setParams", parsedArgs.paramChanges,
                   "Direct parameter input for ad hoc change of the given geometry parameters");
    app.add_option("--verbosity", AppSettings::verbosity, "Restrict console output to different levels");
    app.add_flag("--noProgress", parsedArgs.noProgress, "Log file mode: No percentage updates printed of progress");
    app.add_flag("--loadAutosave", parsedArgs.loadAutosave, "Whether autosave_ file should be used if exists");
    app.add_flag("-r,--reset", parsedArgs.resetOnStart, "Resets simulation status loaded from file");
    app.add_flag("--verbose", verbose, "Verbose console output (all levels)");
    CLI::Option *optOverwrite = app.add_flag("--overwrite", parsedArgs.overwrite,
                                             "Overwrite input file with new results")->excludes(optOfile, optOpath);
    optOfile->excludes(optOverwrite);
    optOpath->excludes(optOverwrite);
    app.set_config("--config");

    try {
        app.parse(argc,argv);
    }
    catch (const CLI::ParseError& e) {
        throw Error(fmt::format("Argument parse error: {}", e.what()));
    }

    if (verbose)
        AppSettings::verbosity = 42;

    //std::cout<<app.config_to_str(true,true);
    for (auto& lim : limits)
        parsedArgs.desLimit.push_back(static_cast<size_t>(lim));

    if (parsedArgs.simDuration == 0 && parsedArgs.desLimit.empty()) {
        throw Error("No end criterion has been set. Use either -t or -d\n");
    }
    return parsedArgs;
}


// \brief Comfort function wrapping around default initialization from command line arguments, taking care of default initializations

SettingsIO::CLIArguments Initializer::initFromArgv(int argc, char **argv, SimulationManager& simManager, std::shared_ptr<MolflowSimulationModel> model) {

// Set local to parse input files the same on all systems
//duplicate, in case we called this function from the test suite and not from main()
#if defined(__APPLE__)
    setlocale(LC_ALL, "en_US.UTF-8");
#else
    std::setlocale(LC_ALL, "en_US.UTF-8");
#endif

    SettingsIO::CLIArguments parsedArgs;
    parsedArgs = parseArguments(argc, argv);

    Log::console_header(1, "Initializing simulation...\n");

    simManager.nbThreads = parsedArgs.nbThreads;
    simManager.noProgress = parsedArgs.noProgress;
    model->otfParams.timeLimit = (double)parsedArgs.simDuration;

    if (simManager.SetUpSimulation()) { //currently only calls CreateCPUHandle()
       throw Error(fmt::format("Error: Setting up simulation units [{} threads]...\n", simManager.nbThreads));
    }

    Log::console_msg_master(1, "Active cores: {}\n", simManager.nbThreads);
    if (parsedArgs.simDuration != 0) {
        Log::console_msg_master(1, "Running simulation for: {} sec\n", parsedArgs.simDuration);
    }
    return parsedArgs;
}

/**
* \brief Initializes the simulation model from a valid input file and handles parameter changes
 */
std::shared_ptr<MolflowSimulationModel>  Initializer::initFromFile(SimulationManager& simManager, std::shared_ptr<MolflowSimulationModel> model,
                              GlobalSimuState *globStatePtr, MolflowUserSettings& userSettings, SettingsIO::CLIArguments& parsedArgs) {
    SettingsIO::prepareIO(parsedArgs);

    std::shared_ptr<MolflowSimulationModel> loadedModel = std::make_shared<MolflowSimulationModel>();
    if (std::filesystem::path(parsedArgs.workFile).extension() == ".xml") {
        loadedModel = CLILoadFromXML(parsedArgs.workFile, !parsedArgs.resetOnStart, model, globStatePtr, userSettings, parsedArgs);
    }
    else {
        throw Error(fmt::format("Invalid file extension for input file detected: {}\n",
                           std::filesystem::path(parsedArgs.workFile).extension().string()));
    }
    if (!parsedArgs.paramFile.empty() || !parsedArgs.paramChanges.empty()) {
        // Apply parameter changes from file
        if (!parsedArgs.paramFile.empty())
            ParameterParser::ParseFile(parsedArgs.paramFile, userSettings.selections);
        if (!parsedArgs.paramChanges.empty())
            ParameterParser::ParseInput(parsedArgs.paramChanges, userSettings.selections);
        ParameterParser::ChangeSimuParams(loadedModel->wp);
        if(ParameterParser::ChangeFacetParams(loadedModel->facets)){
            throw Error("Error in ParameterParser::ChangeFacetParams()");
        }
    }

    // Set desorption limit if used
    if (initDesLimit(loadedModel, *globStatePtr, parsedArgs)) {
        throw Error("Error in Initializer::initDesLimit");
    }

    simManager.simulationChanged = true;
    Log::console_msg_master(2, "Forwarding model to simulation units...\n");
    try {
        simManager.InitSimulation(loadedModel, globStatePtr);
    }
    catch (std::exception &ex) {
        Log::console_error("Failed Initializing simulation units:\n{}\n", ex.what());
        throw ex;
    }
    Log::console_footer(1, "Forwarded successfully.\n");
    return loadedModel;
}

/**
* \brief Initialize simulation from automatically generated test case (prism)
 * \return 0> error code, 0 when ok
 */
/*
void Initializer::initAutoGenerated(SimulationManager *simManager, std::shared_ptr<MolflowSimulationModel> model,
                                   GlobalSimuState *globStatePtr, double ratio, int steps, double angle) {
    if (parsedArgs.prepareIO()) {
        throw Error("Error preparing I/O folders\n");
    }

    loadFromGeneration(model, globStatePtr, ratio, steps, angle);

    // Set desorption limit if used
    if (initDesLimit(model, *globStatePtr)) {
        throw Error("Failed to set desorption limit(s).");
    }

    simManager.simulationChanged = true;
    Log::console_msg_master(2, "Forwarding model to simulation units...\n");
    simManager.InitSimulation(model, globStatePtr);
    Log::console_footer(1, "Forwarded successfully.\n");
}*/

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
* \brief Wrapper for CLI's initFromFile for XML loading with XmlLoader
 */
std::shared_ptr<MolflowSimulationModel> Initializer::CLILoadFromXML(const std::string &fileName, bool loadSimulationState, std::shared_ptr<MolflowSimulationModel> model,
                             GlobalSimuState *globStatePtr, MolflowUserSettings& userSettings, SettingsIO::CLIArguments& parsedArgs) {

    std::shared_ptr<MolflowSimulationModel> loadedModel = std::make_shared<MolflowSimulationModel>();
    GLProgress_CLI prg(fmt::format("Loading geometry from file {}", fileName));
    prg.noProgress = parsedArgs.noProgress; //Suppress percentages if ran in script

    //1. Load Input File (regular XML)
    {
        FlowIO::XmlLoader loader;
        // Easy if XML is split into multiple parts
        // Geometry
        // Settings
        // Previous results
        try {
            loadedModel = loader.LoadGeometry(fileName,
            TimeDependentParameters::GetCatalogParameters(model->tdParams.parameters), prg);
        }
        catch (Error& err) {
            Log::console_error("Error loading file:\n{}", err.what());
        }
        userSettings = loader.userSettings; //persistent for saving
    }

    prg.SetMessage(fmt::format("Loaded geometry ({} Mbytes)", loadedModel->size()/1024/1024));

    //InitializeGeometry();
    loadedModel->InitializeFacets();

    prg.SetMessage("Initializing geometry...",true);
    initSimModel(loadedModel);
    loadedModel->PrepareToRun();

    // 2. Create simulation dataports
    try {

        prg.SetMessage("Resizing global state counters...");
        globStatePtr->Resize(loadedModel);

        // 3. init counters with previous results
        if (loadSimulationState) {
            prg.SetMessage("Loading simulation state...");
            if (parsedArgs.loadAutosave) {
                std::string autosaveFileName = std::filesystem::path(parsedArgs.workFile).filename().string();
                std::string autoSavePrefix = "autosave_";
                autosaveFileName = autoSavePrefix + autosaveFileName;
                
                if (std::filesystem::exists(autosaveFileName)) {
                    prg.SetMessage(fmt::format("Found autosave file {}, loading simulation state...\n",autosaveFileName),true);
                    FlowIO::XmlLoader::LoadSimulationState(autosaveFileName, loadedModel, globStatePtr, prg);
                }
            } else {
                FlowIO::XmlLoader::LoadSimulationState(parsedArgs.workFile, loadedModel, globStatePtr, prg);
            }

            // Update Angle map status
            for(int i = 0; i < loadedModel->facets.size(); i++ ) {
#if defined(MOLFLOW)
                auto f = std::dynamic_pointer_cast<MolflowSimFacet>(loadedModel->facets[i]);
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

    return loadedModel;
}

/**
* \brief Initialize a desorption limit by checking against the given input (vector)
 * \return 0> error code, 0 when ok
 */
int Initializer::initDesLimit(std::shared_ptr<MolflowSimulationModel> model, GlobalSimuState &globState, SettingsIO::CLIArguments& parsedArgs) {
    try { //unti initDesLimit will throw error
        std::lock_guard<std::mutex> lock(model->modelMutex);
    }
    catch (...) {
        return 1;
    }

    model->otfParams.desorptionLimit = 0;

    // Skip desorptions if limit was already reached
    if (!parsedArgs.desLimit.empty()) {
        size_t oldDesNb = globState.globalStats.globalHits.nbDesorbed;
        size_t listSize = parsedArgs.desLimit.size();
        for (size_t l = 0; l < listSize; ++l) {
            model->otfParams.desorptionLimit = parsedArgs.desLimit.front();
            parsedArgs.desLimit.pop_front();

            if (oldDesNb > model->otfParams.desorptionLimit) {
                Log::console_msg_master(1, "Skipping desorption limit: {}\n", model->otfParams.desorptionLimit);
            } else {
                Log::console_msg_master(1, "Starting with desorption limit: {} from {}\n",
                                        model->otfParams.desorptionLimit, oldDesNb);

                
                return 0;
            }
        }
        if (parsedArgs.desLimit.empty()) {
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
int Initializer::initTimeLimit(std::shared_ptr<MolflowSimulationModel> model, SettingsIO::CLIArguments& parsedArgs) {
    try { //unti initTimeLimit will throw error
        std::lock_guard<std::mutex> lock(model->modelMutex);
    }
    catch (...) {
        return 1;
    }

    model->otfParams.timeLimit = parsedArgs.simDuration;
    parsedArgs.simDuration = static_cast<size_t>(parsedArgs.simDuration);

    return 0;
}

// TODO: Combine with loadXML function
/**
* \brief Prepares autosave file that is created in user specified intervals
 * \return output file name
 */
std::string Initializer::getAutosaveFile(SettingsIO::CLIArguments& parsedArgs) {
    // Create copy of input file for autosave
    std::string autoSave;
    if (parsedArgs.autoSaveInterval > 0) {
        autoSave = std::filesystem::path(parsedArgs.workFile).filename().string();

        std::string autoSavePrefix = "autosave_";
        // Check if autosave_ is part of the input filename, if yes, generate a new input file without the prefix
        if (autoSave.size() > autoSavePrefix.size() &&
            std::search(autoSave.begin(), autoSave.begin() + autoSavePrefix.size(), autoSavePrefix.begin(),
                        autoSavePrefix.end()) == autoSave.begin()) {
            // TODO: Revisit wether input/output is acceptable here
            autoSave = std::filesystem::path(parsedArgs.workPath).append(parsedArgs.workFile).filename().string();
            parsedArgs.inputFile = autoSave.substr(autoSavePrefix.size(), autoSave.size() - autoSavePrefix.size());
            Log::console_msg_master(2, "Using autosave file {} for {}\n", autoSave, parsedArgs.inputFile);
        } else {
            // create autosavefile from copy of original
            autoSave = std::filesystem::path(parsedArgs.workPath).append(autoSavePrefix).concat(autoSave).string();
            try {
                if(!parsedArgs.workFile.empty() && std::filesystem::exists(parsedArgs.workFile))
                    std::filesystem::copy_file(parsedArgs.workFile, autoSave,
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