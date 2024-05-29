

#include "Initializer.h"
#include "IO/LoaderXML.h"
#include "ParameterParser.h"
#include <Helper/GLProgress_CLI.hpp>
#include <Simulation/MolflowSimFacet.h>
#include "GLApp/GLTypes.h" //error

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
    long double deslimit_double=0.0;
    SettingsIO::CLIArguments parsedArgs;

    // Define options
    app.add_option("-j,--threads", parsedArgs.nbThreads, "Number of parallel threads to be used for calculation");
    app.add_option("-t,--time", parsedArgs.simDuration, "Finish simulation after this time (in seconds)");
    app.add_option("-d,--ndes", deslimit_double, "Finish simulation after this number of desorbed particles");

    app.add_option("-f,--file", parsedArgs.inputFile, "Input file to load (XML/ZIP only)")
            ->required()
            ->check(CLI::ExistingFile);

    CLI::Option *optOfile = app.add_option("-o,--output", parsedArgs.outputFile,
                                           "Output file name. Example: 'result.zip', default: 'out_{inputFileName}'");
    CLI::Option *optOpath = app.add_option("--outputPath", parsedArgs.outputPath,
                                           "Output path, default: \'Results_{timestamp}\'");
    app.add_option("-s,--statprintInterval", parsedArgs.statprintInterval, fmt::format("Statistics print interval in seconds, 0=disabled, default:{}", parsedArgs.statprintInterval));
    app.add_option("-a,--autosaveInterval", parsedArgs.autoSaveInterval, fmt::format("Autosave interval in seconds, 0=disabled, default:600",parsedArgs.autoSaveInterval));
    app.add_flag("--writeFacetDetails", parsedArgs.outputFacetDetails,
                   "When ready, write a CSV file containing all facet details (incl. physical quantities)");
    app.add_flag("--writeFacetQuantities", parsedArgs.outputFacetQuantities,
                   "When ready, write a CSV file containing all physical quantities for each facet");

    app.add_option("--setParamsByFile", parsedArgs.paramFile,
                   "Parameter file for ad hoc change of the given geometry parameters")
            ->check(CLI::ExistingFile);
    app.add_option("--setParams", parsedArgs.paramChanges,
                   "Direct parameter input for ad hoc change of the given geometry parameters");
    app.add_option("--verbosity", AppSettings::verbosity, fmt::format("Console output verbosity. Levels: 0-4. Default:{}",AppSettings::verbosity));
    app.add_flag("--noProgress", parsedArgs.noProgress, "No percentage updates of progress (useful when output is saved to a log file)");
    app.add_flag("--loadAutosave", parsedArgs.loadAutosave, "Load and continue autosave_ file if exists");
    app.add_flag("-r,--reset", parsedArgs.resetOnStart, "Reset input file simulation (to 0 des.) before starting");
    app.add_flag("--verbose", verbose, "Verbose console output (equivalent to --verbosity 4)");
    CLI::Option *optOverwrite = app.add_flag("--overwrite", parsedArgs.overwrite,
                                             "Overwrite input file with results")->excludes(optOfile, optOpath);
    optOfile->excludes(optOverwrite);
    optOpath->excludes(optOverwrite);
    app.set_config("--config");

    try {
        app.parse(argc,argv);
    }
    catch (const CLI::ParseError& e) {
        app.exit(e,std::cout,std::cerr);
        throw Error("");//Command above printed error
    }

    if (verbose)
        AppSettings::verbosity = 42;

    constexpr long double max_sizet_double = static_cast<long double>(std::numeric_limits<size_t>::max());
    if (deslimit_double < 0.0 || deslimit_double > max_sizet_double) {
        throw Error("Invalid desorption limit {:.4g}, must be larger than 0 and smaller than {:.4g}\n", deslimit_double,max_sizet_double);
    }
    parsedArgs.desLimit = static_cast<size_t>(deslimit_double);
    if (parsedArgs.simDuration == 0 && parsedArgs.desLimit==0) {
        throw Error("No stop condition set. Use either -t (--time) or -d (--ndes)\n");
    }
    return parsedArgs;
}


// \brief Comfort function wrapping around default initialization from command line arguments, taking care of default initializations

SettingsIO::CLIArguments Initializer::initFromArgv(int argc, char **argv, SimulationManager& simManager, const std::shared_ptr<MolflowSimulationModel> model) {

// Set local to parse input files the same on all systems
//duplicate, in case we called this function from the test suite and not from main()
#if defined(__APPLE__)
    setlocale(LC_ALL, "en_US.UTF-8");
#else
    std::setlocale(LC_ALL, "en_US.UTF-8");
#endif

#ifdef _WIN32
    //Set console to UTF-8
    SetConsoleCP(CP_UTF8);
    SetConsoleOutputCP(CP_UTF8);
#endif

    Log::console_header(1, "Parsing and applying arguments...\n");
    SettingsIO::CLIArguments parsedArgs;
    parsedArgs = parseArguments(argc, argv);

    simManager.nbThreads = parsedArgs.nbThreads;
    simManager.noProgress = parsedArgs.noProgress;
    model->otfParams.timeLimit = (double)parsedArgs.simDuration;
    Log::console_msg_master(1, "Creating thread(s)...");
    if (simManager.SetUpSimulation()) { //currently only calls CreateCPUHandle()
       throw Error("Error creating thread(s)");
    }
    Log::console_msg_master(1, "{} thread{} created.\n", simManager.nbThreads, 
        simManager.nbThreads>1 ? "s" : "");
    if (parsedArgs.simDuration != 0) {
        Log::console_msg_master(1, "Running simulation for: {} sec\n", parsedArgs.simDuration);
    }
    return parsedArgs;
}

/**
* \brief Initializes the simulation model from a valid input file and handles parameter changes
 */
std::shared_ptr<MolflowSimulationModel>  Initializer::initFromFile(SimulationManager& simManager, const std::shared_ptr<MolflowSimulationModel> model,
    const std::shared_ptr<GlobalSimuState> globalState, MolflowInterfaceSettings& interfaceSettings, SettingsIO::CLIArguments& parsedArgs) {
    SettingsIO::prepareIO(parsedArgs);

    std::shared_ptr<MolflowSimulationModel> loadedModel = std::make_shared<MolflowSimulationModel>();
    if (std::filesystem::path(parsedArgs.workFile).extension() == ".xml") {
        loadedModel = CLILoadFromXML(parsedArgs.workFile, !parsedArgs.resetOnStart, model, globalState, interfaceSettings, parsedArgs);
    }
    else {
        throw Error("Invalid file extension for input file detected: {}\n",
                           std::filesystem::path(parsedArgs.workFile).extension().string());
    }
    if (!parsedArgs.paramFile.empty() || !parsedArgs.paramChanges.empty()) {
        // Apply parameter changes from file
        if (!parsedArgs.paramFile.empty())
            ParameterParser::ParseFile(parsedArgs.paramFile, interfaceSettings.selections);
        if (!parsedArgs.paramChanges.empty())
            ParameterParser::ParseInput(parsedArgs.paramChanges, interfaceSettings.selections);
        ParameterParser::ChangeSimuParams(loadedModel->sp);
        if(ParameterParser::ChangeFacetParams(loadedModel->facets)){
            throw Error("Error in ParameterParser::ChangeFacetParams()");
        }
    }

    loadedModel->CalcTotalOutgassing(); //Gas mass or facet outgassing / temperature might have changed by parameter changes

    // Verify and set desorption limit (if used)
    if (initDesLimit(loadedModel, globalState, parsedArgs)) {
        throw Error("Error in Initializer::initDesLimit");
    }

    simManager.simulationChanged = true;
    Log::console_msg_master(2, "Initializing simulation... ");
    try {
        simManager.InitSimulation(loadedModel, globalState); //shares model and globalState with simManager
    }
    catch (std::exception &ex) {
        Log::console_error("\nFailed Initializing simulation units:\n{}\n", ex.what());
        throw ex;
    }
    Log::console_footer(1, "done.\n");
    return loadedModel;
}

/**
* \brief Initialize simulation from automatically generated test case (prism)
 * \return 0> error code, 0 when ok
 */
/*
void Initializer::initAutoGenerated(SimulationManager *simManager, const std::shared_ptr<MolflowSimulationModel> model,
                                   const std::shared_ptr<GlobalSimuState> globalState, double ratio, int steps, double angle) {
    if (parsedArgs.prepareIO()) {
        throw Error("Error preparing I/O folders\n");
    }

    loadFromGeneration(model, globalState, ratio, steps, angle);

    // Set desorption limit if used
    if (initDesLimit(model, globalState)) {
        throw Error("Failed to set desorption limit(s).");
    }

    simManager.simulationChanged = true;
    Log::console_msg_master(2, "Forwarding model to simulation units...\n");
    simManager.InitSimulation(model, globalState);
    Log::console_footer(1, "Forwarded successfully.\n");
}*/

/**
* \brief Generate a test case with an oblique prism given an angle
 * \return 0> error code, 0 when ok
 */
/*
int
Initializer::loadFromGeneration(const std::shared_ptr<MolflowSimulationModel> model, const std::shared_ptr<GlobalSimuState> globalState, double ratio,
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
    globalState->Resize(model);
    return 0;
}
*/

/**
* \brief Wrapper for CLI's initFromFile for XML loading with XmlLoader
 */
std::shared_ptr<MolflowSimulationModel> Initializer::CLILoadFromXML(const std::string &fileName, bool loadSimulationState, const std::shared_ptr<MolflowSimulationModel> model,
    const std::shared_ptr<GlobalSimuState> globalState, MolflowInterfaceSettings& interfaceSettings, SettingsIO::CLIArguments& parsedArgs) {

    std::shared_ptr<MolflowSimulationModel> loadedModel;
    GLProgress_CLI prg(fmt::format("Loading geometry from file {}\n", fileName));
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
            loadedModel->otfParams.timeLimit=model->otfParams.timeLimit;
        }
        catch (Error& err) {
            Log::console_error("Error loading file:\n{}", err.what());
            exit(43);
        }
        interfaceSettings = *loader.interfaceSettings; //persistent for saving
    }
    loadedModel->memSizeCache = loadedModel->GetMemSize();
    prg.SetMessage(fmt::format("Loaded geometry ({:.1f} Mbytes)", (double)loadedModel->memSizeCache/1024.0/1024.0));

    //InitializeGeometry();
    loadedModel->InitializeFacets();

    prg.SetMessage("Initializing geometry...",true);
    try {
        initSimModel(loadedModel);
    }
    catch (Error& err) {
        Log::console_error("Couldn't initialize loaded model:\n{}", err.what());
        exit(43);
    }

    // 2. Create simulation dataports
    try {

        prg.SetMessage("Resizing global state counters...");
        globalState->Resize(loadedModel);
        prg.SetProgress(1.0);

        // 3. init counters with previous results
        if (loadSimulationState) {
            prg.SetMessage("Loading simulation state...");
            if (parsedArgs.loadAutosave) {
                std::string autosaveFileName = std::filesystem::path(parsedArgs.workFile).filename().string();
                std::string autoSavePrefix = "autosave_";
                autosaveFileName = autoSavePrefix + autosaveFileName;
                
                if (std::filesystem::exists(autosaveFileName)) {
                    prg.SetMessage(fmt::format("Found autosave file {}, loading simulation state...\n",autosaveFileName),true);
                    FlowIO::XmlLoader::LoadSimulationState(autosaveFileName, loadedModel, globalState, prg);
                }
            } else {
                FlowIO::XmlLoader::LoadSimulationState(parsedArgs.workFile, loadedModel, globalState, prg);
            }
            prg.SetProgress(1.0);

            // Update Angle map status
            for(int i = 0; i < loadedModel->facets.size(); i++ ) {
#if defined(MOLFLOW)
                auto f = std::static_pointer_cast<MolflowSimFacet>(loadedModel->facets[i]);
                if (f->sh.anglemapParams.record) { //Recording, so needs to be updated
                    //Retrieve angle map from hits dp
                    globalState->facetStates[i].recordedAngleMapPdf = f->angleMap.pdf;
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
int Initializer::initDesLimit(const std::shared_ptr<MolflowSimulationModel> model, const std::shared_ptr<GlobalSimuState> globalState, SettingsIO::CLIArguments& parsedArgs) {
    
    try { //unti initDesLimit will throw error
        std::lock_guard<std::mutex> lock(model->modelMutex);
    }
    catch (...) {
        return 1;
    }
    
    model->otfParams.desorptionLimit = parsedArgs.desLimit;
    if (model->otfParams.desorptionLimit > 0) {
        //Check if already reached
        size_t oldDesNb = globalState->globalStats.globalHits.nbDesorbed;
        if (oldDesNb >= model->otfParams.desorptionLimit) {
            Log::console_msg_master(1,
                "Desorption limit ({}) already reached (input file nbDes={}). Consider resetting the simulation (--reset argument)\n",
                model->otfParams.desorptionLimit, oldDesNb);
            return 1;
        }
    }
    return 0; //Success
}

/**
* \brief Initializes time limit for the simulation
 * \return 0> error code, 0 when ok
 */
int Initializer::initTimeLimit(const std::shared_ptr<MolflowSimulationModel> model, SettingsIO::CLIArguments& parsedArgs) {
    try { //unti initTimeLimit will throw error
        std::lock_guard<std::mutex> lock(model->modelMutex);
    }
    catch (...) {
        return 1;
    }

    model->otfParams.timeLimit = static_cast<double>(parsedArgs.simDuration);
    //parsedArgs.simDuration = static_cast<size_t>(parsedArgs.simDuration);

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
*/
void Initializer::initSimModel(const std::shared_ptr<MolflowSimulationModel> model) {

    try { //unti initSimModel will throw error
        std::lock_guard<std::mutex> lock(model->modelMutex);
    }
    catch (...) {
        throw Error("Couldn't lock model mutex");
    }

    model->initialized = false;
    model->structures.resize(model->sh.nbSuper); //Create structures

    for (size_t facIdx = 0; facIdx < model->sh.nbFacet; facIdx++) {
        auto mfFac = std::static_pointer_cast<MolflowSimFacet>(model->facets[facIdx]);

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
            mfFac->InitializeOnLoad(facIdx, model->tdParams);
        }
        catch (const std::exception& err){
            throw Error("Failed to initialize facet (F#{})\n{}\n", facIdx + 1, err.what());
        }

        if ((mfFac->sh.superDest) &&
            ((mfFac->sh.superDest - 1) >= model->sh.nbSuper || mfFac->sh.superDest < 0)) {
            // Geometry error
            throw Error("Invalid structure (wrong link on F#{})\n", facIdx + 1);
        }
    }

    model->PrepareToRun();
}