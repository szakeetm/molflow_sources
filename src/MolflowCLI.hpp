static constexpr const char* molflowCliLogo = R"(
  __  __     _  __ _
 |  \/  |___| |/ _| |_____ __ __
 | |\/| / _ \ |  _| / _ \ V  V /
 |_|  |_\___/_|_| |_\___/\_/\_/
    )"; //Unused, clutters iterative simulations and parameter sweeps

class RuntimeStatPrinter {
    size_t oldHitsNb=0;
    size_t oldDesNb=0;
public:
    RuntimeStatPrinter() = default;
    RuntimeStatPrinter(size_t n_hits, size_t n_des);
    void PrintHeader() const;
    void Print(double elapsedTime, GlobalSimuState& globState, bool printSum = false) const;
};

void ShutdownMPI();
void CLIMainLoop(double& elapsedTime, Chronometer& simTimer, std::shared_ptr<MolflowSimulationModel> model, GlobalSimuState& simuState,
    SimulationManager& simManager, MolflowUserSettings& persistentUserSettings, std::string& autoSave, RuntimeStatPrinter& printer, SettingsIO::CLIArguments& parsedArgs);
void CleanUpMPI();
void WriteResults(std::shared_ptr<MolflowSimulationModel> model, GlobalSimuState& simuState,
    SimulationManager& simManager, MolflowUserSettings& persistentUserSettings, std::string& autoSave, SettingsIO::CLIArguments& parsedArgs);
void HandleIntermediateDesLimit(std::shared_ptr<MolflowSimulationModel> model, GlobalSimuState& simuState,
    SimulationManager& simManager, MolflowUserSettings& persistentUserSettings, bool& endCondition, SettingsIO::CLIArguments& parsedArgs);
void GatherAngleMapRecordings(MolflowSimulationModel& model, GlobalSimuState& globSimState);