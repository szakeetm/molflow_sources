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

// M_PI define
#ifdef _WIN32
#define _USE_MATH_DEFINES // activate defines, e.g. M_PI_2
#endif
#include <cmath>

#include "SimulationManager.h"
#include "SMP.h"
#include "Buffer_shared.h"
#include <fstream>
#include <filesystem>
#include <Parameter.h>
#include <IO/LoaderXML.h>
#include <IO/WriterXML.h>
#include "Simulation/MolflowSimGeom.h"
#include "Initializer.h"
#include "Helper/MathTools.h"
#include <sstream>
#include <Helper/Chronometer.h>
#include <Helper/StringHelper.h>
#include <Helper/ConsoleLogger.h>
#include <ZipLib/ZipFile.h>
#include "versionId.h"
#include "Helper/GLProgress_CLI.hpp"
#include "MolflowCLI.hpp"

//#if defined(MOLFLOW)
#include <SettingsIO.h>
#include <IO/CSVExporter.h>
//#endif

#include "FlowMPI.h"
#include "File.h"

RuntimeStatPrinter::RuntimeStatPrinter(size_t n_hits, size_t n_des) {
    oldHitsNb = n_hits;
    oldDesNb = n_des;
};
void RuntimeStatPrinter::PrintHeader() const{
    // Print Header at the beginning
    Log::console_msg_master(1, "\n");
    if(MFMPI::world_size > 1) {
        Log::console_msg_master(1, "{:<6} ",
                                "Node#");
    }
    Log::console_msg_master(1, "{:<14} {:<20} {:<20} {:<20} {:<20} {:<20} {:<20}\n",
                            "Time",
                            "#Hits (run)", "#Hits (total)","Hit/sec",
                            "#Des (run)", "#Des (total)","Des/sec");
    if(MFMPI::world_size > 1) {
        Log::console_msg_master(1, "{}",std::string(6,'-'));
    }
    Log::console_msg_master(1, "{}\n",std::string(14+20+20+20+20+20+20,'-'));
}
void RuntimeStatPrinter::Print(double elapsedTime, GlobalSimuState& globState, bool printSum) const{
    if(printSum) {
        Log::console_msg_master(1, "{}\n",std::string(6+14+20+20+20+20+20+20,'='));
        Log::console_msg_master(1, "{:<6} ", "x");
    }
    else if(MFMPI::world_size > 1) {
        Log::console_msg(1, "{:<6} ", MFMPI::world_rank);
    }

    Log::console_msg(1,"{:<14.2f} {:<20} {:<20} {:<20.2f} {:<20} {:<20} {:<20.2f}\n",
                        elapsedTime,
                        globState.globalStats.globalHits.nbMCHit - oldHitsNb, globState.globalStats.globalHits.nbMCHit,
                        (double) (globState.globalStats.globalHits.nbMCHit - oldHitsNb) /
                        (elapsedTime),
                        globState.globalStats.globalHits.nbDesorbed - oldDesNb, globState.globalStats.globalHits.nbDesorbed,
                        (double) (globState.globalStats.globalHits.nbDesorbed - oldDesNb) /
                        (elapsedTime));
}


int main(int argc, char** argv) {

// Set local to parse input files the same on all systems
//duplicate, in case we called this function from the test suite and not from main()
#if defined(__APPLE__)
    setlocale(LC_ALL, "en_US.UTF-8");
#else
    std::setlocale(LC_ALL, "en_US.UTF-8");
#endif

#if defined(USE_MPI)
    MFMPI::mpi_initialize();
#endif

    Log::console_msg_master(1, "{} command line mode\n", appTitle);

    // Init necessary components
    SimulationManager simManager{MFMPI::world_rank};
    std::shared_ptr<MolflowSimulationModel> model = std::make_shared<MolflowSimulationModel>();
    GlobalSimuState simuState{};
    MolflowUserSettings persistentUserSettings; //persistent user data that should be written back to a results file when saving
    SettingsIO::CLIArguments parsedArgs;

    // Parse arguments
    try {
        parsedArgs = Initializer::initFromArgv(argc, argv, simManager, model);
    } catch (Error& err) {
        Log::console_error(err.what());
        ShutdownMPI();
        return 41;
    }

#if defined(USE_MPI)
    // Start to transfer simulation data to other nodes via broadcast
    MFMPI::mpi_transfer_simu(parsedArgs);
#endif

    Log::console_msg(1,"Loading parameters from parameter_catalog folder...");
    TimeDependentParameters::LoadParameterCatalog(model->tdParams.parameters);
    Log::console_msg(1, "done.\n");

    try {
        auto loadedModel = Initializer::initFromFile(simManager, model, &simuState, persistentUserSettings, parsedArgs);
        model = loadedModel;
    }
    catch (std::exception& err) {
        Log::console_error("Initializer::initFromFile error:\n{}", err.what());
        ShutdownMPI();
        return 42;
    }

    if(parsedArgs.simDuration == 0 && model->otfParams.desorptionLimit == 0){
        Log::console_error("Neither a time limit nor a desorption limit has been set.\n");
        return 44;
    }
    size_t oldHitsNb = simuState.globalStats.globalHits.nbMCHit;
    size_t oldDesNb = simuState.globalStats.globalHits.nbDesorbed;
    RuntimeStatPrinter printer(oldHitsNb, oldDesNb);
    // Get autosave file name
    std::string autoSave = Initializer::getAutosaveFile(parsedArgs);

    if (parsedArgs.simDuration > 0) {
        Log::console_msg_master(1, "[{}] Preparing simulation for {} seconds from {} desorptions...\n", Util::getTimepointString(), parsedArgs.simDuration, simuState.globalStats.globalHits.nbDesorbed);
    }
    else if (model->otfParams.desorptionLimit > 0) {
        Log::console_msg_master(1, "[{}] Preparing simulation to {} desorptions from {} desorptions...\n", Util::getTimepointString(), model->otfParams.desorptionLimit, simuState.globalStats.globalHits.nbDesorbed);
    }

#if defined(USE_MPI)
    MPI_Barrier(MPI_COMM_WORLD);
    simManager.noProgress = false;
#endif

    // Start async simulation run, check state in following loop
    try {
        simManager.StartSimulation();
    }
    catch (const std::exception& e) {
        Log::console_error("[{}] ERROR: Starting simulation: {}\n",MFMPI::world_rank, e.what());
        Log::console_error("[{}] File folder {} -- {}\n",MFMPI::world_rank,parsedArgs.workPath, parsedArgs.workFile);
        ShutdownMPI();
        return 43;
    }

    Chronometer simTimer;
    simTimer.Start();
    double elapsedTime = 0.0;

    CLIMainLoop(elapsedTime, simTimer, model, simuState,
        simManager, persistentUserSettings, autoSave, printer, parsedArgs);

    simTimer.Stop();
    elapsedTime = simTimer.Elapsed();

    // Terminate simulation
    simManager.StopSimulation();
    simManager.KillSimulation();
    GatherAngleMapRecordings(*model, simuState);
    Log::console_msg(1,"[{}][{}] Simulation finished.\n", MFMPI::world_rank, Util::getTimepointString());

#if defined(USE_MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    printer.PrintHeader();
#if defined(USE_MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    //TODO: Send output to master node for ordered output
    if(elapsedTime > 1e-4) {
        // Global result print --> TODO: ()
        printer.Print(elapsedTime, simuState);
    }

#if defined(USE_MPI)
    MPI_Barrier(MPI_COMM_WORLD);
    MFMPI::mpi_receive_states(model, simuState);
    if (MFMPI::world_rank != 0) {
        // Cleanup all files from nodes tmp path
        if (parsedArgs.outputPath.find("tmp") != std::string::npos) {
            std::filesystem::remove_all(parsedArgs.outputPath);
        }
        if (std::filesystem::exists(autoSave)) {
            std::filesystem::remove(autoSave);
        }
    }
    // Finalize the MPI environment.
    MPI_Finalize();
    if (MFMPI::world_rank != 0) {
        return 0;
    }
#endif //USE_MPI

    //Print sum
    if(elapsedTime > 1e-4) {
        if(MFMPI::world_size > 1)
            printer.Print(elapsedTime, simuState, true);
    }

    if(MFMPI::world_rank == 0){
        if(parsedArgs.outputFacetDetails) {
            FlowIO::Exporter::export_facet_details(&simuState, model.get(), parsedArgs.workPath);

        }
        if(parsedArgs.outputFacetQuantities) {
            FlowIO::Exporter::export_facet_quantities(&simuState, model.get(), parsedArgs.workPath);
        }

        WriteResults(model, simuState, simManager, persistentUserSettings, autoSave, parsedArgs);
    }

    // Cleanup
    SettingsIO::cleanup_files(parsedArgs.outputPath,parsedArgs.workPath);

    return 0; //success
}

void ShutdownMPI() {
#if defined(USE_MPI)
    MPI_Finalize();
#endif
}

void CLIMainLoop(double& elapsedTime, Chronometer& simTimer, std::shared_ptr<MolflowSimulationModel> model, GlobalSimuState& simuState,
    SimulationManager& simManager, MolflowUserSettings& persistentUserSettings, std::string& autoSave, RuntimeStatPrinter& printer, SettingsIO::CLIArguments& parsedArgs) {
    // Simulation runtime loop to check for end conditions and start auto-saving procedures etc.
    bool endCondition = false;
    Log::console_msg_master(1, "[{}] Started simulation.\n", Util::getTimepointString());
    do {
        ProcessSleep(1000);

        elapsedTime = simTimer.Elapsed();
        if (model->otfParams.desorptionLimit != 0)
            endCondition = simuState.globalStats.globalHits.nbDesorbed/* - oldDesNb*/ >= model->otfParams.desorptionLimit;

        if (endCondition) {
            // if there is a next des limit, handle that
            if (!parsedArgs.desLimit.empty()) {
                HandleIntermediateDesLimit(model, simuState, simManager, persistentUserSettings, endCondition, parsedArgs);
            }
        }
        else if (parsedArgs.autoSaveInterval && simTimer.SecondsSinceLastAutosave() > parsedArgs.autoSaveInterval) { // autosave every x seconds
            // Autosave
            GLProgress_CLI prg(fmt::format("[{:.2}s] Creating auto save file {}", elapsedTime, autoSave));
            prg.noProgress = simManager.noProgress;
            FlowIO::XmlWriter writer;
            writer.AppendSimulationStateToFile(autoSave, model, prg, simuState);
            simTimer.UpdateLastAutoSave();
        }

        if (parsedArgs.statprintInterval && simTimer.SecondsSinceLastStatprint() > parsedArgs.statprintInterval) { // autosave every x seconds
            // Print runtime stats
            if ((uint64_t)elapsedTime / parsedArgs.statprintInterval <= 1) {
                printer.PrintHeader();
            }
            printer.Print(elapsedTime, simuState);
            simTimer.UpdateLastStatprintTime();
        }

        // Check for potential time end
        if (parsedArgs.simDuration > 0) {
            endCondition |= (elapsedTime >= (double)parsedArgs.simDuration);
        }
    } while (!endCondition);
}

void WriteResults(std::shared_ptr<MolflowSimulationModel> model, GlobalSimuState& simuState,
    SimulationManager& simManager, MolflowUserSettings& persistentUserSettings, std::string& autoSave, SettingsIO::CLIArguments& parsedArgs) {
    // Export results
        //  a) Use existing autosave as base
        //  b) Create copy of input file
        // update geometry info (in case of param sweep)
        // and simply update simulation results
    bool createZip = std::filesystem::path(parsedArgs.outputFile).extension() == ".zip";
    parsedArgs.outputFile = std::filesystem::path(parsedArgs.outputFile).replace_extension(".xml").string();

    std::string fullOutFile = std::filesystem::path(parsedArgs.outputPath).append(parsedArgs.outputFile).string();
    if (std::filesystem::exists(autoSave)) {
        std::filesystem::rename(autoSave, fullOutFile);
    }
    else if (!parsedArgs.overwrite) {
        // Copy full file description first, in case outputFile is different
        if (!parsedArgs.workFile.empty() && std::filesystem::exists(parsedArgs.workFile)) {
            try {
                std::filesystem::copy_file(parsedArgs.workFile, fullOutFile,
                    std::filesystem::copy_options::overwrite_existing);
            }
            catch (std::filesystem::filesystem_error& e) {
                Log::console_error("Could not copy file to preserve initial file layout: {}\n", e.what());
            }
        }
    }
    GLProgress_CLI prg(fmt::format("Writing file {} ...", fullOutFile));
    prg.noProgress = simManager.noProgress;
    FlowIO::XmlWriter writer(false, true);
    writer.userSettings = persistentUserSettings;
    pugi::xml_document newDoc;
    newDoc.load_file(fullOutFile.c_str());
    writer.SaveGeometry(newDoc, model, prg);
    writer.SaveSimulationState(newDoc, model, prg, simuState);
    writer.WriteXMLToFile(newDoc, fullOutFile);

    if (createZip) {
        prg.SetMessage("Compressing xml to zip...");

        //Zipper library
        std::string fileNameWithZIP = std::filesystem::path(fullOutFile).replace_extension(".zip").string();
        if (std::filesystem::exists(fileNameWithZIP)) { // should be workFile == inputFile
            try {
                std::filesystem::remove(fileNameWithZIP);
            }
            catch (std::exception& e) {
                Log::console_error("Error compressing to \n{}\nMaybe file is in use:\n{}", fileNameWithZIP, e.what());
            }
        }
        try {
            ZipFile::AddFile(fileNameWithZIP, fullOutFile, FileUtils::GetFilename(fullOutFile));
        }
        catch (std::exception& e) {
            Log::console_error("Error compressing to \n{}\nMaybe file is in use:\n{}", fileNameWithZIP, e.what());
        }

        //At this point, if no error was thrown, the compression is successful
        try {
            prg.SetMessage(fmt::format("Successfully compressed to {}, removing {} ...", fileNameWithZIP, fullOutFile));
            std::filesystem::remove(fullOutFile);
        }
        catch (std::exception& e) {
            Log::console_error("Error removing\n{}\nMaybe file is in use:\n{}", fullOutFile, e.what());
        }
    }
}

void HandleIntermediateDesLimit(std::shared_ptr<MolflowSimulationModel> model, GlobalSimuState& simuState,
    SimulationManager& simManager, MolflowUserSettings& persistentUserSettings, bool& endCondition, SettingsIO::CLIArguments& parsedArgs) {
    // First write an intermediate output file
                // 1. Get file name
    std::string outFile = std::filesystem::path(parsedArgs.outputPath)
        .append("desorbed_")
        .concat(std::to_string(model->otfParams.desorptionLimit))
        .concat("_")
        .concat(std::filesystem::path(parsedArgs.outputFile).filename().string()).string();

    try {
        GLProgress_CLI prg(fmt::format("Saving intermediate results... {}", outFile));
        prg.noProgress = simManager.noProgress;
        // 2. Write XML file, use existing file as base or create new file
        FlowIO::XmlWriter writer;
        writer.userSettings = persistentUserSettings; //keep from loaded file
        if (!parsedArgs.workFile.empty() && std::filesystem::exists(parsedArgs.workFile)) {
            try {
                std::filesystem::copy_file(parsedArgs.workFile, outFile,
                    std::filesystem::copy_options::overwrite_existing);
            }
            catch (std::filesystem::filesystem_error& e) {
                Log::console_error("Could not copy file: {}\n", e.what());
            }
        }
        else {
            pugi::xml_document newDoc;
            writer.SaveGeometry(newDoc, model, prg);
            writer.WriteXMLToFile(newDoc, outFile);
        }
        // 3. append updated results
        writer.AppendSimulationStateToFile(outFile, model, prg, simuState);
    }
    catch (std::filesystem::filesystem_error& e) {
        Log::console_error("Warning: Could not create file: {}\n", e.what());
    }
    // Next choose the next desorption limit and start

    model->otfParams.desorptionLimit = parsedArgs.desLimit.front();
    parsedArgs.desLimit.pop_front();
    simManager.ForwardOtfParams(&model->otfParams);
    endCondition = false;
    Log::console_msg_master(1, " Handling next des limit {}\n", model->otfParams.desorptionLimit);

    try {
        ProcessSleep(1000);
        simManager.StartSimulation();
    }
    catch (const std::exception& e) {
        Log::console_error("ERROR: Starting simulation: {}\n", e.what());
        endCondition = true;
    }
}

//Transfers recorded angle maps from "simulation state" to "model"
void GatherAngleMapRecordings(MolflowSimulationModel& model, GlobalSimuState& globSimState) {
    for (int i = 0; i < model.facets.size(); i++) {
#if defined(MOLFLOW)
        auto f = std::dynamic_pointer_cast<MolflowSimFacet>(model.facets[i]);
        if (f->sh.anglemapParams.record) { //Recording, so needs to be updated
            //Retrieve angle map from hits dp
            f->angleMap.pdf = globSimState.facetStates[i].recordedAngleMapPdf;
        }
#endif
    }
}