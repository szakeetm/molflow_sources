//
// Created by pascal on 5/17/21.
//

#include "SettingsIO.h"
#include "SettingsIO_Molflow.h"
#include "CSVExporter.h"
#include <Helper/ConsoleLogger.h>
#include <Helper/StringHelper.h>
#include <filesystem>

// zip
#include <File.h>
#include <ziplib/ZipArchive.h>
#include <ziplib/ZipFile.h>

namespace SettingsIO {
    void export_facet_details(GlobalSimuState* glob, SimulationModel* model){
        //Could use SettingsIO::outputFile instead of fixed name
        //std::string csvFile = std::filesystem::path(SettingsIO::outputFile).replace_extension(".csv").string();
        std::string csvFile = "facet_details.csv";
        csvFile = std::filesystem::path(SettingsIO::workPath).append(csvFile).string();

        if (FlowIO::CSVExporter::ExportAllFacetDetails(csvFile, glob, model)) {
            Log::console_error("Could not write facet details to CSV file %s\n", csvFile.c_str());
        } else {
            Log::console_msg_master(3, "Successfully wrote facet details to CSV file %s\n", csvFile.c_str());
        }
    }

    void export_facet_quantities(GlobalSimuState* glob, SimulationModel* model){
        //Could use SettingsIO::outputFile instead of fixed name
        //std::string csvFile = std::filesystem::path(SettingsIO::outputFile).replace_extension(".csv").string();
        std::string csvFile = "facet_physics.csv";
        csvFile = std::filesystem::path(SettingsIO::workPath).append(csvFile).string();

        if (FlowIO::CSVExporter::ExportPhysicalQuantitiesForFacets(csvFile, glob, model)) {
            Log::console_error("Could not write facet quantities to CSV file %s\n", csvFile.c_str());
        } else {
            Log::console_msg_master(3, "Successfully wrote facet quantities to CSV file %s\n", csvFile.c_str());
        }
    }
} // namespace SettingsIO