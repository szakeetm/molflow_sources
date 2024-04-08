

#pragma once

#include <cstdio>
#include <string>
#include <vector>
#include <memory>

class GlobalSimuState;
class MolflowSimulationModel;

namespace FlowIO {

/**
* \brief Enum relating to individual properties of a facet
*/
    enum class FDetail {
        F_ID,
        F_STICKING,
        F_OPACITY,
        F_STRUCTURE,
        F_LINK,
        F_DESORPTION,
        F_REFLECTION,
        F_TWOSIDED,
        F_VERTEX,
        F_AREA,
        F_TEMP,
        F_2DBOX,
        F_TEXTURE_UV,
        F_MESHSAMPLEPCM,
        F_COUNT,
        F_MEMORY,
        F_PLANARITY,
        F_PROFILE,
        F_IMPINGEMENT,
        F_DENSITY1P,
        F_DENSITYKGP,
        F_PRESSURE,
        F_AVGSPEED,
        F_MCHITS,
        F_EQUIVHITS,
        F_NDESORPTIONS,
        F_EQUIVABS
    };

    struct CSVExporter {
        static std::string FormatCell(FDetail mode, size_t idx, const std::shared_ptr<GlobalSimuState> globalState, const std::shared_ptr<MolflowSimulationModel> model);

        static std::string
        GetLineForFacet(size_t idx, const std::vector<FDetail> &selectedValues, const std::shared_ptr<GlobalSimuState> globalState,
            const std::shared_ptr<MolflowSimulationModel> model);

        static std::string GetFacetDetailsCSV(const std::vector<FDetail> &selectedValues, const std::shared_ptr<GlobalSimuState> globalState,
                                              const std::shared_ptr<MolflowSimulationModel> model);

        static std::string
        GetHeader(const std::vector<FDetail> &selectedValues);

        static int ExportAllFacetDetails(const std::string &fileName, const std::shared_ptr<GlobalSimuState> globalState, const std::shared_ptr<MolflowSimulationModel> model);

        static int
        ExportPhysicalQuantitiesForFacets(const std::string &fileName, const std::shared_ptr<GlobalSimuState> globalState, const std::shared_ptr<MolflowSimulationModel> model);

        static int ValidateCSVFile(const std::string &fileName);
    };

    // export utility functions
    struct Exporter {
        static void export_facet_details(const std::shared_ptr<GlobalSimuState> globalState, const std::shared_ptr<MolflowSimulationModel> mode, std::string& workPath);

        static void export_facet_quantities(const std::shared_ptr<GlobalSimuState> globalState, const std::shared_ptr<MolflowSimulationModel> mode, std::string& workPath);
    };
}
