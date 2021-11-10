//
// Created by Pascal Baehr on 04.11.21.
//

#ifndef MOLFLOW_PROJ_CSVEXPORTER_H
#define MOLFLOW_PROJ_CSVEXPORTER_H

#include <cstdio>
#include <string>
#include <vector>

class GlobalSimuState;
class SimulationModel;

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
        static char *FormatCell(FDetail mode, size_t idx, GlobalSimuState *glob, SimulationModel *model);

        static std::string
        GetLineForFacet(size_t idx, const std::vector<FDetail> &selectedValues, GlobalSimuState *glob,
                        SimulationModel *model);

        static std::string GetFacetDetailsCSV(const std::vector<FDetail> &selectedValues, GlobalSimuState *glob,
                                              SimulationModel *model);

        static std::string
        GetHeader(const std::vector<FDetail> &selectedValues);

        static int ExportAllFacetDetails(const std::string &fileName, GlobalSimuState *glob, SimulationModel *model);

        static int
        ExportPhysicalQuantitiesForFacets(const std::string &fileName, GlobalSimuState *glob, SimulationModel *model);
    };
}

#endif //MOLFLOW_PROJ_CSVEXPORTER_H
