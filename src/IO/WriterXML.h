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

#ifndef MOLFLOW_PROJ_WRITERXML_H
#define MOLFLOW_PROJ_WRITERXML_H

#include <PugiXML/pugixml.hpp>
#include "PugiXML/pugixml.hpp"

#include <string>
#include "Simulation/MolflowSimGeom.h"

class MolflowSimFacet;
namespace FlowIO {

    class Writer {
    protected:

    public:
        //virtual void SaveGeometry(std::string outputFileName, SimulationModel *model) = 0;
    };

    class WriterXML : public Writer {
    protected:
        bool useOldXMLFormat;
        bool update;
    public:
        WriterXML(bool useOldXMLFormat = false, bool update = false);
        //void SaveGeometry(std::string outputFileName, SimulationModel *model) override;
        pugi::xml_node GetRootNode(pugi::xml_document &saveDoc);

        bool SaveXMLToFile(pugi::xml_document &saveDoc, const std::string &outputFileName);
        void SaveGeometry(pugi::xml_document &saveDoc, std::shared_ptr<MolflowSimulationModel> model,
                          const std::vector<size_t> &selection = std::vector<size_t>{});

        bool SaveSimulationState(const std::string &outputFileName, std::shared_ptr<MolflowSimulationModel> model, GlobalSimuState &globState);

        bool SaveSimulationState(pugi::xml_document &saveDoc, std::shared_ptr<MolflowSimulationModel> model, GlobalSimuState &globState);

        void
        SaveFacet(pugi::xml_node facetNode, MolflowSimFacet *facet, size_t nbTotalVertices);

        UserSettings userSettings; //user settings such as selections, facet view settings, parameters and moments, that must be persistent even in CLI
        size_t writeProgress{ 0 };

        void reportWriteStatus(const std::string &slaveStatus) const;
        void reportNewWriteStatus(const std::string &slaveStatus, double newProgress);
        void finishWriteStatus(const std::string &slaveStatus);
        void setWriteProgress(const double newProgress);
        void setWriteProgress(const size_t newProgress);
    };
}

#endif //MOLFLOW_PROJ_WRITERXML_H
