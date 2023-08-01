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
#include <Helper/GLProgress_abstract.hpp>
#include <GLApp/GLFormula.h>

#include <string>
#include "Simulation/MolflowSimGeom.h"

class MolflowSimFacet;
namespace FlowIO {

    class Writer {
    protected:

    public:
        //virtual void SaveGeometry(std::string outputFileName, SimulationModel *model) = 0;
    };

    class XmlWriter : public Writer {
    protected:
        bool useOldXMLFormat;
        bool updateRootNode;
    public:
        XmlWriter(bool useOldXMLFormat = false, bool updateRootNode = false);
        //void SaveGeometry(std::string outputFileName, SimulationModel *model) override;
        pugi::xml_node GetRootNode(pugi::xml_document &saveDoc);

        bool WriteXMLToFile(pugi::xml_document &saveDoc, const std::string &outputFileName); //CLI uses it, prints to console on error
        void SaveGeometry(pugi::xml_document &saveDoc, std::shared_ptr<MolflowSimulationModel> model,
            GLProgress_Abstract& prg, const std::vector<size_t> &selectionToSave = std::vector<size_t>{});

        void WriteConvergenceValues(pugi::xml_document& saveDoc, const std::vector<std::vector<FormulaHistoryDatapoint>>& convergenceData, const std::vector<GLFormula>& appFormulas);

        bool AppendSimulationStateToFile(const std::string &outputFileName, std::shared_ptr<MolflowSimulationModel> model, GLProgress_Abstract& prg, GlobalSimuState &globState);
        bool SaveSimulationState(pugi::xml_document &saveDoc, std::shared_ptr<MolflowSimulationModel> model, GLProgress_Abstract& prg, GlobalSimuState &globState);

        void
        SaveFacet(pugi::xml_node facetNode, MolflowSimFacet *facet, size_t nbTotalVertices);

        MolflowUserSettings userSettings; //user settings such as selections, facet view settings, parameters and moments, that must be persistent even in CLI
    };
}

#endif //MOLFLOW_PROJ_WRITERXML_H
