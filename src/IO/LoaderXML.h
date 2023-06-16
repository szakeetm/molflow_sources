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

#ifndef MOLFLOW_PROJ_LOADERXML_H
#define MOLFLOW_PROJ_LOADERXML_H

#include <set>
#include <map>
#include <GeometryTypes.h>
#include <Formulas.h>
#include "Simulation/MolflowSimGeom.h"
#include "PugiXML/pugixml.hpp"
#include "Simulation/MolflowSimFacet.h"
#include <Helper/GLProgress_abstract.hpp>

namespace FlowIO {
    class Loader {
    protected:
        std::vector<std::vector<IntegratedDesorptionEntry>> IDs;         //integrated distribution function for each time-dependent desorption type
        //std::vector<std::vector<std::pair<double, double>>> CDFs;        //cumulative distribution function for each temperature
    public:
        virtual int LoadGeometry(const std::string &inputFileName, std::shared_ptr<MolflowSimulationModel> model, GLProgress_Abstract& prg) = 0;
    };

    class LoaderXML : public Loader {

    protected:
        void LoadFacet(pugi::xml_node facetNode, MolflowSimFacet *facet, size_t nbTotalVertices);
    public:
        int LoadGeometry(const std::string &inputFileName, std::shared_ptr<MolflowSimulationModel> model, GLProgress_Abstract& prg) override;
        static std::vector<SelectionGroup> LoadSelections(const std::string& inputFileName);
        static int LoadSimulationState(const std::string &inputFileName, std::shared_ptr<MolflowSimulationModel> model,
                                       GlobalSimuState *globState, GLProgress_Abstract& prg);
        static int
        LoadConvergenceValues(const std::string &inputFileName, std::vector<ConvergenceData> *convergenceValues, GLProgress_Abstract& prg);
        GeometrySettings geometrySettings;
    };
}

#endif //MOLFLOW_PROJ_LOADERXML_H
