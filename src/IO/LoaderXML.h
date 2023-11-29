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

#pragma once

#include <set>
#include <map>
#include "PugiXML/pugixml.hpp"

class GLProgress_Abstract;
class MolflowSimulationModel;
class MolflowSimFacet;
class Parameter;
struct TimeDependentParameters;
struct Formulas;
class GlobalSimuState;
struct FacetInterfaceSetting;
struct MolflowInterfaceSettings;
struct CameraView;

namespace FlowIO {

    class Loader {
    protected:
    public:
        virtual std::shared_ptr<MolflowSimulationModel> LoadGeometry(const std::string &inputFileName, const std::vector<Parameter>& catalog, GLProgress_Abstract& prg) = 0;
    };

    class XmlLoader : public Loader {

    protected:
        void LoadFacet(pugi::xml_node facetNode, std::shared_ptr<MolflowSimFacet> facet, FacetInterfaceSetting& fis, size_t nbTotalVertices, const TimeDependentParameters& tdParams);
    public:
        XmlLoader();
        std::shared_ptr<MolflowSimulationModel> LoadGeometry(const std::string& inputFileName, const std::vector<Parameter>& catalog, GLProgress_Abstract& prg) override;
        static int LoadSimulationState(const std::string& inputFileName, const std::shared_ptr<MolflowSimulationModel> model,
            const std::shared_ptr<GlobalSimuState> globalState, GLProgress_Abstract& prg);
        static int
            LoadConvergenceValues(const std::string& inputFileName, const std::shared_ptr<Formulas> appFormulas, GLProgress_Abstract& prg);
        std::unique_ptr<MolflowInterfaceSettings> interfaceSettings; //user settings such as selections, facet view settings, parameters and moments, that must be persistent even in CLI };
        std::unique_ptr<CameraView> XmlToCameraView(const pugi::xml_node& viewNode); //using unique pointer to avoid circular headers
    };
}
