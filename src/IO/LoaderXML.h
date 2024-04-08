

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
