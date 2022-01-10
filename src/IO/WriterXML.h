//
// Created by Pascal Baehr on 30.07.20.
//

#ifndef MOLFLOW_PROJ_WRITERXML_H
#define MOLFLOW_PROJ_WRITERXML_H

#include <PugiXML/pugixml.hpp>
#include "PugiXML/pugixml.hpp"

#include <string>
#include "GeometrySimu.h"

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

        void SaveGeometry(pugi::xml_document &saveDoc, const std::shared_ptr<SimulationModel> &model,
                          const std::vector<size_t> &selection = std::vector<size_t>{});

        bool SaveSimulationState(const std::string &outputFileName, std::shared_ptr<SimulationModel> model, GlobalSimuState &globState);

        bool SaveSimulationState(pugi::xml_document &saveDoc, std::shared_ptr<SimulationModel> model, GlobalSimuState &globState);

        void
        SaveFacet(pugi::xml_node facetNode, SubprocessFacet *facet, size_t nbTotalVertices);

        UserInput uInput;
        double writeProgress{0.0};

        void reportWriteStatus(const std::string &statusString) const;
        void reportNewWriteStatus(const std::string &statusString, double newProgress);
        void finishWriteStatus(const std::string &statusString);
        void setWriteProgress(double newProgress);
    };
}

#endif //MOLFLOW_PROJ_WRITERXML_H
