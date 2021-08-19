//
// Created by Pascal Baehr on 30.07.20.
//

#ifndef MOLFLOW_PROJ_WRITERXML_H
#define MOLFLOW_PROJ_WRITERXML_H

#include <string>
#include <PugiXML/pugixml.hpp>
#include "PugiXML/pugixml.hpp"

#include "GeometrySimu.h"
#include "LoaderXML.h"
#include "WriterXML.h"

namespace FlowIO {

    class Writer {
    protected:

    public:
        //virtual void SaveGeometry(std::string outputFileName, SimulationModel *model) = 0;
    };

    class WriterXML : public Writer {
    protected:
    public:
        //void SaveGeometry(std::string outputFileName, SimulationModel *model) override;
        void SaveGeometry(pugi::xml_document &saveDoc, std::shared_ptr<SimulationModel> model, bool useOldXMLFormat, bool update);

        bool SaveSimulationState(const std::string &outputFileName, std::shared_ptr<SimulationModel> model, GlobalSimuState &globState);

        bool SaveSimulationState(pugi::xml_document &saveDoc, std::shared_ptr<SimulationModel> model, GlobalSimuState &globState);

        void
        SaveFacet(pugi::xml_node facetNode, SubprocessFacet *facet, size_t nbTotalVertices, GeomPrimitive *prim);

        UserInput uInput;
        double writeProgress{0.0};

        void reportWriteStatus(const std::string &statusString) const;
        void reportNewWriteStatus(const std::string &statusString, double newProgress);
        void finishWriteStatus(const std::string &statusString);
        void setWriteProgress(double newProgress);
    };
}

#endif //MOLFLOW_PROJ_WRITERXML_H
