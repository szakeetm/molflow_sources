//
// Created by Pascal Baehr on 30.07.20.
//

#ifndef MOLFLOW_PROJ_WRITERXML_H
#define MOLFLOW_PROJ_WRITERXML_H

#include <string>
#include "GeometrySimu.h"

namespace FlowIO {

    class Writer {
    protected:

    public:
        //virtual void SaveGeometry(std::string outputFileName, SimulationModel *model) = 0;
        virtual int SaveSimulationState(std::string outputFileName, SimulationModel *model, BYTE* buffer) = 0;
    };

    class WriterXML : public Writer {
        protected:
            //void SaveFacet(pugi::xml_node facetNode, SubprocessFacet *facet, size_t nbTotalVertices);
        public:
            //void SaveGeometry(std::string outputFileName, SimulationModel *model) override;
            static void SaveGeometry(const std::string& outputFileName, SimulationModel *model);
            int SaveSimulationState(std::string outputFileName, SimulationModel *model, BYTE* buffer) override;
            static int SaveSimulationState(const std::string& outputFileName, SimulationModel *model, GlobalSimuState& globState);
    };
}

#endif //MOLFLOW_PROJ_WRITERXML_H
