//
// Created by pbahr on 22/02/2021.
//

#ifndef MOLFLOW_PROJ_INTERFACEXML_H
#define MOLFLOW_PROJ_INTERFACEXML_H

#include "LoaderXML.h"
#include "WriterXML.h"

class MolFlow;

namespace FlowIO{
        class LoaderInterfaceXML :
        public LoaderXML {
        public:
            LoaderInterfaceXML() : LoaderXML() {};
            void LoadInterface(pugi::xml_node interfaceNode, MolFlow *mApp);
        };

        class WriterInterfaceXML :
        public WriterXML {
        public:
            void WriteInterface(pugi::xml_node saveDoc, MolFlow *mApp, bool saveSelected);
        };
}
#endif //MOLFLOW_PROJ_INTERFACEXML_H
