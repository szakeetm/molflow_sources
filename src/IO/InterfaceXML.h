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
            static void LoadInterface(pugi::xml_node interfaceNode, MolFlow *mApp);
        };

        class WriterInterfaceXML :
        public WriterXML {
        public:
            WriterInterfaceXML(bool useOldXMLFormat, bool update) : WriterXML(useOldXMLFormat, update){};
            static void WriteInterface(pugi::xml_document &saveDoc, MolFlow *mApp, bool saveSelected);
        };
}
#endif //MOLFLOW_PROJ_INTERFACEXML_H
