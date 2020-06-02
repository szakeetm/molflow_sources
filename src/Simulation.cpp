#include "Simulation.h"
#include "IntersectAABB_shared.h"
#include "Parameter.h"
#include <cstring>
#include <sstream>
#include <cereal/archives/binary.hpp>

/*SuperStructure::SuperStructure()
{
	aabbTree = NULL;
}

SuperStructure::~SuperStructure()
{
	SAFE_DELETE(aabbTree);
}*/

Simulation::Simulation()
{
	totalDesorbed = 0;

    textTotalSize =
    profTotalSize =
    dirTotalSize =
    angleMapTotalSize =
    histogramTotalSize = 0;

    lastLogUpdateOK = true;

    currentParticle = CurrentParticleStatus();
	currentParticle.lastHitFacet = nullptr;

	hasVolatile = false;

	memset(&tmpGlobalResult, 0, sizeof(GlobalHitBuffer));

	sh.nbSuper = 0;

	// AC Mode
	acDensity =
		acMatrix =
		acDensity =
		acDesorb =
		acAbsorb =
		acTArea =
		acRho = 
		acTMatrix =
		acTDensity = acArea =
		nullptr;

	calcACTime = 0.0;
    nbAC = 0;
    prgAC = 0;
    nbACT = 0;

    acLines =
		acTLines = nullptr;
}

Simulation::~Simulation()= default;

int Simulation::ReinitializeParticleLog() {
    tmpParticleLog.clear();
    tmpParticleLog.shrink_to_fit();
    if (ontheflyParams.enableLogging)
        tmpParticleLog.reserve(ontheflyParams.logLimit / ontheflyParams.nbProcess);

    return 0;
}

bool Simulation::UpdateOntheflySimuParams(Dataport *loader) {
    // Connect the dataport


    if (!AccessDataportTimed(loader, 2000)) {
        //SetErrorSub("Failed to connect to loader DP");
        std::cerr << "Failed to connect to loader DP" << std::endl;
        return false;
    }
    std::string inputString(loader->size,'\0');
    BYTE* buffer = (BYTE*)loader->buff;
    std::copy(buffer, buffer + loader->size, inputString.begin());
    std::stringstream inputStream;
    inputStream << inputString;
    cereal::BinaryInputArchive inputArchive(inputStream);

    inputArchive(ontheflyParams);

    ReleaseDataport(loader);

    return true;
}