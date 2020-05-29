#include "Simulation.h"
#include "IntersectAABB_shared.h"
#include "Parameter.h"
#include <cstring>

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

	currentParticle.lastHitFacet = nullptr;

	hasVolatile = false;

	memset(&tmpGlobalResult, 0, sizeof(GlobalHitBuffer));

	sh.nbSuper = 0;
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
		
		acLines =
		acTLines = nullptr;
}

Simulation::~Simulation(){

}

int Simulation::ReinitializeParticleLog() {
    tmpParticleLog.clear();
    tmpParticleLog.shrink_to_fit();
    if (ontheflyParams.enableLogging)
        tmpParticleLog.reserve(ontheflyParams.logLimit / ontheflyParams.nbProcess);

    return 0;
}