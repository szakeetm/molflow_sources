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

Simulation::Simulation(std::string appName , std::string dpName, size_t parentPID, size_t procIdx) : SimulationController(appName, dpName, parentPID, procIdx)
{
	totalDesorbed = 0;

	loadOK = false;
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
