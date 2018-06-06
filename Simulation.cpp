#include "Simulation.h"
#include "IntersectAABB_shared.h"

SuperStructure::SuperStructure()
{
	aabbTree = NULL;
}

SuperStructure::~SuperStructure()
{
	SAFE_DELETE(aabbTree);
}

Simulation::Simulation()
{
	
	hitCacheSize = 0;
	nbLeakSinceUpdate = 0;
	leakCacheSize = 0;
	totalDesorbed = 0;

	loadOK = false;
	sMode = MC_MODE;
	wp.globalHistogramParams.record = false;
	lastHitFacet = NULL;

	hasVolatile = false;

	memset(&tmpGlobalCount, 0, sizeof(tmpGlobalCount));

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
		NULL;
		
		acLines =
		acTLines = NULL;
}
