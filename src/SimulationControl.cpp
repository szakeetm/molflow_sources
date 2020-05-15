/*
Program:     MolFlow+ / Synrad+
Description: Monte Carlo simulator for ultra-high vacuum and synchrotron radiation
Authors:     Jean-Luc PONS / Roberto KERSEVAN / Marton ADY
Copyright:   E.S.R.F / CERN
Website:     https://cern.ch/molflow

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

Full license text: https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html
*/
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
#define NOMINMAX
//#include <windows.h> // For GetTickCount()
#include <process.h> // For _getpid()
#else
//#include <time.h>
#include <sys/time.h>
#include <cstring>
#endif

#include <cmath>
#include <cstdio>
#include <sstream>
#include <cereal/types/utility.hpp>
#include <cereal/archives/binary.hpp>

#include "Simulation.h"
#include "IntersectAABB_shared.h"
#include "Parameter.h"




// Global handles
extern Simulation* sHandle; //Declared at molflowSub.cpp



void InitSimulation() {

	// Global handle allocation
	sHandle = new Simulation();
	InitTick();
}

void ClearSimulation() {

	delete sHandle;
	sHandle = new Simulation;

}

bool LoadSimulation(Dataport *loader) {
	double t1;
	double t0 = GetTick();

	sHandle->loadOK = false;

	SetState(PROCESS_STARTING, "Clearing previous simulation");
	ClearSimulation();

	SetState(PROCESS_STARTING, "Loading simulation");

	sHandle->textTotalSize =
		sHandle->profTotalSize =
		sHandle->dirTotalSize =
		sHandle->angleMapTotalSize =
		sHandle->histogramTotalSize = 0;

	{
		
		std::string inputString(loader->size,'\0');
		BYTE* buffer = (BYTE*)loader->buff;
		std::copy(buffer, buffer + loader->size, inputString.begin());
		std::stringstream inputStream;
		inputStream << inputString;
		cereal::BinaryInputArchive inputArchive(inputStream);

		//Worker params
		inputArchive(sHandle->wp);
		inputArchive(sHandle->ontheflyParams);
		inputArchive(sHandle->CDFs);
		inputArchive(sHandle->IDs);
		inputArchive(sHandle->parameters);
		inputArchive(sHandle->temperatures);
		inputArchive(sHandle->moments);
		inputArchive(sHandle->desorptionParameterIDs);

		//Geometry
		inputArchive(sHandle->sh);
		inputArchive(sHandle->vertices3);

		sHandle->structures.resize(sHandle->sh.nbSuper); //Create structures

		//Facets
		for (size_t i = 0; i < sHandle->sh.nbFacet; i++) { //Necessary because facets is not (yet) a vector in the interface
			SubprocessFacet f;
			inputArchive(
				f.sh,
				f.indices,
				f.vertices2,
				f.outgassingMap,
				f.angleMap.pdf,
				f.textureCellIncrements
			);

			//Some initialization
			if (!f.InitializeOnLoad(i, sHandle->moments.size(), sHandle->histogramTotalSize)) return false;
			// Increase size counters
			//sHandle->histogramTotalSize += 0;
			sHandle->angleMapTotalSize += f.angleMapSize;
			sHandle->dirTotalSize += f.directionSize* (1 + sHandle->moments.size());
			sHandle->profTotalSize += f.profileSize* (1 + sHandle->moments.size());
			sHandle->textTotalSize += f.textureSize* (1 + sHandle->moments.size());

            sHandle->hasVolatile |= f.sh.isVolatile;
            if ((f.sh.superDest || f.sh.isVolatile) && ((f.sh.superDest - 1) >= sHandle->sh.nbSuper || f.sh.superDest < 0)) {
                // Geometry error
                ClearSimulation();
                //ReleaseDataport(loader);
                std::ostringstream err;
                err << "Invalid structure (wrong link on F#" << i + 1 << ")";
                SetErrorSub(err.str().c_str());
                return false;
            }

            if (f.sh.superIdx == -1) { //Facet in all structures
				for (auto& s : sHandle->structures) {
					s.facets.push_back(f);
				}
			}
			else {
				sHandle->structures[f.sh.superIdx].facets.push_back(f); //Assign to structure
			}
		}
	}//inputarchive goes out of scope, file released

	//Initialize global histogram
	FacetHistogramBuffer hist;
    hist.Resize(sHandle->wp.globalHistogramParams);
	sHandle->tmpGlobalHistograms = std::vector<FacetHistogramBuffer>(1 + sHandle->moments.size(), hist);

	//Reserve particle log
	if (sHandle->ontheflyParams.enableLogging)
	    sHandle->tmpParticleLog.reserve(sHandle->ontheflyParams.logLimit / sHandle->ontheflyParams.nbProcess);

	// Build all AABBTrees
	size_t maxDepth=0;
	for (auto& s : sHandle->structures) {
		std::vector<SubprocessFacet*> facetPointers; facetPointers.reserve(s.facets.size());
		for (auto& f : s.facets) {
			facetPointers.push_back(&f);
		}
		s.aabbTree = BuildAABBTree(facetPointers, 0, maxDepth);
	}

	// Initialise simulation

	sHandle->loadOK = true;
	t1 = GetTick();
	printf("  Load %s successful\n", sHandle->sh.name.c_str());
	printf("  Geometry: %zd vertex %zd facets\n", sHandle->vertices3.size(), sHandle->sh.nbFacet);

	printf("  Geom size: %d bytes\n", /*(size_t)(buffer - bufferStart)*/0);
	printf("  Number of stucture: %zd\n", sHandle->sh.nbSuper);
	printf("  Global Hit: %zd bytes\n", sizeof(GlobalHitBuffer));
	printf("  Facet Hit : %zd bytes\n", sHandle->sh.nbFacet * sizeof(FacetHitBuffer));
	printf("  Texture   : %zd bytes\n", sHandle->textTotalSize);
	printf("  Profile   : %zd bytes\n", sHandle->profTotalSize);
	printf("  Direction : %zd bytes\n", sHandle->dirTotalSize);

	printf("  Total     : %zd bytes\n", GetHitsSize());
	printf("  Seed: %lu\n", sHandle->randomGenerator.GetSeed());
	printf("  Loading time: %.3f ms\n", (t1 - t0)*1000.0);
	return true;

}

bool UpdateOntheflySimuParams(Dataport *loader) {
	// Connect the dataport
	

	if (!AccessDataportTimed(loader, 2000)) {
		SetErrorSub("Failed to connect to loader DP");
		return false;
	}
    std::string inputString(loader->size,'\0');
    BYTE* buffer = (BYTE*)loader->buff;
    std::copy(buffer, buffer + loader->size, inputString.begin());
    std::stringstream inputStream;
    inputStream << inputString;
    cereal::BinaryInputArchive inputArchive(inputStream);

    inputArchive(sHandle->ontheflyParams);

	ReleaseDataport(loader);

	return true;
}

void UpdateHits(Dataport *dpHit, Dataport* dpLog,int prIdx, DWORD timeout) {
	switch (sHandle->wp.sMode) {
        case MC_MODE: {
            UpdateMCHits(dpHit, prIdx, sHandle->moments.size(), timeout);
            if (dpLog) UpdateLog(dpLog, timeout);
            break;
        }
        case AC_MODE:
        {
            UpdateACHits(dpHit, prIdx, timeout);
            break;
        }
    }

}

size_t GetHitsSize() {
	return sizeof(GlobalHitBuffer) + sHandle->wp.globalHistogramParams.GetDataSize() +
		sHandle->textTotalSize + sHandle->profTotalSize + sHandle->dirTotalSize + sHandle->angleMapTotalSize + sHandle->histogramTotalSize
		+ sHandle->sh.nbFacet * sizeof(FacetHitBuffer) * (1+sHandle->moments.size());
}

void ResetTmpCounters() {
	SetState(0, "Resetting local cache...", false, true);

	memset(&sHandle->tmpGlobalResult, 0, sizeof(GlobalHitBuffer));
	
	//Reset global histograms
	for (auto& h : sHandle->tmpGlobalHistograms) {
		h.Reset();
	}

	for (int j = 0; j < sHandle->sh.nbSuper; j++) {
		for (auto& f : sHandle->structures[j].facets) {
			f.ResetCounter();
			f.isHit = false;

			//Reset facet histograms
			
				for (auto& t : f.tmpHistograms) {
					t.Reset();
				}
			/*std::vector<TextureCell>(f.texture.size()).swap(f.texture);
			std::vector<ProfileSlice>(f.profile.size()).swap(f.profile);
			std::vector<DirectionCell>(f.direction.size()).swap(f.direction);*/

			for (auto& t : f.texture) {
				std::fill(t.begin(), t.end(), TextureCell());
			}


			for (auto& p : f.profile) {
				std::fill(p.begin(), p.end(), ProfileSlice());
			}

			
			for (auto& d : f.direction) {
				std::fill(d.begin(), d.end(), DirectionCell());
			}

			if (f.sh.anglemapParams.record) {
				ZEROVECTOR(f.angleMap.pdf);
			}
		}
	}

}

void ResetSimulation() {
	sHandle->currentParticle.lastHitFacet = NULL;
	sHandle->totalDesorbed = 0;
	ResetTmpCounters();
	sHandle->tmpParticleLog.clear();
	if (sHandle->acDensity) memset(sHandle->acDensity, 0, sHandle->nbAC * sizeof(ACFLOAT));

}

bool StartSimulation(size_t sMode) {
	sHandle->wp.sMode = sMode;
	switch (sMode) {
	case MC_MODE:
		if (!sHandle->currentParticle.lastHitFacet) StartFromSource();
		return (sHandle->currentParticle.lastHitFacet != NULL);
	case AC_MODE:
		if (sHandle->prgAC != 100) {
			SetErrorSub("AC matrix not calculated");
			return false;
		}
		else {
			sHandle->stepPerSec = 0.0;
			return true;
		}
	}

	SetErrorSub("Unknown simulation mode");
	return false;
}

void RecordHit(const int &type) {
	if (sHandle->tmpGlobalResult.hitCacheSize < HITCACHESIZE) {
		sHandle->tmpGlobalResult.hitCache[sHandle->tmpGlobalResult.hitCacheSize].pos = sHandle->currentParticle.position;
		sHandle->tmpGlobalResult.hitCache[sHandle->tmpGlobalResult.hitCacheSize].type = type;
		sHandle->tmpGlobalResult.hitCacheSize++;
	}
}

void RecordLeakPos() {
	// Source region check performed when calling this routine 
	// Record leak for debugging
	RecordHit(HIT_REF);
	RecordHit(HIT_LAST);
	if (sHandle->tmpGlobalResult.leakCacheSize < LEAKCACHESIZE) {
		sHandle->tmpGlobalResult.leakCache[sHandle->tmpGlobalResult.leakCacheSize].pos = sHandle->currentParticle.position;
		sHandle->tmpGlobalResult.leakCache[sHandle->tmpGlobalResult.leakCacheSize].dir = sHandle->currentParticle.direction;
		sHandle->tmpGlobalResult.leakCacheSize++;
	}
}

bool SimulationRun() {

	// 1s step
	double t0, t1;
	size_t    nbStep = 1;
	bool   goOn;

	if (sHandle->stepPerSec <= 0.0) {
		switch (sHandle->wp.sMode) {
		case MC_MODE:
			nbStep = 250;
			break;
		case AC_MODE:
			nbStep = 1;
			break;
		}

	}
	else {
        nbStep = std::ceil(sHandle->stepPerSec + 0.5);
    }

	t0 = GetTick();
	switch (sHandle->wp.sMode) {
	case MC_MODE:

		goOn = SimulationMCStep(nbStep);
		break;
	case AC_MODE:
		goOn = SimulationACStep(nbStep);
		break;
	}

	t1 = GetTick();
	if(goOn) // don't update on end, this will give a false ratio (SimMCStep could return actual steps instead of plain "false"
	    sHandle->stepPerSec = (1.0 * nbStep) / (t1 - t0); // every 1.0 second

#if defined(_DEBUG)
	printf("Running: stepPerSec = %lf\n", sHandle->stepPerSec);
#endif

	return !goOn;

}