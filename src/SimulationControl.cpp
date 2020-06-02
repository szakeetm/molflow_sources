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
//extern Simulation* sHandle; //Declared at molflowSub.cpp



/*void InitSimulation() {

	// Global handle allocation
	sHandle = new Simulation();
	InitTick();
}*/

int Simulation::SanityCheckGeom() {

    return 0; // all ok
}

void Simulation::ClearSimulation() {

    //loadOK = false;

    textTotalSize =
    profTotalSize =
    dirTotalSize =
    angleMapTotalSize =
    histogramTotalSize = 0;

    this->currentParticle = CurrentParticleStatus();

    this->structures.clear();
    this->CDFs.clear();
    this->IDs.clear();
    this->moments.clear();
    this->parameters.clear();
    this->temperatures.clear();
    this->vertices3.clear();
}

bool Simulation::LoadSimulation(Dataport *loader) {
	double t0 = GetTick();

	//SetState(PROCESS_STARTING, "Clearing previous simulation");
	ClearSimulation();

	//SetState(PROCESS_STARTING, "Loading simulation");

	{
		
		std::string inputString(loader->size,'\0');
		BYTE* buffer = (BYTE*)loader->buff;
		std::copy(buffer, buffer + loader->size, inputString.begin());
		std::stringstream inputStream;
		inputStream << inputString;
		cereal::BinaryInputArchive inputArchive(inputStream);

		//Worker params
		inputArchive(wp);
		inputArchive(ontheflyParams);
		inputArchive(CDFs);
		inputArchive(IDs);
		inputArchive(parameters);
		inputArchive(temperatures);
		inputArchive(moments);
		inputArchive(desorptionParameterIDs);

		//Geometry
		inputArchive(sh);
		inputArchive(vertices3);

		structures.resize(sh.nbSuper); //Create structures

		//Facets
		for (size_t i = 0; i < sh.nbFacet; i++) { //Necessary because facets is not (yet) a vector in the interface
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
			if (!f.InitializeOnLoad(i, moments.size(), histogramTotalSize)) return false;
			// Increase size counters
			//histogramTotalSize += 0;
			angleMapTotalSize += f.angleMapSize;
			dirTotalSize += f.directionSize* (1 + moments.size());
			profTotalSize += f.profileSize* (1 + moments.size());
			textTotalSize += f.textureSize* (1 + moments.size());

            hasVolatile |= f.sh.isVolatile;
            if ((f.sh.superDest || f.sh.isVolatile) && ((f.sh.superDest - 1) >= sh.nbSuper || f.sh.superDest < 0)) {
                // Geometry error
                //ClearSimulation();
                //ReleaseDataport(loader);
                std::ostringstream err;
                err << "Invalid structure (wrong link on F#" << i + 1 << ")";
                //SetErrorSub(err.str().c_str());
                std::cerr << err.str() << std::endl;
                return false;
            }

            if (f.sh.superIdx == -1) { //Facet in all structures
				for (auto& s : structures) {
					s.facets.push_back(f);
				}
			}
            else {
				structures[f.sh.superIdx].facets.push_back(f); //Assign to structure
			}
		}
	}//inputarchive goes out of scope, file released

	//Initialize global histogram
	FacetHistogramBuffer hist;
    hist.Resize(wp.globalHistogramParams);
	tmpGlobalHistograms = std::vector<FacetHistogramBuffer>(1 + moments.size(), hist);

	//Reserve particle log
	if (ontheflyParams.enableLogging)
	    tmpParticleLog.reserve(ontheflyParams.logLimit / ontheflyParams.nbProcess);

	// Build all AABBTrees
	size_t maxDepth=0;
	for (auto& s : structures) {
		std::vector<SubprocessFacet*> facetPointers; facetPointers.reserve(s.facets.size());
		for (auto& f : s.facets) {
			facetPointers.push_back(&f);
		}
		s.aabbTree = BuildAABBTree(facetPointers, 0, maxDepth);
	}

	// Initialise simulation

	//if(!sh.name.empty())
	    //loadOK = true;
	double t1 = GetTick();
	printf("  Load %s successful\n", sh.name.c_str());
	printf("  Geometry: %zd vertex %zd facets\n", vertices3.size(), sh.nbFacet);

	printf("  Geom size: %d bytes\n", /*(size_t)(buffer - bufferStart)*/0);
	printf("  Number of stucture: %zd\n", sh.nbSuper);
	printf("  Global Hit: %zd bytes\n", sizeof(GlobalHitBuffer));
	printf("  Facet Hit : %zd bytes\n", sh.nbFacet * sizeof(FacetHitBuffer));
	printf("  Texture   : %zd bytes\n", textTotalSize);
	printf("  Profile   : %zd bytes\n", profTotalSize);
	printf("  Direction : %zd bytes\n", dirTotalSize);

	printf("  Total     : %zd bytes\n", GetHitsSize());
	printf("  Seed: %lu\n", randomGenerator.GetSeed());
	printf("  Loading time: %.3f ms\n", (t1 - t0)*1000.0);
	return true;

}

void Simulation::UpdateHits(Dataport *dpHit, Dataport* dpLog,int prIdx, DWORD timeout) {
    UpdateMCHits(dpHit, prIdx, moments.size(), timeout);
    if (dpLog) UpdateLog(dpLog, timeout);
}

size_t Simulation::GetHitsSize() {
	return sizeof(GlobalHitBuffer) + wp.globalHistogramParams.GetDataSize() +
		textTotalSize + profTotalSize + dirTotalSize + angleMapTotalSize + histogramTotalSize
		+ sh.nbFacet * sizeof(FacetHitBuffer) * (1+moments.size());
}

void Simulation::ResetTmpCounters() {
	//SetState(0, "Resetting local cache...", false, true);

	memset(&tmpGlobalResult, 0, sizeof(GlobalHitBuffer));
	
	//Reset global histograms
	for (auto& h : tmpGlobalHistograms) {
		h.Reset();
	}

	for (int j = 0; j < sh.nbSuper; j++) {
		for (auto& f : structures[j].facets) {
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

void Simulation::ResetSimulation() {
	currentParticle.lastHitFacet = nullptr;
	totalDesorbed = 0;
	ResetTmpCounters();
	tmpParticleLog.clear();
}

bool Simulation::StartSimulation() {
    if (!currentParticle.lastHitFacet) StartFromSource();
    return (currentParticle.lastHitFacet != nullptr);
}

void Simulation::RecordHit(const int &type) {
	if (tmpGlobalResult.hitCacheSize < HITCACHESIZE) {
		tmpGlobalResult.hitCache[tmpGlobalResult.hitCacheSize].pos = currentParticle.position;
		tmpGlobalResult.hitCache[tmpGlobalResult.hitCacheSize].type = type;
		tmpGlobalResult.hitCacheSize++;
	}
}

void Simulation::RecordLeakPos() {
	// Source region check performed when calling this routine 
	// Record leak for debugging
	RecordHit(HIT_REF);
	RecordHit(HIT_LAST);
	if (tmpGlobalResult.leakCacheSize < LEAKCACHESIZE) {
		tmpGlobalResult.leakCache[tmpGlobalResult.leakCacheSize].pos = currentParticle.position;
		tmpGlobalResult.leakCache[tmpGlobalResult.leakCacheSize].dir = currentParticle.direction;
		tmpGlobalResult.leakCacheSize++;
	}
}