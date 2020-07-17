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
#pragma once

#include <tuple>
#include <vector>
#include <SimulationUnit.h>

#include "Buffer_shared.h" //Facetproperties
#include "SMP.h"
#include "Vector.h"
#include "Random.h"
#include "ProcessControl.h"
#include "GeometrySimu.h"
#include "SimulationController.h"

class Parameter;

// Local simulation structure

class CurrentParticleStatus {
public:
	Vector3d position;    // Position
	Vector3d direction;    // Direction
	double oriRatio; //Represented ratio of desorbed, used for low flux mode

	//Recordings for histogram
	size_t   nbBounces; // Number of hit (current particle) since desorption
    size_t   lastMomentIndex; // Speedup for binary search
	double   distanceTraveled;
	double   flightTime;

	double   velocity;
	double   expectedDecayMoment; //for radioactive gases
	size_t   structureId;        // Current structure
	int      teleportedFrom;   // We memorize where the particle came from: we can teleport back
	SubprocessFacet *lastHitFacet;     // Last hitted facet
	std::vector<SubprocessFacet*> transparentHitBuffer; //Storing this buffer simulation-wide is cheaper than recreating it at every Intersect() call
};

class Simulation : public SimulationUnit {
public:

    Simulation();
    ~Simulation();

    //int controlledLoop();
    int SanityCheckGeom() override;
    void ClearSimulation() override;
    bool LoadSimulation(Dataport *loader) override;
    void ResetSimulation() override;
    bool StartSimulation();
    bool SimulationRun();

    void RecordHit(const int &type);
    void RecordLeakPos();
    size_t GetHitsSize() override;

    int ReinitializeParticleLog() override;
    bool UpdateOntheflySimuParams(Dataport *loader) override;

    void UpdateHits(Dataport *dpHit, Dataport* dpLog,int prIdx, DWORD timeout);
    bool UpdateMCHits(Dataport *dpHit, int prIdx, size_t nbMoments, DWORD timeout);
    void UpdateLog(Dataport *dpLog, DWORD timeout);

    void PerformTeleport(SubprocessFacet *iFacet);
    bool SimulationMCStep(size_t nbStep);
    void IncreaseDistanceCounters(double distanceIncrement);
    bool StartFromSource();
    void PerformBounce(SubprocessFacet *iFacet);
    void PerformTransparentPass(SubprocessFacet *iFacet); //disabled
    void RecordAbsorb(SubprocessFacet *iFacet);
    void RecordHistograms(SubprocessFacet *iFacet);
    void RecordHitOnTexture(SubprocessFacet *f, double time, bool countHit, double velocity_factor, double ortSpeedFactor);
    void RecordDirectionVector(SubprocessFacet *f, double time);
    void ProfileFacet(SubprocessFacet *f, double time, bool countHit, double velocity_factor, double ortSpeedFactor);
    void LogHit(SubprocessFacet *f);
    void RecordAngleMap(SubprocessFacet *collidedFacet);
    void UpdateVelocity(SubprocessFacet *collidedFacet);
    double GenerateRandomVelocity(int CDFId);
    double GenerateDesorptionTime(SubprocessFacet *src);
    double GetStickingAt(SubprocessFacet *f, double time);
    double GetOpacityAt(SubprocessFacet *f, double time);
    void TreatMovingFacet();
    void IncreaseFacetCounter(SubprocessFacet *f, double time, size_t hit, size_t desorb, size_t absorb, double sum_1_per_v,
                              double sum_v_ort);
    void RegisterTransparentPass(SubprocessFacet *facet);

        void ResetTmpCounters();

    /*void GetState();
    size_t GetLocalState();
    void SetErrorSub(const char *message);
    void SetState(size_t state,const char *status,bool changeState = true, bool changeStatus = true);
    void SetStatus(char *status);
    char *GetSimuStatus();
    void SetReady();*/


    //bool endState;

    //Dataport *dpControl;
    //Dataport *dpHit;
    //Dataport *dpLog;

    MersenneTwister randomGenerator;

    GlobalHitBuffer tmpGlobalResult; //Global results since last UpdateMCHits
	std::vector<FacetHistogramBuffer> tmpGlobalHistograms; //Recorded histogram since last UpdateMCHits, 1+nbMoment copies
	std::vector<ParticleLoggerItem> tmpParticleLog; //Recorded particle log since last UpdateMCHits

	//size_t totalDesorbed;           // Total number of desorptions (for this process, not reset on UpdateMCHits)

	std::vector<std::vector<std::pair<double, double>>> CDFs; //cumulative distribution function for each temperature
	std::vector<IntegratedDesorption> IDs; //integrated distribution function for each time-dependent desorption type
	//std::vector<double> temperatures; //keeping track of all temperatures that have a CDF already generated
	std::vector<Moment> moments;      //time values (seconds) when a simulation state is measured
	//std::vector<size_t> desorptionParameterIDs; //time-dependent parameters which are used as desorptions, therefore need to be integrated
	std::vector<Parameter> parameters; //Time-dependent parameters 


	// Geometry
	//GeomProperties sh;
	WorkerParams wp;
	//OntheflySimulationParams ontheflyParams;

	std::vector<Vector3d>   vertices3;        // Vertices
	std::vector<SuperStructure> structures; //They contain the facets  

	//double stepPerSec=0.0;  // Avg number of step per sec
	//Facet size counters
	size_t textTotalSize;  // Texture total size
	size_t profTotalSize;  // Profile total size
	size_t dirTotalSize;   // Direction field total size
	size_t angleMapTotalSize;
	size_t histogramTotalSize;
	//bool loadOK;        // Load OK flag
	//bool lastHitUpdateOK;  // Last hit update timeout
	bool lastLogUpdateOK; // Last log update timeout
	bool hasVolatile;   // Contains volatile facet

	// Particle coordinates (MC)
	CurrentParticleStatus currentParticle;

};
// -- Methods ---------------------------------------------------

/*void RecordHitOnTexture(SubprocessFacet *f, double time, bool countHit, double velocity_factor, double ortSpeedFactor);
void RecordDirectionVector(SubprocessFacet *f, double time);
void ProfileFacet(SubprocessFacet *f, double time, bool countHit, double velocity_factor, double ortSpeedFactor);
void LogHit(SubprocessFacet *f);
void RecordAngleMap(SubprocessFacet* collidedFacet);
void InitSimulation();
void ClearSimulation();
//void SetState(size_t state, const char *status, bool changeState = true, bool changeStatus = true);
void SetErrorSub(const char *msg);
void ClearACMatrix();
bool LoadSimulation(Dataport *loader);
bool UpdateOntheflySimuParams(Dataport *loader);
bool StartSimulation(size_t sMode);
void ResetSimulation();
bool SimulationRun();
bool SimulationMCStep(size_t nbStep);
void IncreaseDistanceCounters(double d);
bool SimulationACStep(int nbStep);
void RecordHit(const int& type);
void RecordLeakPos();
bool StartFromSource();
void PerformBounce(SubprocessFacet *iFacet);
void RecordAbsorb(SubprocessFacet *iFacet);
void RecordHistograms(SubprocessFacet * iFacet);
void PerformTeleport(SubprocessFacet *iFacet);
void PerformTransparentPass(SubprocessFacet *iFacet);
void UpdateHits(Dataport *dpHit, Dataport *dpLog, int prIdx, DWORD timeout);
void UpdateLog(Dataport *dpLog, DWORD timeout);
void UpdateMCHits(Dataport *dpHit, int prIdx, size_t nbMoments, DWORD timeout);
void UpdateACHits(Dataport *dpHit, int prIdx, DWORD timeout);
void ResetTmpCounters();

double GetTick();
size_t   GetHitsSize();
bool ComputeACMatrix(SHELEM_OLD *mesh);

int GetIDId(int paramId);

void   UpdateVelocity(SubprocessFacet *collidedFacet);
double GenerateRandomVelocity(int CDFId);
double GenerateDesorptionTime(SubprocessFacet* src);
double GetStickingAt(SubprocessFacet *src, double time);
double GetOpacityAt(SubprocessFacet *src, double time);
void   IncreaseFacetCounter(SubprocessFacet *f, double time, size_t hit, size_t desorb, size_t absorb, double sum_1_per_v, double sum_v_ort);
void   RegisterTransparentPass(SubprocessFacet *facet);
void   TreatMovingFacet();*/
