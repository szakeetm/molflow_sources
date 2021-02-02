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
#include <mutex>

#include "SimulationUnit.h"
#include "Buffer_shared.h" //Facetproperties
#include "SMP.h"
#include "Vector.h"
#include "Random.h"
#include "ProcessControl.h"
#include "SimulationController.h"
#include <../src/GeometrySimu.h>

class Parameter;

// Local simulation structure

class CurrentParticleStatus {
public:
    double GenerateRandomVelocity(int CDFId, const double rndVal);
    double GenerateDesorptionTime(const SubprocessFacet *src, const double rndVal);
    void IncreaseDistanceCounters(double distanceIncrement);
    bool SimulationMCStep(size_t nbStep, size_t threadNum, size_t remainingDes);
    bool StartFromSource();
    bool UpdateMCHits(GlobalSimuState &globSimuState, size_t nbMoments, DWORD timeout);
    void RecordHitOnTexture(const SubprocessFacet *f, double time, bool countHit, double velocity_factor, double ortSpeedFactor);
    void ProfileFacet(const SubprocessFacet *f, double time, bool countHit, double velocity_factor, double ortSpeedFactor);
    void RecordHit(const int &type);
    void RecordLeakPos();
    void TreatMovingFacet();
    void IncreaseFacetCounter(const SubprocessFacet *f, double time, size_t hit, size_t desorb, size_t absorb,
                              double sum_1_per_v, double sum_v_ort);
    void UpdateVelocity(const SubprocessFacet *collidedFacet);
    void LogHit(SubprocessFacet *f, std::vector<ParticleLoggerItem> &tmpParticleLog);
    void RecordDirectionVector(const SubprocessFacet *f, double time);
    void RecordAngleMap(const SubprocessFacet *collidedFacet);
    void PerformTeleport(SubprocessFacet *iFacet);
    void RegisterTransparentPass(SubprocessFacet *facet);
    void RecordAbsorb(SubprocessFacet *iFacet);
    void PerformBounce(SubprocessFacet *iFacet);
    void RecordHistograms(SubprocessFacet *iFacet);
    bool UpdateHits(GlobalSimuState* globState, DWORD timeout);

    void Reset();
	Vector3d position;    // Position
	Vector3d direction;    // Direction
	double oriRatio; //Represented ratio of desorbed, used for low flux mode

	//Recordings for histogram
	uint64_t totalDesorbed;
	size_t   nbBounces; // Number of hit (current particle) since desorption
    size_t   lastMomentIndex; // Speedup for binary search
    size_t particleId;
	double   distanceTraveled;
    double   generationTime; //Time it was created, constant
    double   particleTime; //Actual time, incremented after every hit. (Flight time = actual time - generation time)
    int      teleportedFrom;   // We memorize where the particle came from: we can teleport back

    double   velocity;
	double   expectedDecayMoment; //for radioactive gases
	size_t   structureId;        // Current structure
	GlobalSimuState tmpState;
	SubprocessFacet *lastHitFacet;     // Last hitted facet
	MersenneTwister randomGenerator;
	SimulationModel *model;
    std::vector<SubprocessFacet*> transparentHitBuffer; //Storing this buffer simulation-wide is cheaper than recreating it at every Intersect() call
    std::vector<SubProcessFacetTempVar> tmpFacetVars; //One per subprocessfacet, for intersect routine
};

class Simulation : public SimulationUnit {
public:

    Simulation();
    Simulation(Simulation&& o) noexcept ;
    virtual ~Simulation() = default;

    int SanityCheckGeom() override;
    void ClearSimulation() override;
    size_t LoadSimulation(SimulationModel *simModel, char *loadStatus) override;
    void ResetSimulation() override;

    size_t GetHitsSize() override;

    int ReinitializeParticleLog() override;
    CurrentParticleStatus* GetParticle() override {return &currentParticle;} ;
    CurrentParticleStatus* GetParticle(size_t i) override {if(i < particles.size()) return &particles.at(i); else return nullptr;} ;
    void SetNParticle(size_t n) override {particles.clear(); particles.resize(n);};
    //bool UpdateHits(int prIdx, DWORD timeout) override;
    //bool UpdateMCHits(GlobalSimuState &globSimuState, size_t nbMoments, DWORD timeout);

    //void PerformTeleport(SubprocessFacet *iFacet, const SimulationModel &model);
    //bool SimulationMCStep(size_t nbStep, size_t threadNum) override;
    //void IncreaseDistanceCounters(double distanceIncrement, CurrentParticleStatus &currentParticle);
    //bool StartFromSource(CurrentParticleStatus &currentParticle);
    //void PerformBounce(SubprocessFacet *iFacet, const SimulationModel &model);
    //void PerformTransparentPass(SubprocessFacet *iFacet); //disabled
    //void RecordAbsorb(SubprocessFacet *iFacet, const SimulationModel &model);
    //void RecordHistograms(SubprocessFacet *iFacet, const SimulationModel &model);
    //void RecordHitOnTexture(SubprocessFacet *f, double time, bool countHit, double velocity_factor,double ortSpeedFactor);
    //void RecordDirectionVector(SubprocessFacet *f, double time, CurrentParticleStatus &currentParticle);
    //void ProfileFacet(SubprocessFacet *f, double time, bool countHit, double velocity_factor, double ortSpeedFactor,const SimulationModel &model);
    //void LogHit(SubprocessFacet *f, const SimulationModel &model, std::vector<ParticleLoggerItem> &tmpParticleLog);
    //void RecordAngleMap(SubprocessFacet *collidedFacet);
    //void UpdateVelocity(SubprocessFacet *collidedFacet, const SimulationModel &model);
    //static double GenerateRandomVelocity(int CDFId, const double rndVal, const SimulationModel &model);
    //static double GenerateDesorptionTime(const SubprocessFacet *src, const double rndVal, const SimulationModel &model);
    //static double GetStickingAt(SubprocessFacet *f, double time, const SimulationModel &model);
    //double GetOpacityAt(SubprocessFacet *f, double time);
    //void TreatMovingFacet(const SimulationModel &model);
    //void IncreaseFacetCounter(SubprocessFacet *f, double time, size_t hit, size_t desorb, size_t absorb, double sum_1_per_v, double sum_v_ort, CurrentParticleStatus &currentParticle);
    //void RegisterTransparentPass(SubprocessFacet *facet, const SimulationModel &model);


	//size_t totalDesorbed;           // Total desorption number (for this process, not reset on UpdateMCHits)

	// Geometry

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
	//std::vector<CurrentParticleStatus> currentParticles;
    //std::vector<FacetHistogramBuffer> tmpGlobalHistograms; //Recorded histogram since last UpdateMCHits, 1+nbMoment copies
    std::vector<ParticleLoggerItem> tmpParticleLog; //Recorded particle log since last UpdateMCHits
    CurrentParticleStatus currentParticle;
    std::vector<CurrentParticleStatus> particles;
    mutable std::timed_mutex tMutex;
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
