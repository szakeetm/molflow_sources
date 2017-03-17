 /*
  File:        Simulation.h
  Description: Monte-Carlo Simulation for UHV
  Program:     MolFlow
  Author:      R. KERSEVAN / J-L PONS / M ADY
  Copyright:   E.S.R.F / CERN

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#define MAX_STRUCT 512

#include "Shared.h"
#include "SMP.h"
#include <vector>
#include "Vector.h"
#include "Parameter.h"

#ifndef _SIMULATIONH_
#define _SIMULATIONH_

// Local facet structure

class FACET {
public:
  SHFACET sh;

  int      *indices;          // Indices (Reference to geometry vertex)
  Vector2d *vertices2;        // Vertices (2D plane space, UV coordinates)
  AHIT     **hits;            // Texture hit recording (taking area, temperature, mass into account)
  double   *inc;              // Texure increment
  BOOL     *largeEnough;      // cells that are NOT too small for autoscaling
  double   fullSizeInc;       // Texture increment of a full texture element
  VHIT     **direction;       // Direction field recording (average)
  //BOOL     *fullElem;         // Direction field recording (only on full element)
  APROFILE **profile;         // Distribution and hit recording

  // Temporary var (used in Intersect for collision)
  double colDist;
  double colU;
  double colV;
  double rw;
  double rh;
  double iw;
  double ih;

  // Temporary var (used in FillHit for hit recording)
  BOOL   hitted;
  BOOL   ready;         // Volatile state
  int    textureSize;   // Texture size (in bytes)
  int    profileSize;   // profile size (in bytes)
  int    directionSize; // direction field size (in bytes)

  /*int CDFid; //Which probability distribution it belongs to (one CDF per temperature)
  int IDid;  //If time-dependent desorption, which is its ID*/
  int globalId; //Global index (to identify when superstructures are present)
  double *outgassingMap; //outgassing map when desorption is based on imported file

  // Facet hit counters
  std::vector<SHHITS> counter;
  void  ResetCounter();
  void	ResizeCounter(size_t nbMoments);
};

// Temporary transparent hit
#define MAX_THIT    16384
extern  FACET     **THits;

// AABBTree node

struct AABBNODE {

  AABB             bb;
  struct AABBNODE *left;
  struct AABBNODE *right;
  FACET          **list;
  int              nbFacet;

};

// Local simulation structure

typedef struct {

  int              nbFacet;  // Number of facet
  FACET          **facets;   // Facet handles
  struct AABBNODE *aabbTree; // Structure AABB tree

} SUPERSTRUCT;

typedef struct {

  SHHITS tmpCount;            // Temporary number of hits (between 2 updates)
  llong  desorptionLimit;       // Maximum number of desorption

  int    hitCacheSize;              // Last hits  
  HIT    hitCache[HITCACHESIZE];       // Last hit history

  size_t    nbLeakSinceUpdate;   // Leaks since last UpdateMC
  size_t	leakCacheSize;		// Leaks from regions with displayed photons since last UpdateMC (registered in cache)
  LEAK		leakCache[LEAKCACHESIZE];      // Leak cache since last UpdateMC

  llong totalDesorbed;           // Total number of desorptions (for this process)


  std::vector<std::vector<std::pair<double,double>>> CDFs; //cumulative distribution function for each temperature
  std::vector<std::vector<std::pair<double,double>>> IDs; //integrated distribution function for each time-dependent desorption type
  std::vector<double> temperatures; //keeping track of all temperatures that have a CDF already generated
  std::vector<double> moments;      //time values (seconds) when a simulation state is measured
  std::vector<size_t> desorptionParameterIDs; //time-dependent parameters which are used as desorptions, therefore need to be integrated
  double latestMoment;

  double totalDesorbedMolecules;  // Number of desorbed molecules from t=0 to latest_moment
  double finalOutgassingRate;     // Outgassing rate at latest_moment (for const. flow calculation)
  double gasMass;
  BOOL   enableDecay;
  double halfLife;
  double timeWindowSize;
  BOOL useMaxwellDistribution; //TRUE: Maxwell-Boltzmann distribution, FALSE: All molecules have the same (V_avg) speed
  BOOL calcConstantFlow;

  std::vector<Parameter> parameters; //Time-dependent parameters

  // Geometry
  char        name[64];         // Global name
  size_t         nbVertex;         // Number of vertex
  size_t         totalFacet;       // Total number of facet
  Vector3d   *vertices3;        // Vertices
  int         nbSuper;          // Number of super structure
  int         curStruct;        // Current structure
  int         teleportedFrom;   // We memorize where the particle came from: we can teleport back
  size_t      nbMoments;        // Number of time moments
  SUPERSTRUCT str[MAX_STRUCT];

  FACET *lastHit;     // Last hitted facet
  double stepPerSec;  // Avg number of step per sec
  size_t textTotalSize;  // Texture total size
  size_t profTotalSize;  // Profile total size
  size_t dirTotalSize;   // Direction field total size
  BOOL loadOK;        // Load OK flag
  BOOL lastUpdateOK;  // Last hit update timeout
  BOOL hasVolatile;   // Contains volatile facet
  BOOL hasDirection;  // Contains direction field
  int  sMode;         // Simulation mode (MC_MODE or AC_MODE)
  double calcACTime;  // AC matrix calculation time

  //double totalOutgassing;
  
  //int    nonIsothermal;
  //HANDLE molflowHandle;

  // Particle coordinates (MC)
  Vector3d pPos;    // Position
  Vector3d pDir;    // Direction
  //int      nbPHit;  // Number of hit (current particle) //Uncommented as it had no function
  //double   distTraveledCurrentParticle; //Distance traveled by particle before absorption
  double   distTraveledSinceUpdate_total; //includes "half" hits, i.e. when particle decays mid-air
  double   distTraveledSinceUpdate_fullHitsOnly; //partial distances not included (for MFP calculation)
  double   velocityCurrentParticle;
  double   flightTimeCurrentParticle;
  double   particleDecayMoment; //for radioactive gases
  //double   temperature;  //Temeperature of the particle (=temp. of last facet hit)

  int motionType;
  Vector3d motionVector1; //base point for rotation
  Vector3d motionVector2; //rotation vector or velocity vector

  // Angular coefficient (opaque facets)
  int     nbAC;
  ACFLOAT *acMatrix;
  ACFLOAT *acDensity;
  ACFLOAT *acDesorb;
  ACFLOAT *acAbsorb;
  ACFLOAT *acRho;
  ACFLOAT *acArea;
  double  *acLines;
  size_t     prgAC;

  // Angular coefficient (transparent facets)
  int     nbACT; 
  ACFLOAT *acTMatrix;
  ACFLOAT *acTDensity;
  ACFLOAT *acTArea;
  double  *acTLines;

  /*
  //Test cube
  double testCubeTime;
  double testSystemDist;
  double testCubeTemp;
  double testCubeDist;
  double testCubeEnterMoment;
  double testCubeEnterDist;
  double testCubeVelocity;
  double testSystemTime;
  llong  testCubeCount;*/
  
#ifdef JACOBI_ITERATION
  ACFLOAT *acDensityTmp;
#endif

} SIMULATION;

// Handle to simulation object
extern SIMULATION *sHandle;

//Just for AC matrix calculation in Molflow, old mesh structure:
typedef struct {

	float   area;     // Area of element
	float   uCenter;  // Center coordinates
	float   vCenter;  // Center coordinates
	int     elemId;   // Element index (MESH array)
	BOOL    full;     // Element is full

} SHELEM;

// -- Macros ---------------------------------------------------




// -- Methods ---------------------------------------------------

void RecordHitOnTexture(FACET *f, double time, BOOL countHit, double velocity_factor, double ortSpeedFactor);
void RecordDirectionVector(FACET *f, double time);
void ProfileFacet(FACET *f, double time, BOOL countHit, double velocity_factor, double ortSpeedFactor);
void InitSimulation();
void ClearSimulation();
void SetState(int state, const char *status, BOOL changeState=TRUE, BOOL changeStatus=TRUE);
void SetErrorSub(const char *msg);
void ClearACMatrix();
BOOL LoadSimulation(Dataport *loader);
BOOL StartSimulation(size_t mode);
void ResetSimulation();
BOOL SimulationRun();
BOOL SimulationMCStep(int nbStep);
BOOL SimulationACStep(int nbStep);
void RecordHit(const int& type);
void RecordLeakPos();
BOOL StartFromSource();
void PerformBounce(FACET *iFacet);
void PerformAbsorb(FACET *iFacet);
void PerformTeleport(FACET *iFacet);
void PerformTransparentPass(FACET *iFacet);
void PolarToCartesian(FACET *iFacet,double theta,double phi,BOOL reverse);
void CartesianToPolar(FACET *iFacet,double *theta,double *phi);
void UpdateHits(Dataport *dpHit,int prIdx,DWORD timeout);
void UpdateMCHits(Dataport *dpHit,int prIdx,size_t nbMoments,DWORD timeout);
void UpdateACHits(Dataport *dpHit,int prIdx,DWORD timeout);
void ResetTmpCounters();
struct AABBNODE *BuildAABBTree(FACET **list,int nbFacet,int depth);
int FindBestCuttingPlane(struct AABBNODE *node,int *left,int *right);
void ComputeBB(struct AABBNODE *node);
void DestroyAABB(struct AABBNODE *node);
void IntersectTree(struct AABBNODE *node);
BOOL Intersect(Vector3d *rayPos,Vector3d *rayDir,double *dist,FACET **iFact,FACET *last);
BOOL Visible(Vector3d *c1,Vector3d *c2,FACET *f1,FACET *f2);
BOOL IsInFacet(FACET *f,const double &u,const double &v);
double GetTick();
size_t   GetHitsSize();
BOOL ComputeACMatrix(SHELEM *mesh);

int GetIDId(int paramId);

void   UpdateVelocity(FACET *collidedFacet);
double GenerateRandomVelocity(int CDFId);
double GenerateDesorptionTime(FACET* src);
double GetStickingAt(FACET *src,double time);
double GetOpacityAt(FACET *src,double time);
void   IncreaseFacetCounter(FACET *f, double time, size_t hit,size_t desorb, size_t absorb, double sum_1_per_v, double sum_v_ort);
void   TreatMovingFacet();


//Static variables shared between IntersectAABB and IntersectAABB_shared routines, declared here
// Minimum number of facet inside a BB
#define MINBB    8

// Maximum AABB tree depth
#define MAXDEPTH 5

#endif /* _SIMULATIONH_ */
