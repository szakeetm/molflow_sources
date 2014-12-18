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



#define PI 3.14159265358979323846
typedef int BOOL;
#define MAX(x,y) (((x)<(y))?(y):(x))
#define MIN(x,y) (((x)<(y))?(x):(y))
#define TRUE  1
#define FALSE 0
#define MY_INFINITY 1.e100
#define SAFE_FREE(x) if(x) { free(x);x=NULL; }
#define SATURATE(x,min,max) {if(x<(min)) x=(min); if(x>(max)) x=(max);}
#define MAX_STRUCT 512
#define TIMELIMIT 0.01

#include "Shared.h"
#include "smp/SMP.h"
#include "Distributions.h"
#include <vector>
#include "Parameter.h"

#ifndef _SIMULATIONH_
#define _SIMULATIONH_

// Local facet structure

typedef struct {

  SHFACET sh;

  int      *indices;          // Indices (Reference to geometry vertex)
  VERTEX2D *vertices2;        // Vertices (2D plane space, UV coordinates)
  AHIT     **hits;            // Texture hit recording (taking area, temperature, mass into account)
  double   *inc;              // Texure increment
  BOOL     *largeEnough;      // cells that are NOT too small for autoscaling
  double   fullSizeInc;       // Texture increment of a full texture element
  VHIT     **direction;       // Direction field recording (average)
  BOOL     *fullElem;         // Direction field recording (only on full element)
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

} FACET;

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
  llong nbDesorbed;           // Total number of desorptions (for this process)             
  llong nbLeakTotal;          // Total number of unexpected leak (simulation error)
  int    nbLastLeak;          // Last leaks
  int    nbHHit;              // Last hits
  llong  maxDesorption;       // Maximum number of desorption
  HIT    pHits[NBHHIT];       // Last hit history
  LEAK   pLeak[NBHLEAK];      // Leak history
  //llong  wallHits[BOUNCEMAX]; // 'Wall collision count before absoprtion' density histogram

  std::vector<std::vector<std::pair<double,double>>> CDFs; //cumulative distribution function for each temperature
  std::vector<std::vector<std::pair<double,double>>> IDs; //integrated distribution function for each time-dependent desorption type
  std::vector<double> temperatures; //keeping track of all temperatures that have a CDF already generated
  std::vector<double> moments;      //time values (seconds) when a simulation state is measured
  std::vector<size_t> desorptionParameterIDs; //time-dependent parameters which are used as desorptions, therefore need to be integrated
  double latestMoment;

  double totalDesorbedMolecules;  // Number of desorbed molecules from t=0 to latest_moment
  double finalOutgassingRate;     // Outgassing rate at latest_moment (for const. flow calculation)
  double gasMass;
  double halfLife;
  double timeWindowSize;
  BOOL useMaxwellDistribution; //TRUE: Maxwell-Boltzmann distribution, FALSE: All molecules have the same (V_avg) speed
  BOOL calcConstantFlow;

  std::vector<Parameter> parameters; //Time-dependent parameters

  // Geometry
  char        name[64];         // Global name
  int         nbVertex;         // Number of vertex
  int         totalFacet;       // Total number of facet
  VERTEX3D   *vertices3;        // Vertices
  int         nbSuper;          // Number of super structure
  int         curStruct;        // Current structure
  size_t      nbMoments;        // Number of time moments
  SUPERSTRUCT str[MAX_STRUCT];

  FACET *lastHit;     // Last hitted facet
  double stepPerSec;  // Avg number of step per sec
  int textTotalSize;  // Texture total size
  int profTotalSize;  // Profile total size
  int dirTotalSize;   // Direction field total size
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
  VERTEX3D pPos;    // Position
  VERTEX3D pDir;    // Direction
  int      nbPHit;  // Number of hit (current particle)
  double   distTraveledCurrentParticle; //Distance traveled by particle before absorption
  double   distTraveledSinceUpdate;
  double   velocityCurrentParticle;
  double   flightTimeCurrentParticle;
  double   particleDecayMoment; //for radioactive gases
  //double   temperature;  //Temeperature of the particle (=temp. of last facet hit)

  // Angular coefficient (opaque facets)
  int     nbAC;
  ACFLOAT *acMatrix;
  ACFLOAT *acDensity;
  ACFLOAT *acDesorb;
  ACFLOAT *acAbsorb;
  ACFLOAT *acRho;
  ACFLOAT *acArea;
  double  *acLines;
  int     prgAC;

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

// -- Macros ---------------------------------------------------




// -- Methods ---------------------------------------------------

void AHIT_FACET(FACET *f, double time, BOOL countHit, double velocity_factor, double ortSpeedFactor);
void DHIT_FACET(FACET *f, double time);
void ProfileFacet(FACET *f, double time, BOOL countHit, double velocity_factor, double ortSpeedFactor);
void InitSimulation();
void ClearSimulation();
void ClearACMatrix();
BOOL LoadSimulation(Dataport *loader);
BOOL StartSimulation(int mode);
void ResetSimulation();
BOOL SimulationRun();
BOOL SimulationMCStep(int nbStep);
BOOL SimulationACStep(int nbStep);
void RecordHit(int type);
void RecordLeakPos();
BOOL StartFromSource();
void PerformBounce(FACET *iFacet);
void PerformAbsorb(FACET *iFacet);
void PerformTeleport(FACET *iFacet);
void PolarToCartesian(FACET *iFacet,double theta,double phi,BOOL reverse);
void CartesianToPolar(FACET *iFacet,double *theta,double *phi);
void UpdateHits(Dataport *dpHit,int prIdx,DWORD timeout);
void UpdateMCHits(Dataport *dpHit,int prIdx,int nbMoments,DWORD timeout);
void UpdateACHits(Dataport *dpHit,int prIdx,DWORD timeout);
void ResetCounter();
struct AABBNODE *BuildAABBTree(FACET **list,int nbFacet,int depth);
int FindBestCuttingPlane(struct AABBNODE *node,int *left,int *right);
void ComputeBB(struct AABBNODE *node);
void DestroyAABB(struct AABBNODE *node);
void IntersectTree(struct AABBNODE *node);
BOOL Intersect(VERTEX3D *rayPos,VERTEX3D *rayDir,double *dist,FACET **iFact,FACET *last);
BOOL Visible(VERTEX3D *c1,VERTEX3D *c2,FACET *f1,FACET *f2);
BOOL IsInFacet(FACET *f,double u,double v);
double GetTick();
long   GetHitsSize();
BOOL ComputeACMatrix(SHELEM *mesh);

int GetIDId(int paramId);

void UpdateVelocity(FACET *collidedFacet);
double GenerateRandomVelocity(int CDFId);
double GenerateDesorptionTime(FACET* src);
double GetStickingAt(FACET *src,double time);
double GetOpacityAt(FACET *src,double time);



#endif /* _SIMULATIONH_ */
