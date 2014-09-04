/*
  File:        Shared.h
  Description: Shared memory structure (inter process communication)
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

#include "Types.h"
#include <Windows.h>
#include <vector>

#ifndef SHAREDH

#define SHAREDH

#define NBHLEAK     2048  // Leak history max length
#define NBHHIT      2048  // Max. displayed number of lines and Porto (OPO)hits.
#define BOUNCEMAX   8192  // 'Wall collision count before absoprtion' histogram
#define MAX_PROCESS 16    // Maximum number of process

#define MC_MODE 0         // Monte Carlo simulation mode
#define AC_MODE 1         // Angular coefficient simulation mode

#define TESTCUBE_STRUCTURE_NO 2 //Special superstructure number for test cube

// -----------------------------------------------------------------
// Hit exchange shared structure       (name: MFLWHITS[masterPID])
//
// SHGHITS
// SHHITSF1
// if(profileF1)   PROFILEF1
// if(textureF1)   TEXTUREF1
// if(directionF1) DIRECTIONF1
// SHHITSF2
// if(profileF2)   PROFILEF2
// if(textureF2)   TEXTUREF2
// if(directionF2) DIRECTIONF2
// ...
// -----------------------------------------------------------------

typedef union {
  
  struct {
    // Counts
    llong nbDesorbed;          // Number of desorbed molec
    llong nbHit;               // Number of hits
	llong nbAbsorbed;          // Number of absorbed molec
	double sum_1_per_speed;    // sum of reciprocials of velocities, used to determine the average speed of molecules in the gas (and not of the hits)
	double sum_v_ort;          // sum of orthogonal speeds of incident velocities, used to determine the pressure
  } hit;

  struct {
    // density
    double desorbed;
    double value;
    double absorbed;
  } density;

} SHHITS;

typedef struct {

  SHHITS total;               // Global counts
  int    mode;                // Simu mode (MC_MODE or AC_MODE)
  llong  nbLeakTotal;         // Total leaks
  int    nbLastLeaks;         // Last leaks
  int    nbHHit;              // Last hits
  HIT    pHits[NBHHIT];       // Hit cache
  LEAK   pLeak[NBHLEAK];      // Leak cache
  TEXTURE_MIN_MAX texture_limits[3]; //Min-max on texture
 /* AHIT   minHit;              // Minimum on texture
  AHIT   maxHit;              // Maximum on texture
  AHIT   minHitMomentsOnly;   // Minimum, not counting constant flow
  AHIT   maxHitMomentsOnly;   // Maximum, not counting constant flow*/
  //llong  wallHits[BOUNCEMAX]; // 'Wall collision count before absoprtion' density histogram
  double distTraveledTotal;
  
} SHGHITS;



// -----------------------------------------------------------------
// Master control shared memory block  (name: MFLWCTRL[masterPID])
// 
// -----------------------------------------------------------------
#define PROCESS_STARTING 0   // Loading state
#define PROCESS_RUN      1   // Running state
#define PROCESS_READY    2   // Waiting state
#define PROCESS_KILLED   3   // Process killed
#define PROCESS_ERROR    4   // Process in error
#define PROCESS_DONE     5   // Simulation ended
#define PROCESS_RUNAC    6   // Computing AC matrix

#define COMMAND_NONE     10  // No change
#define COMMAND_LOAD     11  // Load geometry
#define COMMAND_START    12  // Start simu
#define COMMAND_PAUSE     13  // Pause simu
#define COMMAND_RESET    14  // Reset simu
#define COMMAND_EXIT     15  // Exit
#define COMMAND_CLOSE    16  // Release handles
#define COMMAND_LOADAC   17  // Load mesh and compute AC matrix
#define COMMAND_STEPAC   18  // Perform single iteration step (AC)

static const char *prStates[] = {

"Not started",
"Running",
"Waiting",
"Killed",
"Error",
"Done",
"Computing AC",
"",
"",
"",
"No command",
"Load",
"Start",
"Stop",
"Reset",
"Exit",
"Close",
"ComputeAC",
"Step"

};

typedef struct {
  
  // Process control
  int    states[MAX_PROCESS];        // Process states/commands
  int    cmdParam[MAX_PROCESS];      // Command param 1
  llong  cmdParam2[MAX_PROCESS];     // Command param 2
  char   statusStr[MAX_PROCESS][64]; // Status message
  //float  heartBeat; //a changing int number showing that molflow.exe is alive
} SHMASTER;

// -----------------------------------------------------------------
// Geometry shared structure        (name: MFLWLOAD[masterPID])
//
//  SHGEOM
//  VERTEX3D
//  SHFACET1
//  INDEXF1
//  VERTEX2DF1
//  SHFACET2
//  INDEXF2
//  VERTEX2DF2
//  ...
//  AREA ELEMENTS
// -----------------------------------------------------------------

typedef struct {

  int        nbFacet;   // Number of facets (total)
  int        nbVertex;  // Number of 3D vertices
  int        nbSuper;   // Number of superstructures
  char       name[64];  // (Short file name)
  
  size_t nbMoments; //To pass in advance for memory reservation
  double latestMoment;  
  double totalDesorbedMolecules; //Number of molecules desorbed between t=0 and latest_moment
  double finalOutgassingRate; //Number of outgassing molecules / second at latest_moment (constant flow)
  double gasMass;
  double halfLife;
  double timeWindowSize;
  BOOL useMaxwellDistribution; //TRUE: Maxwell-Boltzmann distribution, FALSE: All molecules have the same (V_avg) speed
  BOOL calcConstantFlow;

  /* //Vectors must be serialized
  std::vector<std::vector<std::pair<double,double>>> CDFs; //cumulative distribution function for each temperature
  std::vector<std::vector<std::pair<double,double>>> IDs; //integrated distribution function for each time-dependent desorption type
  std::vector<Parameter> parameters; //all parameters which are time-dependent
  std::vector<double> temperatures; //keeping track of all temperatures that have a CDF already generated
  std::vector<double> moments;             //moments when a time-dependent simulation state is recorded
  std::vector<size_t> desorptionParameterIDs; //time-dependent parameters which are used as desorptions, therefore need to be integrated
  */

} SHGEOM;

typedef struct {

  // Facet parameters
  double sticking;       // Sticking (0=>reflection  , 1=>absorption)   - can be overridden by time-dependent parameter
  double opacity;        // opacity  (0=>transparent , 1=>opaque)       - can be overridden by time-dependent parameter
  double temperature;    // Facet temperature (Kelvin)                  - can be overridden by time-dependent parameter
  double flow;           // Desorbed flow (in unit *m^3/s) (outgassing) - can be overridden by time-dependent parameter

  int sticking_paramId;    // -1 if use constant value, 0 or more if referencing time-dependent parameter
  int opacity_paramId;     // -1 if use constant value, 0 or more if referencing time-dependent parameter
  int outgassing_paramId;  // -1 if use constant value, 0 or more if referencing time-dependent parameter
  
  int CDFid; //Which probability distribution it belongs to (one CDF per temperature)
  int IDid;  //If time-dependent desorption, which is its ID

  double mass;           // Molecule mass of desorbed flow (in u,u=1.660538782E-27 kg) [CURRENTLY UNUSED, gas mass is a global setting]
  double area;           // Facet area (m^2)
  int    desorbType;     // Desorption type
  double desorbTypeN;    // Exponent in Cos^N desorption type
  int    reflectType;    // Reflection type
  int    profileType;    // Profile type
  int    superIdx;       // Super structure index (Indexed from 0)
  int    superDest;      // Super structure destination index (Indexed from 1, 0=>current)
  int	 teleportDest;   // Teleport destination facet id (for periodic boundary condition) (Indexed from 1, 0=>none)
  BOOL   countDes;       // Count desoprtion (MC texture)
  BOOL   countAbs;       // Count absoprtion (MC texture)
  BOOL   countRefl;      // Count reflection (MC texture)
  BOOL   countTrans;     // Count transparent (MC texture)
  BOOL   countACD;       // Angular coefficient (AC texture)
  BOOL   countDirection; // Record avergare direction (MC texture)
  double maxSpeed;       // Max expected particle velocity (for velocity histogram)
  double accomodationFactor; // Thermal accomodation factor [0..1]

  // Flags
  BOOL   is2sided;     // 2 sided
  BOOL   isProfile;    // Profile facet
  //BOOL   isOpaque;     // Opacity != 0
  BOOL   isTextured;   // texture
  BOOL   isVolatile;   // Volatile facet (absorbtion facet which does not affect particule trajectory)

  //These don't have to be shared
  //BOOL   isVolumeVisible;    //Do we paint the volume on this facet?
  //BOOL   isTextureVisible;	 //Do we paint the texture on this facet?

  // Global hit counters
  SHHITS counter;

  // Normal vector
  VERTEX3D    N;    // normalized
  VERTEX3D    Nuv;  // normal to (u,v) not normlized

  // Axis Aligned Bounding Box (AABB)
  AABB       bb;
  VERTEX3D   center;

  // Geometry
  int    nbIndex;   // Number of index/vertex
  double sign;      // Facet vertex rotation (see Facet::DetectOrientation())

  // Plane basis (O,U,V) (See Geometry::InitializeGeometry() for info)
  VERTEX3D   O;  // Origin
  VERTEX3D   U;  // U vector
  VERTEX3D   V;  // V vector
  VERTEX3D   nU; // Normalized U
  VERTEX3D   nV; // Normalized V

  // Hit/Abs/Des/Density recording on 2D texture map
  int    texWidth;    // Rounded texture resolution (U)
  int    texHeight;   // Rounded texture resolution (V)
  double texWidthD;   // Actual texture resolution (U)
  double texHeightD;  // Actual texture resolution (V)

  int   hitOffset;      // Hit address offset for this facet
  
  //Outgassing map
  BOOL   useOutgassingFile; //has desorption file for cell elements
  double outgassingFileRatio; //desorption file's sample/unit ratio
  int   outgassingMapWidth;
  int   outgassingMapHeight;

} SHFACET;

// -----------------------------------------------------------------
// Mesh shared structure  (name: MFLWLOAD[masterPID])
//
//  SHELEM

typedef struct {

  float   area;     // Area of element
  float   uCenter;  // Center coordinates
  float   vCenter;  // Center coordinates
  int     elemId;   // Element index (MESH array)
  BOOL    full;     // Element is full

} SHELEM;

#endif /* SHAREDH */
