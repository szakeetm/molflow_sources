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

#ifndef SHAREDH
#define SHAREDH

#include "MolflowTypes.h"
#include <Windows.h>
#include <vector>
#include "Vector.h"

#define PROFILE_SIZE  (size_t)100 // Size of profile
#define LEAKCACHESIZE     (size_t)2048  // Leak history max length
#define HITCACHESIZE      (size_t)2048  // Max. displayed number of lines and Porto (OPO)hits.
#define BOUNCEMAX   (size_t)8192  // 'Wall collision count before absoprtion' histogram
#define MAX_PROCESS (size_t)32    // Maximum number of process

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

typedef float ACFLOAT;

typedef union {

	struct {
		// Counts
		llong nbDesorbed;          // Number of desorbed molec
		llong nbHit;               // Number of hits
		llong nbAbsorbed;          // Number of absorbed molec
		double sum_1_per_ort_velocity;    // sum of reciprocials of orthogonal velocity components, used to determine the density, regardless of facet orientation
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

	Vector3d pos;
	int      type;
} HIT;

// Velocity field
typedef struct {
	Vector3d dir;
	llong count;
} VHIT;

typedef struct {

	Vector3d pos;
	Vector3d dir;

} LEAK;

typedef struct {

	int    mode;                // Simu mode (MC_MODE or AC_MODE)

	SHHITS total;               // Global counts
	size_t hitCacheSize;              // Number of valid hits in cache
	size_t lastHitIndex;					//Index of last recorded hit in gHits (turns over when reaches HITCACHESIZE)
	HIT    hitCache[HITCACHESIZE];       // Hit history

	size_t  lastLeakIndex;		  //Index of last recorded leak in gHits (turns over when reaches LEAKCACHESIZE)
	size_t  leakCacheSize;        //Number of valid leaks in the cache
	size_t  nbLeakTotal;         // Total leaks
	LEAK   leakCache[LEAKCACHESIZE];      // Leak history

	TEXTURE_MIN_MAX texture_limits[3]; //Min-max on texture
   /* AHIT   minHit;              // Minimum on texture
	AHIT   maxHit;              // Maximum on texture
	AHIT   minHitMomentsOnly;   // Minimum, not counting constant flow
	AHIT   maxHitMomentsOnly;   // Maximum, not counting constant flow*/
	//llong  wallHits[BOUNCEMAX]; // 'Wall collision count before absoprtion' density histogram
	double distTraveledTotal_total;
	double distTraveledTotal_fullHitsOnly;

} SHGHITS;

struct AnglemapParams {
	bool   record; // Record incident angle 2-dim distribution
	bool hasRecorded;
	size_t phiWidth; //resolution between -PI and +PI
	double thetaLimit; //angle map can have a different resolution under and over the limit. Must be between 0 and PI/2
	size_t thetaLowerRes; //resolution between 0 and angleMapThetaLimit
	size_t thetaHigherRes; //resolution between angleMapThetaLimit and PI/2
} ;

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
"Loading",
"Starting",
"Stopping",
"Resetting",
"Exiting",
"Closing",
"ComputeAC",
"Step"

};

struct Reflection {
	double diffusePart;
	double specularPart;
};

typedef struct {

	// Process control
	int		states[MAX_PROCESS];        // Process states/commands
	size_t    cmdParam[MAX_PROCESS];      // Command param 1
	llong		cmdParam2[MAX_PROCESS];     // Command param 2
	char		statusStr[MAX_PROCESS][128]; // Status message
} SHCONTROL;

// -----------------------------------------------------------------
// Geometry shared structure        (name: MFLWLOAD[masterPID])
//
//  SHGEOM
//  Vector3d
//  SHFACET1
//  INDEXF1
//  Vector2dF1
//  SHFACET2
//  INDEXF2
//  Vector2dF2
//  ...
//  AREA ELEMENTS
// -----------------------------------------------------------------

typedef struct {

	size_t        nbFacet;   // Number of facets (total)
	size_t        nbVertex;  // Number of 3D vertices
	size_t        nbSuper;   // Number of superstructures
	char       name[64];  // (Short file name)

	size_t nbMoments; //To pass in advance for memory reservation
	double latestMoment;
	double totalDesorbedMolecules; //Number of molecules desorbed between t=0 and latest_moment
	double finalOutgassingRate; //Number of outgassing molecules / second at latest_moment (constant flow)
	double gasMass;
	bool enableDecay;
	double halfLife;
	double timeWindowSize;
	bool useMaxwellDistribution; //true: Maxwell-Boltzmann distribution, false: All molecules have the same (V_avg) speed
	bool calcConstantFlow;

	int motionType;
	Vector3d motionVector1; //base point for rotation
	Vector3d motionVector2; //rotation vector or velocity vector

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
	double outgassing;           // (in unit *m^3/s)                      - can be overridden by time-dependent parameter

	int sticking_paramId;    // -1 if use constant value, 0 or more if referencing time-dependent parameter
	int opacity_paramId;     // -1 if use constant value, 0 or more if referencing time-dependent parameter
	int outgassing_paramId;  // -1 if use constant value, 0 or more if referencing time-dependent parameter

	int CDFid; //Which probability distribution it belongs to (one CDF per temperature)
	int IDid;  //If time-dependent desorption, which is its ID

	double mass;           // Molecule mass of desorbed flow (in u,u=1.660538782E-27 kg) [CURRENTLY UNUSED, gas mass is a global setting]
	double area;           // Facet area (m^2)
	int    desorbType;     // Desorption type
	double desorbTypeN;    // Exponent in Cos^N desorption type
	Reflection reflection;
	int    profileType;    // Profile type
	size_t    superIdx;       // Super structure index (Indexed from 0)
	size_t    superDest;      // Super structure destination index (Indexed from 1, 0=>current)
	int	 teleportDest;   // Teleport destination facet id (for periodic boundary condition) (Indexed from 1, 0=>none, -1=>teleport to where it came from)
	bool   countDes;       // Count desoprtion (MC texture)
	bool   countAbs;       // Count absoprtion (MC texture)
	bool   countRefl;      // Count reflection (MC texture)
	bool   countTrans;     // Count transparent (MC texture)
	bool   countACD;       // Angular coefficient (AC texture)
	bool   countDirection; // Record avergare direction (MC texture)
	double maxSpeed;       // Max expected particle velocity (for velocity histogram)
	double accomodationFactor; // Thermal accomodation factor [0..1]
	bool   enableSojournTime;
	double sojournFreq, sojournE;

	// Flags
	bool   is2sided;     // 2 sided
	bool   isProfile;    // Profile facet
	bool   isTextured;   // texture
	bool   isVolatile;   // Volatile facet (absorbtion facet which does not affect particule trajectory)

	// Facet hit counters
	// SHHITS counter; - removed as now it's time-dependent and part of the hits buffer

	// Normal vector
	Vector3d    N;    // normalized
	Vector3d    Nuv;  // normal to (u,v) not normlized

	// Axis Aligned Bounding Box (AABB)
	AABB       bb;
	Vector3d   center;

	// Moving facets
	bool isMoving;

	// Geometry
	size_t    nbIndex;   // Number of index/vertex
	double sign;      // Facet vertex rotation (see Facet::DetectOrientation())

	// Plane basis (O,U,V) (See Geometry::InitializeGeometry() for info)
	Vector3d   O;  // Origin
	Vector3d   U;  // U vector
	Vector3d   V;  // V vector
	Vector3d   nU; // Normalized U
	Vector3d   nV; // Normalized V

	// Hit/Abs/Des/Density recording on 2D texture map
	size_t    texWidth;    // Rounded texture resolution (U)
	size_t    texHeight;   // Rounded texture resolution (V)
	double texWidthD;   // Actual texture resolution (U)
	double texHeightD;  // Actual texture resolution (V)

	size_t   hitOffset;      // Hit address offset for this facet

	//Outgassing map
	bool   useOutgassingFile;   //has desorption file for cell elements
	double outgassingFileRatio; //desorption file's sample/unit ratio
	size_t   outgassingMapWidth; //rounded up outgassing file map width
	size_t   outgassingMapHeight; //rounded up outgassing file map height

	double totalOutgassing; //total outgassing for the given facet

	AnglemapParams anglemapParams;//Incident angle map

} SHFACET;

// -----------------------------------------------------------------
// Mesh shared structure  (name: MFLWLOAD[masterPID])
//
//  SHELEM

//Just for AC matrix calculation in Molflow, old mesh structure:
typedef struct {

	float   area;     // Area of element
	float   uCenter;  // Center coordinates
	float   vCenter;  // Center coordinates
	int     elemId;   // Element index (MESH array)
	bool    full;     // Element is full

} SHELEM;

#endif /* SHAREDH */
