/*
  File:        Types.h
  Description: Various type/macro definitions and util functions
  Program:     MolFlow
  Author:      R. KERSEVAN / J-L PONS / M ADY
  Copyright:   E.S.R.F / CERN / CERN

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#ifndef TYPESH
#define TYPESH

// 64 bit integer declaration
typedef long long llong;

// AC matrix floating type

typedef float ACFLOAT;

// Desorption type

#define DES_NONE    0   // No desorption
#define DES_UNIFORM 1   // Uniform
#define DES_COSINE  2   // cos(theta)
/*#define DES_COSINE2 3   // cos(theta)^2
#define DES_COSINE3 4   // cos(theta)^3
#define DES_COSINE4 5   // cos(theta)^4*/
#define DES_COSINE_N 3 // cos(theta)^N

// Reflection type
#define REF_DIFFUSE 0   // Diffuse (cosine law)
#define REF_MIRROR  1   // Mirror
#define REF_UNIFORM 2   // Uniform (for testing)



  // Profile type

#define REC_NONE       0  // No recording
#define REC_PRESSUREU  1  // Pressure and density profile (U direction)
#define REC_PRESSUREV  2  // Pressure and density profile (V direction)
#define REC_ANGULAR    3  // Angular profile
#define REC_VELOCITY   4 //Velocity distribution
#define REC_ORT_VELOCITY 5 //Orthogonal velocity component

#define PROFILE_SIZE  100 // Size of profile



// Hit type

#define HIT_DES   1
#define HIT_ABS   2
#define HIT_REF   3
#define HIT_TRANS 4
#define HIT_TELEPORT 5
#define HIT_MOVING 6
#define LASTHIT 7

// Geometry structure definitions
class VERTEX3D {
public:


  double x;
  double y;
  double z;
  int selected;
};

struct VERTEX2D {


	double u;
	double v;
};

typedef struct {

  VERTEX3D min;
  VERTEX3D max;

} AABB;

typedef struct {

  VERTEX3D pos;
  int      type;





} HIT;

typedef struct {

  VERTEX3D pos;
  VERTEX3D dir;

} LEAK;

typedef struct {

  int nbPts;
  VERTEX2D *pts;

} MESH;

// Density/Hit field stuff
#define HITMAX 1E38

typedef struct {
	llong count;
	double sum_v_ort;
	double sum_1_per_ort_velocity;
} APROFILE;



typedef struct {
	llong count;
	double sum_v_ort_per_area;
	double sum_1_per_ort_velocity;
} AHIT;

// Velocity field
typedef struct {
	VERTEX3D sumDir;
	llong count;
} VHIT;

//Texture limit types
typedef struct {
	double all;
	double moments_only;
} TEXTURE_MOMENT_TYPE;

typedef struct {
	TEXTURE_MOMENT_TYPE min;
	TEXTURE_MOMENT_TYPE max;
} TEXTURE_MIN_MAX;

typedef struct {
	TEXTURE_MIN_MAX manual;
	TEXTURE_MIN_MAX autoscale;
} TEXTURE_SCALE_TYPE;

/*typedef struct {
	TEXTURE_SCALE_TYPE pressure;
	TEXTURE_SCALE_TYPE impingement;
	TEXTURE_SCALE_TYPE density;
} TEXTURE_LIMITS;*/

struct ProjectedPoint {
	VERTEX2D vertex2d;
	size_t globalId;
};

#define IS_ZERO(x) (fabs((x))<1e-10)

#define DOT2(x1,y1,x2,y2) ((x1)*(x2) + (y1)*(y2))
#define DOT3(x1,y1,z1,x2,y2,z2) ((x1)*(x2) + (y1)*(y2) + (z1)*(z2))

#define DET22(_11,_12,_21,_22) ( (_11)*(_22) - (_21)*(_12) )
#define DET33(_11,_12,_13,_21,_22,_23,_31,_32,_33)  \
  ((_11)*( (_22)*(_33) - (_32)*(_23) ) +            \
   (_12)*( (_23)*(_31) - (_33)*(_21) ) +            \
   (_13)*( (_21)*(_32) - (_31)*(_22) ))


#endif /* TYPESH */
