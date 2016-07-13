/*
  File:        Utils.h
  Description: Various util functions
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

#ifndef _UTILSH_
#define _UTILSH_

#include "Types.h"

typedef struct {

  VERTEX2D  *pts;   // Array of 2D vertex
  int        nbPts; // Number of vertex
  double     sign;  // Polygon orientation

} POLYGON;

typedef struct  {

  VERTEX2D  p;       // Vertex coordinates
  int       mark;    // Cycle detection (0=>not processed, 1=>processed)
  int       isStart; // Possible starting point

  int       nbOut;  // Number of outgoing arc
  int       nbIn;   // Number of incoming arc
  int       VI[2];  // Tangent point detection
  int       VO[2];  // Tangent point detection

} POLYVERTEX;

typedef struct {

  int i1;  // Node 1 index
  int i2;  // Node 2 index
  int s;   // Source polygon (tangent point detection)

} POLYARC;

typedef struct {

  int         nbNode;  // Number of node
  POLYVERTEX *nodes;   // Nodes
  int         nbArc;   // Number of arc
  POLYARC    *arcs;    // Arcs

} POLYGRAPH;

#define IDX(i,nb) (((i)<0)?nb+(i):(i)%(nb))


VERTEX3D CrossProduct(const VERTEX3D &v1, const VERTEX3D &v2);
VERTEX3D operator+ (const VERTEX3D &v1, const VERTEX3D& v2);
VERTEX3D operator-(const VERTEX3D &v1, const VERTEX3D& v2);
double Dot(const VERTEX3D &v1, const VERTEX3D &v2);
VERTEX3D operator*(const VERTEX3D &v1, const double& mult);
VERTEX3D operator*(const double& mult, const VERTEX3D &v1);

double Dot(VERTEX3D *v1,VERTEX3D *v2);
double Dot(const VERTEX2D &v1, const VERTEX2D &v2);
void   Cross(VERTEX3D *result,VERTEX3D *v1,VERTEX3D *v2);
double Norme(VERTEX3D *v);
double Norme(const VERTEX2D &v);
double Norme(VERTEX2D *v1,VERTEX2D *v2);
void   Normalize(VERTEX3D *v);
void   Normalize(VERTEX2D *v);
void   ScalarMult(VERTEX3D *result,double r);
void   ScalarMult(VERTEX2D *result, double r);
void   Sub(VERTEX3D *result,VERTEX3D *v1,VERTEX3D *v2);
void   Sub(VERTEX2D *result, VERTEX2D *v1, VERTEX2D *v2);
void   Add(VERTEX3D *result,VERTEX3D *v1,VERTEX3D *v2);
void   Add(VERTEX2D *result, VERTEX2D *v1, VERTEX2D *v2);
void   ProjectVertex(VERTEX3D *v,VERTEX2D *projected,VERTEX3D *U,VERTEX3D *V,VERTEX3D *origin);
void   Mirror(VERTEX3D *P, VERTEX3D P0, VERTEX3D N);
void   Rotate(VERTEX3D *P, VERTEX3D AXIS_P0, VERTEX3D AXIS_DIR,double theta);
int   IsEqual(const double &a, const double &b, double tolerance=1E-8);
double GetOrientedAngle(VERTEX2D *v1, VERTEX2D *v2);

int    SolveIASM(double *u ,double *v,double *w,
                 VERTEX3D *nuv,VERTEX3D *U,VERTEX3D *V,VERTEX3D *W,VERTEX3D *Z);
int    SolveISSE2(double *u ,double *v,double *w,
                  VERTEX3D *nuv,VERTEX3D *U,VERTEX3D *V,VERTEX3D *W,VERTEX3D *Z);
int    GetPower2(int n);
int		Remainder(int param, int bound);
char  *FormatMemory(size_t size);
char  *FormatMemoryLL(llong size);
int   IsInsideTri(VERTEX2D *p,VERTEX2D *p1,VERTEX2D *p2,VERTEX2D *p3);
int   IsConvex(POLYGON *p,int idx);
int   ContainsConcave(POLYGON *p,int i1,int i2,int i3);
int   EmptyTriangle(POLYGON *p,int i1,int i2,int i3,VERTEX2D *center);
int   IsInPoly(double u,double v,VERTEX2D *pts,int nbPts);
int   IntersectPoly(POLYGON *p1,POLYGON *p2,int *visible2,POLYGON **result);
double GetInterArea(POLYGON *inP1,POLYGON *inP2,int *edgeVisible,float *uC,float *vC,int *nbV,double **lList);
double GetInterAreaBF(POLYGON *inP1,double u0,double v0,double u1,double v1,float *uC,float *vC);
double RoundAngle(double a);

#endif /* _UTILSH_ */
