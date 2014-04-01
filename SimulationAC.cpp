/*
  File:        SimulationAC.c
  Description: Angular Coefficient Simulation for UHV (Physics related routines)
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

#ifdef WIN32
#include <windows.h> // For GetTickCount()
#include <Process.h> // For _getpid()
#else
#include <time.h>
#include <sys/time.h>
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "Simulation.h"
#include "Random.h"

extern void SetErrorSub(char *message);
extern void SetState(int state,char *status);
extern char *GetSimuStatus();
extern int GetLocalState();
//extern void GetState(int sleepTime);
extern void GetState();

// -------------------------------------------------------
// Global handles
// -------------------------------------------------------

extern SIMULATION *sHandle;

// -------------------------------------------------------

// Loop over all opaque element
#define LOOP(_k,_f,_i,_j,_idx)                                      \
for(_k=0;_k<sHandle->str[0].nbFacet;_k++) {                         \
  (_f) = sHandle->str[0].facets[_k];                                \
  for(_j=0;_j<(_f)->sh.texHeight && (_f)->sh.opacity==1.0;_j++) {   \
    for(_i=0;_i<(_f)->sh.texWidth;_i++) {

#define END_LOOP(_f,_idx)                                                  \
  _idx++;}}                                                                \
if((_f)->sh.opacity!=1.0) _idx+=((_f)->sh.texWidth*(_f)->sh.texHeight); }


// -------------------------------------------------------

void GetCenter(FACET *f,SHELEM *mesh,int idx,VERTEX3D *c) {

  c->x = f->sh.O.x + f->sh.U.x*mesh[idx].uCenter + f->sh.V.x*mesh[idx].vCenter;
  c->y = f->sh.O.y + f->sh.U.y*mesh[idx].uCenter + f->sh.V.y*mesh[idx].vCenter;
  c->z = f->sh.O.z + f->sh.U.z*mesh[idx].uCenter + f->sh.V.z*mesh[idx].vCenter;

}

// -------------------------------------------------------

void ClearACMatrix() {

  SAFE_FREE(sHandle->acMatrix);
  SAFE_FREE(sHandle->acDensity);
  SAFE_FREE(sHandle->acDesorb);
  SAFE_FREE(sHandle->acAbsorb);
  SAFE_FREE(sHandle->acRho);
  SAFE_FREE(sHandle->acArea);
  SAFE_FREE(sHandle->acLines);
  SAFE_FREE(sHandle->acTMatrix);
  SAFE_FREE(sHandle->acTDensity);
  SAFE_FREE(sHandle->acTArea);
  SAFE_FREE(sHandle->acTLines);
  sHandle->nbAC = 0;
  sHandle->nbACT = 0;
  sHandle->prgAC = 0;

#ifdef JACOBI_ITERATION
  SAFE_FREE(sHandle->acDensityTmp);
#endif

}

// -------------------------------------------------------

BOOL ComputeACMatrix(SHELEM *mesh) {

  int      idx,idx1,idx2,i1,i2,j1,j2,k1,k2,nbElem;
  FACET   *f1,*f2;
  int      sz,nbO=0,nbB=0,nbE=0,nbV;
  VERTEX3D c1,c2;
  double   r2,cos1,cos2,vf,pv;
  double   t0,t1;
  int      p;

  t0 = GetTick();

  idx1 = 0;
  idx2 = 0;

  SetState(PROCESS_RUNAC,GetSimuStatus());

  // Count number of elements
  sHandle->nbAC = 0;
  for(k1=0;k1<sHandle->str[0].nbFacet;k1++) {
    f1 = sHandle->str[0].facets[k1];
    if(f1->sh.opacity==1.0f) {
      sHandle->nbAC += f1->sh.texHeight * f1->sh.texWidth;
    } else if(f1->sh.opacity==0.0f) {
      sHandle->nbACT += f1->sh.texHeight * f1->sh.texWidth;
    } else {
      // partial transparent facet not supported with AC
      SetErrorSub("AC does not handle partial opacity");
      return FALSE;
    }
  }

  // Allocate memory for angular coefficient 
  // (we keep only the stricly lower triangular part)
  nbElem = (sHandle->nbAC * (sHandle->nbAC-1))/2;
  sz = sizeof(ACFLOAT) * nbElem;
  sHandle->acMatrix = (ACFLOAT *)malloc(sz);
  if( !sHandle->acMatrix ) {
    SetErrorSub("Not enough memory for AC matrix");
    return FALSE;
  }
  memset(sHandle->acMatrix,0,sz);

  // Allocate memory for various vectors
  sz = sizeof(ACFLOAT) * sHandle->nbAC;

  sHandle->acArea = (ACFLOAT *)malloc(sz);
  if( !sHandle->acArea ) {
    SetErrorSub("Not enough memory for AC matrix");
    return FALSE;
  }

  sHandle->acDensity = (ACFLOAT *)malloc(sz);
  if( !sHandle->acDensity ) {
    SetErrorSub("Not enough memory for AC matrix");
    return FALSE;
  }

#ifdef JACOBI_ITERATION
  sHandle->acDensityTmp = (ACFLOAT *)malloc(sz);
  if( !sHandle->acDensityTmp ) {
    SetErrorSub("Not enough memory for AC matrix");
    return FALSE;
  }
#endif

  sHandle->acDesorb = (ACFLOAT *)malloc(sz);
  if( !sHandle->acDesorb ) {
    SetErrorSub("Not enough memory for AC matrix");
    return FALSE;
  }

  sHandle->acAbsorb = (ACFLOAT *)malloc(sz);
  if( !sHandle->acAbsorb ) {
    SetErrorSub("Not enough memory for AC matrix");
    return FALSE;
  }

  sHandle->acRho = (ACFLOAT *)malloc(sz);
  if( !sHandle->acRho ) {
    SetErrorSub("Not enough memory for AC matrix");
    return FALSE;
  }

  sHandle->acLines = (double *)malloc(sizeof(double) * sHandle->nbAC);
  if( !sHandle->acLines ) {
    SetErrorSub("Not enough memory for AC matrix");
    return FALSE;
  }

  SetState(PROCESS_RUNAC,GetSimuStatus());

  // Init vectors
  idx = 0;
  memset(sHandle->acDensity,0,sz);
  memset(sHandle->acDesorb,0,sz);
  memset(sHandle->acAbsorb,0,sz);
  memset(sHandle->acLines,0,sizeof(double) * sHandle->nbAC);
  LOOP(k1,f1,i1,j1,idx1)
    sHandle->acArea[idx] = mesh[idx1].area;
    // Normalize density for Q=1,M=28,T=20degC
    // TODO: Update for random gas mass/temperature
    if(f1->sh.desorbType) sHandle->acDesorb[idx] = (ACFLOAT)(1.0/(11.7*sHandle->sourceArea));
    sHandle->acRho[idx] = (ACFLOAT)(1.0 - f1->sh.sticking);
    idx++;
  END_LOOP(f1,idx1)

  // Compute AC matrix
  idx1 = 0;
  idx = 0;
  LOOP(k1,f1,i1,j1,idx1)

    GetCenter(f1,mesh,idx1,&c1);

    idx2=0;
    LOOP(k2,f2,i2,j2,idx2)

      if( idx2<idx1 ) {

        // Progress
        p = (int)( ((double)idx * 99.0 /(double)nbElem) + 0.5);
        if( sHandle->prgAC!=p ) {
          sHandle->prgAC = p;
          GetState();
          if(GetLocalState()==COMMAND_PAUSE) {
            sHandle->prgAC=0;
            return FALSE;
          }
          SetState(PROCESS_RUNAC,GetSimuStatus());
        }

        if( mesh[idx2].area>0.0f && mesh[idx1].area>0.0f ) {

          // Compute view factor
          // (Area of surface element is included whithin the iteration)
          GetCenter(f2,mesh,idx2,&c2);

          r2 = DOT3(c1.x-c2.x,c1.y-c2.y,c1.z-c2.z,
                    c1.x-c2.x,c1.y-c2.y,c1.z-c2.z);
          // cos1 = cos(theta1) * r
          cos1 = DOT3(f1->sh.N.x  , f1->sh.N.y , f1->sh.N.z ,
                      c2.x - c1.x , c2.y - c1.y , c2.z - c1.z);
          // cos2 = cos(theta2) * r
          cos2 = DOT3(f2->sh.N.x  , f2->sh.N.y , f2->sh.N.z ,
                      c1.x - c2.x , c1.y - c2.y , c1.z - c2.z);

          if(cos1>0.0 && cos2>0.0 && r2>0.0) {
            if( Visible(&c1,&c2,f1,f2) ) {
              vf = (cos1 * cos2) / (PI * r2 * r2);
              sHandle->acMatrix[idx] = (ACFLOAT)vf;
            } else {
              nbO++; // Obstacle
            }
          } else {
            nbB++; // Not visible (Back to back)
          }         
        } else {
          nbB++; // Not visible (outside geometry)
        }
        idx++;
      }

    END_LOOP(f2,idx2); // End loop k2

  END_LOOP(f1,idx1); // End loop k1


  // Avoid divergence by renormalizing lines

  for(i1=0;i1<sHandle->nbAC;i1++) {

      double sum = 0.0;
      int inc;

      // AC matrix format (strictly lower triangular part, diagonal not included)
      // 0 
      // 1 2
      // 3 4 5
      // ...

      idx = (i1*i1 - i1)/2;
      for(j1=0;j1<i1;j1++,idx++)
        sum +=  sHandle->acMatrix[idx] * sHandle->acArea[j1];
      idx = (i1*i1 + 3*i1)/2;
      for(j1=i1+1,inc=i1;j1<sHandle->nbAC;j1++,idx+=inc) {
        sum +=  sHandle->acMatrix[idx] * sHandle->acArea[j1];
        inc++;
      }

      if(sum>1e-10) sHandle->acLines[i1] = 1.0 / sum;

  }

  t1 = GetTick();
  printf("AC matrix calculation succesful (%dx%d)\n",sHandle->nbAC,sHandle->nbAC);
  nbV = nbElem-nbO-nbB;
  pv = (double)nbV * 100.0 / (double)nbElem;
  printf("Obstacle:%d Nvisible:%d Not null:%d (%.2f%%)\n",nbO,nbB,nbV,pv);
  sHandle->calcACTime = (t1-t0);
  printf("Calculation time: %.3f s\n",sHandle->calcACTime);
  sHandle->prgAC = 100; // AC matrix calculation done

  return TRUE;

}

BOOL SimulationACStep(int nbStep) {

  int      i,inc,j,idx,step;

  if( sHandle->prgAC!=100 ) {
    return FALSE;
  }

  step = 0;
  // Run iterations
  while( (sHandle->maxDesorption==0 || sHandle->nbDesorbed<sHandle->maxDesorption) && step<nbStep ) {

    // Perform iteration
    // density[i] = rho[i] * Sum{j=0,nbAC-1}{ AC(i,j)*area(j)*density(j) } + desorb[i]
    // rho = 1-sticking
    idx = 0;
    for(i=0;i<sHandle->nbAC;i++) {

      double sum = 0.0;
      ACFLOAT fSum;

      // AC matrix format (strictly lower triangular part, diagonal not included)
      // 0 
      // 1 2
      // 3 4 5
      // ...

      if( sHandle->acLines[i]>0.0 ) {
      
        idx = (i*i - i)/2;
        for(j=0;j<i;j++,idx++)
          sum +=  (sHandle->acMatrix[idx] * sHandle->acDensity[j] * sHandle->acArea[j]);
        idx = (i*i + 3*i)/2;
        for(j=i+1,inc=i;j<sHandle->nbAC;j++,idx+=inc) {
          sum +=  (sHandle->acMatrix[idx] * sHandle->acDensity[j] * sHandle->acArea[j]);
          inc++;
        }

        fSum = (ACFLOAT)(sum * sHandle->acLines[i]);

#ifdef JACOBI_ITERATION
        sHandle->acDensityTmp[i] = sHandle->acRho[i] * fSum + sHandle->acDesorb[i];
#else
        sHandle->acDensity[i] = sHandle->acRho[i] * fSum + sHandle->acDesorb[i];
#endif
        sHandle->acAbsorb[i] = (1.0f - sHandle->acRho[i]) * fSum;

      } else {

#ifdef JACOBI_ITERATION
        sHandle->acDensityTmp[i] = 0.0;
#endif
        sHandle->acDensity[i] = 0.0;
        sHandle->acAbsorb[i] = 0.0;

      }

    }

#ifdef JACOBI_ITERATION
    memcpy(sHandle->acDensity , sHandle->acDensityTmp, sizeof(ACFLOAT)*sHandle->nbAC );
#endif
    sHandle->nbDesorbed++;
    step++;
  
  }

  return sHandle->maxDesorption==0 || sHandle->nbDesorbed<sHandle->maxDesorption;

}

void UpdateACHits(Dataport *dpHit,int prIdx,DWORD timeout) {

  SHGHITS *gHits;
  FACET   *f;
  AHIT    *shTexture;
  SHHITS  *fHits;
  int      i,j,k,idx,nbE;
  double  sumVal=0.0;
  double  sumAbs=0.0;
  double  sumDes=0.0;

  if( !AccessDataportTimed(dpHit,timeout) ) return;

  gHits = (SHGHITS *)dpHit->buff;
  gHits->texture_limits[0].max.all = 0.0;
  gHits->texture_limits[0].min.all = 0.0;
  gHits->mode = AC_MODE;
  gHits->total.hit.nbDesorbed = sHandle->nbDesorbed;

  // Update texture
  idx = 0;
  for(k=0;k<sHandle->str[0].nbFacet;k++) {

    f = sHandle->str[0].facets[k];
    sumDes = 0.0;
    sumAbs = 0.0;
    sumVal = 0.0;
    nbE = 0;
    shTexture = (AHIT *)((char *)dpHit->buff + (f->sh.hitOffset + sizeof(SHHITS) + f->profileSize));
    for(j=0;j<f->sh.texHeight && f->sh.opacity==1.0f;j++) {
      for(i=0;i<f->sh.texWidth;i++) {
        ACFLOAT val = (ACFLOAT)(sHandle->acDensity[idx]);
        shTexture[i + j*f->sh.texWidth].count = val;
        if( val > gHits->texture_limits[0].max.all ) gHits->texture_limits[0].max.all = val;
        // Normalize over facet area
        sumVal += (double)(sHandle->acDensity[idx]*sHandle->acArea[idx]/f->sh.area);
        sumAbs += (double)(sHandle->acAbsorb[idx]*sHandle->acArea[idx]/f->sh.area);
        sumDes += (double)(sHandle->acDesorb[idx]*sHandle->acArea[idx]/f->sh.area);
        idx++;
      }
    }
    fHits = (SHHITS  *)((char *)dpHit->buff + f->sh.hitOffset);
    fHits->density.value    = sumVal;
    fHits->density.desorbed = sumDes;
    fHits->density.absorbed = sumAbs;

  }


  ReleaseDataport(dpHit);

}