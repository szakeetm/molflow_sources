/*
File:        SimulationMC.c
Description: Monte-Carlo Simulation for UHV (Physics related routines) 
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



#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "Simulation.h"
#include "Random.h"

extern SIMULATION *sHandle;
extern void SetErrorSub(char *message);

// -------------------------------------------------------------
// Compute area of all the desorption facet
// -------------------------------------------------------------

void ComputeSourceArea() {
	int i,j,k,l,tSize;
	FACET *f;
	double scale;
	float scale_precomputed;

	// Compute the outgassing of all source facet
	sHandle->sourceArea = 0.0;
	//sHandle->totalOutgassing = 0.0;  //calculated locally
	for(j=0;j<sHandle->nbSuper;j++) {
		for(i=0;i<sHandle->str[j].nbFacet;i++) {
			f = sHandle->str[j].facets[i];
			if (f->sh.useOutgassingFile) {
				for (l=0;l<(f->sh.outgassingMapWidth*f->sh.outgassingMapHeight);l++)
					sHandle->sourceArea+=f->outgassingMap[l];
			}
			else if( f->sh.desorbType!=DES_NONE ) sHandle->sourceArea+=f->sh.flow;  //instead of area, we count with outgassing rate (sh.flow)
		}
	}

	// Update texture increment for MC
	scale_precomputed=(float)(40.0/(sqrt(8.0*8.31/(PI*sHandle->gasMass*0.001))));
	for(j=0;j<sHandle->nbSuper;j++) {
		for(i=0;i<sHandle->str[j].nbFacet;i++) {
			f = sHandle->str[j].facets[i];
			if(f->inc) {
				tSize = f->sh.texWidth*f->sh.texHeight;
				scale = scale_precomputed;
				//if (!sHandle->nonIsothermal) //isothermal system
				scale=scale/sqrt((float)f->sh.temperature);
				if(f->sh.opacity>0.0 && f->sh.opacity!=1.0) scale = scale/f->sh.opacity;
				if(f->sh.is2sided) scale = scale / 2.0f;
				for(k=0;k<tSize;k++) f->inc[k] *= (AHIT)scale;
				f->fullSizeInc*=(float)scale;
			}
		}
	}

}





// -------------------------------------------------------

void PolarToCartesian(FACET *iFacet,double theta,double phi,BOOL reverse) {

	VERTEX3D U,V,N;
	double u,v,n;

	// Polar in (nU,nV,N) to Cartesian(x,y,z) transformation  ( nU = U/|U| , nV = V/|V| )
	// tetha is the angle to the normal of the facet N, phi to U
	// ! See Geometry::InitializeGeometry() for further informations on the (U,V,N) basis !
	// (nU,nV,N) and (x,y,z) are both left handed

	//This should be a speed-up routine, but I didn't experience any speed difference so I commented it out. Marton
	/*#ifdef WIN32
	_asm {                    // FPU stack
	fld qword ptr [theta]
	fsincos                 // cos(t)        sin(t)
	fld qword ptr [phi]
	fsincos                 // cos(p)        sin(p) cos(t) sin(t)
	fmul st(0),st(3)        // cos(p)*sin(t) sin(p) cos(t) sin(t)
	fstp qword ptr [u]      // sin(p)        cos(t) sin(t)
	fmul st(0),st(2)        // sin(p)*sin(t) cos(t) sin(t)
	fstp qword ptr [v]      // cos(t) sin(t)
	fstp qword ptr [n]      // sin(t)
	fstp qword ptr [dummy]  // Flush the sin(t)
	}
	#else*/
	u = sin(theta)*cos(phi);
	v = sin(theta)*sin(phi);
	n = cos(theta);
	//#endif

	// Get the (nU,nV,N) orthonormal basis of the facet
	U = iFacet->sh.nU;
	V = iFacet->sh.nV;
	N = iFacet->sh.N;
	if (reverse) {
		N.x=N.x*(-1.0);
		N.y=N.y*(-1.0);
		N.z=N.z*(-1.0);
	}

	// Basis change (nU,nV,N) -> (x,y,z)
	sHandle->pDir.x = u*U.x + v*V.x + n*N.x;
	sHandle->pDir.y = u*U.y + v*V.y + n*N.y;
	sHandle->pDir.z = u*U.z + v*V.z + n*N.z;

}

// -------------------------------------------------------

void CartesianToPolar(FACET *iFacet,double *theta,double *phi) {

	// Get polar coordinates of the incoming particule direction in the (U,V,N) facet space.
	// Note: The facet is parallel to (U,V), we use its (nU,nV,N) orthonormal basis here.
	// (nU,nV,N) and (x,y,z) are both left handed

	// Cartesian(x,y,z) to polar in (nU,nV,N) transformation

	// Basis change (x,y,z) -> (nU,nV,N)
	// We use the fact that (nU,nV,N) belongs to SO(3)
	double u = DOT3(sHandle->pDir.x,sHandle->pDir.y,sHandle->pDir.z,
		iFacet->sh.nU.x,iFacet->sh.nU.y,iFacet->sh.nU.z);
	double v = DOT3(sHandle->pDir.x,sHandle->pDir.y,sHandle->pDir.z,
		iFacet->sh.nV.x,iFacet->sh.nV.y,iFacet->sh.nV.z);
	double n = DOT3(sHandle->pDir.x,sHandle->pDir.y,sHandle->pDir.z,
		iFacet->sh.N.x,iFacet->sh.N.y,iFacet->sh.N.z);

	// (u,v,n) -> (theta,phi)
	double rho = sqrt( v*v + u*u );
	*theta = acos(n);              // Angle to normal (PI/2 => PI)
	*phi = asin(v/rho);
	if( u<0.0 ) *phi = PI - *phi;  // Angle to U

}

// -------------------------------------------------------

void UpdateMCHits(Dataport *dpHit,int prIdx,int nbMoments,DWORD timeout) {

	BYTE *buffer;
	SHGHITS *gHits;
	AHIT minHitsOld,minHitsMomentsOnlyOld;
	AHIT maxHitsOld,maxHitsMomentsOnlyOld;
	int i,j,s,x,y;
	int nb;
#ifdef _DEBUG
	double t0,t1;
	t0 = GetTick();
#endif
	sHandle->lastUpdateOK = AccessDataportTimed(dpHit,timeout);
	if( !sHandle->lastUpdateOK ) return;

	buffer = (BYTE*)dpHit->buff;
	gHits = (SHGHITS *)buffer;

	// Global hits and leaks
	gHits->total.hit.nbHit      += sHandle->tmpCount.hit.nbHit;
	gHits->total.hit.nbAbsorbed += sHandle->tmpCount.hit.nbAbsorbed;
	gHits->total.hit.nbDesorbed += sHandle->tmpCount.hit.nbDesorbed;
	gHits->distTraveledTotal    += sHandle->distTraveledSinceUpdate;
	minHitsOld=gHits->minHit;
	maxHitsOld=gHits->maxHit;
	minHitsMomentsOnlyOld=gHits->minHitMomentsOnly;
	maxHitsMomentsOnlyOld=gHits->maxHitMomentsOnly;
	gHits->minHit=gHits->minHitMomentsOnly=HITMAX;
	gHits->maxHit=gHits->maxHitMomentsOnly=0;
	gHits->mode = MC_MODE;
	//for(i=0;i<BOUNCEMAX;i++) gHits->wallHits[i] += sHandle->wallHits[i];

	// Leak
	nb = gHits->nbLastLeaks;
	for(i=0;i<sHandle->nbLastLeak && i<NBHLEAK;i++)
		gHits->pLeak[(i+nb) % NBHLEAK] = sHandle->pLeak[i];
	gHits->nbLeakTotal += sHandle->nbLeakTotal;
	gHits->nbLastLeaks += sHandle->nbLastLeak;

	// HHit (Only prIdx 0)
	if( prIdx==0 ) {
		gHits->nbHHit = sHandle->nbHHit;
		memcpy(gHits->pHits,sHandle->pHits,NBHHIT*sizeof(HIT));
	}

	// Facets
	for(s=0;s<sHandle->nbSuper;s++) {
		for(i=0;i<sHandle->str[s].nbFacet;i++) {

			FACET *f = sHandle->str[s].facets[i];
			if( f->hitted ) {

				SHHITS *fFit = (SHHITS *)(buffer + f->sh.hitOffset);
				fFit->hit.nbAbsorbed += f->sh.counter.hit.nbAbsorbed;
				//printf("\n%d %g + %g",i,(double)fFit->hit.nbHit,(double)f->sh.counter.hit.nbHit);
				fFit->hit.nbDesorbed += f->sh.counter.hit.nbDesorbed;
				fFit->hit.nbHit += f->sh.counter.hit.nbHit;

				if( f->sh.isProfile ) {
					for (int m=0;m<(1+nbMoments);m++) {
						llong *shProfile = (llong *)(buffer + (f->sh.hitOffset + sizeof(SHHITS) + m*f->profileSize));
						for(j=0;j<PROFILE_SIZE;j++)
							shProfile[j] += f->profile[m][j];
					}
				}

				if( f->sh.isTextured ) {
					for (int m=0;m<(1+nbMoments);m++) {
						AHIT *shTexture = (AHIT *)(buffer + (f->sh.hitOffset + sizeof(SHHITS) + f->profileSize*(1+nbMoments) + m*f->textureSize));
						AHIT correctionFactor = ((m==0)?1.0f:(float)((sHandle->desorptionStopTime-sHandle->desorptionStartTime)
							/sHandle->timeWindowSize));
						for(y=0;y<f->sh.texHeight;y++) {
							for(x=0;x<f->sh.texWidth;x++) {
								int add = x + y*f->sh.texWidth;
								AHIT val = shTexture[add] + f->hits[m][add];

								//Global autoscale
								if((val * correctionFactor)>gHits->maxHit && f->largeEnough[add])
									gHits->maxHit=val * correctionFactor;

								if (val>0.0f && (val * correctionFactor)<gHits->minHit && f->largeEnough[add])
									gHits->minHit=val * correctionFactor;

								//Autoscale ignoring constant flow
								if (m!=0) {
									if((val * correctionFactor)>gHits->maxHitMomentsOnly && f->largeEnough[add])
										gHits->maxHitMomentsOnly=val * correctionFactor;

									if (val>0.0f && (val * correctionFactor)<gHits->minHitMomentsOnly && f->largeEnough[add])
										gHits->minHitMomentsOnly=val * correctionFactor;
								}

								shTexture[add] = val;

							}
						}
					}
				}


				if( f->sh.countDirection ) {
					for (int m=0;m<(1+nbMoments);m++) {
						VHIT *shDir = (VHIT *)(buffer + (f->sh.hitOffset + sizeof(SHHITS) + f->profileSize*(1+nbMoments) + f->textureSize*(1+nbMoments) + f->directionSize*m));
						for(y=0;y<f->sh.texHeight;y++) {
							for(x=0;x<f->sh.texWidth;x++) {
								int add = x + y*f->sh.texWidth;
								shDir[add].dir.x += f->direction[m][add].dir.x;
								shDir[add].dir.y += f->direction[m][add].dir.y;
								shDir[add].dir.z += f->direction[m][add].dir.z;
								shDir[add].sumSpeed += f->direction[m][add].sumSpeed;
								shDir[add].count += f->direction[m][add].count;
							}
						}
					}
				}

			} // End if(hitted)
		} // End nbFacet
	} // End nbSuper

	//if there were no textures:
	if (gHits->minHit==HITMAX) {
		gHits->minHit=minHitsOld;
		gHits->minHitMomentsOnly=minHitsMomentsOnlyOld;
	}
	if (gHits->maxHit==0) {
		gHits->maxHit=maxHitsOld;
		gHits->maxHitMomentsOnly=maxHitsMomentsOnlyOld;
	}

	ReleaseDataport(dpHit);

	//printf("\nResetCounter called from UpdateMCHits");
	/*printf("\nTest cube pass count:%g"
		"\nAvg. test cube time (total time/#desorbed) (s):%g"
		"\nAvg. test cube dist (cm):%g"
		"\nAvg. system flight time (s):%g"
		"\nRatio (cube time/system flight time):%g"
		"\nAvg. syst distance (m):%g"
		"\nAvg. test cube temperature (K):%g"
		"\nAvg. test cube velocity (m/s):%g\n"
		"\nAvg. cube dist/Avg. velocity:%g\n\n",
		(double)sHandle->testCubeCount,
		sHandle->testCubeTime/(double)gHits->total.hit.nbDesorbed,
		sHandle->testCubeDist/(double)gHits->total.hit.nbDesorbed,
		sHandle->testSystemTime/(double)gHits->total.hit.nbDesorbed,
		sHandle->testCubeTime/sHandle->testSystemTime,
		sHandle->testSystemDist/(double)gHits->total.hit.nbDesorbed,
		sHandle->testCubeTemp/(double)sHandle->testCubeCount,
		sHandle->testCubeVelocity/(double)sHandle->testCubeCount,
		sHandle->testCubeDist/(double)gHits->total.hit.nbDesorbed/sHandle->testCubeVelocity*(double)sHandle->testCubeCount);*/
	ResetCounter();

#ifdef _DEBUG
	t1 = GetTick();
	printf("Update hits: %f us\n",(t1-t0)*1000000.0);
#endif

}






// -------------------------------------------------------------
// Compute particle teleport
// -------------------------------------------------------------

void PerformTeleport(FACET *iFacet) {

	double inPhi,inTheta;
	//Search destination
	FACET *destination;
	BOOL found=FALSE;
	BOOL revert=FALSE;
	int i,j;
	//Look in which superstructure is the destination facet:
	for (i=0;i<sHandle->nbSuper&&(!found);i++) {
		for (j=0;j<sHandle->str[i].nbFacet&&(!found);j++) {
			if ((iFacet->sh.teleportDest-1)==sHandle->str[i].facets[j]->globalId) {
				destination=sHandle->str[i].facets[j];
				sHandle->curStruct = destination->sh.superIdx; //change current superstructure
				found=TRUE;
			}
		}
	}
	if (!found) printf("Global Id %d not found",iFacet->sh.teleportDest-1);

	// Count this hit as a transparent pass
	RecordHit(HIT_TELEPORT);
	if( iFacet->hits && iFacet->sh.countTrans ) AHIT_FACET(iFacet,sHandle->flightTimeCurrentParticle);
	if( iFacet->direction && iFacet->sh.countDirection ) DHIT_FACET(iFacet,sHandle->flightTimeCurrentParticle);

	// Relaunch particle from new facet
	CartesianToPolar(iFacet,&inTheta,&inPhi);
	PolarToCartesian(destination,inTheta,inPhi,FALSE);
	// Move particle to teleport destination point
	sHandle->pPos.x = destination->sh.O.x+iFacet->colU*destination->sh.U.x+iFacet->colV*destination->sh.V.x;
	sHandle->pPos.y = destination->sh.O.y+iFacet->colU*destination->sh.U.y+iFacet->colV*destination->sh.V.y;
	sHandle->pPos.z = destination->sh.O.z+iFacet->colU*destination->sh.U.z+iFacet->colV*destination->sh.V.z;
	RecordHit(HIT_TELEPORT);
	sHandle->lastHit = destination;

}

// -------------------------------------------------------------
// Perform nbStep simulation steps (a step is a bounce)
// -------------------------------------------------------------

BOOL SimulationMCStep(int nbStep) {

	FACET   *collidedFacet;
	double   d;
	BOOL     found;
	int      i;

	// Perform simulation steps
	for(i=0;i<nbStep;i++) {

		found = Intersect(&(sHandle->pPos),&(sHandle->pDir),&d,&collidedFacet,sHandle->lastHit);

		if( found ) {

			// Move particle to intersection point
			sHandle->pPos.x += d*sHandle->pDir.x;
			sHandle->pPos.y += d*sHandle->pDir.y;
			sHandle->pPos.z += d*sHandle->pDir.z;
			sHandle->distTraveledCurrentParticle += d;
			sHandle->flightTimeCurrentParticle += d/100.0/sHandle->velocityCurrentParticle; //conversion from cm to m

			if (!sHandle->calcConstantFlow && sHandle->flightTimeCurrentParticle>sHandle->latestMoment) {
				RecordHit(LASTHIT);
				sHandle->distTraveledSinceUpdate+=sHandle->distTraveledCurrentParticle;
				if( !StartFromSource() )
					// maxDesorption reached
						return FALSE;
				RecordHit(HIT_DES);
			} else {

				if (collidedFacet->sh.teleportDest) {
					PerformTeleport(collidedFacet); 
				}	else if( collidedFacet->sh.sticking>0.0 ) {

					//if( collidedFacet->sh.sticking>0.0 && rnd()<(1000*PI/collidedFacet->sh.area/sHandle->velocityCurrentParticle/100.0) ) { //if absorbed
					if( collidedFacet->sh.sticking==1.0 || rnd()<(collidedFacet->sh.sticking) ) { //if absorbed
						//if( collidedFacet->sh.sticking==1.0 || rnd()<(collidedFacet->sh.sticking*145.469*sqrt(sHandle->temperature/sHandle->gasMass)/sHandle->velocityCurrentParticle) ) { //if absorbed
						collidedFacet->sh.counter.hit.nbAbsorbed++;
						sHandle->tmpCount.hit.nbAbsorbed++;
						//sHandle->counter.hit.nbAbsorbed++;
						sHandle->distTraveledSinceUpdate+=sHandle->distTraveledCurrentParticle;
						/*sHandle->testSystemDist+=sHandle->distTraveledCurrentParticle;
						sHandle->testSystemTime+=sHandle->flightTimeCurrentParticle;*/
						RecordHit(HIT_ABS);
						//TODO: Histogram count routine here
						if( collidedFacet->hits && collidedFacet->sh.countAbs ) AHIT_FACET(collidedFacet,sHandle->flightTimeCurrentParticle);
						if( collidedFacet->direction && collidedFacet->sh.countDirection ) DHIT_FACET(collidedFacet,sHandle->flightTimeCurrentParticle);

						//if( sHandle->nbPHit < BOUNCEMAX && !collidedFacet->sh.desorbType )  
						//  sHandle->wallHits[sHandle->nbPHit]++;

						if( !StartFromSource() )
							// maxDesorption reached
								return FALSE;
						RecordHit(HIT_DES);
					} else {
						PerformBounce(collidedFacet);
					}

				} else {
					PerformBounce(collidedFacet);
				}

				if (!collidedFacet->sh.teleportDest){
					// Hit count
					sHandle->tmpCount.hit.nbHit++;
					//sHandle->counter.hit.nbHit++;
					collidedFacet->sh.counter.hit.nbHit++;
				}
			}
		} else {

			// Leak (simulation error)
			if(sHandle->nbLastLeak<NBHLEAK) {
				// Record leak for debugging
				sHandle->pLeak[sHandle->nbLastLeak].pos = sHandle->pPos;
				sHandle->pLeak[sHandle->nbLastLeak].dir = sHandle->pDir;
				sHandle->nbLastLeak++;
			}
			sHandle->nbLeakTotal++;
			if( !StartFromSource() )
				// maxDesorption reached
					return FALSE;
		}
	}
	return TRUE;
}

// -------------------------------------------------------------
// Launch a ray from a source facet (starting point location with 
// uniform distribution on the whole desoprtion surface !not anymmore! changed to outgassing). The ray 
// direction is choosen according to the desorption type.
// -------------------------------------------------------------

BOOL StartFromSource() {
	BOOL found_side1 = FALSE;
	BOOL found_side2 = FALSE;
	BOOL found;
	BOOL foundInMap=FALSE;
	int mapPositionW,mapPositionH;
	FACET *src=NULL;
	double srcRnd;
	double A = 0.0;
	int i=0,j=0,w,h;
	int nbTry = 0;

	// Check end of simulation
	if( sHandle->maxDesorption>0 ) {
		if( sHandle->nbDesorbed>=sHandle->maxDesorption ) {
			sHandle->lastHit=NULL;
			return FALSE;
		}
	}

	// Select source
	srcRnd = rnd() * sHandle->sourceArea;

	while(!found_side1 && !found_side2 && j<sHandle->nbSuper) {
		i = 0;
		while(!found_side1 && !found_side2 && i<sHandle->str[j].nbFacet) {
			FACET *f = sHandle->str[j].facets[i];
			if (f->sh.useOutgassingFile) { //Using SynRad-generated outgassing map
				for (w=0;w<f->sh.outgassingMapWidth && !found_side1 && !found_side2;w++)
					for (h=0;h<f->sh.outgassingMapHeight && !found_side1 && !found_side2;h++) {
						double flow=f->outgassingMap[h*f->sh.outgassingMapWidth+w];
						if (flow>0.0) {
							foundInMap = found_side1 = (srcRnd>=A) && (srcRnd<(A+flow*((f->sh.is2sided)?0.5:1.0))); //2-sided facets have half outgassing on each side
							if (foundInMap) {mapPositionW=w;mapPositionH=h;}
							A += flow*((f->sh.is2sided)?0.5:1.0);
							if (f->sh.is2sided) { //check the other side
								foundInMap = found_side2 = (srcRnd>=A) && (srcRnd<A+flow*0.5);
								if (foundInMap) {mapPositionW=w;mapPositionH=h;}
								A += flow*0.5;
							}
						}
					}
			} else if( f->sh.desorbType!=DES_NONE ) {
				found_side1 = (srcRnd>=A) && (srcRnd<(A+f->sh.flow*((f->sh.is2sided)?0.5:1.0))); //2-sided facets have half outgassing on each side
				A += f->sh.flow*((f->sh.is2sided)?0.5:1.0);
				if (f->sh.is2sided) { //check the other side
					found_side2 = (srcRnd>=A) && (srcRnd<A+f->sh.flow*0.5);
					A += f->sh.flow*0.5;
				}
			}
			if(!found_side1 && !found_side2) i++;
		}
		if(!found_side1 && !found_side2) j++;
	}
	if( !found_side1 && !found_side2) {
		SetErrorSub("No starting point, aborting");
		return FALSE;
	}
	src = sHandle->str[j].facets[i];

	sHandle->lastHit = src;
	sHandle->distTraveledCurrentParticle=0.0;  //for mean free path calculations
	sHandle->flightTimeCurrentParticle=sHandle->desorptionStartTime+(sHandle->desorptionStopTime-sHandle->desorptionStartTime)*rnd();

	if (sHandle->useMaxwellDistribution) sHandle->velocityCurrentParticle=GenerateRandomVelocity(src->CDFid);
	
	else sHandle->velocityCurrentParticle=145.469*sqrt(src->sh.temperature/sHandle->gasMass);  //sqrt(8*R/PI/1000)=145.47
	sHandle->temperature = src->sh.temperature; //Thermalize particle

	found = FALSE;


	// Choose a starting point
	while(!found && nbTry<1000) {
		double u,v;

		if (foundInMap) {
			double uLength=sqrt(pow(src->sh.U.x,2)+pow(src->sh.U.y,2)+pow(src->sh.U.z,2));
			double vLength=sqrt(pow(src->sh.V.x,2)+pow(src->sh.V.y,2)+pow(src->sh.V.z,2));
			u=((double)mapPositionW+rnd())/src->sh.outgassingFileRatio/uLength;
			v=((double)mapPositionH+rnd())/src->sh.outgassingFileRatio/vLength;
		} else {
			u = rnd();
			v = rnd();
		}
		if( IsInFacet(src,u,v) ) {

			// (U,V) -> (x,y,z)
			sHandle->pPos.x = src->sh.O.x + u*src->sh.U.x + v*src->sh.V.x;
			sHandle->pPos.y = src->sh.O.y + u*src->sh.U.y + v*src->sh.V.y;
			sHandle->pPos.z = src->sh.O.z + u*src->sh.U.z + v*src->sh.V.z;
			src->colU = u;
			src->colV = v;
			found = TRUE;

		}
		nbTry++;
	}

	if( !found ) {
		// Get the center, if the center is not included in the facet, a leak is generated.
		sHandle->pPos = sHandle->str[j].facets[i]->sh.center;
	}

	//See docs/theta_gen.png for further details on angular distribution generation
	switch( src->sh.desorbType ) {
	case DES_UNIFORM:
		PolarToCartesian(src,acos(rnd()),rnd()*2.0*PI,found_side2);
		break;
	case DES_NONE: //for file-based
	case DES_COSINE:
		PolarToCartesian(src,acos(sqrt(rnd())),rnd()*2.0*PI,found_side2);
		break;
		/*case DES_COSINE2:
		PolarToCartesian(src,acos(pow(rnd(),1.0/3.0)),rnd()*2.0*PI,found_side2);
		break;
		case DES_COSINE3:
		PolarToCartesian(src,acos(pow(rnd(),1.0/4.0)),rnd()*2.0*PI,found_side2);
		break;
		case DES_COSINE4:
		PolarToCartesian(src,acos(pow(rnd(),1.0/5.0)),rnd()*2.0*PI,found_side2);
		break;*/
	case DES_COSINE_N:
		PolarToCartesian(src,acos(pow(rnd(),1.0/(src->sh.desorbTypeN+1.0))),rnd()*2.0*PI,found_side2);
		break;
	}

	// Current structure
	sHandle->curStruct = src->sh.superIdx;

	// Count
	src->sh.counter.hit.nbDesorbed++;
	src->hitted = TRUE;
	sHandle->nbDesorbed++;
	sHandle->tmpCount.hit.nbDesorbed++;
	sHandle->nbPHit = 0;
	if( src->hits && src->sh.countDes ) AHIT_FACET(src,sHandle->flightTimeCurrentParticle);
	if( src->direction && src->sh.countDirection ) DHIT_FACET(src,sHandle->flightTimeCurrentParticle);

	// Reset volatile state
	if( sHandle->hasVolatile ) {
		for(j=0;j<sHandle->nbSuper;j++) {
			for(i=0;i<sHandle->str[j].nbFacet;i++) {
				sHandle->str[j].facets[i]->ready = TRUE;
			}
		}
	}

	found_side1=FALSE;
	found_side2=FALSE;
	return TRUE;

}

// -------------------------------------------------------------
// Compute bounce against a facet
// -------------------------------------------------------------

void PerformBounce(FACET *iFacet) {

	double inPhi,inTheta;
	BOOL revert=FALSE;

	// Handle super structure link facet
	if( iFacet->sh.superDest ) {

		/*if (iFacet->sh.superDest==TESTCUBE_STRUCTURE_NO) { //entering test cube
			sHandle->testCubeEnterMoment=sHandle->flightTimeCurrentParticle;
			sHandle->testCubeEnterDist=sHandle->distTraveledCurrentParticle;
			sHandle->testCubeTemp+=sHandle->temperature;
			sHandle->testCubeVelocity+=sHandle->velocityCurrentParticle;
			sHandle->testCubeCount++;
		}

		if (sHandle->curStruct==(TESTCUBE_STRUCTURE_NO-1) && iFacet->sh.superDest!=TESTCUBE_STRUCTURE_NO) { //exiting test cube
			sHandle->testCubeTime+=sHandle->flightTimeCurrentParticle-sHandle->testCubeEnterMoment;
			sHandle->testCubeDist+=sHandle->distTraveledCurrentParticle-sHandle->testCubeEnterDist;
		}*/

		sHandle->curStruct = iFacet->sh.superDest - 1;
		// Count this hit as a transparent pass
		RecordHit(HIT_TRANS);
		if( iFacet->hits && iFacet->sh.countTrans ) AHIT_FACET(iFacet,sHandle->flightTimeCurrentParticle);
		if( iFacet->direction && iFacet->sh.countDirection ) DHIT_FACET(iFacet,sHandle->flightTimeCurrentParticle);
		return;

	}

	// Handle volatile facet
	if( iFacet->sh.isVolatile ) {

		if( iFacet->ready ) {
			iFacet->sh.counter.hit.nbAbsorbed++;
			iFacet->ready = FALSE;
			if( iFacet->hits && iFacet->sh.countAbs ) AHIT_FACET(iFacet,sHandle->flightTimeCurrentParticle);
			if( iFacet->direction && iFacet->sh.countDirection ) DHIT_FACET(iFacet,sHandle->flightTimeCurrentParticle);
		}
		return;

	}

	if( iFacet->sh.is2sided ) {
		// We may need to revert normal in case of 2 sided hit
		revert = DOT3(sHandle->pDir.x,sHandle->pDir.y,sHandle->pDir.z,
			iFacet->sh.N.x,iFacet->sh.N.y,iFacet->sh.N.z)>0.0;
	}

	if( iFacet->hits && iFacet->sh.countRefl ) AHIT_FACET(iFacet,sHandle->flightTimeCurrentParticle);
	if( iFacet->direction && iFacet->sh.countDirection ) DHIT_FACET(iFacet,sHandle->flightTimeCurrentParticle);

	// Relaunch particle
	UpdateVelocity(iFacet);
	sHandle->temperature = iFacet->sh.temperature; //Thermalize particle

	switch( iFacet->sh.reflectType ) {
	case REF_MIRROR:
		CartesianToPolar(iFacet,&inTheta,&inPhi);
		PolarToCartesian(iFacet,PI-inTheta,inPhi,FALSE);
		break;
	case REF_DIFFUSE:
		//See docs/theta_gen.png for further details on angular distribution generation
		PolarToCartesian(iFacet,acos(sqrt(rnd())),rnd()*2.0*PI,FALSE);
		break;
	}

	if( revert ) {
		sHandle->pDir.x = -sHandle->pDir.x;
		sHandle->pDir.y = -sHandle->pDir.y;
		sHandle->pDir.z = -sHandle->pDir.z;
	}

	RecordHit(HIT_REF);
	sHandle->lastHit = iFacet;
	sHandle->nbPHit++;


}

void AHIT_FACET(FACET *_f,double time) {                            
	int tu = (int)((_f)->colU * (_f)->sh.texWidthD);  
	int tv = (int)((_f)->colV * (_f)->sh.texHeightD); 
	int add = tu+tv*((_f)->sh.texWidth);              
	float tempfactor=(sHandle->nonIsothermal)?(float)(1.0/sqrt(sHandle->temperature)):1.0f; //check

	for (int m=0;m<=sHandle->nbMoments;m++)
		if (m==0 || abs(time-sHandle->moments[m-1])<sHandle->timeWindowSize/2.0) {
			//(_f)->hits[m][add]+=tempfactor*(_f)->inc[add]*145.469*sqrt(sHandle->temperature/sHandle->gasMass)/sHandle->velocityCurrentParticle;
			(_f)->hits[m][add]+=tempfactor*(_f)->inc[add];
		}
}

void DHIT_FACET(FACET *f,double time) {                            
	int tu = (int)(f->colU * f->sh.texWidthD);  
	int tv = (int)(f->colV * f->sh.texHeightD); 
	int add = tu+tv*(f->sh.texWidth);              

	for (size_t m=0;m<=sHandle->moments.size();m++) {
		if (m==0 || abs(time-sHandle->moments[m-1])<sHandle->timeWindowSize/2.0) {
			f->direction[m][add].dir.x += sHandle->pDir.x*sHandle->velocityCurrentParticle;
			f->direction[m][add].dir.y += sHandle->pDir.y*sHandle->velocityCurrentParticle;
			f->direction[m][add].dir.z += sHandle->pDir.z*sHandle->velocityCurrentParticle;
			f->direction[m][add].sumSpeed += sHandle->velocityCurrentParticle;
			f->direction[m][add].count++;

		}
	}

}

void UpdateVelocity(FACET *collidedFacet) {
	//thermalize perfectly
	if (sHandle->useMaxwellDistribution) sHandle->velocityCurrentParticle=GenerateRandomVelocity(collidedFacet->CDFid);
	else sHandle->velocityCurrentParticle=145.469*sqrt(collidedFacet->sh.temperature/sHandle->gasMass);
	/*//DEBUG!!!!
	std::vector<llong> count(301,0);
	for (int i=0;i<5000000;i++) {
		size_t index=(int)(GenerateRandomVelocity(0)*0.1);
		if (index>300) {
			index=300;
		}
		count[index]++;
	}*/
}

double GenerateRandomVelocity(int CDFId){
	//return InterpolateY(rnd(),sHandle->CDFs[CDFId],TRUE);
	return FastLookupY(rnd(),sHandle->CDFs[CDFId],FALSE);
}