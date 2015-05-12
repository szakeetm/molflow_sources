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
#include "Utils.h"

extern SIMULATION *sHandle;
extern void SetErrorSub(char *message);

// -------------------------------------------------------------
// Compute area of all the desorption facet
// -------------------------------------------------------------

void CalcTotalOutgassing() { 
	int i, j, k, tSize;
	FACET *f;
	double scale;
	//float scale_precomputed;

	// Update texture increment for MC
	//scale_precomputed=(float)(40.0/(sqrt(8.0*8.31/(PI*sHandle->gasMass*0.001))));
	for (j = 0; j < sHandle->nbSuper; j++) {
		for (i = 0; i < sHandle->str[j].nbFacet; i++) {
			f = sHandle->str[j].facets[i];
			if (f->inc) {
				tSize = f->sh.texWidth*f->sh.texHeight;
				scale = 1.0;
				if (f->sh.is2sided) scale = scale * 0.5;
				for (k = 0; k < tSize; k++) f->inc[k] *= scale;
				f->fullSizeInc *= scale;
			}
		}
	}

}





// -------------------------------------------------------

void PolarToCartesian(FACET *iFacet, double theta, double phi, BOOL reverse) {

	VERTEX3D U, V, N;
	double u, v, n;

	// Polar in (nU,nV,N) to Cartesian(x,y,z) transformation  ( nU = U/|U| , nV = V/|V| )
	// tetha is the angle to the normal of the facet N, phi to U
	// ! See Geometry::InitializeGeometry() for further informations on the (U,V,N) basis !
	// (nU,nV,N) and (x,y,z) are both left handed

	//This should be a speed-up routine, but I didn't experience any speed difference so I commented it out. Marton
	/*#ifdef WIN64
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
		N.x = N.x*(-1.0);
		N.y = N.y*(-1.0);
		N.z = N.z*(-1.0);
	}

	// Basis change (nU,nV,N) -> (x,y,z)
	sHandle->pDir.x = u*U.x + v*V.x + n*N.x;
	sHandle->pDir.y = u*U.y + v*V.y + n*N.y;
	sHandle->pDir.z = u*U.z + v*V.z + n*N.z;

}

// -------------------------------------------------------

void CartesianToPolar(FACET *iFacet, double *theta, double *phi) {

	// Get polar coordinates of the incoming particule direction in the (U,V,N) facet space.
	// Note: The facet is parallel to (U,V), we use its (nU,nV,N) orthonormal basis here.
	// (nU,nV,N) and (x,y,z) are both left handed

	// Cartesian(x,y,z) to polar in (nU,nV,N) transformation

	// Basis change (x,y,z) -> (nU,nV,N)
	// We use the fact that (nU,nV,N) belongs to SO(3)
	double u = DOT3(sHandle->pDir.x, sHandle->pDir.y, sHandle->pDir.z,
		iFacet->sh.nU.x, iFacet->sh.nU.y, iFacet->sh.nU.z);
	double v = DOT3(sHandle->pDir.x, sHandle->pDir.y, sHandle->pDir.z,
		iFacet->sh.nV.x, iFacet->sh.nV.y, iFacet->sh.nV.z);
	double n = DOT3(sHandle->pDir.x, sHandle->pDir.y, sHandle->pDir.z,
		iFacet->sh.N.x, iFacet->sh.N.y, iFacet->sh.N.z);

	// (u,v,n) -> (theta,phi)
	double rho = sqrt(v*v + u*u);
	*theta = acos(n);              // Angle to normal (PI/2 => PI)
	*phi = asin(v / rho);
	if (u < 0.0) *phi = PI - *phi;  // Angle to U

}

// -------------------------------------------------------

void UpdateMCHits(Dataport *dpHit, int prIdx, size_t nbMoments, DWORD timeout) {

	BYTE *buffer;
	SHGHITS *gHits;
	TEXTURE_MIN_MAX texture_limits_old[3];
	int i, j, s, x, y;
	int nb;
#ifdef _DEBUG
	double t0,t1;
	t0 = GetTick();
#endif
	sHandle->lastUpdateOK = AccessDataportTimed(dpHit, timeout);
	if (!sHandle->lastUpdateOK) return;

	buffer = (BYTE*)dpHit->buff;
	gHits = (SHGHITS *)buffer;

	// Global hits and leaks: adding local hits to shared memory
	gHits->total.hit.nbHit += sHandle->tmpCount.hit.nbHit;
	gHits->total.hit.nbAbsorbed += sHandle->tmpCount.hit.nbAbsorbed;
	gHits->total.hit.nbDesorbed += sHandle->tmpCount.hit.nbDesorbed;
	gHits->distTraveledTotal += sHandle->distTraveledSinceUpdate;
	
	//Memorize current limits, then do a min/max search
	for (i = 0; i < 3; i++) {
		texture_limits_old[i] = gHits->texture_limits[i];
		gHits->texture_limits[i].min.all = gHits->texture_limits[i].min.moments_only = HITMAX;
		gHits->texture_limits[i].max.all = gHits->texture_limits[i].max.moments_only = 0;
	}

	gHits->mode = MC_MODE;
	//for(i=0;i<BOUNCEMAX;i++) gHits->wallHits[i] += sHandle->wallHits[i];

	// Leak
	nb = gHits->nbLastLeaks;
	for (i = 0; i < sHandle->nbLastLeak && i < NBHLEAK; i++)
		gHits->pLeak[(i + nb) % NBHLEAK] = sHandle->pLeak[i];
	gHits->nbLeakTotal += sHandle->nbLeakTotal;
	gHits->nbLastLeaks += sHandle->nbLastLeak;

	// HHit (Only prIdx 0)
	if (prIdx == 0) {
		gHits->nbHHit = sHandle->nbHHit;
		memcpy(gHits->pHits, sHandle->pHits, NBHHIT*sizeof(HIT));
	}

	// Facets
	for (s = 0; s < sHandle->nbSuper; s++) {
		for (i = 0; i < sHandle->str[s].nbFacet; i++) {

			FACET *f = sHandle->str[s].facets[i];
			if (f->hitted) {

				SHHITS *fFit = (SHHITS *)(buffer + f->sh.hitOffset);
				fFit->hit.nbAbsorbed += f->sh.counter.hit.nbAbsorbed;
				//printf("\n%d %g + %g",i,(double)fFit->hit.nbHit,(double)f->sh.counter.hit.nbHit);
				fFit->hit.nbDesorbed += f->sh.counter.hit.nbDesorbed;
				fFit->hit.nbHit += f->sh.counter.hit.nbHit;
				fFit->hit.sum_1_per_ort_velocity += f->sh.counter.hit.sum_1_per_ort_velocity;
				fFit->hit.sum_v_ort += f->sh.counter.hit.sum_v_ort;

				if (f->sh.isProfile) {
					for (int m = 0; m < (1 + nbMoments); m++) {
						APROFILE *shProfile = (APROFILE *)(buffer + (f->sh.hitOffset + sizeof(SHHITS)+m*f->profileSize));
						for (j = 0; j < PROFILE_SIZE; j++) {
							shProfile[j].count += f->profile[m][j].count;
							shProfile[j].sum_1_per_ort_velocity += f->profile[m][j].sum_1_per_ort_velocity;
							shProfile[j].sum_v_ort += f->profile[m][j].sum_v_ort;
						}
					}
				}

				if (f->sh.isTextured) {
					for (int m = 0; m < (1 + nbMoments); m++) {
						AHIT *shTexture = (AHIT *)(buffer + (f->sh.hitOffset + sizeof(SHHITS)+f->profileSize*(1 + nbMoments) + m*f->textureSize));
						//double dCoef = gHits->total.hit.nbDesorbed * 1E4 * sHandle->gasMass / 1000 / 6E23 * MAGIC_CORRECTION_FACTOR;  //1E4 is conversion from m2 to cm2
						double timeCorrection = m == 0 ? sHandle->finalOutgassingRate : (sHandle->totalDesorbedMolecules) / sHandle->timeWindowSize;
						
						for (y = 0; y < f->sh.texHeight; y++) {
							for (x = 0; x < f->sh.texWidth; x++) {
								int add = x + y*f->sh.texWidth;

								//Add temporary hit counts
								shTexture[add].count += f->hits[m][add].count;
								shTexture[add].sum_1_per_ort_velocity += f->hits[m][add].sum_1_per_ort_velocity;
								shTexture[add].sum_v_ort_per_area += f->hits[m][add].sum_v_ort_per_area;

								double val[3];  //pre-calculated autoscaling values (Pressure, imp.rate, density)

								val[0] = shTexture[add].sum_v_ort_per_area*timeCorrection; //pressure without dCoef_pressure
								val[1] = shTexture[add].count*f->inc[add] * timeCorrection; //imp.rate without dCoef
								val[2] = f->inc[add] * shTexture[add].sum_1_per_ort_velocity* timeCorrection; //particle density without dCoef

								//Global autoscale
								for (int v = 0; v<3; v++) {
									if (val[v]>gHits->texture_limits[v].max.all && f->largeEnough[add])
										gHits->texture_limits[v].max.all = val[v];

									if (val[v] > 0.0 && val[v]<gHits->texture_limits[v].min.all && f->largeEnough[add])
										gHits->texture_limits[v].min.all = val[v];

									//Autoscale ignoring constant flow (moments only)
									if (m != 0) {
										if (val[v]>gHits->texture_limits[v].max.moments_only && f->largeEnough[add])
											gHits->texture_limits[v].max.moments_only = val[v];

										if (val[v] > 0.0 && val[v] < gHits->texture_limits[v].min.moments_only && f->largeEnough[add])
											gHits->texture_limits[v].min.moments_only = val[v];
									}
								}

							}
						}
					}
				}


				if (f->sh.countDirection) {
					for (int m = 0; m < (1 + nbMoments); m++) {
						VHIT *shDir = (VHIT *)(buffer + (f->sh.hitOffset + sizeof(SHHITS)+f->profileSize*(1 + nbMoments) + f->textureSize*(1 + nbMoments) + f->directionSize*m));
						for (y = 0; y < f->sh.texHeight; y++) {
							for (x = 0; x < f->sh.texWidth; x++) {
								int add = x + y*f->sh.texWidth;
								shDir[add].sumDir.x += f->direction[m][add].sumDir.x;
								shDir[add].sumDir.y += f->direction[m][add].sumDir.y;
								shDir[add].sumDir.z += f->direction[m][add].sumDir.z;
								//shDir[add].sumSpeed += f->direction[m][add].sumSpeed;
								shDir[add].count += f->direction[m][add].count;
							}
						}
					}
				}

			} // End if(hitted)
		} // End nbFacet
	} // End nbSuper

	//if there were no textures:
	for (int v = 0; v < 3; v++) {
		if (gHits->texture_limits[v].min.all == HITMAX) gHits->texture_limits[v].min.all = texture_limits_old[v].min.all;
		if (gHits->texture_limits[v].min.moments_only == HITMAX) gHits->texture_limits[v].min.moments_only = texture_limits_old[v].min.moments_only;
		if (gHits->texture_limits[v].max.all == 0.0) gHits->texture_limits[v].max.all = texture_limits_old[v].max.all;
		if (gHits->texture_limits[v].max.moments_only == 0.0) gHits->texture_limits[v].max.moments_only = texture_limits_old[v].max.moments_only;
	}

	ReleaseDataport(dpHit);
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

	double inPhi, inTheta;
	//Search destination
	FACET *destination;
	BOOL found = FALSE;
	BOOL revert = FALSE;
	int i, j;
	//Look in which superstructure is the destination facet:
	for (i = 0; i < sHandle->nbSuper && (!found); i++) {
		for (j = 0; j < sHandle->str[i].nbFacet && (!found); j++) {
			if ((iFacet->sh.teleportDest - 1) == sHandle->str[i].facets[j]->globalId) {
				destination = sHandle->str[i].facets[j];
				sHandle->curStruct = destination->sh.superIdx; //change current superstructure
				found = TRUE;
			}
		}
	}
	if (!found) {
		char err[128];
		sprintf(err, "Teleport destination of facet %d not found (facet %d does not exist)", iFacet->globalId + 1, iFacet->sh.teleportDest);
		SetErrorSub(err);
		return;
	}
	// Count this hit as a transparent pass
	RecordHit(HIT_TELEPORT);
	if (iFacet->hits && iFacet->sh.countTrans) AHIT_FACET(iFacet, sHandle->flightTimeCurrentParticle, TRUE, 2.0, 2.0);
	if (iFacet->direction && iFacet->sh.countDirection) DHIT_FACET(iFacet, sHandle->flightTimeCurrentParticle);
	ProfileFacet(iFacet, sHandle->flightTimeCurrentParticle, TRUE, 2.0, 2.0);

	// Relaunch particle from new facet
	CartesianToPolar(iFacet, &inTheta, &inPhi);
	PolarToCartesian(destination, inTheta, inPhi, FALSE);
	// Move particle to teleport destination point
	double u = iFacet->colU;
	double v = iFacet->colV;
	sHandle->pPos.x = destination->sh.O.x + u*destination->sh.U.x + v*destination->sh.V.x;
	sHandle->pPos.y = destination->sh.O.y + u*destination->sh.U.y + v*destination->sh.V.y;
	sHandle->pPos.z = destination->sh.O.z + u*destination->sh.U.z + v*destination->sh.V.z;
	RecordHit(HIT_TELEPORT);
	int nbTry = 0;
	if (!IsInFacet(destination, u, v)) { //source and destination facets not the same shape, would generate leak
		// Choose a new starting point
		RecordHit(HIT_ABS);
		BOOL found = FALSE;
		while (!found && nbTry < 1000) {
			u = rnd();
			v = rnd();
			if (IsInFacet(destination, u, v)) {
				found = TRUE;
				sHandle->pPos.x = destination->sh.O.x + u*destination->sh.U.x + v*destination->sh.V.x;
				sHandle->pPos.y = destination->sh.O.y + u*destination->sh.U.y + v*destination->sh.V.y;
				sHandle->pPos.z = destination->sh.O.z + u*destination->sh.U.z + v*destination->sh.V.z;
				RecordHit(HIT_DES);
			}
		}
		nbTry++;
	}
	
	sHandle->lastHit = destination;

	//Count hits on teleport facets
	/*iFacet->sh.counter.hit.nbAbsorbed++;
	destination->sh.counter.hit.nbDesorbed++;*/

	double ortVelocity = sHandle->velocityCurrentParticle*abs(DOT3(
		sHandle->pDir.x, sHandle->pDir.y, sHandle->pDir.z,
		iFacet->sh.N.x, iFacet->sh.N.y, iFacet->sh.N.z));
	iFacet->sh.counter.hit.nbHit++;
	iFacet->sh.counter.hit.sum_1_per_ort_velocity += 2.0 / ortVelocity;
	iFacet->sh.counter.hit.sum_v_ort += 2.0*(sHandle->useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity;
	iFacet->hitted = TRUE;
	/*destination->sh.counter.hit.sum_1_per_ort_velocity += 2.0 / sHandle->velocityCurrentParticle;
	destination->sh.counter.hit.sum_v_ort += sHandle->velocityCurrentParticle*abs(DOT3(
	sHandle->pDir.x, sHandle->pDir.y, sHandle->pDir.z,
	destination->sh.N.x, destination->sh.N.y, destination->sh.N.z));*/
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
	for (i = 0; i < nbStep; i++) {

		found = Intersect(&(sHandle->pPos), &(sHandle->pDir), &d, &collidedFacet, sHandle->lastHit);

		if (found) {

			// Move particle to intersection point
			sHandle->pPos.x += d*sHandle->pDir.x;
			sHandle->pPos.y += d*sHandle->pDir.y;
			sHandle->pPos.z += d*sHandle->pDir.z;
			sHandle->distTraveledCurrentParticle += d;
			sHandle->flightTimeCurrentParticle += d / 100.0 / sHandle->velocityCurrentParticle; //conversion from cm to m

			if ((!sHandle->calcConstantFlow && sHandle->flightTimeCurrentParticle > sHandle->latestMoment)
				|| sHandle->particleDecayMoment<sHandle->flightTimeCurrentParticle) {
				//hit time over the measured period - we create a new particle
				//OR particle has decayed
				RecordHit(LASTHIT);
				sHandle->distTraveledSinceUpdate += sHandle->distTraveledCurrentParticle;
				if (!StartFromSource())
					// maxDesorption reached
					return FALSE;
			}
			else { //hit within measured time, particle still alive
				if (collidedFacet->sh.teleportDest) {
					PerformTeleport(collidedFacet);
				}
				else if (GetStickingAt(collidedFacet,sHandle->flightTimeCurrentParticle) == 1.0 || ((GetStickingAt(collidedFacet,sHandle->flightTimeCurrentParticle) > 0.0) && (rnd() < (GetStickingAt(collidedFacet,sHandle->flightTimeCurrentParticle))))) {
					//absorbed
					PerformAbsorb(collidedFacet);
					sHandle->distTraveledSinceUpdate += sHandle->distTraveledCurrentParticle;
					if (!StartFromSource())
						// maxDesorption reached
						return FALSE;
				}
				else {
					//reflected
					PerformBounce(collidedFacet);
				}
			} //end hit within measured time
		} //end intersection found
		else {
			// Leak (simulation error)
			RecordLeakPos();
			sHandle->nbLeakTotal++;
			if (!StartFromSource())
				// maxDesorption reached
				return FALSE;
		}
	}
	return TRUE;
}

// -------------------------------------------------------------
// Launch a ray from a source facet. The ray 
// direction is choosen according to the desorption type.
// -------------------------------------------------------------

BOOL StartFromSource() {
	BOOL found_side1 = FALSE;
	BOOL found_side2 = FALSE;
	BOOL found;
	BOOL foundInMap = FALSE;
	int mapPositionW, mapPositionH;
	FACET *src = NULL;
	double srcRnd;
	double A = 0.0;
	int i = 0, j = 0, w, h;
	int nbTry = 0;

	// Check end of simulation
	if (sHandle->maxDesorption > 0) {
		if (sHandle->nbDesorbed >= sHandle->maxDesorption) {
			sHandle->lastHit = NULL;
			return FALSE;
		}
	}

	// Select source
	srcRnd = rnd() * sHandle->totalDesorbedMolecules;

	while (!found_side1 && !found_side2 && j < sHandle->nbSuper) {
		i = 0;
		while (!found_side1 && !found_side2 && i < sHandle->str[j].nbFacet) {
			FACET *f = sHandle->str[j].facets[i];
			if (f->sh.desorbType != DES_NONE) { //there is some kind of outgassing
				if (f->sh.useOutgassingFile) { //Using SynRad-generated outgassing map
					for (w = 0; w < f->sh.outgassingMapWidth && !found_side1 && !found_side2; w++)
						for (h = 0; h<f->sh.outgassingMapHeight && !found_side1 && !found_side2; h++) {
							double outgassing = sHandle->latestMoment * f->outgassingMap[h*f->sh.outgassingMapWidth + w] / (1.38E-23*f->sh.temperature);
							if (outgassing>0.0) {
								foundInMap = found_side1 = (srcRnd >= A) && (srcRnd < (A + outgassing*((f->sh.is2sided) ? 0.5 : 1.0))); //2-sided facets have half outgassing on each side
								if (foundInMap) { mapPositionW = w; mapPositionH = h; }
								A += outgassing*((f->sh.is2sided) ? 0.5 : 1.0);
								if (f->sh.is2sided) { //check the other side
									foundInMap = found_side2 = (srcRnd >= A) && (srcRnd < A + outgassing*0.5);
									if (foundInMap) { mapPositionW = w; mapPositionH = h; }
									A += outgassing*0.5;
								}
							}
						}
				}
				else  { //constant or time-dependent outgassing
					double outgassing=
						(f->sh.outgassing_paramId>=0)
						?sHandle->IDs[f->sh.IDid].back().second/ (1.38E-23*f->sh.temperature)
						:sHandle->latestMoment*f->sh.flow / (1.38E-23*f->sh.temperature);
					found_side1 = (srcRnd >= A) && (srcRnd < (A + outgassing *((f->sh.is2sided) ? 0.5 : 1.0))); //2-sided facets have half outgassing on each side
					A += outgassing *((f->sh.is2sided) ? 0.5 : 1.0);
					if (f->sh.is2sided) { //check the other side
						found_side2 = (srcRnd >= A) && (srcRnd < A + outgassing*0.5);
						A += outgassing*0.5;
					}
				}
			}
			if (!found_side1 && !found_side2) i++;
		}
		if (!found_side1 && !found_side2) j++;
	}
	if (!found_side1 && !found_side2) {
		SetErrorSub("No starting point, aborting");
		return FALSE;
	}
	src = sHandle->str[j].facets[i];

	sHandle->lastHit = src;
	sHandle->distTraveledCurrentParticle = 0.0;  //for mean free path calculations
	//sHandle->flightTimeCurrentParticle = sHandle->desorptionStartTime + (sHandle->desorptionStopTime - sHandle->desorptionStartTime)*rnd();
	sHandle->flightTimeCurrentParticle=GenerateDesorptionTime(src);
	if (sHandle->useMaxwellDistribution) sHandle->velocityCurrentParticle = GenerateRandomVelocity(src->sh.CDFid);
	else sHandle->velocityCurrentParticle = 145.469*sqrt(src->sh.temperature / sHandle->gasMass);  //sqrt(8*R/PI/1000)=145.47
	if (sHandle->halfLife < 9e99) { //decaying gas
		sHandle->particleDecayMoment = sHandle->flightTimeCurrentParticle + sHandle->halfLife*1.44269*-log(rnd()); //1.44269=1/ln2
		//Exponential distribution PDF: probability of 't' life = 1/TAU*exp(-t/TAU) where TAU = half_life/ln2
		//Exponential distribution CDF: probability of life shorter than 't" = 1-exp(-t/TAU)
		//Equation: rnd()=1-exp(-t/TAU)
		//Solution: t=TAU*-log(1-rnd()) and 1-rnd()=rnd() therefore t=half_life/ln2*-log(rnd())
	}
	else {
		sHandle->particleDecayMoment = 1e100;
	}
	//sHandle->temperature = src->sh.temperature; //Thermalize particle

	found = FALSE;


	// Choose a starting point
	while (!found && nbTry < 1000) {
		double u, v;

		if (foundInMap) {
			double uLength = sqrt(pow(src->sh.U.x, 2) + pow(src->sh.U.y, 2) + pow(src->sh.U.z, 2));
			double vLength = sqrt(pow(src->sh.V.x, 2) + pow(src->sh.V.y, 2) + pow(src->sh.V.z, 2));
			u = ((double)mapPositionW + rnd()) / src->sh.outgassingFileRatio / uLength;
			v = ((double)mapPositionH + rnd()) / src->sh.outgassingFileRatio / vLength;
		}
		else {
			u = rnd();
			v = rnd();
		}
		if (IsInFacet(src, u, v)) {

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

	if (!found) {
		// Get the center, if the center is not included in the facet, a leak is generated.
		sHandle->pPos = sHandle->str[j].facets[i]->sh.center;
	}

	if (src->sh.isMoving && sHandle->motionType) RecordHit(HIT_MOVING);
	else RecordHit(HIT_DES); //create blue hit point for created particle

	//See docs/theta_gen.png for further details on angular distribution generation
	switch (src->sh.desorbType) {
	case DES_UNIFORM:
		PolarToCartesian(src, acos(rnd()), rnd()*2.0*PI, found_side2);
		break;
	case DES_NONE: //for file-based
	case DES_COSINE:
		PolarToCartesian(src, acos(sqrt(rnd())), rnd()*2.0*PI, found_side2);
		break;
	case DES_COSINE_N:
		PolarToCartesian(src, acos(pow(rnd(), 1.0 / (src->sh.desorbTypeN + 1.0))), rnd()*2.0*PI, found_side2);
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

	if (src->sh.isMoving) {
		TreatMovingFacet();
	}

	double ortVelocity = sHandle->velocityCurrentParticle*abs(DOT3(
		sHandle->pDir.x, sHandle->pDir.y, sHandle->pDir.z,
		src->sh.N.x, src->sh.N.y, src->sh.N.z));
	src->sh.counter.hit.sum_1_per_ort_velocity += 2.0 / ortVelocity;
	src->sh.counter.hit.sum_v_ort += (sHandle->useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity;
	//Desorption doesn't contribute to angular profiles
	ProfileFacet(src, sHandle->flightTimeCurrentParticle, FALSE, 2.0, 1.0);
	if (src->hits && src->sh.countDes) AHIT_FACET(src, sHandle->flightTimeCurrentParticle, TRUE, 2.0, 1.0);
	//if (src->direction && src->sh.countDirection) DHIT_FACET(src, sHandle->flightTimeCurrentParticle);

	// Reset volatile state
	if (sHandle->hasVolatile) {
		for (j = 0; j < sHandle->nbSuper; j++) {
			for (i = 0; i < sHandle->str[j].nbFacet; i++) {
				sHandle->str[j].facets[i]->ready = TRUE;
			}
		}
	}

	found_side1 = FALSE;
	found_side2 = FALSE;
	return TRUE;
}

// -------------------------------------------------------------
// Compute bounce against a facet
// -------------------------------------------------------------

void PerformBounce(FACET *iFacet) {

	double inPhi, inTheta;
	BOOL revert = FALSE;

	// Handle super structure link facet
	if (iFacet->sh.superDest) {
		sHandle->curStruct = iFacet->sh.superDest - 1;
		// Count this hit as a transparent pass
		RecordHit(HIT_TRANS);
		ProfileFacet(iFacet, sHandle->flightTimeCurrentParticle, TRUE, 2.0, 2.0);
		if (iFacet->hits && iFacet->sh.countTrans) AHIT_FACET(iFacet, sHandle->flightTimeCurrentParticle, TRUE, 2.0, 2.0);
		if (iFacet->direction && iFacet->sh.countDirection) DHIT_FACET(iFacet, sHandle->flightTimeCurrentParticle);
		return;

	}

	// Handle volatile facet
	if (iFacet->sh.isVolatile) {

		if (iFacet->ready) {
			iFacet->sh.counter.hit.nbAbsorbed++;
			iFacet->ready = FALSE;
			ProfileFacet(iFacet, sHandle->flightTimeCurrentParticle, TRUE, 2.0, 1.0);
			if (iFacet->hits && iFacet->sh.countAbs) AHIT_FACET(iFacet, sHandle->flightTimeCurrentParticle, TRUE, 2.0, 1.0);
			if (iFacet->direction && iFacet->sh.countDirection) DHIT_FACET(iFacet, sHandle->flightTimeCurrentParticle);
		}
		return;

	}

	if (iFacet->sh.is2sided) {
		// We may need to revert normal in case of 2 sided hit
		revert = DOT3(sHandle->pDir.x, sHandle->pDir.y, sHandle->pDir.z,
			iFacet->sh.N.x, iFacet->sh.N.y, iFacet->sh.N.z) > 0.0;
	}

	//Texture/Profile incoming hit
	sHandle->tmpCount.hit.nbHit++; //global
	iFacet->sh.counter.hit.nbHit++; //hit facet
	//Register (orthogonal) velocity
	double ortVelocity = sHandle->velocityCurrentParticle*abs(DOT3(
		sHandle->pDir.x, sHandle->pDir.y, sHandle->pDir.z,
		iFacet->sh.N.x, iFacet->sh.N.y, iFacet->sh.N.z));
	iFacet->sh.counter.hit.sum_1_per_ort_velocity += 1.0 / ortVelocity;
	iFacet->sh.counter.hit.sum_v_ort += (sHandle->useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity;
	if (iFacet->hits && iFacet->sh.countRefl) AHIT_FACET(iFacet, sHandle->flightTimeCurrentParticle, TRUE, 1.0, 1.0);
	if (iFacet->direction && iFacet->sh.countDirection) DHIT_FACET(iFacet, sHandle->flightTimeCurrentParticle);
	ProfileFacet(iFacet, sHandle->flightTimeCurrentParticle, TRUE, 1.0, 1.0);

	// Relaunch particle
	UpdateVelocity(iFacet);
	//sHandle->temperature = iFacet->sh.temperature; //Thermalize particle


	switch (iFacet->sh.reflectType) {
	case REF_MIRROR:
		CartesianToPolar(iFacet, &inTheta, &inPhi);
		PolarToCartesian(iFacet, PI - inTheta, inPhi, FALSE);
		break;
	case REF_DIFFUSE:
		//See docs/theta_gen.png for further details on angular distribution generation
		PolarToCartesian(iFacet, acos(sqrt(rnd())), rnd()*2.0*PI, FALSE);
		break;
	case REF_UNIFORM:
		PolarToCartesian(iFacet, acos(rnd()), rnd()*2.0*PI, FALSE);
		break;
	}

	if (revert) {
		sHandle->pDir.x = -sHandle->pDir.x;
		sHandle->pDir.y = -sHandle->pDir.y;
		sHandle->pDir.z = -sHandle->pDir.z;
	}

	if (iFacet->sh.isMoving) {
		TreatMovingFacet();
	}

	//Texture/Profile outgoing particle
	//Register outgoing velocity
	ortVelocity = sHandle->velocityCurrentParticle*abs(DOT3(
		sHandle->pDir.x, sHandle->pDir.y, sHandle->pDir.z,
		iFacet->sh.N.x, iFacet->sh.N.y, iFacet->sh.N.z));
	iFacet->sh.counter.hit.sum_1_per_ort_velocity += 1.0 / ortVelocity;
	iFacet->sh.counter.hit.sum_v_ort += (sHandle->useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity;
	if (iFacet->hits && iFacet->sh.countRefl) AHIT_FACET(iFacet, sHandle->flightTimeCurrentParticle, FALSE, 1.0, 1.0); //count again for outward velocity
	ProfileFacet(iFacet, sHandle->flightTimeCurrentParticle, FALSE, 1.0, 1.0);
	//no direction count on outgoing

	if (iFacet->sh.isMoving && sHandle->motionType) RecordHit(HIT_MOVING);
	else RecordHit(HIT_REF);
	sHandle->lastHit = iFacet;
	sHandle->nbPHit++;


}

void PerformAbsorb(FACET *iFacet) {
	sHandle->tmpCount.hit.nbHit++; //global	
	sHandle->tmpCount.hit.nbAbsorbed++;
	iFacet->sh.counter.hit.nbHit++;
	iFacet->sh.counter.hit.nbAbsorbed++;
	RecordHit(HIT_ABS);
	double ortVelocity = sHandle->velocityCurrentParticle*abs(DOT3(
		sHandle->pDir.x, sHandle->pDir.y, sHandle->pDir.z,
		iFacet->sh.N.x, iFacet->sh.N.y, iFacet->sh.N.z));
	iFacet->sh.counter.hit.sum_1_per_ort_velocity += 2.0 / ortVelocity;
	iFacet->sh.counter.hit.sum_v_ort += (sHandle->useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity;
	ProfileFacet(iFacet, sHandle->flightTimeCurrentParticle, TRUE, 2.0, 1.0);
	if (iFacet->hits && iFacet->sh.countAbs) AHIT_FACET(iFacet, sHandle->flightTimeCurrentParticle, TRUE, 2.0, 1.0);
	if (iFacet->direction && iFacet->sh.countDirection) DHIT_FACET(iFacet, sHandle->flightTimeCurrentParticle);
}

void AHIT_FACET(FACET *f, double time, BOOL countHit, double velocity_factor, double ortSpeedFactor) {

	int tu = (int)(f->colU * f->sh.texWidthD);
	int tv = (int)(f->colV * f->sh.texHeightD);
	int add = tu + tv*(f->sh.texWidth);
	double ortVelocity = (sHandle->useMaxwellDistribution ? 1.0 : 1.1781)*sHandle->velocityCurrentParticle*abs(DOT3(
		sHandle->pDir.x, sHandle->pDir.y, sHandle->pDir.z,
		f->sh.N.x, f->sh.N.y, f->sh.N.z)); //surface-orthogonal velocity component

	for (size_t m = 0; m <= sHandle->nbMoments; m++)
		if (m == 0 || abs(time - sHandle->moments[m - 1]) < sHandle->timeWindowSize / 2.0) {
			if (countHit) f->hits[m][add].count++;
			f->hits[m][add].sum_1_per_ort_velocity += velocity_factor / ortVelocity;
			f->hits[m][add].sum_v_ort_per_area += ortSpeedFactor*ortVelocity*f->inc[add]; // sum ortho_velocity[m/s] / cell_area[cm2]
		}
}

void DHIT_FACET(FACET *f, double time) {
	int tu = (int)(f->colU * f->sh.texWidthD);
	int tv = (int)(f->colV * f->sh.texHeightD);
	int add = tu + tv*(f->sh.texWidth);

	for (size_t m = 0; m <= sHandle->moments.size(); m++) {
		if (m == 0 || abs(time - sHandle->moments[m - 1]) < sHandle->timeWindowSize / 2.0) {
			f->direction[m][add].sumDir.x += sHandle->pDir.x*sHandle->velocityCurrentParticle;
			f->direction[m][add].sumDir.y += sHandle->pDir.y*sHandle->velocityCurrentParticle;
			f->direction[m][add].sumDir.z += sHandle->pDir.z*sHandle->velocityCurrentParticle;
			f->direction[m][add].count++;
		}
	}

}

void ProfileFacet(FACET *f, double time, BOOL countHit, double velocity_factor, double ortSpeedFactor) {

	int pos;
	int nbMoments = (int)sHandle->moments.size();
	double dot = 1.0;

	switch (f->sh.profileType) {

	case REC_ANGULAR: {
		if (countHit) {
			dot = DOT3(f->sh.N.x, f->sh.N.y, f->sh.N.z, sHandle->pDir.x, sHandle->pDir.y, sHandle->pDir.z);
			double theta = acos(abs(dot));              // Angle to normal (PI/2 => PI)
			int pos = (int)(theta / (PI / 2)*((double)PROFILE_SIZE)); // To Grad
			//int pos = (int)(cos(theta)*(double)PROFILE_SIZE); // COS(theta)
			SATURATE(pos, 0, PROFILE_SIZE - 1);
			for (int m = 0; m <= nbMoments; m++) {
				if (m == 0 || abs(time - sHandle->moments[m - 1]) < sHandle->timeWindowSize / 2.0) {
					f->profile[m][pos].count++;
				}
			}
		}
	}
		break;

	case REC_PRESSUREU:
	case REC_PRESSUREV:

		pos = (int)((f->sh.profileType == REC_PRESSUREU?f->colU:f->colV)*(double)PROFILE_SIZE);
		SATURATE(pos, 0, PROFILE_SIZE - 1);
		for (int m = 0; m <= nbMoments; m++) {
			if (m == 0 || abs(time - sHandle->moments[m - 1]) < sHandle->timeWindowSize / 2.0) {
				if (countHit) f->profile[m][pos].count++;
				double ortVelocity = sHandle->velocityCurrentParticle*abs(DOT3(f->sh.N.x, f->sh.N.y, f->sh.N.z,
					sHandle->pDir.x, sHandle->pDir.y, sHandle->pDir.z));
				f->profile[m][pos].sum_1_per_ort_velocity += velocity_factor / ortVelocity;
				f->profile[m][pos].sum_v_ort += ortSpeedFactor*(sHandle->useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity;
			}
		}
		break;

	case REC_ORT_VELOCITY:
		if (countHit) dot = abs(DOT3(f->sh.N.x, f->sh.N.y, f->sh.N.z, sHandle->pDir.x, sHandle->pDir.y, sHandle->pDir.z));  //cos(theta)
	case REC_VELOCITY:
		if (countHit) {
			pos = (int)(dot*sHandle->velocityCurrentParticle / f->sh.maxSpeed*(double)PROFILE_SIZE);
			SATURATE(pos, 0, PROFILE_SIZE - 1);
			for (int m = 0; m <= nbMoments; m++) {
				if (m == 0 || abs(time - sHandle->moments[m - 1]) < sHandle->timeWindowSize / 2.0) {
					f->profile[m][pos].count++;
				}
			}
		}
		break;
	}
}

void UpdateVelocity(FACET *collidedFacet) {
	if (collidedFacet->sh.accomodationFactor>0.9999) { //speedup for the most common case: perfect thermalization
		if (sHandle->useMaxwellDistribution) sHandle->velocityCurrentParticle = GenerateRandomVelocity(collidedFacet->sh.CDFid);
		else sHandle->velocityCurrentParticle = 145.469*sqrt(collidedFacet->sh.temperature / sHandle->gasMass);
	}
	else {
		double oldSpeed2 = pow(sHandle->velocityCurrentParticle, 2);
		double newSpeed2;
		if (sHandle->useMaxwellDistribution) newSpeed2 = pow(GenerateRandomVelocity(collidedFacet->sh.CDFid), 2);
		else newSpeed2 = /*145.469*/ 29369.939*(collidedFacet->sh.temperature / sHandle->gasMass);
		//sqrt(29369)=171.3766= sqrt(8*R*1000/PI)*3PI/8, that is, the constant part of the v_avg=sqrt(8RT/PI/m/0.001)) found in literature, multiplied by
		//the corrective factor of 3PI/8 that accounts for moving from volumetric speed distribution to wall collision speed distribution
		sHandle->velocityCurrentParticle = sqrt(oldSpeed2 + (newSpeed2 - oldSpeed2)*collidedFacet->sh.accomodationFactor);
	}
}

double GenerateRandomVelocity(int CDFId){
	//return FastLookupY(rnd(),sHandle->CDFs[CDFId],FALSE);
	return InterpolateX(rnd(), sHandle->CDFs[CDFId], TRUE);
}

double GenerateDesorptionTime(FACET *src){
	if (src->sh.outgassing_paramId>=0) { //time-dependent desorption
		return InterpolateX(rnd()*sHandle->IDs[src->sh.IDid].back().second,sHandle->IDs[src->sh.IDid],TRUE);
	} else {
		return rnd()*sHandle->latestMoment; //continous desorption between 0 and latestMoment
	}
}

double GetStickingAt(FACET *f,double time) {
	if (f->sh.sticking_paramId==-1) //constant sticking
		return f->sh.sticking;
	else return InterpolateY(time,sHandle->parameters[f->sh.sticking_paramId].values,TRUE);
}

double GetOpacityAt(FACET *f,double time) {
	if (f->sh.opacity_paramId==-1) //constant sticking
		return f->sh.opacity;
	else return InterpolateY(time,sHandle->parameters[f->sh.opacity_paramId].values,TRUE);
}

void TreatMovingFacet() {
	VERTEX3D localVelocityToAdd;
	if (sHandle->motionType == 1) {
		localVelocityToAdd = sHandle->motionVector2;
	}
	else if (sHandle->motionType == 2) {
		VERTEX3D distanceVector;
		Sub(&distanceVector, &sHandle->pPos, &sHandle->motionVector1); //distance from base
		ScalarMult(&distanceVector, 0.01); //cm->m
		Cross(&localVelocityToAdd, &sHandle->motionVector2, &distanceVector);
	}
	VERTEX3D oldVelocity, newVelocity;
	oldVelocity = sHandle->pDir;
	ScalarMult(&oldVelocity, sHandle->velocityCurrentParticle);
	Add(&newVelocity, &oldVelocity, &localVelocityToAdd);
	sHandle->pDir = newVelocity;
	Normalize(&sHandle->pDir);
	sHandle->velocityCurrentParticle = Norme(&newVelocity);
}