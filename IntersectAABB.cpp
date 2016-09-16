/*
File:        IntersectAABB.cpp
Description: Ray geometry intersection (Using AABB tree optimisation)
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
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include "Simulation.h"
#include "Random.h"

// -----------------------------------------------------------
// AABB tree stuff
// -----------------------------------------------------------

// Temporary for intersection
extern  double    intMinLgth;
extern  BOOL      intFound;
extern  VERTEX3D  intD;
extern  VERTEX3D  intZ;
extern  int       intNbTHits;
extern  double    iRx;
extern  double    iRy;
extern  double    iRz;
extern  BOOL      nullRx;
extern  BOOL      nullRy;
extern  BOOL      nullRz;
extern  VERTEX3D *rayPos;
extern  VERTEX3D *rayDir;
extern  FACET   **iFacet;
extern  FACET    *fLast;
extern  double    tNear;
extern  double    tFar;
extern  double    it1, it2;
extern  BOOL      AABBHit;

BOOL Intersect(VERTEX3D *rPos, VERTEX3D *rDir,  // Source ray (rayDir vector must be normalized)
	double *dist,                   // Distance to collision point
	FACET **iFact, FACET *last) {    // Collided facet, previous collision


	intMinLgth = 1e100;
	intFound = FALSE;
	intNbTHits = 0;
	rayPos = rPos;
	rayDir = rDir;
	intD.x = -rayDir->x;
	intD.y = -rayDir->y;
	intD.z = -rayDir->z;
	nullRx = (rayDir->x == 0.0);
	nullRy = (rayDir->y == 0.0);
	nullRz = (rayDir->z == 0.0);
	if (!nullRx) iRx = 1.0 / rayDir->x;
	if (!nullRy) iRy = 1.0 / rayDir->y;
	if (!nullRz) iRz = 1.0 / rayDir->z;
	iFacet = iFact;
	fLast = last;

	IntersectTree(sHandle->str[sHandle->curStruct].aabbTree);

	if (intFound) {

		FACET *f = *iFacet;
		*dist = intMinLgth;


		//ProfileFacet(f,sHandle->flightTimeCurrentParticle+*dist/100.0/sHandle->velocityCurrentParticle);
		f->hitted = TRUE;

		// Second pass for transparent hits
		for (int i = 0; i<intNbTHits; i++) {

			f = THits[i];
			if (f->colDist < intMinLgth) {
				double directionFactor = abs(DOT3(
					sHandle->pDir.x, sHandle->pDir.y, sHandle->pDir.z,
					f->sh.N.x, f->sh.N.y, f->sh.N.z));
				IncreaseFacetCounter(f, sHandle->flightTimeCurrentParticle + f->colDist / 100.0 / sHandle->velocityCurrentParticle, 1, 0, 0, 2.0 / (sHandle->velocityCurrentParticle*directionFactor), 2.0*(sHandle->useMaxwellDistribution ? 1.0 : 1.1781)*sHandle->velocityCurrentParticle*directionFactor);

				f->hitted = TRUE;
				if (f->hits && f->sh.countTrans) AHIT_FACET(f, sHandle->flightTimeCurrentParticle + f->colDist / 100.0 / sHandle->velocityCurrentParticle,
					TRUE, 2.0, 2.0);
				if (f->direction && f->sh.countDirection) DHIT_FACET(f, sHandle->flightTimeCurrentParticle + f->colDist / 100.0 / sHandle->velocityCurrentParticle);
				ProfileFacet(f, sHandle->flightTimeCurrentParticle + f->colDist / 100.0 / sHandle->velocityCurrentParticle,
					TRUE, 2.0, 2.0);
			}
		}
	}

	return intFound;

}

