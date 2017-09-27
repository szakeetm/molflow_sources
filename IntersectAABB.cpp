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
//#include <malloc.h>
#include <string.h>
#include <math.h>
#include "IntersectAABB_shared.h"
#include "Random.h"
#include "Simulation.h"

std::tuple<bool,FACET*,double> Intersect(const Vector3d& rayPos, const Vector3d& rayDir, /*FACET* lastHitFacet,*/ FACET**& THitCache) {
	// Source ray (rayDir vector must be normalized)
	// lastHit is to avoid detecting twice the same collision
	// returns bool found (is there a collision), pointer to collided facet, double d (distance to collision)

	bool nullRx = (rayDir.x == 0.0);
	bool nullRy = (rayDir.y == 0.0);
	bool nullRz = (rayDir.z == 0.0);
	Vector3d inverseRayDir;
	if (!nullRx) inverseRayDir.x = 1.0 / rayDir.x;
	if (!nullRy) inverseRayDir.y = 1.0 / rayDir.y;
	if (!nullRz) inverseRayDir.z = 1.0 / rayDir.z;
	
	//Global variables, easier for recursion:
	size_t intNbTHits = 0;

	//Output values
	bool found=false;
	FACET *collidedFacet;
	double minLength=1e100;

	IntersectTree(sHandle->str[sHandle->curStruct].aabbTree,rayPos, -1.0*rayDir,sHandle->lastHit,
		nullRx,nullRy,nullRz,inverseRayDir,
		intNbTHits,THitCache,found,collidedFacet,minLength); //output params

	if (found) {

		//ProfileFacet(f,sHandle->flightTimeCurrentParticle+*dist/100.0/sHandle->velocityCurrentParticle);
		collidedFacet->hitted = true;

		// Second pass for transparent hits
		for (int i = 0; i<intNbTHits; i++) {

			FACET* f = THitCache[i];
			if (collidedFacet->colDist < minLength) {
				double directionFactor = abs(Dot(sHandle->pDir,collidedFacet->sh.N));
				IncreaseFacetCounter(f, sHandle->flightTimeCurrentParticle + collidedFacet->colDist / 100.0 / sHandle->velocityCurrentParticle, 1, 0, 0, 2.0 / (sHandle->velocityCurrentParticle*directionFactor), 2.0*(sHandle->useMaxwellDistribution ? 1.0 : 1.1781)*sHandle->velocityCurrentParticle*directionFactor);

				collidedFacet->hitted = true;
				if (collidedFacet->hits && collidedFacet->sh.countTrans) {
					RecordHitOnTexture(f, sHandle->flightTimeCurrentParticle + collidedFacet->colDist / 100.0 / sHandle->velocityCurrentParticle,
						true, 2.0, 2.0);
				}
				if (collidedFacet->direction && collidedFacet->sh.countDirection) {
					RecordDirectionVector(f, sHandle->flightTimeCurrentParticle + collidedFacet->colDist / 100.0 / sHandle->velocityCurrentParticle);
				}
				ProfileFacet(f, sHandle->flightTimeCurrentParticle + collidedFacet->colDist / 100.0 / sHandle->velocityCurrentParticle,
					true, 2.0, 2.0);
				if (collidedFacet->sh.anglemapParams.record) RecordAngleMap(f);
			}
		}
	}
	return std::make_tuple(found,collidedFacet,minLength);

}

bool Visible(Vector3d *c1, Vector3d *c2, FACET *f1, FACET *f2, FACET** THitCache) {
	//For AC matrix calculation

	Vector3d rayPos = *c1;
	Vector3d rayDir = *c2 - *c1;

	bool nullRx = (rayDir.x == 0.0);
	bool nullRy = (rayDir.y == 0.0);
	bool nullRz = (rayDir.z == 0.0);
	Vector3d inverseRayDir;
	if (!nullRx) inverseRayDir.x = 1.0 / rayDir.x;
	if (!nullRy) inverseRayDir.y = 1.0 / rayDir.y;
	if (!nullRz) inverseRayDir.z = 1.0 / rayDir.z;

	//Global variables, easier for recursion:
	size_t intNbTHits = 0;

	//Output values
	bool found;
	FACET *collidedFacet;
	double minLength;

	IntersectTree(sHandle->str[0].aabbTree, rayPos,-1.0*rayDir, 
		f1, nullRx, nullRy, nullRz, inverseRayDir, intNbTHits, THitCache, found,collidedFacet,minLength);

	if (found) {
		if (collidedFacet != f2) {
			// Obstacle found
			return false;
		}
	}

	return true;
}