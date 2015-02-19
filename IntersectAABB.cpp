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

// Minimum number of facet inside a BB
#define MINBB    8

// Maximum AABB tree depth
#define MAXDEPTH 5

// Temporary for intersection
static  double    intMinLgth;
static  BOOL      intFound;
static  VERTEX3D  intD;
static  VERTEX3D  intZ;
static  int       intNbTHits;
static  double    iRx;
static  double    iRy;
static  double    iRz;
static  BOOL      nullRx;
static  BOOL      nullRy;
static  BOOL      nullRz;
static  VERTEX3D *rayPos;
static  VERTEX3D *rayDir;
static  FACET   **iFacet;
static  FACET    *fLast;
static  double    tNear;
static  double    tFar;
static  double    it1,it2;
static  BOOL      AABBHit;

// -----------------------------------------------------------

int FindBestCuttingPlane(struct AABBNODE *node,int *left,int *right) {

	// AABB tree balancing

	double mX = (node->bb.min.x + node->bb.max.x)/2.0;
	double mY = (node->bb.min.y + node->bb.max.y)/2.0;
	double mZ = (node->bb.min.z + node->bb.max.z)/2.0;
	int Lxy=0;
	int Lxz=0;
	int Lyz=0;
	int ret;
	int i;
	double best = 1e100;
	double Cx,Cy,Cz;

	for(i=0;i<node->nbFacet;i++) {
		FACET *f = node->list[i];
		if(f->sh.center.x > mX) Lyz++;
		if(f->sh.center.y > mY) Lxz++;
		if(f->sh.center.z > mZ) Lxy++;
	}

	Cx = fabs((double)Lyz-(double)(node->nbFacet)/2.0);
	if(  Cx < best ) {
		best = Cx;
		*left = node->nbFacet - Lyz;
		*right = Lyz;
		ret = 1;
	}
	Cy = fabs((double)Lxz-(double)(node->nbFacet)/2.0);
	if(  Cy < best ) {
		best = Cy;
		*left = node->nbFacet - Lxz;
		*right = Lxz;
		ret = 2;
	}
	Cz = fabs((double)Lxy-(double)(node->nbFacet)/2.0);
	if(  Cz < best ) {
		best = Cz;
		*left = node->nbFacet - Lxy;
		*right = Lxy;
		ret = 3;
	}

	return ret;

}

// -----------------------------------------------------------

void ComputeBB(struct AABBNODE *node) {

	int i;

	node->bb.max.x = -1e100;
	node->bb.max.y = -1e100;
	node->bb.max.z = -1e100;
	node->bb.min.x =  1e100;
	node->bb.min.y =  1e100;
	node->bb.min.z =  1e100;

	for(i=0;i<node->nbFacet;i++) {
		FACET *f = node->list[i];
		if( f->sh.bb.min.x < node->bb.min.x ) node->bb.min.x = f->sh.bb.min.x;
		if( f->sh.bb.min.y < node->bb.min.y ) node->bb.min.y = f->sh.bb.min.y;
		if( f->sh.bb.min.z < node->bb.min.z ) node->bb.min.z = f->sh.bb.min.z;
		if( f->sh.bb.max.x > node->bb.max.x ) node->bb.max.x = f->sh.bb.max.x;
		if( f->sh.bb.max.y > node->bb.max.y ) node->bb.max.y = f->sh.bb.max.y;
		if( f->sh.bb.max.z > node->bb.max.z ) node->bb.max.z = f->sh.bb.max.z;
	}

}

// -----------------------------------------------------------

struct AABBNODE *BuildAABBTree(FACET **list,int nb,int depth) {

	int    i,l,r,cut,nbl=0,nbr=0;
	double m;
	struct AABBNODE *newNode;

	if(depth>=MAXDEPTH) return NULL;

	newNode = (struct AABBNODE *)malloc(sizeof(struct AABBNODE));
	memset(newNode,0,sizeof(struct AABBNODE));
	newNode->nbFacet = nb;
	newNode->list = list;
	ComputeBB(newNode);
	cut = FindBestCuttingPlane(newNode,&l,&r);

	if( l>=MINBB && r>=MINBB ) {

		// We can cut
		FACET **lList = (FACET **)malloc(sizeof(FACET *)*l);
		FACET **rList = (FACET **)malloc(sizeof(FACET *)*r);
		switch(cut) {

		case 1: // yz
			m = (newNode->bb.min.x + newNode->bb.max.x)/2.0;
			for(i=0;i<newNode->nbFacet;i++) {
				FACET *f = newNode->list[i];
				if(f->sh.center.x > m) rList[nbr++] = f;
				else                   lList[nbl++] = f;
			}
			break;

		case 2: // xz
			m = (newNode->bb.min.y + newNode->bb.max.y)/2.0;
			for(i=0;i<newNode->nbFacet;i++) {
				FACET *f = newNode->list[i];
				if(f->sh.center.y > m) rList[nbr++] = f;
				else                   lList[nbl++] = f;
			}
			break;

		case 3: // xy
			m = (newNode->bb.min.z + newNode->bb.max.z)/2.0;
			for(i=0;i<newNode->nbFacet;i++) {
				FACET *f = newNode->list[i];
				if(f->sh.center.z > m) rList[nbr++] = f;
				else                   lList[nbl++] = f;
			}
			break;

		}
		newNode->left  = BuildAABBTree(lList,nbl,depth+1);
		newNode->right = BuildAABBTree(rList,nbr,depth+1);

	}

	return newNode;

}

// -----------------------------------------------------------

void DestroyAABB(struct AABBNODE *node) {

	if( node!=NULL ) {
		DestroyAABB(node->left);
		DestroyAABB(node->right);
		free(node->list);
		free(node);
	}

}

// -----------------------------------------------------------
// Ray AABB intersection check (slabs method)

#define IntersectBB(n,lab)                                                   \
	\
	AABBHit = FALSE;                                                           \
	\
	if( nullRx ) {                                                             \
	if( rayPos->x < (n)->bb.min.x || rayPos->x > (n)->bb.max.x ) goto lab;   \
	tNear =  -1e100;                                                         \
	tFar  =   1e100;                                                         \
	} else {                                                                   \
	it1 = ((n)->bb.min.x - rayPos->x)*iRx;                                   \
	it2 = ((n)->bb.max.x - rayPos->x)*iRx;                                   \
	if( it2>it1 ) {                                                          \
	tFar  = it2;                                                           \
	tNear = it1;                                                           \
	} else {                                                                 \
	tFar  = it1;                                                           \
	tNear = it2;                                                           \
	}                                                                        \
	}                                                                          \
	if( tFar<0.0 ) goto lab;                                                   \
	\
	if( nullRy ) {                                                             \
	if( rayPos->y < (n)->bb.min.y || rayPos->y > (n)->bb.max.y ) goto lab;   \
	} else {                                                                   \
	it1 = ((n)->bb.min.y - rayPos->y)*iRy;                                   \
	it2 = ((n)->bb.max.y - rayPos->y)*iRy;                                   \
	if( it2>it1 ) {                                                          \
	if( it2<tFar  ) tFar  = it2;                                           \
	if( it1>tNear ) tNear = it1;                                           \
	} else {                                                                 \
	if( it1<tFar  ) tFar  = it1;                                           \
	if( it2>tNear ) tNear = it2;                                           \
	}                                                                        \
	}                                                                          \
	if( tNear>tFar || tFar<0.0 ) goto lab;                                     \
	\
	if( nullRz ) {                                                             \
	if( rayPos->z < (n)->bb.min.z || rayPos->z > (n)->bb.max.z ) goto lab;   \
	} else {                                                                   \
	it1 = ((n)->bb.min.z - rayPos->z)*iRz;                                   \
	it2 = ((n)->bb.max.z - rayPos->z)*iRz;                                   \
	if( it2>it1 ) {                                                          \
	if( it2<tFar  ) tFar  = it2;                                           \
	if( it1>tNear ) tNear = it1;                                           \
	} else {                                                                 \
	if( it1<tFar  ) tFar  = it1;                                           \
	if( it2>tNear ) tNear = it2;                                           \
	}                                                                        \
	}                                                                          \
	if( tNear>tFar || tFar<0.0 ) goto lab;                                     \
	AABBHit = TRUE;                                                            \
lab:


BOOL Visible(VERTEX3D *c1,VERTEX3D *c2,FACET *f1,FACET *f2) {

	VERTEX3D r;
	FACET *iF;

	// Check if there is an obstacle between c1 and c2
	intMinLgth = 1e100;
	intFound = FALSE;
	intNbTHits=0;
	rayPos = c1;
	r.x = (c2->x - c1->x);
	r.y = (c2->y - c1->y);
	r.z = (c2->z - c1->z);
	rayDir = &r;
	intD.x = -rayDir->x;
	intD.y = -rayDir->y;
	intD.z = -rayDir->z;
	nullRx = (rayDir->x == 0.0);
	nullRy = (rayDir->y == 0.0);
	nullRz = (rayDir->z == 0.0);
	if( !nullRx ) iRx = 1.0 / rayDir->x;
	if( !nullRy ) iRy = 1.0 / rayDir->y;
	if( !nullRz ) iRz = 1.0 / rayDir->z;
	iFacet = &iF;
	fLast = f1;

	IntersectTree(sHandle->str[0].aabbTree);

	if(intFound) {
		if( iF != f2 ) {
			// Obstacle found
			return FALSE;
		}
	}

	return TRUE;

}

// -----------------------------------------------------------
BOOL RaySphereIntersect(VERTEX3D *center,double radius,VERTEX3D *rPos,VERTEX3D *rDir,double *dist) {

	// Perform ray-sphere intersection
	double B,C,D;
	VERTEX3D s;
	s.x = (rPos->x - center->x);
	s.y = (rPos->y - center->y);
	s.z = (rPos->z - center->z);

	//|rDir|=1 => A=1
	//A = DOT3(rDir->x,rDir->y,rDir->z,rDir->x,rDir->y,rDir->z);
	B = 2.0 * DOT3(rDir->x,rDir->y,rDir->z,s.x,s.y,s.z);
	C = DOT3(s.x,s.y,s.z,s.x,s.y,s.z) - radius*radius;
	D = B*B - 4*C;

	if(D>=0.0) {
		double rD = sqrt(D);
		double t0 = (-B-rD)/2.0;
		double t1 = (-B+rD)/2.0;
		*dist = MIN(t0,t1);
		return (*dist>=0.0);
	}

	return FALSE;

}

// -----------------------------------------------------------

BOOL Intersect(VERTEX3D *rPos,VERTEX3D *rDir,  // Source ray (rayDir vector must be normalized)
			   double *dist,                   // Distance to collision point
			   FACET **iFact,FACET *last) {    // Collided facet, previous collision

				   int i;
				   intMinLgth = 1e100;
				   intFound = FALSE;
				   intNbTHits=0;
				   rayPos = rPos;
				   rayDir = rDir;
				   intD.x = -rayDir->x;
				   intD.y = -rayDir->y;
				   intD.z = -rayDir->z;
				   nullRx = (rayDir->x == 0.0);
				   nullRy = (rayDir->y == 0.0);
				   nullRz = (rayDir->z == 0.0);
				   if( !nullRx ) iRx = 1.0 / rayDir->x;
				   if( !nullRy ) iRy = 1.0 / rayDir->y;
				   if( !nullRz ) iRz = 1.0 / rayDir->z;
				   iFacet = iFact;
				   fLast = last;

				   IntersectTree(sHandle->str[sHandle->curStruct].aabbTree);

				   if(intFound) {

					   FACET *f = *iFacet;
					   *dist = intMinLgth;


					   //ProfileFacet(f,sHandle->flightTimeCurrentParticle+*dist/100.0/sHandle->velocityCurrentParticle);
					   f->hitted = TRUE;

					   // Second pass for transparent hits
					   for(i=0;i<intNbTHits;i++) {

						   f = THits[i];
						   if( f->colDist < intMinLgth ) {
							   double directionFactor = abs(DOT3(
								   sHandle->pDir.x, sHandle->pDir.y, sHandle->pDir.z,
								   f->sh.N.x, f->sh.N.y, f->sh.N.z));
							   f->sh.counter.hit.nbHit++;
							   f->sh.counter.hit.sum_1_per_ort_velocity += 2.0 / (sHandle->velocityCurrentParticle*directionFactor);
							   f->sh.counter.hit.sum_v_ort += 2.0*(sHandle->useMaxwellDistribution ? 1.0 : 1.1781)*sHandle->velocityCurrentParticle*directionFactor;
							   f->hitted = TRUE;
							   if( f->hits && f->sh.countTrans ) AHIT_FACET(f,sHandle->flightTimeCurrentParticle+f->colDist/100.0/sHandle->velocityCurrentParticle,
								   TRUE,2.0,2.0);
							   if (f->direction && f->sh.countDirection) DHIT_FACET(f, sHandle->flightTimeCurrentParticle + f->colDist / 100.0 / sHandle->velocityCurrentParticle);
							   ProfileFacet(f, sHandle->flightTimeCurrentParticle + f->colDist / 100.0 / sHandle->velocityCurrentParticle,
								   TRUE, 2.0, 2.0);
						   }
					   }
				   }

				   return intFound;

}

// -----------------------------------------------------------

void IntersectTree(struct AABBNODE *node) {

	// Returns the intersection between an oriented ray and
	// the geometry. If several intersection occurs, the
	// closest to rayPos is returned.
	// Returns TRUE is a collision occurs, FALSE otherwise.


	// Method: 3x3 Sytem solving for ray/rectangle intersection. 
	// Solve the vector equation u*U + v*V + d*D = Z (using Cramer's rule)
	// nuv = u^v (for faster calculation)
	int i;

	if( node->left==NULL || node->right==NULL ) {

		// Leaf
		for(i=0;i<node->nbFacet;i++) {

			FACET *f = node->list[i];

			// Do no check last collided facet
			if( f==fLast )
				continue;

			// Eliminate "back facet"
			if( (f->sh.is2sided) || (DOT3(f->sh.N.x,f->sh.N.y,f->sh.N.z,rayDir->x,rayDir->y,rayDir->z) < 0.0) ) {

				double u,v,d;

				// ---------------------------------------------------------------------------------
				// Ray/rectange instersection. Find (u,v,dist) and check 0<=u<=1, 0<=v<=1, dist>=0
				// ---------------------------------------------------------------------------------

				double det = DOT3(f->sh.Nuv.x,f->sh.Nuv.y,f->sh.Nuv.z,intD.x,intD.y,intD.z);

				if(det!=0.0) {

					double iDet = 1.0 / det;
					intZ.x = rayPos->x - f->sh.O.x;
					intZ.y = rayPos->y - f->sh.O.y;
					intZ.z = rayPos->z - f->sh.O.z;

					u = iDet * DET33( intZ.x , f->sh.V.x , intD.x ,
						intZ.y , f->sh.V.y , intD.y ,
						intZ.z , f->sh.V.z , intD.z );

					if( u>=0.0 && u<=1.0 ) {

						v = iDet * DET33( f->sh.U.x , intZ.x , intD.x ,
							f->sh.U.y , intZ.y , intD.y ,
							f->sh.U.z , intZ.z , intD.z );

						if( v>=0.0 && v<=1.0 ) {

							d = iDet * DOT3(f->sh.Nuv.x,f->sh.Nuv.y,f->sh.Nuv.z,intZ.x,intZ.y,intZ.z);

							if(d>0.0) {

								// Now check intersection with the facet polygon (in the u,v space)
								// This check could be avoided on rectangular facet.
								if( IsInFacet(f,u,v) ) {

									double time = sHandle->flightTimeCurrentParticle + d / 100.0 / sHandle->velocityCurrentParticle;
									if( (GetOpacityAt(f,time) == 1.0) || (rnd()<GetOpacityAt(f,time)) ) {

										// Hard hit
										if( d < intMinLgth ) {
											*iFacet = f;
											intFound = TRUE;
											intMinLgth = d;
											f->colU = u;
											f->colV = v;
										}

									} else {

										// Pass on partial transparent facet
										if( f->sh.isProfile || f->hits ) {
											f->colDist = d;
											f->colU = u;
											f->colV = v;
											if(intNbTHits<MAX_THIT)
												THits[intNbTHits++]=f;
										}
									}
								} // IsInFacet
							} // d range
						} // u range
					} // v range
				} // det==0
			} // dot<0
		} // end for

	} /* end is Leaf */ else {
		IntersectBB(node->left,jmp1);
		if( AABBHit ) IntersectTree(node->left);
		IntersectBB(node->right,jmp2);
		if( AABBHit ) IntersectTree(node->right);
	}

}

BOOL IsInFacet(FACET *f,double u,double v) {

	// 2D polygon "is inside" solving
	// Using the "Jordan curve theorem" (we intersect in v direction here)

	int n_updown,n_found,j;
	double x1,x2,y1,y2,a,minx,maxx;

	n_updown=0;
	n_found=0;

	for (j = 0; j < f->sh.nbIndex-1; j++) {

		x1 = f->vertices2[j].u;
		y1 = f->vertices2[j].v;
		x2 = f->vertices2[j+1].u;
		y2 = f->vertices2[j+1].v;

		if( x2>x1 ) { minx=x1;maxx=x2; } 
		else        { minx=x2;maxx=x1; }

		if (u > minx && u <= maxx) {
			a = (y2 - y1) / (x2 - x1);
			if ((a*(u-x1) + y1) < v) {
				n_updown = n_updown + 1;
			} else {
				n_updown = n_updown - 1;
			}
			n_found++;
		}

	}

	// Last point
	x1 = f->vertices2[j].u;
	y1 = f->vertices2[j].v;
	x2 = f->vertices2[0].u;
	y2 = f->vertices2[0].v;

	if( x2>x1 ) { minx=x1;maxx=x2; } 
	else        { minx=x2;maxx=x1; }

	if (u > minx && u <= maxx) {
		a = (y2 - y1) / (x2 - x1);
		if ((a*(u-x1) + y1) < v) {
			n_updown = n_updown + 1;
		} else {
			n_updown = n_updown - 1;
		}
		n_found++;
	}

	if (n_updown<0) n_updown=-n_updown;
	return (((n_found/2)&1) ^ ((n_updown/2)&1));

}
