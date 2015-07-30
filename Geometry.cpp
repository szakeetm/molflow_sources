/*
File:        Geometry.cpp
Description: Main geometry class (Handles sets of facets)
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

#include "Geometry.h"
#include "ASELoader.h"
#include "Utils.h"
#include <malloc.h>
#include <string.h>
#include <math.h>
#include "GLApp/GLMatrix.h"
#include "GLApp\GLMessageBox.h"
#include "Molflow.h"
#include "GLApp\GLWindowManager.h"
#include "Distributions.h" //InterpolateY

#define WRITEVAL(_value,_type) *((_type *)buffer)=_value;buffer += sizeof(_type)

/*
//Leak detection
#ifdef _DEBUG
#define _CRTDBG_MAP_ALLOC
#define DEBUG_NEW new(_NORMAL_BLOCK, __FILE__, __LINE__)
#define new DEBUG_NEW
#endif
*/

using namespace pugi;
extern MolFlow *mApp;

Geometry::Geometry() {

	facets = NULL;
	vertices3 = NULL;
	polyList = 0;
	selectList = 0;
	selectList2 = 0;
	selectList3 = 0;
	arrowList = 0;
	sphereList = 0;
	texture_limits[0].autoscale.min.all = texture_limits[0].autoscale.min.moments_only =
		texture_limits[1].autoscale.min.all = texture_limits[1].autoscale.min.moments_only =
		texture_limits[2].autoscale.min.all = texture_limits[2].autoscale.min.moments_only =
		texture_limits[0].manual.min.all = texture_limits[0].manual.min.moments_only =
		texture_limits[1].manual.min.all = texture_limits[1].manual.min.moments_only =
		texture_limits[2].manual.min.all = texture_limits[2].manual.min.moments_only = 0.0;
	texture_limits[0].autoscale.max.all = texture_limits[0].autoscale.max.moments_only =
		texture_limits[1].autoscale.max.all = texture_limits[1].autoscale.max.moments_only =
		texture_limits[2].autoscale.max.all = texture_limits[2].autoscale.max.moments_only =
		texture_limits[0].manual.max.all = texture_limits[0].manual.max.moments_only =
		texture_limits[1].manual.max.all = texture_limits[1].manual.max.moments_only =
		texture_limits[2].manual.max.all = texture_limits[2].manual.max.moments_only = 1.0;
	textureMode = 0;
	autoNorme = TRUE;
	centerNorme = TRUE;
	normeRatio = 1.0f;
	texAutoScale = TRUE;
	texAutoScaleIncludeConstantFlow = TRUE;
	texLogScale = FALSE;

	sh.nbSuper = 0;
	viewStruct = -1;
	strcpy(strPath, "");
	Clear();

	//Debug memory check
	//_ASSERTE (!_CrtDumpMemoryLeaks());;
	_ASSERTE(_CrtCheckMemory());
}

Geometry::~Geometry() {
	Clear();
}

void Geometry::Clear() {

	// Free memory
	if (facets) {
		for (int i = 0; i < sh.nbFacet; i++)
			SAFE_DELETE(facets[i]);
		free(facets);
	}
	if (vertices3) free(vertices3);
	for (int i = 0; i < sh.nbSuper; i++) {
		SAFE_FREE(strName[i]);
		SAFE_FREE(strFileName[i]);
	}
	memset(strName, 0, MAX_SUPERSTR * sizeof(char *));
	memset(strFileName, 0, MAX_SUPERSTR * sizeof(char *));
	DeleteGLLists(TRUE, TRUE);

	// Init default
	facets = NULL;         // Facets array
	vertices3 = NULL;      // Facets vertices in (x,y,z) space
	sh.nbFacet = 0;        // Number of facets
	sh.nbVertex = 0;       // Number of vertex
	isLoaded = FALSE;      // isLoaded flag
	texture_limits[0].autoscale.min.all = texture_limits[0].autoscale.min.moments_only =
		texture_limits[1].autoscale.min.all = texture_limits[1].autoscale.min.moments_only =
		texture_limits[2].autoscale.min.all = texture_limits[2].autoscale.min.moments_only =
		texture_limits[0].manual.min.all = texture_limits[0].manual.min.moments_only =
		texture_limits[1].manual.min.all = texture_limits[1].manual.min.moments_only =
		texture_limits[2].manual.min.all = texture_limits[2].manual.min.moments_only = 0.0;
	texture_limits[0].autoscale.max.all = texture_limits[0].autoscale.max.moments_only =
		texture_limits[1].autoscale.max.all = texture_limits[1].autoscale.max.moments_only =
		texture_limits[2].autoscale.max.all = texture_limits[2].autoscale.max.moments_only =
		texture_limits[0].manual.max.all = texture_limits[0].manual.max.moments_only =
		texture_limits[1].manual.max.all = texture_limits[1].manual.max.moments_only =
		texture_limits[2].manual.max.all = texture_limits[2].manual.max.moments_only = 1.0;
	sh.nbSuper = 0;          // Structure number
	nbSelected = 0;          // Number of selected facet

	memset(lineList, 0, sizeof(lineList));
	memset(strName, 0, sizeof(strName));
	memset(strFileName, 0, sizeof(strFileName));

	// Init OpenGL material
	memset(&whiteMaterial, 0, sizeof(GLMATERIAL));
	whiteMaterial.Diffuse.r = 0.9f;
	whiteMaterial.Diffuse.g = 0.9f;
	whiteMaterial.Diffuse.b = 0.9f;
	whiteMaterial.Ambient.r = 0.9f;
	whiteMaterial.Ambient.g = 0.9f;
	whiteMaterial.Ambient.b = 0.9f;

	memset(&fillMaterial, 0, sizeof(GLMATERIAL));
	fillMaterial.Diffuse.r = 0.6f;
	fillMaterial.Diffuse.g = 0.65f;
	fillMaterial.Diffuse.b = 0.65f;
	fillMaterial.Ambient.r = 0.45f;
	fillMaterial.Ambient.g = 0.41f;
	fillMaterial.Ambient.b = 0.41f;

	memset(&arrowMaterial, 0, sizeof(GLMATERIAL));
	arrowMaterial.Diffuse.r = 0.4f;
	arrowMaterial.Diffuse.g = 0.2f;
	arrowMaterial.Diffuse.b = 0.0f;
	arrowMaterial.Ambient.r = 0.6f;
	arrowMaterial.Ambient.g = 0.3f;
	arrowMaterial.Ambient.b = 0.0f;

	nbSelectedHist = 0;    // Selection history
	nbSelectedHistVertex = 0;

	//Debug memory check
	//_ASSERTE (!_CrtDumpMemoryLeaks());;
	_ASSERTE(_CrtCheckMemory());
}

void Geometry::CalcTotalOutGassing() {

	/*totalOutgassing = 0.0;
	totalInFlux = 0.0;
	for (int i = 0; i<sh.nbFacet; i++) {
	if (facets[i]->sh.desorbType>0) { //has desorption
	if (facets[i]->sh.useOutgassingFile) { //uses outgassing file values
	for (int l = 0; l < (facets[i]->sh.outgassingMapWidth*facets[i]->sh.outgassingMapHeight); l++) {
	totalOutgassing += facets[i]->outgassingMap[l];
	totalInFlux += facets[i]->outgassingMap[l] / (1.38E-23*facets[i]->sh.temperature);
	}
	} else { //uses absolute outgassing values
	totalOutgassing += facets[i]->sh.flow;
	totalInFlux += facets[i]->sh.flow / (1.38E-23*facets[i]->sh.temperature);
	}
	}
	}
	//totalOutgassing*=0.001; //conversion from "unit*l/s" to "unit*m3/sec"

	if (mApp->globalSettings) mApp->globalSettings->UpdateOutgassing();*/
}

void Geometry::CalculateFacetParam(int facet) {

	Facet *f = facets[facet];

	// Calculate facet normal
	VERTEX3D p0 = vertices3[facets[facet]->indices[0]];
	VERTEX3D v1;
	VERTEX3D v2;
	BOOL consecutive = TRUE;
	int ind = 2;

	// TODO: Handle possible collinear consequtive vectors
	int i0 = facets[facet]->indices[0];
	int i1 = facets[facet]->indices[1];
	while (ind < f->sh.nbIndex && consecutive) {
		int i2 = facets[facet]->indices[ind++];

		Sub(&v1, vertices3 + i1, vertices3 + i0); // v1 = P0P1
		Sub(&v2, vertices3 + i2, vertices3 + i1); // v2 = P1P2
		Cross(&(f->sh.N), &v1, &v2);              // Cross product
		consecutive = (Norme(&(f->sh.N)) < 1e-11);
	}
	f->collinear = consecutive; //mark for later that this facet was on a line
	Normalize(&(f->sh.N));                  // Normalize

	// Calculate Axis Aligned Bounding Box
	f->sh.bb.min.x = 1e100;
	f->sh.bb.min.y = 1e100;
	f->sh.bb.min.z = 1e100;
	f->sh.bb.max.x = -1e100;
	f->sh.bb.max.y = -1e100;
	f->sh.bb.max.z = -1e100;

	for (int i = 0; i < f->sh.nbIndex; i++) {
		VERTEX3D p = vertices3[f->indices[i]];
		if (p.x < f->sh.bb.min.x) f->sh.bb.min.x = p.x;
		if (p.y < f->sh.bb.min.y) f->sh.bb.min.y = p.y;
		if (p.z < f->sh.bb.min.z) f->sh.bb.min.z = p.z;
		if (p.x > f->sh.bb.max.x) f->sh.bb.max.x = p.x;
		if (p.y > f->sh.bb.max.y) f->sh.bb.max.y = p.y;
		if (p.z > f->sh.bb.max.z) f->sh.bb.max.z = p.z;
	}

	// Facet center (AABB center)
	f->sh.center.x = (f->sh.bb.max.x + f->sh.bb.min.x) / 2.0;
	f->sh.center.y = (f->sh.bb.max.y + f->sh.bb.min.y) / 2.0;
	f->sh.center.z = (f->sh.bb.max.z + f->sh.bb.min.z) / 2.0;

	// Plane equation
	double A = f->sh.N.x;
	double B = f->sh.N.y;
	double C = f->sh.N.z;
	double D = -Dot(&(f->sh.N), &p0);

	// Facet planeity
	int nb = f->sh.nbIndex;
	double max = 0.0;
	for (int i = 1; i<nb; i++) {
		VERTEX3D p = vertices3[f->indices[i]];
		double d = A * p.x + B * p.y + C * p.z + D;
		if (fabs(d)>fabs(max)) max = d;
	}

	// Plane coef
	f->a = A;
	f->b = B;
	f->c = C;
	f->d = D;
	f->err = max;

	f->sh.maxSpeed = 4.0 * sqrt(2.0*8.31*f->sh.temperature / 0.001 / mApp->worker.gasMass);
}

void Geometry::InitializeGeometry(int facet_number) {


	// Perform various precalculation here for a faster simulation.

	//GLProgress *initGeoPrg = new GLProgress("Initializing geometry...","Please wait");
	//initGeoPrg->SetProgress(0.0);
	//initGeoPrg->SetVisible(TRUE);
	if (facet_number == -1) { //bounding box for all vertices
		bb.min.x = 1e100;
		bb.min.y = 1e100;
		bb.min.z = 1e100;
		bb.max.x = -1e100;
		bb.max.y = -1e100;
		bb.max.z = -1e100;

		// Axis Aligned Bounding Box
		for (int i = 0; i < sh.nbVertex; i++) {
			VERTEX3D p = vertices3[i];
			if (!(vertices3[i].selected == FALSE || vertices3[i].selected == TRUE)) vertices3[i].selected = FALSE; //initialize selection
			if (p.x < bb.min.x) bb.min.x = p.x;
			if (p.y < bb.min.y) bb.min.y = p.y;
			if (p.z < bb.min.z) bb.min.z = p.z;
			if (p.x > bb.max.x) bb.max.x = p.x;
			if (p.y > bb.max.y) bb.max.y = p.y;
			if (p.z > bb.max.z) bb.max.z = p.z;
		}
	}
	else { //bounding box only for the changed facet
		for (int i = 0; i < facets[facet_number]->sh.nbIndex; i++) {
			VERTEX3D p = vertices3[facets[facet_number]->indices[i]];
			//if(!(vertices3[i].selected==FALSE || vertices3[i].selected==TRUE)) vertices3[i].selected=FALSE; //initialize selection
			if (p.x < bb.min.x) bb.min.x = p.x;
			if (p.y < bb.min.y) bb.min.y = p.y;
			if (p.z < bb.min.z) bb.min.z = p.z;
			if (p.x > bb.max.x) bb.max.x = p.x;
			if (p.y > bb.max.y) bb.max.y = p.y;
			if (p.z > bb.max.z) bb.max.z = p.z;
		}
	}

	center.x = (bb.max.x + bb.min.x) / 2.0;
	center.y = (bb.max.y + bb.min.y) / 2.0;
	center.z = (bb.max.z + bb.min.z) / 2.0;

	// Choose an orthogonal (U,V) 2D basis for each facet. (This method can be 
	// improved. stub). The algorithm chooses the longest vedge for the U vector.
	// then it computes V (orthogonal to U and N). Afterwards, U and V are rescaled 
	// so each facet vertex are included in the rectangle defined by uU + vV (0<=u<=1 
	// and 0<=v<=1) of origin f->sh.O, U and V are always orthogonal and (U,V,N) 
	// form a 3D left handed orthogonal basis (not necessary orthonormal).
	// This coordinates system allows to prevent from possible "almost degenerated"
	// basis on fine geometry. It also greatly eases the facet/ray instersection routine 
	// and ref/abs/des hit recording and visualization. In order to ease some calculations, 
	// nU et nV (normalized U et V) are also stored in the Facet structure.
	// The local coordinates of facet vertex are stored in (U,V) coordinates (vertices2).

	int fOffset = sizeof(SHGHITS);
	for (int i = 0; i < sh.nbFacet; i++) {
		//initGeoPrg->SetProgress((double)i/(double)sh.nbFacet);
		if ((facet_number == -1) || (i == facet_number)) { //permits to initialize only one facet
			// Main facet params
			CalculateFacetParam(i);

			// Current facet
			Facet *f = facets[i];

			/*
			// Search longest edge
			double dMax = 0.0;
			int    i1Max,i2Max;
			for(int j=0;j<f->sh.nbIndex;j++) {
			int j1 = IDX(j+1,f->sh.nbIndex);
			int i1 = f->indices[j];
			int i2 = f->indices[j1];
			double xl = (vertices3[i1].x - vertices3[i2].x);
			double yl = (vertices3[i1].y - vertices3[i2].y);
			double zl = (vertices3[i1].z - vertices3[i2].z);
			double l = DOT3( xl,yl,zl,xl,yl,zl ); //distance of vertices at i1 and i2
			if( l>dMax + 1e-8 ) { //disregard differences below double precision
			dMax = l;
			i1Max = i1;
			i2Max = i2;
			}
			}

			// First vertex
			int i0 = f->indices[0];
			VERTEX3D p0 = vertices3[i0];
			*/

			VERTEX3D p0 = vertices3[f->indices[0]];
			VERTEX3D p1 = vertices3[f->indices[1]];
			VERTEX3D U, V;

			/* OLD fashion (no longueur used)
			// Intersection line (directed by U and including p0)
			U.x = f->c;
			U.y = 0.0;
			U.z = -f->a;
			double nU = Norme(&U);

			if( IS_ZERO(nU) ) {
			// Plane parallel to (O,x,z)
			U.x = 1.0;
			U.y = 0.0;
			U.z = 0.0;
			} else {
			ScalarMult(&U,1.0/nU); // Normalize U
			}
			*/
			/*
			U.x = vertices3[i2Max].x - vertices3[i1Max].x;
			U.y = vertices3[i2Max].y - vertices3[i1Max].y;
			U.z = vertices3[i2Max].z - vertices3[i1Max].z;
			*/
			U.x = p1.x - p0.x;
			U.y = p1.y - p0.y;
			U.z = p1.z - p0.z;

			double nU = Norme(&U);
			ScalarMult(&U, 1.0 / nU); // Normalize U

			// Construct a normal vector V
			Cross(&V, &(f->sh.N), &U); // |U|=1 and |N|=1 => |V|=1

			// u,v vertices (we start with p0 at 0,0)
			f->vertices2[0].u = 0.0;
			f->vertices2[0].v = 0.0;
			VERTEX2D min; min.u = 0.0; min.v = 0.0;
			VERTEX2D max; max.u = 0.0; max.v = 0.0;

			for (int j = 1; j < f->sh.nbIndex; j++) {

				VERTEX3D p = vertices3[f->indices[j]];
				VERTEX3D v;
				Sub(&v, &p, &p0); // v = p0p
				f->vertices2[j].u = Dot(&U, &v);  // Project p on U along the V direction
				f->vertices2[j].v = Dot(&V, &v);  // Project p on V along the U direction

				// Bounds
				if (f->vertices2[j].u > max.u) max.u = f->vertices2[j].u;
				if (f->vertices2[j].v > max.v) max.v = f->vertices2[j].v;
				if (f->vertices2[j].u < min.u) min.u = f->vertices2[j].u;
				if (f->vertices2[j].v < min.v) min.v = f->vertices2[j].v;

			}

			// Calculate facet area (Meister/Gauss formula)
			double A = 0.0;
			for (int j = 0; j < f->sh.nbIndex; j++) {
				int j1 = IDX(j + 1, f->sh.nbIndex);
				A += f->vertices2[j].u*f->vertices2[j1].v - f->vertices2[j1].u*f->vertices2[j].v;
			}
			f->sh.area = fabs(0.5 * A);

			// Compute the 2D basis (O,U,V)
			double uD = (max.u - min.u);
			double vD = (max.v - min.v);

			// Origin
			f->sh.O.x = min.u * U.x + min.v * V.x + p0.x;
			f->sh.O.y = min.u * U.y + min.v * V.y + p0.y;
			f->sh.O.z = min.u * U.z + min.v * V.z + p0.z;

			// Rescale U and V vector
			f->sh.nU = U;
			ScalarMult(&U, uD);
			f->sh.U = U;

			f->sh.nV = V;
			ScalarMult(&V, vD);
			f->sh.V = V;

			Cross(&(f->sh.Nuv), &U, &V);

			// Rescale u,v coordinates
			for (int j = 0; j < f->sh.nbIndex; j++) {

				VERTEX2D p = f->vertices2[j];
				f->vertices2[j].u = (p.u - min.u) / uD;
				f->vertices2[j].v = (p.v - min.v) / vD;

			}

			// Detect non visible edge
			f->InitVisibleEdge();

			// Detect orientation
			f->DetectOrientation();

			// Hit address
			f->sh.hitOffset = fOffset;
			fOffset += f->GetHitsSize((int)mApp->worker.moments.size());
		}
	}


	if (facet_number == -1) {
		BuildGLList();
		mApp->UpdateModelParams();
		mApp->UpdateFacetParams();
	}
	//initGeoPrg->SetVisible(FALSE);
	//SAFE_DELETE(initGeoPrg);

	//Debug memory check
	//_ASSERTE (!_CrtDumpMemoryLeaks());
	_ASSERTE(_CrtCheckMemory());
}

void Geometry::RebuildLists() {
	BuildGLList();
}

DWORD Geometry::GetGeometrySize() {

	Worker  *work = &mApp->worker;
	// Compute number of bytes allocated
	DWORD memoryUsage = 0;
	memoryUsage += sizeof(SHGEOM);
	memoryUsage += sh.nbVertex * sizeof(VERTEX3D);
	for (int i = 0; i < sh.nbFacet; i++)
		memoryUsage += facets[i]->GetGeometrySize();

	//CDFs
	memoryUsage += sizeof(size_t); //number of CDFs
	for (auto i : work->CDFs) {
		memoryUsage += sizeof(size_t); //CDF size
		memoryUsage += i.size() * 2 * sizeof(double);
	}

	//IDs
	memoryUsage += sizeof(size_t); //number of IDs
	for (auto i : work->IDs) {
		memoryUsage += sizeof(size_t); //ID size
		memoryUsage += i.size() * 2 * sizeof(double);
	}

	//Parameters
	memoryUsage += sizeof(size_t); //number of parameters
	for (auto i : work->parameters) {
		memoryUsage += sizeof(size_t); //parameter size
		memoryUsage += i.values.size() * 2 * sizeof(double);
	}
	memoryUsage += sizeof(size_t); //number of temperatures
	memoryUsage += sizeof(double)*(int)(work->temperatures).size(); //temperatures

	//moments size already passed
	memoryUsage += sizeof(double)*(int)(work->moments).size(); //moments

	memoryUsage += sizeof(size_t); //number of desparamIDs
	memoryUsage += sizeof(size_t)*(int)(work->desorptionParameterIDs).size(); //desparamIDs
	return memoryUsage;
}

void Geometry::CopyGeometryBuffer(BYTE *buffer) {

	// Build shared buffer for geometry (see Shared.h)
	// Basically we serialize all data and let the subprocesses read them

	/*
	Memory map:

	-->bufferStart
	SHGEOM (nbFacets, time-dep parameters, gas mass, etc.)
	vertices3 (nbVertex times VERTEX3D struct)
	FOR EACH FACET
	SHFACET
	indices (nbIndex times int)
	vertices2 (nbIndex times VERTEX2D struct)
	[outgassingMap (height*width*double)]
	-->incBuff
	[inc Map: for each facet with texture, height*width*double]
	CDFs.size()
	CDFs
	IDs.size()
	[IDs]
	parameters.size()
	[parameters]
	temperatures.size()
	temperatures

	*/


	Worker *w = &mApp->worker;
	int fOffset = sizeof(SHGHITS);
	SHGEOM *shGeom = (SHGEOM *)buffer;

	sh.nbMoments = w->moments.size();
	sh.latestMoment = w->latestMoment;
	sh.totalDesorbedMolecules = w->totalDesorbedMolecules;
	sh.finalOutgassingRate = w->finalOutgassingRate;
	sh.gasMass = w->gasMass;
	sh.halfLife = w->halfLife;
	sh.timeWindowSize = w->timeWindowSize;
	sh.useMaxwellDistribution = w->useMaxwellDistribution;
	sh.calcConstantFlow = w->calcConstantFlow;

	sh.motionType = w->motionType;
	sh.motionVector1 = w->motionVector1;
	sh.motionVector2 = w->motionVector2;

	memcpy(shGeom, &(this->sh), sizeof(SHGEOM));
	buffer += sizeof(SHGEOM);

	memcpy(buffer, vertices3, sizeof(VERTEX3D)*sh.nbVertex);
	buffer += sizeof(VERTEX3D)*sh.nbVertex;

	for (int i = 0; i < sh.nbFacet; i++) {
		Facet *f = facets[i];
		f->sh.hitOffset = fOffset;
		fOffset += f->GetHitsSize((int)mApp->worker.moments.size());
		memcpy(buffer, &(f->sh), sizeof(SHFACET));
		buffer += sizeof(SHFACET);
		memcpy(buffer, f->indices, sizeof(int)*f->sh.nbIndex);
		buffer += sizeof(int)*f->sh.nbIndex;
		memcpy(buffer, f->vertices2, sizeof(VERTEX2D)*f->sh.nbIndex);
		buffer += sizeof(VERTEX2D)*f->sh.nbIndex;
		if (f->sh.useOutgassingFile) {
			memcpy(buffer, f->outgassingMap, sizeof(double)*f->sh.outgassingMapWidth*f->sh.outgassingMapHeight);
			buffer += sizeof(double)*f->sh.outgassingMapWidth*f->sh.outgassingMapHeight;
		}
	}

	// Add surface elements area (reciprocal)
	for (int k = 0; k < sh.nbFacet; k++) {
		Facet *f = facets[k];
		DWORD add = 0;
		if (f->sh.isTextured) {

			if (f->mesh) {

				for (int j = 0; j < f->sh.texHeight; j++) {
					for (int i = 0; i < f->sh.texWidth; i++) {
						double area = f->mesh[add].area*(f->sh.is2sided ? 2.0 : 1.0);
						if (area > 0.0) {
							// Use the sign bit to store isFull flag
							if (f->mesh[add].full)
							{
								WRITEVAL(-1.0 / area, double);
							}
							else
							{
								WRITEVAL(1.0 / area, double);
							}
						}
						else {
							WRITEVAL(0.0, double);
						}
						add++;
					}
				}

			}
			else {

				double rw = Norme(&(f->sh.U)) / (double)(f->sh.texWidthD);
				double rh = Norme(&(f->sh.V)) / (double)(f->sh.texHeightD);
				double area = rw*rh;

				for (int j = 0; j < f->sh.texHeight; j++) {
					for (int i = 0; i < f->sh.texWidth; i++) {
						if (area > 0.0) {
							WRITEVAL(1.0 / area, double);
						}
						else {
							WRITEVAL(0.0, double);
						}
					}
				}

			}
		}
	}


	//CDFs
	WRITEVAL(w->CDFs.size(), size_t);
	for (size_t i = 0; i < w->CDFs.size(); i++) {
		WRITEVAL(w->CDFs[i].size(), size_t);
		for (size_t j = 0; j < w->CDFs[i].size(); j++) {
			WRITEVAL(w->CDFs[i][j].first, double);
			WRITEVAL(w->CDFs[i][j].second, double);
		}
	}

	//IDs
	WRITEVAL(w->IDs.size(), size_t);
	for (size_t i = 0; i < w->IDs.size(); i++) {
		WRITEVAL(w->IDs[i].size(), size_t);
		for (size_t j = 0; j < w->IDs[i].size(); j++) {
			WRITEVAL(w->IDs[i][j].first, double);
			WRITEVAL(w->IDs[i][j].second, double);
		}
	}

	//Parameters
	/*WRITEVAL(w->parameters.size(), size_t);
	  for (size_t i = 0; i < w->parameters.size(); i++) {
	  WRITEVAL(w->parameters[i].values.size(), size_t);
	  for (size_t j=0;j<w->parameters[i].values.size();j++) {
	  WRITEVAL(w->parameters[i].values[j].first, double);
	  WRITEVAL(w->parameters[i].values[j].second, double);
	  }
	  }*/
	WRITEVAL(w->parameters.size(), size_t);
	for (auto i : w->parameters) {
		WRITEVAL(i.values.size(), size_t);
		for (auto j : i.values) {
			WRITEVAL(j.first, double);
			WRITEVAL(j.second, double);
		}
	}

	//Temperatures
	WRITEVAL(w->temperatures.size(), size_t);
	for (size_t i = 0; i < w->temperatures.size(); i++) {
		WRITEVAL(w->temperatures[i], double);
	}

	//Time moments
	//WRITEVAL(w->moments.size(), size_t); //nbMoments already passed
	for (size_t i = 0; i < w->moments.size(); i++) {
		WRITEVAL(w->moments[i], double);
	}

	//Desorption parameter IDs
	WRITEVAL(w->desorptionParameterIDs.size(), size_t);
	for (size_t i = 0; i < w->desorptionParameterIDs.size(); i++) {
		WRITEVAL(w->desorptionParameterIDs[i], size_t);
	}

}

void Geometry::SetAutoNorme(BOOL enable) {
	autoNorme = enable;
}

BOOL Geometry::GetAutoNorme() {
	return autoNorme;
}

void Geometry::SetCenterNorme(BOOL enable) {
	centerNorme = enable;
}

BOOL Geometry::GetCenterNorme() {
	return centerNorme;
}

void Geometry::SetNormeRatio(float r) {
	normeRatio = r;
}

float Geometry::GetNormeRatio() {
	return normeRatio;
}

DWORD Geometry::GetHitsSize(std::vector<double> *moments) {

	// Compute number of bytes allocated
	DWORD memoryUsage = 0;
	memoryUsage += sizeof(SHGHITS);
	for (int i = 0; i < sh.nbFacet; i++) {
		memoryUsage += facets[i]->GetHitsSize((int)mApp->worker.moments.size());
	}

	return memoryUsage;
}

DWORD Geometry::GetMaxElemNumber() {

	int nbElem = 0;
	for (int i = 0; i < sh.nbFacet; i++) {
		Facet *f = facets[i];
		if (f->mesh) nbElem += f->sh.texWidth*f->sh.texHeight;
		else          return 0;
	}
	return nbElem;

}


void Geometry::CopyElemBuffer(BYTE *buffer) {

	int idx = 0;
	for (int i = 0; i < sh.nbFacet; i++) {
		Facet *f = facets[i];
		int sz = f->sh.texWidth * f->sh.texHeight * sizeof(SHELEM);
		memcpy(buffer + idx, f->mesh, sz);
		idx += sz;
	}

}


void Geometry::BuildShapeList() {

	// Shapes used for direction field rendering

	// 3D arrow (direction field)
	int nbDiv = 10;
	double alpha = 2.0*PI / (double)nbDiv;

	arrowList = glGenLists(1);
	glNewList(arrowList, GL_COMPILE);

	glEnable(GL_CULL_FACE);
	glCullFace(GL_FRONT);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glEnable(GL_LIGHTING);
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_BLEND);

	glBegin(GL_TRIANGLES);

	// Arrow
	for (int i = 0; i < nbDiv; i++) {

		double y1 = sin(alpha*(double)i);
		double z1 = cos(alpha*(double)i);
		double y2 = sin(alpha*(double)((i + 1) % nbDiv));
		double z2 = cos(alpha*(double)((i + 1) % nbDiv));

		glNormal3d(0.0, y1, z1);
		glVertex3d(-0.5, 0.5*y1, 0.5*z1);
		glNormal3d(1.0, 0.0, 0.0);
		glVertex3d(0.5, 0.0, 0.0);
		glNormal3d(0.0, y2, z2);
		glVertex3d(-0.5, 0.5*y2, 0.5*z2);

	}

	// Cap facets
	for (int i = 0; i < nbDiv; i++) {

		double y1 = sin(alpha*(double)i);
		double z1 = cos(alpha*(double)i);
		double y2 = sin(alpha*(double)((i + 1) % nbDiv));
		double z2 = cos(alpha*(double)((i + 1) % nbDiv));

		glNormal3d(-1.0, 0.0, 0.0);
		glVertex3d(-0.5, 0.5*y1, 0.5*z1);
		glNormal3d(-1.0, 0.0, 0.0);
		glVertex3d(-0.5, 0.5*y2, 0.5*z2);
		glNormal3d(-1.0, 0.0, 0.0);
		glVertex3d(-0.5, 0.0, 0.0);

	}

	glEnd();
	glEndList();

	// Shpere list (isotropic case)
	int nbPhi = 16;
	int nbTetha = 7;
	double dphi = 2.0*PI / (double)(nbPhi);
	double dtetha = PI / (double)(nbTetha + 1);

	sphereList = glGenLists(1);
	glNewList(sphereList, GL_COMPILE);

	glEnable(GL_CULL_FACE);
	glCullFace(GL_FRONT);
	glDisable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glEnable(GL_LIGHTING);
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_BLEND);

	glBegin(GL_TRIANGLES);

	for (int i = 0; i <= nbTetha; i++) {
		for (int j = 0; j < nbPhi; j++) {

			VERTEX3D v1, v2, v3, v4;

			v1.x = sin(dtetha*(double)i)*cos(dphi*(double)j);
			v1.y = sin(dtetha*(double)i)*sin(dphi*(double)j);
			v1.z = cos(dtetha*(double)i);

			v2.x = sin(dtetha*(double)(i + 1))*cos(dphi*(double)j);
			v2.y = sin(dtetha*(double)(i + 1))*sin(dphi*(double)j);
			v2.z = cos(dtetha*(double)(i + 1));

			v3.x = sin(dtetha*(double)(i + 1))*cos(dphi*(double)(j + 1));
			v3.y = sin(dtetha*(double)(i + 1))*sin(dphi*(double)(j + 1));
			v3.z = cos(dtetha*(double)(i + 1));

			v4.x = sin(dtetha*(double)i)*cos(dphi*(double)(j + 1));
			v4.y = sin(dtetha*(double)i)*sin(dphi*(double)(j + 1));
			v4.z = cos(dtetha*(double)i);

			if (i < nbTetha) {
				glNormal3d(v1.x, v1.y, v1.z);
				glVertex3d(v1.x, v1.y, v1.z);
				glNormal3d(v2.x, v2.y, v2.z);
				glVertex3d(v2.x, v2.y, v2.z);
				glNormal3d(v3.x, v3.y, v3.z);
				glVertex3d(v3.x, v3.y, v3.z);
			}

			if (i > 0) {
				glNormal3d(v1.x, v1.y, v1.z);
				glVertex3d(v1.x, v1.y, v1.z);
				glNormal3d(v3.x, v3.y, v3.z);
				glVertex3d(v3.x, v3.y, v3.z);
				glNormal3d(v4.x, v4.y, v4.z);
				glVertex3d(v4.x, v4.y, v4.z);
			}

		}
	}

	glEnd();
	glEndList();


}

void Geometry::BuildSelectList() {

	nbSelected = 0;

	selectList = glGenLists(1);
	glNewList(selectList, GL_COMPILE);
	/*
	glDisable(GL_CULL_FACE);
	glDisable(GL_LIGHTING);
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_BLEND);

	if (antiAliasing){
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);	//glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
	//glBlendFunc(GL_ONE,GL_ZERO);
	}
	glLineWidth(2.0f);


	for(int i=0;i<sh.nbFacet;i++ ) {
	Facet *f = facets[i];
	if( f->selected ) {
	//DrawFacet(f,FALSE);
	DrawFacet(f,1,1,1);
	nbSelected++;
	}
	}
	glLineWidth(1.0f);
	if (antiAliasing) {
	glDisable(GL_BLEND);
	glDisable(GL_LINE_SMOOTH);
	}*/
	glEndList();



	// Second list for usage with POLYGON_OFFSET
	selectList2 = glGenLists(1);
	glNewList(selectList2, GL_COMPILE);
	/*
	glDisable(GL_CULL_FACE);
	glDisable(GL_LIGHTING);
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_BLEND);

	if (antiAliasing){
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);	//glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
	}
	glLineWidth(2.0f);

	for(int i=0;i<sh.nbFacet;i++ ) {
	Facet *f = facets[i];
	if( f->selected )
	{
	//DrawFacet(f,TRUE,FALSE,TRUE);
	DrawFacet(f,1,1,1);
	}
	}
	glLineWidth(1.0f);
	if (antiAliasing) {
	glDisable(GL_BLEND);
	glDisable(GL_LINE_SMOOTH);
	}*/
	glEndList();


	// Third list with hidden (hole join) edge visible
	selectList3 = glGenLists(1);
	glNewList(selectList3, GL_COMPILE);

	glDisable(GL_CULL_FACE);
	glDisable(GL_LIGHTING);
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_BLEND);
	if (mApp->antiAliasing){
		glEnable(GL_LINE_SMOOTH);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);	//glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
	}
	glLineWidth(2.0f);

	for (int i = 0; i < sh.nbFacet; i++) {
		Facet *f = facets[i];
		if (f->selected) {
			//DrawFacet(f,FALSE,TRUE,TRUE);
			DrawFacet(f, 1, 1, 1);
			nbSelected++;
		}
	}
	glLineWidth(1.0f);
	if (mApp->antiAliasing) {
		glDisable(GL_BLEND);
		glDisable(GL_LINE_SMOOTH);
	}
	glEndList();


}

void Geometry::UpdateSelection() {

	DeleteGLLists();
	BuildSelectList();

}

void Geometry::BuildGLList() {

	// Compile geometry for OpenGL
	for (int j = 0; j < sh.nbSuper; j++) {
		lineList[j] = glGenLists(1);
		glNewList(lineList[j], GL_COMPILE);
		for (int i = 0; i < sh.nbFacet; i++) {
			if (facets[i]->sh.superIdx == j)
				DrawFacet(facets[i]);
		}
		glEndList();
	}

	polyList = glGenLists(1);
	glNewList(polyList, GL_COMPILE);
	DrawPolys();
	glEndList();

	BuildSelectList();

}

// -----------------------------------------------------------
// Collapse stuff
// -----------------------------------------------------------

Facet *Geometry::MergeFacet(Facet *f1, Facet *f2) {
	mApp->changedSinceSave = TRUE;
	// Merge 2 facets into 1 when possible and create a new facet
	// otherwise return NULL.
	int  c1;
	int  c2;
	int  l;
	BOOL end = FALSE;
	Facet *nF = NULL;

	if (GetCommonEdges(f1, f2, &c1, &c2, &l)) {
		int commonNo = f1->sh.nbIndex + f2->sh.nbIndex - 2 * l;
		if (commonNo == 0) { //two identical facets, so return a copy of f1
			nF = new Facet(f1->sh.nbIndex);
			nF->Copy(f1);
			for (int i = 0; i < f1->sh.nbIndex; i++)
				nF->indices[i] = f1->GetIndex(i);
			return nF;
		}
		int nbI = 0;
		nF = new Facet(commonNo);
		// Copy params from f1
		//nF->Copy(f1);
		nF->Copy(f1);

		if (l == f1->sh.nbIndex) {

			// f1 absorbed, copy indices from f2
			for (int i = 0; i < f2->sh.nbIndex - l; i++)
				nF->indices[nbI++] = f2->GetIndex(c2 + 2 + i);

		}
		else if (l == f2->sh.nbIndex) {

			// f2 absorbed, copy indices from f1
			for (int i = 0; i < f1->sh.nbIndex - l; i++)
				nF->indices[nbI++] = f1->GetIndex(c1 + l + i);

		}
		else {

			// Copy indices from f1
			for (int i = 0; i < f1->sh.nbIndex - (l - 1); i++)
				nF->indices[nbI++] = f1->GetIndex(c1 + l + i);
			// Copy indices from f2
			for (int i = 0; i < f2->sh.nbIndex - (l + 1); i++)
				nF->indices[nbI++] = f2->GetIndex(c2 + 2 + i);

		}

	}

	return nF;

}

void Geometry::Collapse(double vT, double fT, double lT, BOOL doSelectedOnly, GLProgress *prg) {
	mApp->changedSinceSave = TRUE;
	Facet *fi, *fj;
	Facet *merged;


	vThreshold = vT;
	double totalWork = (1.0 + (double)(fT > 0.0) + (double)(lT > 0.0)); //for progress indicator
	// Collapse vertex
	if (vT > 0.0) {
		CollapseVertex(prg, totalWork);
		RemoveCollinear();
		RemoveNullFacet();
		InitializeGeometry();
	}


	if (fT > 0.0) {

		// Collapse facets
		int i = 0;
		prg->SetMessage("Collapsing facets...");
		while (i < sh.nbFacet) {
			prg->SetProgress((1.0 + ((double)i / (double)sh.nbFacet)) / totalWork);
			fi = facets[i];
			// Search a coplanar facet
			int j = i + 1;
			while ((!doSelectedOnly || fi->selected) && j < sh.nbFacet) {
				fj = facets[j];
				merged = NULL;
				if ((!doSelectedOnly || fj->selected) && fi->IsCoplanar(fj, fT)) {
					// Collapse
					merged = MergeFacet(fi, fj);
					if (merged) {
						// Replace the old 2 facets by the new one
						SAFE_DELETE(fi);
						SAFE_DELETE(fj);
						for (int k = j; k < sh.nbFacet - 1; k++)
							facets[k] = facets[k + 1];
						sh.nbFacet--;
						facets[i] = merged;
						//InitializeGeometry(i);
						//SetFacetTexture(i,facets[i]->tRatio,facets[i]->hasMesh);  //rebuild mesh
						fi = facets[i];
						mApp->RenumberSelections(j);
						mApp->RenumberFormulas(j);
						j = i + 1;
					}
				}
				if (!merged) j++;
			}
			i++;
		}
	}
	//Collapse collinear sides. Takes some time, so only if threshold>0
	prg->SetMessage("Collapsing collinear sides...");
	if (lT > 0.0) {
		for (int i = 0; i < sh.nbFacet; i++) {
			prg->SetProgress((1.0 + (double)(fT > 0.0) + ((double)i / (double)sh.nbFacet)) / totalWork);
			if (!doSelectedOnly || facets[i]->selected)
				MergecollinearSides(facets[i], lT);
		}
	}
	prg->SetMessage("Rebuilding geometry...");
	for (int i = 0; i < sh.nbFacet; i++) {

		Facet *f = facets[i];

		// Revert orientation if normal has been swapped
		// This happens when the second vertex is no longer convex
		VERTEX3D n, v1, v2;
		double   d;
		int i0 = facets[i]->indices[0];
		int i1 = facets[i]->indices[1];
		int i2 = facets[i]->indices[2];

		Sub(&v1, vertices3 + i1, vertices3 + i0); // v1 = P0P1
		Sub(&v2, vertices3 + i2, vertices3 + i1); // v2 = P1P2
		Cross(&n, &v1, &v2);                      // Cross product
		d = Dot(&n, &(f->sh.N));
		if (d < 0.0) f->SwapNormal();

	}



	// Delete old resources
	for (int i = 0; i < sh.nbSuper; i++)
		DeleteGLLists(TRUE, TRUE);

	// Reinitialise geom
	InitializeGeometry();

}

void Geometry::MergecollinearSides(Facet *f, double lT) {
	mApp->changedSinceSave = TRUE;
	BOOL collinear;
	double linTreshold = cos(lT*PI / 180);
	// Merge collinear sides
	for (int k = 0; (k < f->sh.nbIndex&&f->sh.nbIndex>3); k++){
		k = k;
		do {
			collinear = FALSE;
			int p0 = f->indices[k];
			int p1 = f->indices[(k + 1) % f->sh.nbIndex];
			int p2 = f->indices[(k + 2) % f->sh.nbIndex]; //to compare last side with first too
			VERTEX3D p0p1;
			VERTEX3D p0p2;
			Sub(&p0p1, &vertices3[p1], &vertices3[p0]);
			Sub(&p0p2, &vertices3[p2], &vertices3[p1]);
			Normalize(&p0p1);
			Normalize(&p0p2);
			collinear = (Dot(&p0p1, &p0p2) >= linTreshold);
			if (collinear&&f->sh.nbIndex > 3) { //collinear
				for (int l = (k + 1) % f->sh.nbIndex; l < f->sh.nbIndex - 1; l++){
					f->indices[l] = f->indices[l + 1];
				}
				f->sh.nbIndex--;
			}
		} while (collinear&&f->sh.nbIndex > 3);
	}
}

void Geometry::Rebuild() {

	// Rebuild internal structure on geometry change

	// Remove texture (improvement TODO)
	for (int i = 0; i < sh.nbFacet; i++)
		if (facets[i]->sh.isTextured)
			facets[i]->SetTexture(0.0, 0.0, FALSE);

	// Delete old resources
	DeleteGLLists(TRUE, TRUE);

	// Reinitialise geom
	InitializeGeometry();

}

int Geometry::InvalidateDeviceObjects() {

	DeleteGLLists(TRUE, TRUE);
	DELETE_LIST(arrowList);
	DELETE_LIST(sphereList);
	for (int i = 0; i < sh.nbFacet; i++)
		facets[i]->InvalidateDeviceObjects();

	return GL_OK;

}

int Geometry::RestoreDeviceObjects() {

	if (!IsLoaded()) return GL_OK;

	for (int i = 0; i < sh.nbFacet; i++) {
		Facet *f = facets[i];
		f->RestoreDeviceObjects();
		BuildFacetList(f);
	}

	BuildGLList();

	return GL_OK;

}

void Geometry::BuildFacetList(Facet *f) {

	// Rebuild OpenGL geomtetry with texture
	if (f->sh.isTextured) {

		// Facet geometry
		glNewList(f->glList, GL_COMPILE);
		if (f->sh.nbIndex == 3) {
			glBegin(GL_TRIANGLES);
			FillFacet(f, TRUE);
			glEnd();
		}
		else if (f->sh.nbIndex == 4) {
			glBegin(GL_QUADS);
			FillFacet(f, TRUE);
			glEnd();
		}
		else {
			glBegin(GL_TRIANGLES);
			Triangulate(f, TRUE);
			glEnd();
		}
		glEndList();

	}

}

void Geometry::SetFacetTexture(int facet, double ratio, BOOL mesh) {

	Facet *f = facets[facet];
	double nU = Norme(&(f->sh.U));
	double nV = Norme(&(f->sh.V));
	if (!f->SetTexture(nU*ratio, nV*ratio, mesh)) {
		char errMsg[512];
		sprintf(errMsg, "Not enough memory to build mesh on Facet %d. ", facet + 1);
		throw Error(errMsg);
	}
	f->tRatio = ratio;
	BuildFacetList(f);

}

// -----------------------------------------------------------
// Testing purpose function, construct a PIPE
// -----------------------------------------------------------
void  Geometry::BuildPipe(double L, double R, double s, int step) {
	Clear();


	//mApp->ClearAllSelections();
	//mApp->ClearAllViews();
	sprintf(sh.name, "PIPE%g", L / R);

	int nbDecade = 0;
	int nbTF = 9 * nbDecade;
	int nbTV = 4 * nbTF;

	sh.nbVertex = 2 * step + nbTV;
	if (!(vertices3 = (VERTEX3D *)malloc(sh.nbVertex * sizeof(VERTEX3D))))
		throw Error("Couldn't allocate memory for vertices");
	memset(vertices3, 0, sh.nbVertex * sizeof(VERTEX3D));

	sh.nbFacet = step + 2 + nbTF;
	sh.nbSuper = 1;
	if (!(facets = (Facet **)malloc(sh.nbFacet * sizeof(Facet *))))
		throw Error("Couldn't allocate memory for facets");
	memset(facets, 0, sh.nbFacet * sizeof(Facet *));
	// Vertices
	for (int i = 0; i < step; i++) {
		double angle = (double)i / (double)step * 2 * PI;
		vertices3[2 * i + nbTV].x = R*cos(angle);
		vertices3[2 * i + nbTV].y = R*sin(angle);
		vertices3[2 * i + nbTV].z = 0.0;
		vertices3[2 * i + 1 + nbTV].x = R*cos(angle);
		vertices3[2 * i + 1 + nbTV].y = R*sin(angle);
		vertices3[2 * i + 1 + nbTV].z = L;
	}

	try {
		// Cap facet
		facets[0 + nbTF] = new Facet(step);
		facets[0 + nbTF]->sh.sticking = 1.0;
		facets[0 + nbTF]->sh.desorbType = DES_COSINE;
		facets[0 + nbTF]->sh.flow = 1.0;
		for (int i = 0; i < step; i++)
			facets[0 + nbTF]->indices[i] = 2 * i + nbTV;

		facets[1 + nbTF] = new Facet(step);
		facets[1 + nbTF]->sh.sticking = 1.0;
		facets[1 + nbTF]->sh.desorbType = DES_NONE;
		for (int i = 0; i < step; i++)
			facets[1 + nbTF]->indices[step - i - 1] = 2 * i + 1 + nbTV;

		// Wall facet
		for (int i = 0; i < step; i++) {
			facets[i + 2 + nbTF] = new Facet(4);
			facets[i + 2 + nbTF]->sh.reflectType = REF_DIFFUSE;
			facets[i + 2 + nbTF]->sh.sticking = s;
			facets[i + 2 + nbTF]->indices[0] = 2 * i + nbTV;
			facets[i + 2 + nbTF]->indices[1] = 2 * i + 1 + nbTV;
			if (i < step - 1) {
				facets[i + 2 + nbTF]->indices[2] = 2 * (i + 1) + 1 + nbTV;
				facets[i + 2 + nbTF]->indices[3] = 2 * (i + 1) + nbTV;
			}
			else {
				facets[i + 2 + nbTF]->indices[2] = 1 + nbTV;
				facets[i + 2 + nbTF]->indices[3] = 0 + nbTV;
			}
		}

		// Volatile facet
		for (int d = 0; d < nbDecade; d++) {
			for (int i = 0; i < 9; i++) {

				double z = (double)(i + 1) * pow(10, (double)d);
				int idx = d * 36 + i * 4;

				vertices3[idx + 0].x = -R;
				vertices3[idx + 0].y = R;
				vertices3[idx + 0].z = z;
				vertices3[idx + 1].x = R;
				vertices3[idx + 1].y = R;
				vertices3[idx + 1].z = z;
				vertices3[idx + 2].x = R;
				vertices3[idx + 2].y = -R;
				vertices3[idx + 2].z = z;
				vertices3[idx + 3].x = -R;
				vertices3[idx + 3].y = -R;
				vertices3[idx + 3].z = z;

				facets[9 * d + i] = new Facet(4);
				facets[9 * d + i]->sh.sticking = 0.0;
				facets[9 * d + i]->sh.opacity = 0.0;
				facets[9 * d + i]->sh.isVolatile = TRUE;
				facets[9 * d + i]->indices[0] = idx + 0;
				facets[9 * d + i]->indices[1] = idx + 1;
				facets[9 * d + i]->indices[2] = idx + 2;
				facets[9 * d + i]->indices[3] = idx + 3;

			}
		}
	}
	catch (std::bad_alloc) {
		Clear();
		throw Error("Couldn't reserve memory for the facets");
	}
	catch (...) {
		throw Error("Unspecified Error while building pipe");
	}
	InitializeGeometry();
	CalcTotalOutGassing();
	isLoaded = TRUE;
	strName[0] = _strdup("Pipe");
	strFileName[0] = _strdup("pipe.txt");

}

// -----------------------------------------------------------
// File handling
// -----------------------------------------------------------

void Geometry::UpdateName(FileReader *file) {

	char *p = strrchr(file->GetName(), '\\');
	if (!p) p = strrchr(file->GetName(), '/');

	if (!p)
		strcpy(sh.name, file->GetName());
	else
		strcpy(sh.name, p + 1);
}

void Geometry::UpdateName(char *fileName) {

	char *p = strrchr(fileName, '\\');
	if (!p) p = strrchr(fileName, '/');

	if (!p)
		strcpy(sh.name, fileName);
	else
		strcpy(sh.name, p + 1);
}

void Geometry::AdjustProfile() {

	// Backward compatibily with TXT profile (To be improved)
	for (int i = 0; i < sh.nbFacet; i++) {
		Facet *f = facets[i];
		if (f->sh.profileType == REC_PRESSUREU) {
			VERTEX3D v0;
			Sub(&v0, vertices3 + (f->indices[1]), vertices3 + (f->indices[0])); // v0 = P0P1
			double n0 = Norme(&v0);
			double nU = Norme(&(f->sh.U));
			if (IS_ZERO(n0 - nU)) f->sh.profileType = REC_PRESSUREU; // Select U
			else f->sh.profileType = REC_PRESSUREV; // Select V
		}
	}

}

void Geometry::LoadASE(FileReader *file, GLProgress *prg) {

	Clear();

	//mApp->ClearAllSelections();
	//mApp->ClearAllViews();
	ASELoader ase(file);
	ase.Load();

	// Compute total of facet
	sh.nbFacet = 0;
	for (int i = 0; i < ase.nbObj; i++) sh.nbFacet += ase.OBJ[i].nb_face;

	// Allocate mem
	sh.nbVertex = 3 * sh.nbFacet;
	facets = (Facet **)malloc(sh.nbFacet * sizeof(Facet *));
	memset(facets, 0, sh.nbFacet * sizeof(Facet *));
	vertices3 = (VERTEX3D *)malloc(sh.nbVertex * sizeof(VERTEX3D));
	memset(vertices3, 0, sh.nbVertex * sizeof(VERTEX3D));

	// Fill 
	int nb = 0;
	for (int i = 0; i < ase.nbObj; i++) {

		for (int j = 0; j < ase.OBJ[i].nb_face; j++) {
			vertices3[3 * nb + 0] = ase.OBJ[i].pts[ase.OBJ[i].face[j].v1];
			vertices3[3 * nb + 1] = ase.OBJ[i].pts[ase.OBJ[i].face[j].v2];
			vertices3[3 * nb + 2] = ase.OBJ[i].pts[ase.OBJ[i].face[j].v3];
			facets[nb] = new Facet(3);
			facets[nb]->indices[0] = 3 * nb + 0;
			facets[nb]->indices[1] = 3 * nb + 1;
			facets[nb]->indices[2] = 3 * nb + 2;
			nb++;
		}

	}

	UpdateName(file);
	sh.nbSuper = 1;
	strName[0] = _strdup(sh.name);
	strFileName[0] = _strdup(file->GetName());
	char *e = strrchr(strName[0], '.');
	if (e) *e = 0;
	InitializeGeometry();
	isLoaded = TRUE;

}

void Geometry::LoadSTR(FileReader *file, GLProgress *prg) {

	char nPath[512];
	char fPath[512];
	char fName[512];
	char sName[512];
	int nF, nV;
	Facet **F;
	VERTEX3D *V;
	FileReader *fr;

	Clear();

	//mApp->ClearAllSelections();
	//mApp->ClearAllViews();
	// Load multiple structure file
	sh.nbSuper = file->ReadInt();

	strcpy(fPath, file->ReadLine());
	strcpy(nPath, FileUtils::GetPath(file->GetName()));

	for (int n = 0; n < sh.nbSuper; n++) {

		int i1 = file->ReadInt();
		int i2 = file->ReadInt();
		fr = NULL;
		strcpy(sName, file->ReadWord());
		strName[n] = _strdup(sName);
		char *e = strrchr(strName[n], '.');
		if (e) *e = 0;

		sprintf(fName, "%s\\%s", nPath, sName);
		if (FileUtils::Exist(fName)) {
			fr = new FileReader(fName);
			strcpy(strPath, nPath);
		}
		else {
			sprintf(fName, "%s\\%s", fPath, sName);
			if (FileUtils::Exist(fName)) {
				fr = new FileReader(fName);
				strcpy(strPath, fPath);
			}
		}

		if (!fr) {
			char errMsg[512];
			sprintf(errMsg, "Cannot find %s", sName);
			throw Error(errMsg);
		}

		strFileName[n] = _strdup(sName);
		LoadTXTGeom(fr, &nV, &nF, &V, &F, n);
		Merge(nV, nF, V, F);
		SAFE_FREE(V);
		SAFE_FREE(F);
		delete fr;

	}

	UpdateName(file);
	InitializeGeometry();
	AdjustProfile();
	isLoaded = TRUE;

}

void Geometry::LoadSTL(FileReader *file, GLProgress *prg, double scaleFactor) {


	//mApp->ClearAllSelections();
	//mApp->ClearAllViews();
	char *w;

	prg->SetMessage("Clearing current geometry...");
	Clear();

	// First pass
	prg->SetMessage("Counting facets in STL file...");
	//file->ReadKeyword("solid");
	file->ReadLine(); // solid name
	w = file->ReadWord();
	while (strcmp(w, "facet") == 0) {
		sh.nbFacet++;
		file->JumpSection("endfacet");
		w = file->ReadWord();
	}
	if (strcmp(w, "endsolid") != 0) throw Error("Unexpected or not supported STL keyword, 'endsolid' required");

	// Allocate mem
	sh.nbVertex = 3 * sh.nbFacet;
	facets = (Facet **)malloc(sh.nbFacet * sizeof(Facet *));
	memset(facets, 0, sh.nbFacet * sizeof(Facet *));
	vertices3 = (VERTEX3D *)malloc(sh.nbVertex * sizeof(VERTEX3D));
	memset(vertices3, 0, sh.nbVertex * sizeof(VERTEX3D));

	// Second pass
	prg->SetMessage("Reading facets...");
	file->SeekStart();
	//file->ReadKeyword("solid");
	file->ReadLine();
	for (int i = 0; i < sh.nbFacet; i++) {

		double p = (double)i / (double)(sh.nbFacet);
		prg->SetProgress(p);

		file->ReadKeyword("facet");
		file->ReadKeyword("normal");
		file->ReadDouble();
		file->ReadDouble();
		file->ReadDouble();
		file->ReadKeyword("outer");
		file->ReadKeyword("loop");

		file->ReadKeyword("vertex");
		vertices3[3 * i + 0].x = file->ReadDouble()*scaleFactor;
		vertices3[3 * i + 0].y = file->ReadDouble()*scaleFactor;
		vertices3[3 * i + 0].z = file->ReadDouble()*scaleFactor;

		file->ReadKeyword("vertex");
		vertices3[3 * i + 1].x = file->ReadDouble()*scaleFactor;
		vertices3[3 * i + 1].y = file->ReadDouble()*scaleFactor;
		vertices3[3 * i + 1].z = file->ReadDouble()*scaleFactor;

		file->ReadKeyword("vertex");
		vertices3[3 * i + 2].x = file->ReadDouble()*scaleFactor;
		vertices3[3 * i + 2].y = file->ReadDouble()*scaleFactor;
		vertices3[3 * i + 2].z = file->ReadDouble()*scaleFactor;

		file->ReadKeyword("endloop");
		file->ReadKeyword("endfacet");

		try {
			facets[i] = new Facet(3);
		}
		catch (std::bad_alloc& badalloc) {

			throw Error("Out of memory");

		}
		facets[i]->indices[0] = 3 * i + 0;
		facets[i]->indices[1] = 3 * i + 2;
		facets[i]->indices[2] = 3 * i + 1;

	}

	sh.nbSuper = 1;
	UpdateName(file);
	strName[0] = _strdup(sh.name);
	strFileName[0] = _strdup(file->GetName());
	char *e = strrchr(strName[0], '.');
	if (e) *e = 0;
	prg->SetMessage("Initializing geometry...");
	InitializeGeometry();
	isLoaded = TRUE;

}

void Geometry::LoadTXT(FileReader *file, GLProgress *prg) {


	//mApp->ClearAllSelections();
	//mApp->ClearAllViews();
	Clear();
	UpdateName(file);
	sh.nbSuper = 1;
	strName[0] = _strdup(sh.name);
	strFileName[0] = _strdup(file->GetName());
	LoadTXTGeom(file, &(sh.nbVertex), &(sh.nbFacet), &vertices3, &facets);
	char *e = strrchr(strName[0], '.');
	if (e) *e = 0;
	InitializeGeometry();
	AdjustProfile();
	isLoaded = TRUE;

}

void Geometry::InsertTXT(FileReader *file, GLProgress *prg, BOOL newStr) {

	//Clear();
	int structId = viewStruct;
	if (structId == -1) structId = 0;
	InsertTXTGeom(file, &(sh.nbVertex), &(sh.nbFacet), &vertices3, &facets, structId, newStr);
	//UpdateName(file);
	//sh.nbSuper = 1;
	//strName[0] = _strdup(sh.name);
	//strFileName[0] = _strdup(file->GetName());
	char *e = strrchr(strName[0], '.');
	if (e) *e = 0;
	InitializeGeometry();
	AdjustProfile();
	isLoaded = TRUE;

}

void Geometry::InsertSTL(FileReader *file, GLProgress *prg, double scaleFactor, BOOL newStr) {

	//Clear();
	int structId = viewStruct;
	if (structId == -1) structId = 0;
	InsertSTLGeom(file, &(sh.nbVertex), &(sh.nbFacet), &vertices3, &facets, structId, scaleFactor, newStr);
	//UpdateName(file);
	//sh.nbSuper = 1;
	//strName[0] = _strdup(sh.name);
	//strFileName[0] = _strdup(file->GetName());
	char *e = strrchr(strName[0], '.');
	if (e) *e = 0;
	InitializeGeometry();
	//AdjustProfile();
	isLoaded = TRUE;

}

void Geometry::InsertGEO(FileReader *file, GLProgress *prg, BOOL newStr) {

	//Clear();
	int structId = viewStruct;
	if (structId == -1) structId = 0;
	InsertGEOGeom(file, &(sh.nbVertex), &(sh.nbFacet), &vertices3, &facets, structId, newStr);
	//UpdateName(file);
	//sh.nbSuper = 1;
	//strName[0] = _strdup(sh.name);
	//strFileName[0] = _strdup(file->GetName());
	char *e = strrchr(strName[0], '.');
	if (e) *e = 0;
	InitializeGeometry();
	//AdjustProfile();
	isLoaded = TRUE;

}

void Geometry::InsertSYN(FileReader *file, GLProgress *prg, BOOL newStr) {

	int structId = viewStruct;
	if (structId == -1) structId = 0;
	InsertSYNGeom(file, &(sh.nbVertex), &(sh.nbFacet), &vertices3, &facets, structId, newStr);
	char *e = strrchr(strName[0], '.');
	if (e) *e = 0;
	InitializeGeometry();
	//AdjustProfile();
}

void Geometry::LoadTXTGeom(FileReader *file, int *nbV, int *nbF, VERTEX3D **V, Facet ***F, int strIdx) {

	file->ReadInt(); // Unused
	tNbHit = file->ReadLLong();
	tNbLeak = file->ReadLLong();
	tNbDesorption = file->ReadLLong();
	tNbDesorptionMax = file->ReadLLong();

	int nV = file->ReadInt();
	int nF = file->ReadInt();

	// Allocate memory
	Facet   **f = (Facet **)malloc(nF * sizeof(Facet *));
	memset(f, 0, nF * sizeof(Facet *));
	VERTEX3D *v = (VERTEX3D *)malloc(nV * sizeof(VERTEX3D));

	// Read geometry vertices
	for (int i = 0; i < nV; i++) {
		v[i].x = file->ReadDouble();
		v[i].y = file->ReadDouble();
		v[i].z = file->ReadDouble();
	}

	// Read geometry facets (indexed from 1)
	for (int i = 0; i < nF; i++) {
		int nb = file->ReadInt();
		f[i] = new Facet(nb);
		for (int j = 0; j < nb; j++)
			f[i]->indices[j] = file->ReadInt() - 1;
	}

	// Read facets params
	for (int i = 0; i < nF; i++) {
		f[i]->LoadTXT(file);
		while ((f[i]->sh.superDest) > sh.nbSuper) { //If facet refers to a structure that doesn't exist, create it
				AddStruct("TXT linked");
		}
		f[i]->sh.superIdx = strIdx;
	}

	SAFE_FREE(*V);
	SAFE_FREE(*F);

	*nbV = nV;
	*nbF = nF;
	*V = v;
	*F = f;

}

// -----------------------------------------------------------
// -----------------------------------------------------------

void Geometry::InsertTXTGeom(FileReader *file, int *nbVertex, int *nbFacet, VERTEX3D **vertices3, Facet ***facets, int strIdx, BOOL newStruct) {

	UnSelectAll();

	//tNbHit = file->ReadLLong();
	//tNbLeak = file->ReadInt();
	//tNbDesorption = file->ReadLLong();
	//tNbDesorptionMax = file->ReadLLong(); 
	for (int i = 0; i < 5; i++) file->ReadInt(); //leading lines

	int nbNewVertex = file->ReadInt();
	int nbNewFacets = file->ReadInt();

	// Allocate memory
	*facets = (Facet **)realloc(*facets, (nbNewFacets + *nbFacet) * sizeof(Facet **));
	memset(*facets + *nbFacet, 0, nbNewFacets * sizeof(Facet *));
	//*vertices3 = (VERTEX3D*)realloc(*vertices3,(nbNewVertex+*nbVertex) * sizeof(VERTEX3D));
	VERTEX3D *tmp_vertices3 = (VERTEX3D *)malloc((nbNewVertex + *nbVertex) * sizeof(VERTEX3D));
	memmove(tmp_vertices3, *vertices3, (*nbVertex)*sizeof(VERTEX3D));
	memset(tmp_vertices3 + *nbVertex, 0, nbNewVertex * sizeof(VERTEX3D));
	SAFE_FREE(*vertices3);
	*vertices3 = tmp_vertices3;

	// Read geometry vertices
	for (int i = *nbVertex; i < (*nbVertex + nbNewVertex); i++) {
		(*vertices3 + i)->x = file->ReadDouble();
		(*vertices3 + i)->y = file->ReadDouble();
		(*vertices3 + i)->z = file->ReadDouble();
		(*vertices3 + i)->selected = FALSE;
	}

	// Read geometry facets (indexed from 1)
	for (int i = *nbFacet; i < (*nbFacet + nbNewFacets); i++) {
		int nb = file->ReadInt();
		*(*facets + i) = new Facet(nb);
		(*facets)[i]->selected = TRUE;
		for (int j = 0; j < nb; j++)
			(*facets)[i]->indices[j] = file->ReadInt() - 1 + *nbVertex;
	}

	// Read facets params
	for (int i = *nbFacet; i < (*nbFacet + nbNewFacets); i++) {
		(*facets)[i]->LoadTXT(file);
		while (((*facets)[i]->sh.superDest) > sh.nbSuper) { //If facet refers to a structure that doesn't exist, create it
			AddStruct("TXT linked");
		}
		if (newStruct) {
			(*facets)[i]->sh.superIdx = sh.nbSuper;
		}
		else {
			(*facets)[i]->sh.superIdx = strIdx;
		}
	}

	*nbVertex += nbNewVertex;
	*nbFacet += nbNewFacets;
	if (newStruct) AddStruct("Inserted TXT file");

}

void Geometry::InsertGEOGeom(FileReader *file, int *nbVertex, int *nbFacet, VERTEX3D **vertices3, Facet ***facets, int strIdx, BOOL newStruct) {


	UnSelectAll();
	char tmp[512];

	file->ReadKeyword("version"); file->ReadKeyword(":");
	int version2;
	version2 = file->ReadInt();
	if (version2 > GEOVERSION) {
		char errMsg[512];
		sprintf(errMsg, "Unsupported GEO version V%d", version2);
		throw Error(errMsg);
	}

	file->ReadKeyword("totalHit"); file->ReadKeyword(":");
	file->ReadLLong();
	file->ReadKeyword("totalDes"); file->ReadKeyword(":");
	file->ReadLLong();
	file->ReadKeyword("totalLeak"); file->ReadKeyword(":");
	file->ReadLLong();
	if (version2 >= 12) {
		file->ReadKeyword("totalAbs"); file->ReadKeyword(":");
		file->ReadLLong();
		if (version2 >= 15) {
			file->ReadKeyword("totalDist_total");
		}
		else { //between versions 12 and 15
			file->ReadKeyword("totalDist");
		}
		file->ReadKeyword(":");
		file->ReadDouble();
		if (version2 >= 15) {
			file->ReadKeyword("totalDist_fullHitsOnly"); file->ReadKeyword(":");
			file->ReadDouble();
		}
	}
	file->ReadKeyword("maxDes"); file->ReadKeyword(":");
	file->ReadLLong();
	file->ReadKeyword("nbVertex"); file->ReadKeyword(":");
	int nbNewVertex = file->ReadInt();
	file->ReadKeyword("nbFacet"); file->ReadKeyword(":");
	int nbNewFacets = file->ReadInt();
	file->ReadKeyword("nbSuper"); file->ReadKeyword(":");
	int nbNewSuper = file->ReadInt();
	int nbF = 0;
	int nbV = 0;
	if (version2 >= 2) {
		file->ReadKeyword("nbFormula"); file->ReadKeyword(":");
		nbF = file->ReadInt();
		file->ReadKeyword("nbView"); file->ReadKeyword(":");
		nbV = file->ReadInt();
	}
	int nbS = 0;
	if (version2 >= 8) {
		file->ReadKeyword("nbSelection"); file->ReadKeyword(":");
		nbS = file->ReadInt();
	}
	if (version2 >= 7) {
		file->ReadKeyword("gasMass"); file->ReadKeyword(":");
		/*gasMass = */file->ReadDouble();
	}
	if (version2 >= 10) { //time-dependent version
		file->ReadKeyword("userMoments"); file->ReadKeyword("{");
		file->ReadKeyword("nb"); file->ReadKeyword(":");
		int nb = file->ReadInt();

		for (int i = 0; i < nb; i++)
			file->ReadString();
		file->ReadKeyword("}");
	}
	if (version2 >= 11) { //gas pulse parameters
		file->ReadKeyword("desorptionStart"); file->ReadKeyword(":");
		file->ReadDouble();
		file->ReadKeyword("desorptionStop"); file->ReadKeyword(":");
		file->ReadDouble();
		file->ReadKeyword("timeWindow"); file->ReadKeyword(":");
		file->ReadDouble();
		file->ReadKeyword("useMaxwellian"); file->ReadKeyword(":");
		file->ReadInt();
	}

	if (version2 >= 12) { //2013.aug.22
		file->ReadKeyword("calcConstantFlow"); file->ReadKeyword(":");
		file->ReadInt();
	}

	if (version2 >= 2) {
		file->ReadKeyword("formulas"); file->ReadKeyword("{");
		for (int i = 0; i < nbF; i++) {
			char tmpName[256];
			char tmpExpr[512];
			strcpy(tmpName, file->ReadString());
			strcpy(tmpExpr, file->ReadString());
			mApp->OffsetFormula(tmpExpr, sh.nbFacet);
			mApp->AddFormula(tmpName, tmpExpr);
		}
		file->ReadKeyword("}");

		file->ReadKeyword("views"); file->ReadKeyword("{");
		for (int i = 0; i < nbV; i++) {
			char tmpName[256];
			AVIEW v;
			strcpy(tmpName, file->ReadString());
			v.projMode = file->ReadInt();
			v.camAngleOx = file->ReadDouble();
			v.camAngleOy = file->ReadDouble();
			v.camDist = file->ReadDouble();
			v.camOffset.x = file->ReadDouble();
			v.camOffset.y = file->ReadDouble();
			v.camOffset.z = file->ReadDouble();
			v.performXY = file->ReadInt();

			v.vLeft = file->ReadDouble();
			v.vRight = file->ReadDouble();
			v.vTop = file->ReadDouble();
			v.vBottom = file->ReadDouble();
			mApp->AddView(tmpName, v);
		}
		file->ReadKeyword("}");
	}

	if (version2 >= 8) {
		file->ReadKeyword("selections"); file->ReadKeyword("{");
		for (int i = 0; i < nbS; i++) {
			ASELECTION s;
			char tmpName[256];
			strcpy(tmpName, file->ReadString());
			s.name = _strdup(tmpName);
			s.nbSel = file->ReadInt();
			s.selection = (int *)malloc((s.nbSel)*sizeof(int));

			for (int j = 0; j < s.nbSel; j++) {
				s.selection[j] = file->ReadInt() + sh.nbFacet; //offset facet number by current number of facets
			}
			mApp->AddSelection(s.name, s);
		}
		file->ReadKeyword("}");
	}

	file->ReadKeyword("structures"); file->ReadKeyword("{");
	for (int i = 0; i < nbNewSuper; i++) {
		strName[i] = _strdup(file->ReadString());
		// For backward compatibilty with STR
		/* //Commented out for GEO
		sprintf(tmp, "%s.txt", strName[i]);
		strFileName[i] = _strdup(tmp);
		*/
	}
	file->ReadKeyword("}");

	// Reallocate memory
	*facets = (Facet **)realloc(*facets, (nbNewFacets + *nbFacet) * sizeof(Facet **));
	memset(*facets + *nbFacet, 0, nbNewFacets * sizeof(Facet *));
	//*vertices3 = (VERTEX3D*)realloc(*vertices3,(nbNewVertex+*nbVertex) * sizeof(VERTEX3D));
	VERTEX3D *tmp_vertices3 = (VERTEX3D *)malloc((nbNewVertex + *nbVertex) * sizeof(VERTEX3D));
	memmove(tmp_vertices3, *vertices3, (*nbVertex)*sizeof(VERTEX3D));
	memset(tmp_vertices3 + *nbVertex, 0, nbNewVertex * sizeof(VERTEX3D));
	SAFE_FREE(*vertices3);
	*vertices3 = tmp_vertices3;

	// Read geometry vertices
	file->ReadKeyword("vertices"); file->ReadKeyword("{");
	for (int i = *nbVertex; i < (*nbVertex + nbNewVertex); i++) {
		// Check idx
		int idx = file->ReadInt();
		if (idx != i - *nbVertex + 1) throw Error(file->MakeError("Wrong vertex index !"));
		(*vertices3 + i)->x = file->ReadDouble();
		(*vertices3 + i)->y = file->ReadDouble();
		(*vertices3 + i)->z = file->ReadDouble();
		(*vertices3 + i)->selected = FALSE;
	}
	file->ReadKeyword("}");

	if (version2 >= 6) {
		// Read leaks
		file->ReadKeyword("leaks"); file->ReadKeyword("{");
		file->ReadKeyword("nbLeak"); file->ReadKeyword(":");
		int nbleak2 = file->ReadInt();
		for (int i = 0; i < nbleak2; i++) {
			int idx = file->ReadInt();
			//if( idx != i ) throw Error(file->MakeError("Wrong leak index !"));
			file->ReadDouble();
			file->ReadDouble();
			file->ReadDouble();

			file->ReadDouble();
			file->ReadDouble();
			file->ReadDouble();
		}
		file->ReadKeyword("}");

		// Read hit cache
		file->ReadKeyword("hits"); file->ReadKeyword("{");
		file->ReadKeyword("nbHHit"); file->ReadKeyword(":");
		int nbHHit2 = file->ReadInt();
		for (int i = 0; i < nbHHit2; i++) {
			int idx = file->ReadInt();
			//if( idx != i ) throw Error(file->MakeError("Wrong hit cache index !"));
			file->ReadDouble();
			file->ReadDouble();
			file->ReadDouble();

			file->ReadInt();
		}
		file->ReadKeyword("}");
	}

	// Read geometry facets (indexed from 1)
	for (int i = *nbFacet; i < (*nbFacet + nbNewFacets); i++) {
		file->ReadKeyword("facet");
		// Check idx
		int idx = file->ReadInt();
		if (idx != i + 1 - *nbFacet) throw Error(file->MakeError("Wrong facet index !"));
		file->ReadKeyword("{");
		file->ReadKeyword("nbIndex");
		file->ReadKeyword(":");
		int nb = file->ReadInt();

		if (nb < 3) {
			char errMsg[512];
			sprintf(errMsg, "Facet %d has only %d vertices. ", i, nb);
			throw Error(errMsg);
		}

		*(*facets + i) = new Facet(nb);
		(*facets)[i]->LoadGEO(file, version2, nbNewVertex);
		(*facets)[i]->selected = TRUE;
		for (int j = 0; j < nb; j++)
			(*facets)[i]->indices[j] += *nbVertex;
		file->ReadKeyword("}");
		if (newStruct) {
			(*facets)[i]->sh.superIdx += sh.nbSuper;
			if ((*facets)[i]->sh.superDest > 0) (*facets)[i]->sh.superDest += sh.nbSuper;
		}
		else {
			(*facets)[i]->sh.superIdx += strIdx;
			if ((*facets)[i]->sh.superDest > 0) (*facets)[i]->sh.superDest += strIdx;
		}
	}


	*nbVertex += nbNewVertex;
	*nbFacet += nbNewFacets;
	if (newStruct) sh.nbSuper += nbNewSuper;
	else if (sh.nbSuper < strIdx + nbNewSuper) sh.nbSuper = strIdx + nbNewSuper;

}

void Geometry::InsertSYNGeom(FileReader *file, int *nbVertex, int *nbFacet, VERTEX3D **vertices3, Facet ***facets, int strIdx, BOOL newStruct) {


	UnSelectAll();
	char tmp[512];

	file->ReadKeyword("version"); file->ReadKeyword(":");
	int version2;
	version2 = file->ReadInt();
	if (version2 > SYNVERSION) {
		char errMsg[512];
		sprintf(errMsg, "Unsupported SYN version V%d", version2);
		throw Error(errMsg);
	}

	file->ReadKeyword("totalHit"); file->ReadKeyword(":");
	file->ReadLLong();
	file->ReadKeyword("totalDes"); file->ReadKeyword(":");
	file->ReadLLong();
	if (version2 >= 6) {
		file->ReadKeyword("no_scans"); file->ReadKeyword(":");
		/*loaded_no_scans = */file->ReadDouble();
	}
	file->ReadKeyword("totalLeak"); file->ReadKeyword(":");
	file->ReadLLong();
	if (version2 > 2) {
		file->ReadKeyword("totalFlux"); file->ReadKeyword(":");
		file->ReadDouble();
		file->ReadKeyword("totalPower"); file->ReadKeyword(":");
		file->ReadDouble();
	}
	file->ReadKeyword("maxDes"); file->ReadKeyword(":");
	file->ReadLLong();
	file->ReadKeyword("nbVertex"); file->ReadKeyword(":");
	int nbNewVertex = file->ReadInt();
	file->ReadKeyword("nbFacet"); file->ReadKeyword(":");
	int nbNewFacets = file->ReadInt();
	file->ReadKeyword("nbSuper"); file->ReadKeyword(":");
	int nbNewSuper = file->ReadInt();
	int nbF = 0;
	int nbV = 0;

	file->ReadKeyword("nbFormula"); file->ReadKeyword(":");
	nbF = file->ReadInt();
	file->ReadKeyword("nbView"); file->ReadKeyword(":");
	nbV = file->ReadInt();
	int nbS = 0;

	file->ReadKeyword("nbSelection"); file->ReadKeyword(":");
	nbS = file->ReadInt();
	if (version2 > 1) {
		file->ReadKeyword("nbRegions"); file->ReadKeyword(":");
		int nbR = file->ReadInt();
		//result=PARfileList(nbR);

		file->ReadKeyword("PARfiles"); file->ReadKeyword("{");
		for (int i = 0; i < nbR; i++) {
			/*char tmp[512];
			strcpy(tmp,file->ReadString());
			result.fileNames[i]=_strdup(tmp);*/
			file->ReadString();
		}
		file->ReadKeyword("}");
	}

	file->ReadKeyword("formulas"); file->ReadKeyword("{");
	for (int i = 0; i < nbF; i++) {
		char tmpName[256];
		char tmpExpr[512];
		strcpy(tmpName, file->ReadString());
		strcpy(tmpExpr, file->ReadString());
		//mApp->AddFormula(tmpName,tmpExpr); //we don't add SynRad formulas to MolFlow

	}
	file->ReadKeyword("}");

	file->ReadKeyword("views"); file->ReadKeyword("{");
	for (int i = 0; i < nbV; i++) {
		char tmpName[256];
		AVIEW v;
		strcpy(tmpName, file->ReadString());
		v.projMode = file->ReadInt();
		v.camAngleOx = file->ReadDouble();
		v.camAngleOy = file->ReadDouble();
		v.camDist = file->ReadDouble();
		v.camOffset.x = file->ReadDouble();
		v.camOffset.y = file->ReadDouble();
		v.camOffset.z = file->ReadDouble();
		v.performXY = file->ReadInt();

		v.vLeft = file->ReadDouble();
		v.vRight = file->ReadDouble();
		v.vTop = file->ReadDouble();
		v.vBottom = file->ReadDouble();
		mApp->AddView(tmpName, v);
	}
	file->ReadKeyword("}");

	file->ReadKeyword("selections"); file->ReadKeyword("{");
	for (int i = 0; i < nbS; i++) {
		ASELECTION s;
		char tmpName[256];
		strcpy(tmpName, file->ReadString());
		s.name = _strdup(tmpName);
		s.nbSel = file->ReadInt();
		s.selection = (int *)malloc((s.nbSel)*sizeof(int));

		for (int j = 0; j < s.nbSel; j++) {
			s.selection[j] = file->ReadInt() + sh.nbFacet;
		}
		mApp->AddSelection(s.name, s);
	}
	file->ReadKeyword("}");

	file->ReadKeyword("structures"); file->ReadKeyword("{");
	for (int i = 0; i < nbNewSuper; i++) {
		strName[i] = _strdup(file->ReadString());
		// For backward compatibilty with STR
		/*sprintf(tmp, "%s.txt", strName[i]);
		strFileName[i] = _strdup(tmp);*/
	}
	file->ReadKeyword("}");

	// Reallocate memory
	*facets = (Facet **)realloc(*facets, (nbNewFacets + *nbFacet) * sizeof(Facet **));
	memset(*facets + *nbFacet, 0, nbNewFacets * sizeof(Facet *));
	//*vertices3 = (VERTEX3D*)realloc(*vertices3,(nbNewVertex+*nbVertex) * sizeof(VERTEX3D));
	VERTEX3D *tmp_vertices3 = (VERTEX3D *)malloc((nbNewVertex + *nbVertex) * sizeof(VERTEX3D));
	memmove(tmp_vertices3, *vertices3, (*nbVertex)*sizeof(VERTEX3D));
	memset(tmp_vertices3 + *nbVertex, 0, nbNewVertex * sizeof(VERTEX3D));
	SAFE_FREE(*vertices3);
	*vertices3 = tmp_vertices3;

	// Read geometry vertices
	file->ReadKeyword("vertices"); file->ReadKeyword("{");
	for (int i = *nbVertex; i < (*nbVertex + nbNewVertex); i++) {
		// Check idx
		int idx = file->ReadInt();
		if (idx != i - *nbVertex + 1) throw Error(file->MakeError("Wrong vertex index !"));
		(*vertices3 + i)->x = file->ReadDouble();
		(*vertices3 + i)->y = file->ReadDouble();
		(*vertices3 + i)->z = file->ReadDouble();
		(*vertices3 + i)->selected = FALSE;
	}
	file->ReadKeyword("}");


	// Read leaks
	file->ReadKeyword("leaks"); file->ReadKeyword("{");
	file->ReadKeyword("nbLeak"); file->ReadKeyword(":");
	int nbleak_local = file->ReadInt();
	for (int i = 0; i < nbleak_local; i++) {
		int idx = file->ReadInt();
		//if( idx != i ) throw Error(file->MakeError("Wrong leak index !"));
		file->ReadDouble();
		file->ReadDouble();
		file->ReadDouble();

		file->ReadDouble();
		file->ReadDouble();
		file->ReadDouble();
	}
	file->ReadKeyword("}");

	// Read hit cache
	file->ReadKeyword("hits"); file->ReadKeyword("{");
	file->ReadKeyword("nbHHit"); file->ReadKeyword(":");
	int nbHHit_local = file->ReadInt();
	for (int i = 0; i < nbHHit_local; i++) {
		int idx = file->ReadInt();
		//if( idx != i ) throw Error(file->MakeError("Wrong hit cache index !"));
		file->ReadDouble();
		file->ReadDouble();
		file->ReadDouble();
		file->ReadDouble();
		file->ReadDouble();
		file->ReadInt();
	}
	file->ReadKeyword("}");


	// Read geometry facets (indexed from 1)
	for (int i = *nbFacet; i < (*nbFacet + nbNewFacets); i++) {
		file->ReadKeyword("facet");
		// Check idx
		int idx = file->ReadInt();
		if (idx != i + 1 - *nbFacet) throw Error(file->MakeError("Wrong facet index !"));
		file->ReadKeyword("{");
		file->ReadKeyword("nbIndex");
		file->ReadKeyword(":");
		int nb = file->ReadInt();

		if (nb < 3) {
			char errMsg[512];
			sprintf(errMsg, "Facet %d has only %d vertices. ", i, nb);
			throw Error(errMsg);
		}

		*(*facets + i) = new Facet(nb);
		(*facets)[i]->LoadSYN(file, version2, nbNewVertex);
		(*facets)[i]->selected = TRUE;
		for (int j = 0; j < nb; j++)
			(*facets)[i]->indices[j] += *nbVertex;
		file->ReadKeyword("}");
		if (newStruct) {
			(*facets)[i]->sh.superIdx += sh.nbSuper;
		}
		else {
			(*facets)[i]->sh.superIdx = strIdx;
		}
	}


	*nbVertex += nbNewVertex;
	*nbFacet += nbNewFacets;
	if (newStruct) sh.nbSuper += nbNewSuper;
	//return result;
}

void Geometry::InsertSTLGeom(FileReader *file, int *nbVertex, int *nbFacet, VERTEX3D **vertices3, Facet ***facets, int strIdx, double scaleFactor, BOOL newStruct) {

	UnSelectAll();
	char *w;

	int nbNewFacets = 0;
	// First pass
	file->ReadKeyword("solid");
	file->ReadLine(); // solid name
	w = file->ReadWord();
	while (strcmp(w, "facet") == 0) {
		nbNewFacets++;
		file->JumpSection("endfacet");
		w = file->ReadWord();
	}
	if (strcmp(w, "endsolid") != 0) throw Error("Unexpected or not supported STL keyword, 'endsolid' required");

	// Allocate memory
	int nbNewVertex = 3 * nbNewFacets;
	*facets = (Facet **)realloc(*facets, (nbNewFacets + *nbFacet) * sizeof(Facet **));
	memset(*facets + *nbFacet, 0, nbNewFacets * sizeof(Facet *));
	//*vertices3 = (VERTEX3D*)realloc(*vertices3,(nbNewVertex+*nbVertex) * sizeof(VERTEX3D));
	VERTEX3D *tmp_vertices3 = (VERTEX3D *)malloc((nbNewVertex + *nbVertex) * sizeof(VERTEX3D));
	memmove(tmp_vertices3, *vertices3, (*nbVertex)*sizeof(VERTEX3D));
	memset(tmp_vertices3 + *nbVertex, 0, nbNewVertex * sizeof(VERTEX3D));
	SAFE_FREE(*vertices3);
	*vertices3 = tmp_vertices3;

	// Second pass
	file->SeekStart();
	file->ReadKeyword("solid");
	file->ReadLine();
	for (int i = 0; i < nbNewFacets; i++) {

		file->ReadKeyword("facet");
		file->ReadKeyword("normal");
		file->ReadDouble();
		file->ReadDouble();
		file->ReadDouble();
		file->ReadKeyword("outer");
		file->ReadKeyword("loop");

		file->ReadKeyword("vertex");
		(*vertices3)[*nbVertex + 3 * i + 0].x = file->ReadDouble()*scaleFactor;
		(*vertices3)[*nbVertex + 3 * i + 0].y = file->ReadDouble()*scaleFactor;
		(*vertices3)[*nbVertex + 3 * i + 0].z = file->ReadDouble()*scaleFactor;

		file->ReadKeyword("vertex");
		(*vertices3)[*nbVertex + 3 * i + 1].x = file->ReadDouble()*scaleFactor;
		(*vertices3)[*nbVertex + 3 * i + 1].y = file->ReadDouble()*scaleFactor;
		(*vertices3)[*nbVertex + 3 * i + 1].z = file->ReadDouble()*scaleFactor;

		file->ReadKeyword("vertex");
		(*vertices3)[*nbVertex + 3 * i + 2].x = file->ReadDouble()*scaleFactor;
		(*vertices3)[*nbVertex + 3 * i + 2].y = file->ReadDouble()*scaleFactor;
		(*vertices3)[*nbVertex + 3 * i + 2].z = file->ReadDouble()*scaleFactor;

		file->ReadKeyword("endloop");
		file->ReadKeyword("endfacet");

		*(*facets + i + *nbFacet) = new Facet(3);
		(*facets)[i + *nbFacet]->selected = TRUE;
		(*facets)[i + *nbFacet]->indices[0] = *nbVertex + 3 * i + 0;
		(*facets)[i + *nbFacet]->indices[1] = *nbVertex + 3 * i + 1;
		(*facets)[i + *nbFacet]->indices[2] = *nbVertex + 3 * i + 2;

		if (newStruct) {
			(*facets)[i + *nbFacet]->sh.superIdx = sh.nbSuper;
		}
		else {
			(*facets)[i + *nbFacet]->sh.superIdx = strIdx;
		}
	}

	*nbVertex += nbNewVertex;
	*nbFacet += nbNewFacets;
	if (newStruct) AddStruct("Inserted STL file");

}

void Geometry::SaveProfileTXT(FileWriter *file) {
	// Profiles
	for (int j = 0; j < PROFILE_SIZE; j++)
		file->Write("\n");
}

void Geometry::SaveProfileGEO(FileWriter *file, Dataport *dpHit, int super, BOOL saveSelected, BOOL crashSave) {

	BYTE *buffer;
	if (!crashSave && !saveSelected) buffer = (BYTE *)dpHit->buff;
	file->Write("profiles {\n");
	// Profiles
	int nbProfile = 0;
	int *profileFacet = (int *)malloc((sh.nbFacet)*sizeof(int));
	for (int i = 0; i < sh.nbFacet; i++)
		if ((!saveSelected && !crashSave) && facets[i]->sh.isProfile)
			profileFacet[nbProfile++] = i;

	file->Write(" number: "); file->WriteInt(nbProfile, "\n");
	file->Write(" facets: ");
	for (int i = 0; i < nbProfile; i++) //doesn't execute when crashSave or saveSelected...
		file->WriteInt(profileFacet[i], "\t");

	file->Write("\n");
	for (size_t m = 0; (m <= mApp->worker.moments.size()) || (m == 0); m++){
		char tmp[128];
		sprintf(tmp, " moment %d {\n", m);
		file->Write(tmp);
		for (int j = 0; j < PROFILE_SIZE; j++) {
			for (int i = 0; i < nbProfile; i++) { //doesn't execute when crashSave or saveSelected...
				Facet *f = GetFacet(profileFacet[i]);
				APROFILE *profilePtr = (APROFILE *)(buffer + f->sh.hitOffset + sizeof(SHHITS) + m*sizeof(APROFILE)*PROFILE_SIZE);
				//char tmp2[128];
				file->WriteLLong(profilePtr[j].count, "\t");
				file->WriteDouble(profilePtr[j].sum_1_per_ort_velocity, "\t");
				file->WriteDouble(profilePtr[j].sum_v_ort);
				file->Write("\t");
			}
			if (nbProfile>0) file->Write("\n");
		}
		file->Write(" }\n");
	}
	file->Write("}\n");
	SAFE_FREE(profileFacet);
}

void Geometry::LoadProfile(FileReader *file, Dataport *dpHit, int version) {

	AccessDataport(dpHit);
	BYTE *buffer = (BYTE *)dpHit->buff;
	file->ReadKeyword("profiles"); file->ReadKeyword("{");
	// Profiles
	int nbProfile;
	file->ReadKeyword("number"); file->ReadKeyword(":"); nbProfile = file->ReadInt();
	int *profileFacet = (int *)malloc((nbProfile)*sizeof(int));
	file->ReadKeyword("facets"); file->ReadKeyword(":");
	for (int i = 0; i < nbProfile; i++)
		profileFacet[i] = file->ReadInt();
	for (size_t m = 0; m <= mApp->worker.moments.size() || (version < 10 && m == 0); m++){
		if (version >= 10) {
			file->ReadKeyword("moment");
			if (m != file->ReadInt()) {
				char errMsg[512];
				sprintf(errMsg, "Unexpected profile moment");
				throw Error(errMsg);
				break;
			}
			file->ReadKeyword("{");
		}

		for (int j = 0; j < PROFILE_SIZE; j++) {
			for (int i = 0; i < nbProfile; i++) {
				Facet *f = GetFacet(profileFacet[i]);
				APROFILE *profilePtr = (APROFILE *)(buffer + f->sh.hitOffset + sizeof(SHHITS) + m*PROFILE_SIZE*sizeof(APROFILE));
				profilePtr[j].count = file->ReadLLong();
				if (version >= 13) profilePtr[j].sum_1_per_ort_velocity = file->ReadDouble();
				if (version >= 13) profilePtr[j].sum_v_ort = file->ReadDouble();
			}
		}
		if (version >= 10) file->ReadKeyword("}");
	}
	ReleaseDataport(dpHit);
	SAFE_FREE(profileFacet);
}

void Geometry::LoadGEO(FileReader *file, GLProgress *prg, LEAK *pleak, int *nbleak, HIT *pHits, int *nbHHit, int *version, Worker *worker) {

	//mApp->ClearAllSelections();
	//mApp->ClearAllViews();
	prg->SetMessage("Clearing current geometry...");
	Clear();
	//mApp->ClearFormula();

	// Globals
	char tmp[512];
	prg->SetMessage("Reading GEO file header...");
	file->ReadKeyword("version"); file->ReadKeyword(":");
	*version = file->ReadInt();
	if (*version > GEOVERSION) {
		char errMsg[512];
		sprintf(errMsg, "Unsupported GEO version V%d", *version);
		throw Error(errMsg);
	}

	file->ReadKeyword("totalHit"); file->ReadKeyword(":");
	tNbHit = file->ReadLLong();
	file->ReadKeyword("totalDes"); file->ReadKeyword(":");
	tNbDesorption = file->ReadLLong();
	file->ReadKeyword("totalLeak"); file->ReadKeyword(":");
	tNbLeak = file->ReadLLong();
	if (*version >= 12) {
		file->ReadKeyword("totalAbs"); file->ReadKeyword(":");
		tNbAbsorption = file->ReadLLong();
		if (*version >= 15) {
			file->ReadKeyword("totalDist_total");
		}
		else { //between versions 12 and 15
			file->ReadKeyword("totalDist");
		}
		file->ReadKeyword(":");
		distTraveledTotal_total = file->ReadDouble();
		if (*version >= 15) {
			file->ReadKeyword("totalDist_fullHitsOnly"); file->ReadKeyword(":");
			distTraveledTotal_fullHitsOnly = file->ReadDouble();
		}
	}
	else {
		tNbAbsorption = 0;
		distTraveledTotal_total = 0.0;
		distTraveledTotal_fullHitsOnly = 0.0;
	}
	file->ReadKeyword("maxDes"); file->ReadKeyword(":");
	tNbDesorptionMax = file->ReadLLong();
	file->ReadKeyword("nbVertex"); file->ReadKeyword(":");
	sh.nbVertex = file->ReadInt();
	file->ReadKeyword("nbFacet"); file->ReadKeyword(":");
	sh.nbFacet = file->ReadInt();
	file->ReadKeyword("nbSuper"); file->ReadKeyword(":");
	sh.nbSuper = file->ReadInt();
	int nbF = 0; std::vector<std::vector<string>> loadFormulas;
	int nbV = 0;
	if (*version >= 2) {
		file->ReadKeyword("nbFormula"); file->ReadKeyword(":");
		nbF = file->ReadInt(); loadFormulas.reserve(nbF);
		file->ReadKeyword("nbView"); file->ReadKeyword(":");
		nbV = file->ReadInt();
	}
	int nbS = 0;
	if (*version >= 8) {
		file->ReadKeyword("nbSelection"); file->ReadKeyword(":");
		nbS = file->ReadInt();
	}
	if (*version >= 7) {
		file->ReadKeyword("gasMass"); file->ReadKeyword(":");
		worker->gasMass = file->ReadDouble();
	}
	if (*version >= 10) { //time-dependent version
		file->ReadKeyword("userMoments"); file->ReadKeyword("{");
		file->ReadKeyword("nb"); file->ReadKeyword(":");
		int nb = file->ReadInt();

		for (int i = 0; i < nb; i++) {
			char tmpExpr[512];
			strcpy(tmpExpr, file->ReadString());
			worker->userMoments.push_back(tmpExpr);
			worker->AddMoment(mApp->worker.ParseMoment(tmpExpr));
		}
		file->ReadKeyword("}");
	}
	if (*version >= 11) { //pulse version
		file->ReadKeyword("desorptionStart"); file->ReadKeyword(":");
		/*worker->desorptionStartTime =*/ file->ReadDouble();
		file->ReadKeyword("desorptionStop"); file->ReadKeyword(":");
		/*worker->desorptionStopTime =*/ file->ReadDouble();
		file->ReadKeyword("timeWindow"); file->ReadKeyword(":");
		worker->timeWindowSize = file->ReadDouble();
		file->ReadKeyword("useMaxwellian"); file->ReadKeyword(":");
		worker->useMaxwellDistribution = file->ReadInt();
	}
	if (*version >= 12) { //2013.aug.22
		file->ReadKeyword("calcConstantFlow"); file->ReadKeyword(":");
		worker->calcConstantFlow = file->ReadInt();
	}
	if (*version >= 2) {
		file->ReadKeyword("formulas"); file->ReadKeyword("{");
		for (int i = 0; i < nbF; i++) {
			char tmpName[256];
			char tmpExpr[512];
			strcpy(tmpName, file->ReadString());
			strcpy(tmpExpr, file->ReadString());
			//mApp->AddFormula(tmpName, tmpExpr); //parse after selection groups are loaded
			std::vector<string> newFormula;
			newFormula.push_back(tmpName);
			newFormula.push_back(tmpExpr);
			loadFormulas.push_back(newFormula);
		}
		file->ReadKeyword("}");

		file->ReadKeyword("views"); file->ReadKeyword("{");
		for (int i = 0; i < nbV; i++) {
			char tmpName[256];
			AVIEW v;
			strcpy(tmpName, file->ReadString());
			v.projMode = file->ReadInt();
			v.camAngleOx = file->ReadDouble();
			v.camAngleOy = file->ReadDouble();
			v.camDist = file->ReadDouble();
			v.camOffset.x = file->ReadDouble();
			v.camOffset.y = file->ReadDouble();
			v.camOffset.z = file->ReadDouble();
			v.performXY = file->ReadInt();

			v.vLeft = file->ReadDouble();
			v.vRight = file->ReadDouble();
			v.vTop = file->ReadDouble();
			v.vBottom = file->ReadDouble();
			mApp->AddView(tmpName, v);
		}
		file->ReadKeyword("}");
	}

	if (*version >= 8) {
		file->ReadKeyword("selections"); file->ReadKeyword("{");
		for (int i = 0; i < nbS; i++) {
			ASELECTION s;
			char tmpName[256];
			strcpy(tmpName, file->ReadString());
			s.name = _strdup(tmpName);
			s.nbSel = file->ReadInt();
			s.selection = (int *)malloc((s.nbSel)*sizeof(int));

			for (int j = 0; j < s.nbSel; j++) {
				s.selection[j] = file->ReadInt();
			}
			mApp->AddSelection(s.name, s);
		}
		file->ReadKeyword("}");
	}

	for (int i = 0; i < nbF; i++) { //parse formulas now that selection groups are loaded
		mApp->AddFormula(loadFormulas[i][0].c_str(), loadFormulas[i][1].c_str());
	}

	file->ReadKeyword("structures"); file->ReadKeyword("{");
	for (int i = 0; i < sh.nbSuper; i++) {
		strName[i] = _strdup(file->ReadString());
		// For backward compatibilty with STR
		sprintf(tmp, "%s.txt", strName[i]);
		strFileName[i] = _strdup(tmp);
	}
	file->ReadKeyword("}");

	// Allocate memory
	facets = (Facet **)malloc(sh.nbFacet * sizeof(Facet *));
	memset(facets, 0, sh.nbFacet * sizeof(Facet *));
	vertices3 = (VERTEX3D *)malloc(sh.nbVertex * sizeof(VERTEX3D));

	// Read vertices
	prg->SetMessage("Reading vertices...");
	file->ReadKeyword("vertices"); file->ReadKeyword("{");
	for (int i = 0; i < sh.nbVertex; i++) {
		// Check idx
		int idx = file->ReadInt();
		if (idx != i + 1) throw Error(file->MakeError("Wrong vertex index !"));
		vertices3[i].x = file->ReadDouble();
		vertices3[i].y = file->ReadDouble();
		vertices3[i].z = file->ReadDouble();
		vertices3[i].selected = FALSE;
	}
	file->ReadKeyword("}");

	if (*version >= 6) {
		prg->SetMessage("Reading leaks and hits...");
		// Read leaks
		file->ReadKeyword("leaks"); file->ReadKeyword("{");
		file->ReadKeyword("nbLeak"); file->ReadKeyword(":");
		*nbleak = file->ReadInt();
		for (int i = 0; i < *nbleak; i++) {
			int idx = file->ReadInt();
			if (idx != i) throw Error(file->MakeError("Wrong leak index !"));
			(pleak + i)->pos.x = file->ReadDouble();
			(pleak + i)->pos.y = file->ReadDouble();
			(pleak + i)->pos.z = file->ReadDouble();

			(pleak + i)->dir.x = file->ReadDouble();
			(pleak + i)->dir.y = file->ReadDouble();
			(pleak + i)->dir.z = file->ReadDouble();
		}
		file->ReadKeyword("}");

		// Read hit cache
		file->ReadKeyword("hits"); file->ReadKeyword("{");
		file->ReadKeyword("nbHHit"); file->ReadKeyword(":");
		*nbHHit = file->ReadInt();
		for (int i = 0; i < *nbHHit; i++) {
			int idx = file->ReadInt();
			if (idx != i) throw Error(file->MakeError("Wrong hit cache index !"));
			(pHits + i)->pos.x = file->ReadDouble();
			(pHits + i)->pos.y = file->ReadDouble();
			(pHits + i)->pos.z = file->ReadDouble();

			(pHits + i)->type = file->ReadInt();
		}
		file->ReadKeyword("}");
	}

	// Read facets
	prg->SetMessage("Reading facets...");
	for (int i = 0; i < sh.nbFacet; i++) {
		file->ReadKeyword("facet");
		// Check idx
		int idx = file->ReadInt();
		if (idx != i + 1) throw Error(file->MakeError("Wrong facet index !"));
		file->ReadKeyword("{");
		file->ReadKeyword("nbIndex");
		file->ReadKeyword(":");
		int nbI = file->ReadInt();
		if (nbI < 3) {
			char errMsg[512];
			sprintf(errMsg, "Facet %d has only %d vertices. ", i + 1, nbI);
			throw Error(errMsg);
		}
		prg->SetProgress((float)i / sh.nbFacet);
		facets[i] = new Facet(nbI);
		facets[i]->LoadGEO(file, *version, sh.nbVertex);
		file->ReadKeyword("}");
	}

	InitializeGeometry();
	//AdjustProfile();
	isLoaded = TRUE;
	UpdateName(file);

	// Update mesh
	prg->SetMessage("Building mesh...");
	for (int i = 0; i < sh.nbFacet; i++) {
		double p = (double)i / (double)sh.nbFacet;
		prg->SetProgress(p);
		Facet *f = facets[i];
		if (!f->SetTexture(f->sh.texWidthD, f->sh.texHeightD, f->hasMesh)) {
			char errMsg[512];
			sprintf(errMsg, "Not enough memory to build mesh on Facet %d. ", i + 1);
			throw Error(errMsg);
		}
		BuildFacetList(f);
		double nU = Norme(&(f->sh.U));
		f->tRatio = f->sh.texWidthD / nU;
	}

}

void Geometry::LoadSYN(FileReader *file, GLProgress *prg, LEAK *pleak, int *nbleak, HIT *pHits, int *nbHHit, int *version) {



	//mApp->ClearAllSelections();
	//mApp->ClearAllViews();
	prg->SetMessage("Clearing current geometry...");
	Clear();
	//mApp->ClearFormula();

	// Globals
	char tmp[512];
	prg->SetMessage("Reading SYN file header...");
	file->ReadKeyword("version"); file->ReadKeyword(":");
	*version = file->ReadInt();
	if (*version > SYNVERSION) {
		char errMsg[512];
		sprintf(errMsg, "Unsupported SYN version V%d", *version);
		throw Error(errMsg);
	}
	file->ReadKeyword("totalHit"); file->ReadKeyword(":");
	tNbHit = 0; file->ReadLLong();
	file->ReadKeyword("totalDes"); file->ReadKeyword(":");
	tNbDesorption = 0; file->ReadLLong();
	if (*version >= 6) {
		file->ReadKeyword("no_scans"); file->ReadKeyword(":");
		/*loaded_no_scans = */file->ReadDouble();
	}
	file->ReadKeyword("totalLeak"); file->ReadKeyword(":");
	tNbLeak = 0; file->ReadLLong();
	if (*version > 2) {
		file->ReadKeyword("totalFlux"); file->ReadKeyword(":");
		file->ReadDouble();
		file->ReadKeyword("totalPower"); file->ReadKeyword(":");
		file->ReadDouble();
	}
	file->ReadKeyword("maxDes"); file->ReadKeyword(":");
	tNbDesorptionMax = 0; file->ReadLLong();
	file->ReadKeyword("nbVertex"); file->ReadKeyword(":");
	sh.nbVertex = file->ReadInt();
	file->ReadKeyword("nbFacet"); file->ReadKeyword(":");
	sh.nbFacet = file->ReadInt();
	file->ReadKeyword("nbSuper"); file->ReadKeyword(":");
	sh.nbSuper = file->ReadInt();
	int nbF = 0;
	int nbV = 0;

	file->ReadKeyword("nbFormula"); file->ReadKeyword(":");
	nbF = file->ReadInt();
	file->ReadKeyword("nbView"); file->ReadKeyword(":");
	nbV = file->ReadInt();
	int nbS = 0;
	file->ReadKeyword("nbSelection"); file->ReadKeyword(":");
	nbS = file->ReadInt();

	if (*version > 1) {
		file->ReadKeyword("nbRegions"); file->ReadKeyword(":");
		int nbR = file->ReadInt();
		//result=PARfileList(nbR);

		file->ReadKeyword("PARfiles"); file->ReadKeyword("{");
		for (int i = 0; i < nbR; i++) {
			/*char tmp[512];
			strcpy(tmp,file->ReadString());
			result.fileNames[i]=_strdup(tmp);*/
			file->ReadString();
		}
		file->ReadKeyword("}");
	}

	file->ReadKeyword("formulas"); file->ReadKeyword("{");
	for (int i = 0; i < nbF; i++) {
		char tmpName[256];
		char tmpExpr[512];
		strcpy(tmpName, file->ReadString());
		strcpy(tmpExpr, file->ReadString());
		mApp->AddFormula(tmpName, tmpExpr);
	}
	file->ReadKeyword("}");

	file->ReadKeyword("views"); file->ReadKeyword("{");
	for (int i = 0; i < nbV; i++) {
		char tmpName[256];
		AVIEW v;
		strcpy(tmpName, file->ReadString());
		v.projMode = file->ReadInt();
		v.camAngleOx = file->ReadDouble();
		v.camAngleOy = file->ReadDouble();
		v.camDist = file->ReadDouble();
		v.camOffset.x = file->ReadDouble();
		v.camOffset.y = file->ReadDouble();
		v.camOffset.z = file->ReadDouble();
		v.performXY = file->ReadInt();
		v.vLeft = file->ReadDouble();
		v.vRight = file->ReadDouble();
		v.vTop = file->ReadDouble();
		v.vBottom = file->ReadDouble();
		mApp->AddView(tmpName, v);
	}
	file->ReadKeyword("}");

	file->ReadKeyword("selections"); file->ReadKeyword("{");
	for (int i = 0; i < nbS; i++) {
		ASELECTION s;
		char tmpName[256];
		strcpy(tmpName, file->ReadString());
		s.name = _strdup(tmpName);
		s.nbSel = file->ReadInt();
		s.selection = (int *)malloc((s.nbSel)*sizeof(int));

		for (int j = 0; j < s.nbSel; j++) {
			s.selection[j] = file->ReadInt();
		}
		mApp->AddSelection(s.name, s);
	}
	file->ReadKeyword("}");

	file->ReadKeyword("structures"); file->ReadKeyword("{");
	for (int i = 0; i < sh.nbSuper; i++) {
		strName[i] = _strdup(file->ReadString());
		// For backward compatibilty with STR
		sprintf(tmp, "%s.txt", strName[i]);
		strFileName[i] = _strdup(tmp);
	}
	file->ReadKeyword("}");

	// Allocate memory
	facets = (Facet **)malloc(sh.nbFacet * sizeof(Facet *));
	memset(facets, 0, sh.nbFacet * sizeof(Facet *));
	vertices3 = (VERTEX3D *)malloc(sh.nbVertex * sizeof(VERTEX3D));

	// Read vertices
	prg->SetMessage("Reading vertices...");
	file->ReadKeyword("vertices"); file->ReadKeyword("{");
	for (int i = 0; i < sh.nbVertex; i++) {
		// Check idx
		int idx = file->ReadInt();
		if (idx != i + 1) throw Error(file->MakeError("Wrong vertex index !"));
		vertices3[i].x = file->ReadDouble();
		vertices3[i].y = file->ReadDouble();
		vertices3[i].z = file->ReadDouble();
		vertices3[i].selected = FALSE;
	}
	file->ReadKeyword("}");
	prg->SetMessage("Reading leaks and hits...");
	// Read leaks
	file->ReadKeyword("leaks"); file->ReadKeyword("{");
	file->ReadKeyword("nbLeak"); file->ReadKeyword(":");
	*nbleak = 0;
	int nbleak_local = file->ReadInt();
	for (int i = 0; i < nbleak_local; i++) {
		int idx = file->ReadInt();
		if (idx != i) throw Error(file->MakeError("Wrong leak index !"));
		//(pleak+i)->pos.x = 
		file->ReadDouble();
		//(pleak+i)->pos.y = 
		file->ReadDouble();
		//(pleak+i)->pos.z = 
		file->ReadDouble();

		//(pleak+i)->dir.x = 
		file->ReadDouble();
		//(pleak+i)->dir.y = 
		file->ReadDouble();
		//(pleak+i)->dir.z = 
		file->ReadDouble();
	}
	file->ReadKeyword("}");

	// Read hit cache
	file->ReadKeyword("hits"); file->ReadKeyword("{");
	file->ReadKeyword("nbHHit"); file->ReadKeyword(":");
	*nbHHit = 0;
	int nbHHit_local = file->ReadInt();
	for (int i = 0; i < nbHHit_local; i++) {
		int idx = file->ReadInt();
		if (idx != i) throw Error(file->MakeError("Wrong hit cache index !"));
		//(pHits+i)->pos.x = 
		file->ReadDouble();
		//(pHits+i)->pos.y = 
		file->ReadDouble();
		//(pHits+i)->pos.z = 
		file->ReadDouble();
		//(pHits+i)->dF = 
		file->ReadDouble();
		//(pHits+i)->dP = 
		file->ReadDouble();
		//(pHits+i)->type = 
		file->ReadInt();
	}
	file->ReadKeyword("}");
	// Read facets
	prg->SetMessage("Reading facets...");
	for (int i = 0; i < sh.nbFacet; i++) {
		file->ReadKeyword("facet");
		// Check idx
		int idx = file->ReadInt();
		if (idx != i + 1) throw Error(file->MakeError("Wrong facet index !"));
		file->ReadKeyword("{");
		file->ReadKeyword("nbIndex");
		file->ReadKeyword(":");
		int nbI = file->ReadInt();
		if (nbI < 3) {
			char errMsg[512];
			sprintf(errMsg, "Facet %d has only %d vertices. ", i, nbI);
			throw Error(errMsg);
		}
		prg->SetProgress((float)i / sh.nbFacet);
		facets[i] = new Facet(nbI);
		facets[i]->LoadSYN(file, *version, sh.nbVertex);
		file->ReadKeyword("}");
	}

	prg->SetMessage("Initalizing geometry and building mesh...");
	InitializeGeometry();
	//AdjustProfile();
	isLoaded = TRUE;
	UpdateName(file);

	// Update meshgHits->distTraveledTotal
	prg->SetMessage("Drawing textures...");
	for (int i = 0; i < sh.nbFacet; i++) {
		double p = (double)i / (double)sh.nbFacet;
		prg->SetProgress(p);
		Facet *f = facets[i];
		//f->SetTexture(f->sh.texWidthD,f->sh.texHeightD,f->hasMesh);
		BuildFacetList(f);
		//double nU = Norme(&(f->sh.U));
		//f->tRatio = f->sh.texWidthD / nU;
	}
	//return result;
}

bool Geometry::LoadTextures(FileReader *file, GLProgress *prg, Dataport *dpHit, int version) {

	if (file->SeekFor("{textures}")) {
		char tmp[256];
		//versions 3+
		// Block dpHit during the whole disc reading

		AccessDataport(dpHit);
		//AHIT readVal;

		// Globals
		BYTE *buffer = (BYTE *)dpHit->buff;
		SHGHITS *gHits = (SHGHITS *)buffer;

		gHits->total.hit.nbHit = tNbHit;
		gHits->total.hit.nbDesorbed = tNbDesorption;
		gHits->total.hit.nbAbsorbed = tNbAbsorption;
		gHits->nbLeakTotal = tNbLeak;
		gHits->distTraveledTotal_total = distTraveledTotal_total;
		gHits->distTraveledTotal_fullHitsOnly = distTraveledTotal_fullHitsOnly;

		// Read facets
		if (version >= 13) {
			file->ReadKeyword("min_pressure_all"); file->ReadKeyword(":");
			gHits->texture_limits[0].min.all = file->ReadDouble();
			file->ReadKeyword("min_pressure_moments_only"); file->ReadKeyword(":");
			gHits->texture_limits[0].min.moments_only = file->ReadDouble();
			file->ReadKeyword("max_pressure_all"); file->ReadKeyword(":");
			gHits->texture_limits[0].max.all = file->ReadDouble();
			file->ReadKeyword("max_pressure_moments_only"); file->ReadKeyword(":");
			gHits->texture_limits[0].max.moments_only = file->ReadDouble();

			file->ReadKeyword("min_impingement_all"); file->ReadKeyword(":");
			gHits->texture_limits[1].min.all = file->ReadDouble();
			file->ReadKeyword("min_impingement_moments_only"); file->ReadKeyword(":");
			gHits->texture_limits[1].min.moments_only = file->ReadDouble();
			file->ReadKeyword("max_impingement_all"); file->ReadKeyword(":");
			gHits->texture_limits[1].max.all = file->ReadDouble();
			file->ReadKeyword("max_impingement_moments_only"); file->ReadKeyword(":");
			gHits->texture_limits[1].max.moments_only = file->ReadDouble();

			file->ReadKeyword("min_density_all"); file->ReadKeyword(":");
			gHits->texture_limits[2].min.all = file->ReadDouble();
			file->ReadKeyword("min_density_moments_only"); file->ReadKeyword(":");
			gHits->texture_limits[2].min.moments_only = file->ReadDouble();
			file->ReadKeyword("max_density_all"); file->ReadKeyword(":");
			gHits->texture_limits[2].max.all = file->ReadDouble();
			file->ReadKeyword("max_density_moments_only"); file->ReadKeyword(":");
			gHits->texture_limits[2].max.moments_only = file->ReadDouble();

			for (size_t m = 0; m <= mApp->worker.moments.size() || (m == 0 /*&& version<10*/); m++){
				//if (version>=10) {
				file->ReadKeyword("moment");
				if (m != file->ReadInt()) {
					throw Error("Unexpected profile moment");
					break;
				}
				file->ReadKeyword("{");
				//}

				for (int i = 0; i < sh.nbFacet; i++) {
					Facet *f = facets[i];
					if (f->hasMesh/* || version<8*/) {
						prg->SetProgress((double)(i + m*sh.nbFacet) / (double)(mApp->worker.moments.size()*sh.nbFacet)*0.33 + 0.66);
						file->ReadKeyword("texture_facet");
						// Check idx
						int idx = file->ReadInt();

						if (idx != i + 1) {
							sprintf(tmp, "Wrong facet index. Expected %d, read %d.", i + 1, idx);
							throw Error(file->MakeError(tmp));
						}
						file->ReadKeyword("{");

						int ix, iy;

						///Load textures, for GEO file version 3+

						int profSize = (f->sh.isProfile) ? ((1 + (int)mApp->worker.moments.size())*(PROFILE_SIZE*sizeof(APROFILE))) : 0;
						int h = (f->sh.texHeight);
						int w = (f->sh.texWidth);
						AHIT *hits = (AHIT *)((BYTE *)gHits + (f->sh.hitOffset + sizeof(SHHITS) + profSize + m*w*h*sizeof(AHIT)));

						int texWidth_file, texHeight_file;
						//In case of rounding errors, the file might contain different texture dimensions than expected.
						if (version >= 14) {
							file->ReadKeyword("width"); file->ReadKeyword(":"); texWidth_file = file->ReadInt();
							file->ReadKeyword("height"); file->ReadKeyword(":"); texHeight_file = file->ReadInt();
						}
						else {
							texWidth_file = f->sh.texWidth;
							texHeight_file = f->sh.texHeight;
						}


						for (iy = 0; iy<(MIN(f->sh.texHeight, texHeight_file)); iy++) { //MIN: If stored texture is larger, don't read extra cells
							for (ix = 0; ix<(MIN(f->sh.texWidth, texWidth_file)); ix++) { //MIN: If stored texture is larger, don't read extra cells
								hits[iy*f->sh.texWidth + ix].count = file->ReadLLong();
								hits[iy*f->sh.texWidth + ix].sum_1_per_ort_velocity = file->ReadDouble();
								hits[iy*f->sh.texWidth + ix].sum_v_ort_per_area = file->ReadDouble();
							}
							for (int ie = 0; ie < texWidth_file - f->sh.texWidth; ie++) {//Executed if file texture is bigger than expected texture
								//Read extra cells from file without doing anything
								file->ReadLLong();
								file->ReadDouble();
								file->ReadDouble();
							}
						}
						for (int ie = 0; ie < texHeight_file - f->sh.texHeight; ie++) {//Executed if file texture is bigger than expected texture
							//Read extra cells ffrom file without doing anything
							for (int iw = 0; iw < texWidth_file; iw++) {
								file->ReadLLong();
								file->ReadDouble();
								file->ReadDouble();
							}
						}
						file->ReadKeyword("}");
					}
				}
				/*if (version>=10)*/ file->ReadKeyword("}");
			}
		}
		ReleaseDataport(dpHit);
		//Debug memory check
		//_ASSERTE (!_CrtDumpMemoryLeaks());;
		_ASSERTE(_CrtCheckMemory());
		return true;

	}
	else
	{
		//old versions
		return false;
	}

}

void Geometry::SaveGEO(FileWriter *file, GLProgress *prg, Dataport *dpHit, std::vector<std::string> userMoments, Worker *worker,
	BOOL saveSelected, LEAK *pleak, int *nbleakSave, HIT *pHits, int *nbHHitSave, BOOL crashSave) {


	prg->SetMessage("Counting hits...");
	if (!IsLoaded()) throw Error("Nothing to save !");


	// Block dpHit during the whole disc writing
	if (!crashSave && !saveSelected) AccessDataport(dpHit);

	// Globals
	BYTE *buffer;
	if (!crashSave && !saveSelected) buffer = (BYTE *)dpHit->buff;
	SHGHITS *gHits;
	if (!crashSave && !saveSelected) gHits = (SHGHITS *)buffer;

	double dCoef = 1.0;
	int ix, iy;

	/*switch(gHits->mode) {

	case MC_MODE:
	if( gHits->total.hit.nbDesorbed>0 ) {
	dCoef = (float)totalOutgassing / (float)gHits->total.hit.nbDesorbed / 8.31 * gasMass / 100;
	texMinAutoscale = gHits->minHit * dCoef;
	texMaxAutoscale = gHits->maxHit * dCoef;
	} else {
	texMinAutoscale = gHits->minHit;
	texMaxAutoscale = gHits->maxHit;
	}
	break;

	case AC_MODE:
	texMinAutoscale = gHits->minHit;
	texMaxAutoscale = gHits->maxHit;
	break;

	}*/

	prg->SetMessage("Writing geometry details...");
	file->Write("version:"); file->WriteInt(GEOVERSION, "\n");
	file->Write("totalHit:"); file->WriteLLong((!crashSave && !saveSelected) ? gHits->total.hit.nbHit : 0, "\n");
	file->Write("totalDes:"); file->WriteLLong((!crashSave && !saveSelected) ? gHits->total.hit.nbDesorbed : 0, "\n");
	file->Write("totalLeak:"); file->WriteLLong((!crashSave && !saveSelected) ? gHits->nbLeakTotal : 0, "\n");
	file->Write("totalAbs:"); file->WriteLLong((!crashSave && !saveSelected) ? gHits->total.hit.nbAbsorbed : 0, "\n");
	file->Write("totalDist_total:"); file->WriteDouble((!crashSave && !saveSelected) ? gHits->distTraveledTotal_total : 0, "\n");
	file->Write("totalDist_fullHitsOnly:"); file->WriteDouble((!crashSave && !saveSelected) ? gHits->distTraveledTotal_fullHitsOnly : 0, "\n");
	file->Write("maxDes:"); file->WriteLLong((!crashSave && !saveSelected) ? tNbDesorptionMax : 0, "\n");
	file->Write("nbVertex:"); file->WriteInt(sh.nbVertex, "\n");
	file->Write("nbFacet:"); file->WriteInt(saveSelected ? nbSelected : sh.nbFacet, "\n");
	file->Write("nbSuper:"); file->WriteInt(sh.nbSuper, "\n");
	file->Write("nbFormula:"); file->WriteInt((!saveSelected) ? mApp->nbFormula : 0, "\n");
	file->Write("nbView:"); file->WriteInt(mApp->nbView, "\n");
	file->Write("nbSelection:"); file->WriteInt((!saveSelected) ? mApp->nbSelection : 0, "\n");
	file->Write("gasMass:"); file->WriteDouble(worker->gasMass, "\n");

	file->Write("userMoments {\n");
	file->Write(" nb:"); file->WriteInt((int)userMoments.size());
	for (size_t u = 0; u < userMoments.size(); u++) {
		file->Write("\n \"");
		file->Write(userMoments[u].c_str());
		file->Write("\"");
	}
	file->Write("\n}\n");

	file->Write("desorptionStart:"); file->WriteDouble(/*worker->desorptionStartTime*/0.0, "\n");
	file->Write("desorptionStop:"); file->WriteDouble(/*worker->desorptionStopTime*/1.0, "\n");
	file->Write("timeWindow:"); file->WriteDouble(worker->timeWindowSize, "\n");
	file->Write("useMaxwellian:"); file->WriteInt(worker->useMaxwellDistribution, "\n");
	file->Write("calcConstantFlow:"); file->WriteInt(worker->calcConstantFlow, "\n");

	file->Write("formulas {\n");
	if (!saveSelected){
		for (int i = 0; i < mApp->nbFormula; i++) {
			file->Write("  \"");
			file->Write(mApp->formulas[i].parser->GetName());
			file->Write("\" \"");
			file->Write(mApp->formulas[i].parser->GetExpression());
			file->Write("\"\n");
		}
	}
	file->Write("}\n");

	file->Write("views {\n");
	for (int i = 0; i < mApp->nbView; i++) {
		file->Write("  \"");
		file->Write(mApp->views[i].name);
		file->Write("\"\n");
		file->WriteInt(mApp->views[i].projMode, " ");
		file->WriteDouble(mApp->views[i].camAngleOx, " ");
		file->WriteDouble(mApp->views[i].camAngleOy, " ");
		file->WriteDouble(mApp->views[i].camDist, " ");
		file->WriteDouble(mApp->views[i].camOffset.x, " ");
		file->WriteDouble(mApp->views[i].camOffset.y, " ");
		file->WriteDouble(mApp->views[i].camOffset.z, " ");
		file->WriteInt(mApp->views[i].performXY, " ");
		file->WriteDouble(mApp->views[i].vLeft, " ");
		file->WriteDouble(mApp->views[i].vRight, " ");
		file->WriteDouble(mApp->views[i].vTop, " ");
		file->WriteDouble(mApp->views[i].vBottom, "\n");
	}
	file->Write("}\n");

	file->Write("selections {\n");
	for (int i = 0; (i < mApp->nbSelection) && !saveSelected; i++) { //don't save selections when exporting part of the geometry (saveSelected)
		file->Write("  \"");
		file->Write(mApp->selections[i].name);
		file->Write("\"\n ");
		file->WriteInt(mApp->selections[i].nbSel, "\n");
		for (int j = 0; j < mApp->selections[i].nbSel; j++) {
			file->Write("  ");
			file->WriteInt(mApp->selections[i].selection[j], "\n");
		}
		//file->Write("\n");
	}
	file->Write("}\n");

	file->Write("structures {\n");
	for (int i = 0; i < sh.nbSuper; i++) {
		file->Write("  \"");
		file->Write(strName[i]);
		file->Write("\"\n");
	}
	file->Write("}\n");
	//vertices
	prg->SetMessage("Writing vertices...");
	file->Write("vertices {\n");
	for (int i = 0; i < sh.nbVertex; i++) {
		prg->SetProgress(0.33*((double)i / (double)sh.nbVertex));
		file->Write("  ");
		file->WriteInt(i + 1, " ");
		file->WriteDouble(vertices3[i].x, " ");
		file->WriteDouble(vertices3[i].y, " ");
		file->WriteDouble(vertices3[i].z, "\n");
	}
	file->Write("}\n");

	//leaks
	prg->SetMessage("Writing leaks...");
	file->Write("leaks {\n");
	file->Write("  nbLeak:"); file->WriteInt((!crashSave && !saveSelected) ? *nbleakSave : 0, "\n");
	for (int i = 0; (i < *nbleakSave) && (!crashSave && !saveSelected); i++) {

		file->Write("  ");
		file->WriteInt(i, " ");
		file->WriteDouble((pleak + i)->pos.x, " ");
		file->WriteDouble((pleak + i)->pos.y, " ");
		file->WriteDouble((pleak + i)->pos.z, " ");

		file->WriteDouble((pleak + i)->dir.x, " ");
		file->WriteDouble((pleak + i)->dir.y, " ");
		file->WriteDouble((pleak + i)->dir.z, "\n");
	}
	file->Write("}\n");

	//hit cache (lines and dots)
	prg->SetMessage("Writing hit cache...");
	file->Write("hits {\n");
	file->Write("  nbHHit:"); file->WriteInt((!crashSave && !saveSelected) ? *nbHHitSave : 0, "\n");
	for (int i = 0; (i < *nbHHitSave) && (!crashSave && !saveSelected); i++) {

		file->Write("  ");
		file->WriteInt(i, " ");
		file->WriteDouble((pHits + i)->pos.x, " ");
		file->WriteDouble((pHits + i)->pos.y, " ");
		file->WriteDouble((pHits + i)->pos.z, " ");

		file->WriteInt((pHits + i)->type, "\n");
	}
	file->Write("}\n");

	//facets

	prg->SetMessage("Writing facets...");

	for (int i = 0, k = 0; i < sh.nbFacet; i++) {
		prg->SetProgress(0.33 + ((double)i / (double)sh.nbFacet) *0.33);
		if (!saveSelected || facets[i]->selected) facets[i]->SaveGEO(file, k++);
	}

	prg->SetMessage("Writing profiles...");
	SaveProfileGEO(file, dpHit, -1, saveSelected, crashSave);

	///Save textures, for GEO file version 3+
	char tmp[256];
	file->Write("{textures}\n");

	file->Write("min_pressure_all:"); file->WriteDouble(
		(!crashSave && !saveSelected) ? gHits->texture_limits[0].min.all : 0, "\n");
	file->Write("min_pressure_moments_only:"); file->WriteDouble(
		(!crashSave && !saveSelected) ? gHits->texture_limits[0].min.moments_only : 0, "\n");
	file->Write("max_pressure_all:"); file->WriteDouble(
		(!crashSave && !saveSelected) ? gHits->texture_limits[0].max.all : 1, "\n");
	file->Write("max_pressure_moments_only:"); file->WriteDouble(
		(!crashSave && !saveSelected) ? gHits->texture_limits[0].max.moments_only : 1, "\n");

	file->Write("min_impingement_all:"); file->WriteDouble(
		(!crashSave && !saveSelected) ? gHits->texture_limits[1].min.all : 0, "\n");
	file->Write("min_impingement_moments_only:"); file->WriteDouble(
		(!crashSave && !saveSelected) ? gHits->texture_limits[1].min.moments_only : 0, "\n");
	file->Write("max_impingement_all:"); file->WriteDouble(
		(!crashSave && !saveSelected) ? gHits->texture_limits[1].max.all : 1, "\n");
	file->Write("max_impingement_moments_only:"); file->WriteDouble(
		(!crashSave && !saveSelected) ? gHits->texture_limits[1].max.moments_only : 1, "\n");

	file->Write("min_density_all:"); file->WriteDouble(
		(!crashSave && !saveSelected) ? gHits->texture_limits[2].min.all : 0, "\n");
	file->Write("min_density_moments_only:"); file->WriteDouble(
		(!crashSave && !saveSelected) ? gHits->texture_limits[2].min.moments_only : 0, "\n");
	file->Write("max_density_all:"); file->WriteDouble(
		(!crashSave && !saveSelected) ? gHits->texture_limits[2].max.all : 1, "\n");
	file->Write("max_density_moments_only:"); file->WriteDouble(
		(!crashSave && !saveSelected) ? gHits->texture_limits[2].max.moments_only : 1, "\n");

	//Selections
	//SaveSelections();

	prg->SetMessage("Writing textures...");
	for (size_t m = 0; m <= mApp->worker.moments.size(); m++){
		sprintf(tmp, "moment %d {\n", m);
		file->Write(tmp);
		for (int i = 0; i < sh.nbFacet; i++) {
			prg->SetProgress((double)(i + m*sh.nbFacet) / (double)(mApp->worker.moments.size()*sh.nbFacet)*0.33 + 0.66);
			Facet *f = facets[i];
			if (f->hasMesh) {
				int h = (f->sh.texHeight);
				int w = (f->sh.texWidth);
				int profSize = (f->sh.isProfile) ? (PROFILE_SIZE*sizeof(APROFILE)*(1 + (int)mApp->worker.moments.size())) : 0;
				AHIT *hits;
				if (!crashSave && !saveSelected) hits = (AHIT *)((BYTE *)gHits + (f->sh.hitOffset + sizeof(SHHITS) + profSize + m*w*h*sizeof(AHIT)));

				//char tmp[256];
				sprintf(tmp, " texture_facet %d {\n", i + 1);
				file->Write(tmp);
				file->Write("width:"); file->WriteInt(f->sh.texWidth); file->Write(" height:"); file->WriteInt(f->sh.texHeight); file->Write("\n");
				for (iy = 0; iy < h; iy++) {
					for (ix = 0; ix < w; ix++) {
						file->WriteLLong((!crashSave && !saveSelected) ? hits[iy*f->sh.texWidth + ix].count : 0, "\t");
						file->WriteDouble((!crashSave && !saveSelected) ? hits[iy*f->sh.texWidth + ix].sum_1_per_ort_velocity : 0, "\t");
						file->WriteDouble((!crashSave && !saveSelected) ? hits[iy*f->sh.texWidth + ix].sum_v_ort_per_area : 0, "\t");
					}
					file->Write("\n");
				}
				file->Write(" }\n");
			}
		}
		file->Write("}\n");
	}

	if (!crashSave && !saveSelected) ReleaseDataport(dpHit);

}

void Geometry::SaveTXT(FileWriter *file, Dataport *dpHit, BOOL saveSelected) {

	if (!IsLoaded()) throw Error("Nothing to save !");

	// Unused
	file->WriteInt(0, "\n");

	// Block dpHit during the whole disc writing
	AccessDataport(dpHit);

	// Globals
	BYTE *buffer = (BYTE *)dpHit->buff;
	SHGHITS *gHits = (SHGHITS *)buffer;

	// Unused
	file->WriteLLong(gHits->total.hit.nbHit, "\n");
	file->WriteLLong(gHits->nbLeakTotal, "\n");
	file->WriteLLong(gHits->total.hit.nbDesorbed, "\n");
	file->WriteLLong(tNbDesorptionMax, "\n");

	file->WriteInt(sh.nbVertex, "\n");
	file->WriteInt(saveSelected ? nbSelected : sh.nbFacet, "\n");

	// Read geometry vertices
	for (int i = 0; i < sh.nbVertex; i++) {
		file->WriteDouble(vertices3[i].x, " ");
		file->WriteDouble(vertices3[i].y, " ");
		file->WriteDouble(vertices3[i].z, "\n");
	}

	// Facets
	for (int i = 0; i < sh.nbFacet; i++) {
		Facet *f = facets[i];
		int j;
		if (saveSelected) {
			if (f->selected) {
				file->WriteInt(f->sh.nbIndex, " ");
				for (j = 0; j < f->sh.nbIndex - 1; j++)
					file->WriteInt(f->indices[j] + 1, " ");
				file->WriteInt(f->indices[j] + 1, "\n");
			}
		}
		else {
			file->WriteInt(f->sh.nbIndex, " ");
			for (j = 0; j < f->sh.nbIndex - 1; j++)
				file->WriteInt(f->indices[j] + 1, " ");
			file->WriteInt(f->indices[j] + 1, "\n");
		}
	}

	// Params
	for (int i = 0; i < sh.nbFacet; i++) {

		// Update facet hits from shared mem
		Facet *f = facets[i];
		SHHITS *shF = (SHHITS *)(buffer + f->sh.hitOffset);
		memcpy(&(f->sh.counter), shF, sizeof(SHHITS));
		if (saveSelected) {
			if (f->selected) f->SaveTXT(file);
		}
		else {
			f->SaveTXT(file);
		}

	}

	SaveProfileTXT(file);

	ReleaseDataport(dpHit);

}

void Geometry::ExportTextures(FILE *file, int mode, Dataport *dpHit, BOOL saveSelected) {

	//if(!IsLoaded()) throw Error("Nothing to save !");

	// Block dpHit during the whole disc writing
	BYTE *buffer = NULL;
	if (dpHit)
		if (AccessDataport(dpHit))
			buffer = (BYTE *)dpHit->buff;

	// Globals
	//BYTE *buffer = (BYTE *)dpHit->buff;
	//SHGHITS *gHits = (SHGHITS *)buffer;

	for (size_t m = 0; m <= mApp->worker.moments.size(); m++){
		if (m == 0) fprintf(file, " moment 0 (Constant Flow){\n");
		else fprintf(file, " moment %d (%g s){\n", m, mApp->worker.moments[m - 1]);
		// Facets

		for (int i = 0; i < sh.nbFacet; i++) {
			Facet *f = facets[i];


			if (f->selected) {
				AHIT *hits = NULL;
				VHIT *dirs = NULL;
				fprintf(file, "FACET%d\n", i + 1);

				if (f->mesh || f->sh.countDirection) {
					char tmp[256];
					double dCoef = 1.0;
					if (!buffer) return;
					SHGHITS *shGHit = (SHGHITS *)buffer;
					int nbMoments = (int)mApp->worker.moments.size();
					int profSize = (f->sh.isProfile) ? (PROFILE_SIZE*sizeof(APROFILE)*(1 + nbMoments)) : 0;
					int w = f->sh.texWidth;
					int h = f->sh.texHeight;
					int tSize = w*h*sizeof(AHIT);
					int dSize = w*h*sizeof(VHIT);
					if (f->mesh) hits = (AHIT *)((BYTE *)buffer + (f->sh.hitOffset + sizeof(SHHITS) + profSize + m*tSize));
					if (f->sh.countDirection) dirs = (VHIT *)((BYTE *)buffer + (f->sh.hitOffset + sizeof(SHHITS) + profSize*(1 + nbMoments) + tSize*(1 + nbMoments) + m*dSize));


					switch (mode) {

					case 1: // Element area
						for (int j = 0; j < h; j++) {
							for (int i = 0; i < w; i++) {
								sprintf(tmp, "%g", f->mesh[i + j*w].area);
								if (tmp) fprintf(file, "%s", tmp);
								if (j < w - 1)
									fprintf(file, "\t");
							}
							fprintf(file, "\n");
						}
						break;

					case 2: //MC Hits
						for (int j = 0; j < h; j++) {
							for (int i = 0; i < w; i++) {
								fprintf(file, "%g", (double)hits[i + j*w].count);
								fprintf(file, "\t");
							}
							fprintf(file, "\n");
						}
						break;

					case 3: //Impingement rate
						dCoef = /*totalInFlux*/ 1.0 / shGHit->total.hit.nbDesorbed * 1E4; //1E4: conversion m2->cm2
						if (shGHit->mode == MC_MODE) dCoef *= (mApp->worker.displayedMoment == 0)
							? mApp->worker.finalOutgassingRate : (mApp->worker.totalDesorbedMolecules / mApp->worker.timeWindowSize);
						for (int j = 0; j < h; j++) {
							for (int i = 0; i < w; i++) {
								fprintf(file, "%g", (double)hits[i + j*w].count / (f->mesh[i + j*w].area*(f->sh.is2sided ? 2.0 : 1.0))*dCoef);
								fprintf(file, "\t");
							}
							fprintf(file, "\n");
						}
						break;

					case 4: //Particle density
					{
						dCoef = /*totalInFlux*/ 1.0 / shGHit->total.hit.nbDesorbed * 1E4; //1E4: conversion m2->cm2
						if (shGHit->mode == MC_MODE) dCoef *= (mApp->worker.displayedMoment == 0)
							? mApp->worker.finalOutgassingRate : (mApp->worker.totalDesorbedMolecules / mApp->worker.timeWindowSize);
						for (int j = 0; j < h; j++) {
							for (int i = 0; i < w; i++) {
								double v_ort_avg = 2.0*(double)hits[i + j*w].count / hits[i + j*w].sum_1_per_ort_velocity;
								double imp_rate = hits[i + j*w].count / (f->mesh[i + j*w].area*(f->sh.is2sided ? 2.0 : 1.0))*dCoef;
								double rho = 2.0*imp_rate / v_ort_avg;
								fprintf(file, "%g", rho);
								fprintf(file, "\t");
							}
							fprintf(file, "\n");
						}
						break;
					}
					case 5: //Gas density
					{
						dCoef = /*totalInFlux*/ 1.0 / shGHit->total.hit.nbDesorbed * 1E4; //1E4: conversion m2->cm2
						if (shGHit->mode == MC_MODE) dCoef *= (mApp->worker.displayedMoment == 0)
							? mApp->worker.finalOutgassingRate : (mApp->worker.totalDesorbedMolecules / mApp->worker.timeWindowSize);
						for (int j = 0; j < h; j++) {
							for (int i = 0; i < w; i++) {
								double v_ort_avg = 2.0*(double)hits[i + j*w].count / hits[i + j*w].sum_1_per_ort_velocity;
								double imp_rate = hits[i + j*w].count / (f->mesh[i + j*w].area*(f->sh.is2sided ? 2.0 : 1.0))*dCoef;
								double rho = 2.0*imp_rate / v_ort_avg;
								double rho_mass = rho*mApp->worker.gasMass / 1000.0 / 6E23;
								fprintf(file, "%g", rho_mass);
								fprintf(file, "\t");
							}
							fprintf(file, "\n");
						}
						break;
					}
					case 6:  // Pressure [mbar]

						// Lock during update
						dCoef = /*totalInFlux*/ 1.0 / shGHit->total.hit.nbDesorbed * 1E4 * (mApp->worker.gasMass / 1000 / 6E23) *0.0100;  //1E4 is conversion from m2 to cm2, 0.01: Pa->mbar
						if (shGHit->mode == MC_MODE) dCoef *= (mApp->worker.displayedMoment == 0)
							? mApp->worker.finalOutgassingRate : (mApp->worker.totalDesorbedMolecules / mApp->worker.timeWindowSize);
						for (int j = 0; j < h; j++) {
							for (int i = 0; i < w; i++) {
								fprintf(file, "%g", hits[i + j*w].sum_v_ort_per_area*dCoef);
								fprintf(file, "\t");
							}
							fprintf(file, "\n");
						}
						break;


					case 7: // Average velocity
						for (int j = 0; j < h; j++) {
							for (int i = 0; i < w; i++) {
								fprintf(file, "%g", 2.0*(double)hits[i + j*w].count / hits[i + j*w].sum_1_per_ort_velocity);
								fprintf(file, "\t");
							}
							fprintf(file, "\n");
						}
						break;

					case 8: // Velocity vector
						for (int j = 0; j < h; j++) {
							for (int i = 0; i < w; i++) {
								if (f->sh.countDirection) {
									sprintf(tmp, "%g,%g,%g",
										dirs[i + j*w].sumDir.x,
										dirs[i + j*w].sumDir.y,
										dirs[i + j*w].sumDir.z);
								}
								else {
									sprintf(tmp, "Direction not recorded");
								}
								fprintf(file, "%s", tmp);
								fprintf(file, "\t");
							}
							fprintf(file, "\n");
						}
						break;

					case 9: // Velocity vector Count
						for (int j = 0; j < h; j++) {
							for (int i = 0; i < w; i++) {
								if (f->sh.countDirection) {
									sprintf(tmp, "%I64d", dirs[i + j*w].count);
								}
								else {
									sprintf(tmp, "None");
								}
								fprintf(file, "%s", tmp);
								fprintf(file, "\t");
							}
							fprintf(file, "\n");
						}
						break;
					}
				}
				else {
					fprintf(file, "No mesh.\n");
				}
				fprintf(file, "\n"); //Current facet exported.
			}
		}
		fprintf(file, " }\n");
	}
	ReleaseDataport(dpHit);
}

void Geometry::ExportProfiles(FILE *file, int isTXT, Dataport *dpHit, Worker *worker) {

	char sep = isTXT ? '\t' : ',';
	//if(!IsLoaded()) throw Error("Nothing to save !");

	// Block dpHit during the whole disc writing
	BYTE *buffer = NULL;
	if (dpHit)
		if (AccessDataport(dpHit))
			buffer = (BYTE *)dpHit->buff;

	static const char* profType[] = {
		"None",
		"Pressure U [mbar]",
		"Pressure V [mbar]",
		"Incident angle [deg]",
		"Speed [m/s]",
		"Ort. velocity [m/s]" };

	// Globals
	//BYTE *buffer = (BYTE *)dpHit->buff;
	//SHGHITS *gHits = (SHGHITS *)buffer;
	std::ostringstream header;
	header << "Facet number" << sep << "Profile_type" << sep << "O_x" << sep << "O_y" << sep << "O_z" << sep << "U_x" << sep << "U_y" << sep << "U_z" << sep;
	header << "V_x" << sep << "V_y" << sep << "V_z" << sep << "U_length" << sep << "V_length" << sep << "Center_x" << sep << "Center_y" << sep << "Center_z" << sep << "V_max" << sep << "Hits" << sep;
	for (int i = 0; i < PROFILE_SIZE; i++)
		header << i + 1 << sep;
	header << '\n';

	fputs(header.str().c_str(), file);

	for (size_t m = 0; m <= mApp->worker.moments.size(); m++){
		if (m == 0) fputs(" moment 0 (Constant Flow){\n", file);
		else fprintf(file, " moment %d (%g s){\n", m, mApp->worker.moments[m - 1]);
		// Facets

		for (int i = 0; i < sh.nbFacet; i++) {
			Facet *f = facets[i];


			if (f->selected) {
				APROFILE *prof = NULL;
				std::ostringstream line;

				line << i + 1 << sep << profType[f->sh.profileType] << sep << f->sh.O.x << sep << f->sh.O.y << sep << f->sh.O.z << sep << f->sh.U.x << sep << f->sh.U.y << sep << f->sh.U.z << sep;
				line << f->sh.V.x << sep << f->sh.V.y << sep << f->sh.V.z << sep << Norme(&f->sh.U) << sep << Norme(&f->sh.V) << sep << f->sh.center.x << sep << f->sh.center.y << sep << f->sh.center.z << sep << f->sh.maxSpeed << sep << f->sh.counter.hit.nbHit << sep;


				if (f->sh.isProfile) {
					double dCoef = 1.0;
					if (!buffer) return;
					SHGHITS *shGHit = (SHGHITS *)buffer;
					double nbDes = (shGHit->total.hit.nbDesorbed > 0) ? (double)shGHit->total.hit.nbDesorbed : 1.0;
					int profOffset = PROFILE_SIZE*sizeof(APROFILE)*m;
					prof = (APROFILE*)((BYTE *)buffer + (f->sh.hitOffset + sizeof(SHHITS) + profOffset));
					double scaleX, scaleY;
					switch (f->sh.profileType) {
					case REC_PRESSUREU:
					case REC_PRESSUREV:
						scaleY = 1.0 / nbDes / (f->sh.area / (double)PROFILE_SIZE*1E-4)* worker->gasMass / 1000 / 6E23 * 0.0100; //0.01: Pa->mbar
						/*scaleY *= ((worker->displayedMoment == 0) ? 1.0 : ((worker->desorptionStopTime - worker->desorptionStartTime)
						/ worker->timeWindowSize)); //correction for time window length*/
						scaleY *= ((worker->displayedMoment == 0) ? worker->finalOutgassingRate : (worker->totalDesorbedMolecules
							/ worker->timeWindowSize));
						if (f->sh.is2sided) scaleY *= 0.5;
						//if(f->sh.opacity>0.0) scaleY *= f->sh.opacity;
						//if(IS_ZERO(f->sh.opacity)) scaleY*=2; //transparent profiles are profiled only once...

						for (int j = 0; j < PROFILE_SIZE; j++)
							line << prof[j].sum_v_ort*scaleY << sep;
						break;
					case REC_VELOCITY:
					case REC_ORT_VELOCITY:
						scaleX = f->sh.maxSpeed / (double)PROFILE_SIZE;
						for (int j = 0; j < PROFILE_SIZE; j++)
							line << (double)prof[j].count / f->sh.counter.hit.nbHit << sep;
						break;
					case REC_ANGULAR:
						scaleX = 90.0 / (double)PROFILE_SIZE;
						for (int j = 0; j < PROFILE_SIZE; j++)
							line << (double)prof[j].count / f->sh.counter.hit.nbHit << sep;
						break;
					}
				}
				else {
					line << "No profile.";
				}
				line << '\n';
				fputs(line.str().c_str(), file);
			}
		}
		fputs(" }\n", file);
	}
	ReleaseDataport(dpHit);
}

void Geometry::ImportDesorption_DES(FileReader *file) {

	//if(!IsLoaded()) throw Error("Nothing to save !");

	// Block dpHit during the whole disc writing

	for (int i = 0; i < sh.nbFacet; i++) { //clear previous desorption maps
		facets[i]->hasOutgassingFile = FALSE;
		facets[i]->sh.useOutgassingFile = FALSE;
		facets[i]->sh.desorbType = DES_NONE; //clear previously set desorptions
		facets[i]->selected = FALSE;
		facets[i]->UnselectElem();
	}
	// Facets
	while (!file->IsEof()) {
		file->ReadKeyword("facet");
		int facetId = file->ReadInt() - 1;
		file->ReadKeyword("{");
		if (!(facetId >= 0 && facetId < sh.nbFacet)) {
			file->MakeError("Invalid facet Id (loaded desorption file for a different geometry?)");
			return;
		}

		Facet *f = facets[facetId];
		f->hasOutgassingFile = TRUE;
		f->sh.useOutgassingFile = TRUE; //turn on file usage by default
		f->sh.desorbType = DES_COSINE; //auto-set to cosine
		Select(f);
		file->ReadKeyword("cell_size_cm"); file->ReadKeyword(":");
		double ratio = f->sh.outgassingFileRatio = file->ReadDouble();
		if (f->sh.outgassingFileRatio != 0.0) {
			f->sh.outgassingFileRatio = 1.0 / f->sh.outgassingFileRatio; //cell size -> samples per cm
			ratio = f->sh.outgassingFileRatio;
		}
		double nU = Norme(&(f->sh.U));
		double nV = Norme(&(f->sh.V));
		int w = f->sh.outgassingMapWidth = (int)ceil(nU*ratio); //double precision written to file
		int h = f->sh.outgassingMapHeight = (int)ceil(nV*ratio); //double precision written to file
		f->outgassingMap = (double*)malloc(w*h*sizeof(double));
		if (!f->outgassingMap) throw Error("Not enough memory to store outgassing map.");
		for (int i = 0; i < w; i++) {
			for (int j = 0; j < h; j++) {
				f->outgassingMap[i + j*w] = file->ReadDouble();
			}
		}
		file->ReadKeyword("}");
	}
	UpdateSelection();
	//InitializeGeometry();

	//Debug memory check
	//_ASSERTE (!_CrtDumpMemoryLeaks());;
	_ASSERTE(_CrtCheckMemory());
}

void Geometry::SaveSTR(Dataport *dpHit, BOOL saveSelected) {

	if (!IsLoaded()) throw Error("Nothing to save !");
	if (sh.nbSuper < 1) throw Error("Cannot save single structure in STR format");

	// Block dpHit during the whole disc writting
	AccessDataport(dpHit);
	for (int i = 0; i < sh.nbSuper; i++)
		SaveSuper(dpHit, i);
	ReleaseDataport(dpHit);

}

void Geometry::SaveSuper(Dataport *dpHit, int s) {

	char fName[512];
	sprintf(fName, "%s/%s", strPath, strFileName[s]);
	FileWriter *file = new FileWriter(fName);

	// Unused
	file->WriteInt(0, "\n");

	// Globals
	BYTE *buffer = (BYTE *)dpHit->buff;
	SHGHITS *gHits = (SHGHITS *)buffer;

	//Extract data of the specified super structure
	llong totHit = 0;
	llong totAbs = 0;
	llong totDes = 0;
	int *refIdx = (int *)malloc(sh.nbVertex*sizeof(int));
	memset(refIdx, 0xFF, sh.nbVertex*sizeof(int));
	int nbV = 0;
	int nbF = 0;

	for (int i = 0; i < sh.nbFacet; i++) {
		Facet *f = facets[i];
		if (f->sh.superIdx == s) {
			totHit += f->sh.counter.hit.nbHit;
			totAbs += f->sh.counter.hit.nbAbsorbed;
			totDes += f->sh.counter.hit.nbDesorbed;
			for (int j = 0; j < f->sh.nbIndex; j++)
				refIdx[f->indices[j]] = 1;
			nbF++;
		}
	}

	for (int i = 0; i < sh.nbVertex; i++) {
		if (refIdx[i] >= 0) {
			refIdx[i] = nbV;
			nbV++;
		}
	}

	file->WriteLLong(totHit, "\n");
	file->WriteLLong(gHits->nbLeakTotal, "\n");
	file->WriteLLong(totDes, "\n");
	file->WriteLLong(tNbDesorptionMax, "\n");

	file->WriteInt(nbV, "\n");
	file->WriteInt(nbF, "\n");

	// Read geometry vertices
	for (int i = 0; i < sh.nbVertex; i++) {
		if (refIdx[i] >= 0) {
			file->WriteDouble(vertices3[i].x, " ");
			file->WriteDouble(vertices3[i].y, " ");
			file->WriteDouble(vertices3[i].z, "\n");
		}
	}

	// Facets
	for (int i = 0; i < sh.nbFacet; i++) {
		Facet *f = facets[i];
		int j;
		if (f->sh.superIdx == s) {
			file->WriteInt(f->sh.nbIndex, " ");
			for (j = 0; j < f->sh.nbIndex - 1; j++)
				file->WriteInt(refIdx[f->indices[j]] + 1, " ");
			file->WriteInt(refIdx[f->indices[j]] + 1, "\n");
		}
	}

	// Params
	for (int i = 0; i < sh.nbFacet; i++) {

		// Update facet hits from shared mem
		Facet *f = facets[i];
		if (f->sh.superIdx == s) {
			SHHITS *shF = (SHHITS *)(buffer + f->sh.hitOffset);
			memcpy(&(f->sh.counter), shF, sizeof(SHHITS));
			f->SaveTXT(file);
		}

	}

	SaveProfileTXT(file);

	SAFE_DELETE(file);
	free(refIdx);

}

void Geometry::RemoveFromStruct(int numToDel) {
	mApp->changedSinceSave = TRUE;
	int nb = 0;
	for (int i = 0; i < sh.nbFacet; i++)
		if (facets[i]->sh.superIdx == numToDel) nb++;

	if (nb == 0) return;
	/*
	if(sh.nbFacet-nb==0) {
	// Remove all
	Clear();
	return;
	}
	*/

	Facet   **f = (Facet **)malloc((sh.nbFacet - nb) * sizeof(Facet *));

	nb = 0;
	for (int i = 0; i < sh.nbFacet; i++) {
		if (facets[i]->sh.superIdx == numToDel) {
			delete facets[i];
			mApp->RenumberSelections(nb);
			mApp->RenumberFormulas(nb);
		}
		else {
			f[nb++] = facets[i];
		}
	}

	SAFE_FREE(facets);
	facets = f;
	sh.nbFacet = nb;
	CalcTotalOutGassing();

}

/*BOOL AskToReset_Geom(Worker *work) {

if (work->nbHit>0) {
int rep = GLMessageBox::Display("This will reset simulation data.","Geometry change",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONWARNING);
if( rep != GLDLG_OK ) {
return FALSE;
}
}
work->Reset(m_fTime);
mApp->nbDesStart=0;
mApp->nbHitStart = 0;

if(mApp->profilePlotter) mApp->profilePlotter->Update(m_fTime,TRUE);
if(mApp->texturePlotter) mApp->texturePlotter->Update(m_fTime,TRUE);
return TRUE;

}*/

BOOL Geometry::IsLoaded() {
	return isLoaded;
}

void Geometry::ImportDesorption_SYN(
	FileReader *file, const size_t &source, const double &time,
	const size_t &mode, const double &eta0, const double &alpha, const double &cutoffdose,
	const std::vector<std::pair<double, double>> &convDistr,
	GLProgress *prg){


	//UnSelectAll();
	char tmp[512];
	std::vector<double> xdims, ydims;
	double no_scans = 1.0;

	file->ReadKeyword("version"); file->ReadKeyword(":");
	int version2;
	version2 = file->ReadInt();
	if (version2 > SYNVERSION) {
		char errMsg[512];
		sprintf(errMsg, "Unsupported SYN version V%d", version2);
		throw Error(errMsg);
	}

	//now read number of facets
	file->ReadKeyword("totalHit"); file->ReadKeyword(":");
	file->ReadLLong();
	file->ReadKeyword("totalDes"); file->ReadKeyword(":");
	file->ReadLLong();
	if (version2 >= 6) {
		file->ReadKeyword("no_scans"); file->ReadKeyword(":");
		no_scans = file->ReadDouble();
	}
	file->ReadKeyword("totalLeak"); file->ReadKeyword(":");
	file->ReadLLong();
	if (version2 > 2) {
		file->ReadKeyword("totalFlux"); file->ReadKeyword(":");
		file->ReadDouble();
		file->ReadKeyword("totalPower"); file->ReadKeyword(":");
		file->ReadDouble();
	}
	file->ReadKeyword("maxDes"); file->ReadKeyword(":");
	file->ReadLLong();
	file->ReadKeyword("nbVertex"); file->ReadKeyword(":");
	file->ReadInt();
	file->ReadKeyword("nbFacet"); file->ReadKeyword(":");
	int nbNewFacet = file->ReadInt(); //gotcha! :)
	xdims.reserve(nbNewFacet);
	ydims.reserve(nbNewFacet);

	//now go for the facets to get their texture ratio
	for (int i = 0; i < nbNewFacet && i < GetNbFacet(); i++) {
		prg->SetProgress(0.5*(double)i / (double)MIN(nbNewFacet, GetNbFacet()));
		file->JumpSection("facet");
		// Check idx
		int idx = file->ReadInt();
		if (idx != i + 1) throw Error(file->MakeError("Wrong facet index !"));

		file->JumpSection("texDimX"); file->ReadKeyword(":");
		xdims.push_back(file->ReadDouble());
		//if (!IS_ZERO(xdims[i])) Select(GetFacet(i));
		file->ReadKeyword("texDimY"); file->ReadKeyword(":");
		ydims.push_back(file->ReadDouble());
	}

	//now read actual textures
	//read header
	file->SeekFor("{textures}");
	file->ReadKeyword("minHit_MC"); file->ReadKeyword(":");
	file->ReadLLong();
	file->ReadKeyword("maxHit_MC"); file->ReadKeyword(":");
	file->ReadLLong();
	file->ReadKeyword("minHit_flux"); file->ReadKeyword(":");
	file->ReadDouble();
	file->ReadKeyword("maxHit_flux"); file->ReadKeyword(":");
	file->ReadDouble();
	file->ReadKeyword("minHit_power"); file->ReadKeyword(":");
	file->ReadDouble();
	file->ReadKeyword("maxHit_power"); file->ReadKeyword(":");
	file->ReadDouble();

	//read texture values
	for (int i = 0; i < nbNewFacet && i < GetNbFacet(); i++) {
		prg->SetProgress(0.5 + 0.5*(double)i / (double)MIN(nbNewFacet, GetNbFacet()));
		if (!IS_ZERO(xdims[i])) { //has texture
			Facet *f = GetFacet(i);
			if (f->selected) {
				f->hasOutgassingFile = TRUE;
				f->sh.useOutgassingFile = TRUE; //turn on file usage by default
				f->sh.desorbType = DES_COSINE; //auto-set to cosine
			}
			//Select(f);
			file->ReadKeyword("texture_facet");
			// Check idx
			int idx = file->ReadInt();

			if (idx != i + 1) {
				sprintf(tmp, "Wrong facet index. Expected %d, read %d.", i + 1, idx);
				throw Error(file->MakeError(tmp));
			}

			//Now load values
			file->ReadKeyword("{");

			int ix, iy;
			f->sh.outgassingMapWidth = (int)(xdims[i] - 1e-9) + 1;
			f->sh.outgassingMapHeight = (int)(ydims[i] - 1e-9) + 1;

			if (f->selected) {
				f->sh.outgassingFileRatio = xdims[i] / Norme(&(f->sh.U));
				f->outgassingMap = (double*)malloc(f->sh.outgassingMapWidth*f->sh.outgassingMapHeight*sizeof(double));
				if (!f->outgassingMap) throw Error("Not enough memory to store outgassing map.");
				memset(f->outgassingMap, 0, f->sh.outgassingMapWidth*f->sh.outgassingMapHeight*sizeof(double)); //set inital values to zero
				f->totalDose = f->sh.totalOutgassing = f->totalFlux = 0.0;
			}

			int texWidth_file, texHeight_file;
			//In case of rounding errors, the file might contain different texture dimensions than expected.
			if (version2 >= 8) {
				file->ReadKeyword("width"); file->ReadKeyword(":"); texWidth_file = file->ReadInt();
				file->ReadKeyword("height"); file->ReadKeyword(":"); texHeight_file = file->ReadInt();
			}
			else {
				texWidth_file = f->sh.outgassingMapWidth;
				texHeight_file = f->sh.outgassingMapHeight;
			}

			for (iy = 0; iy<(MIN(f->sh.outgassingMapHeight, texHeight_file)); iy++) { //MIN: If stored texture is larger, don't read extra cells
				for (ix = 0; ix<(MIN(f->sh.outgassingMapWidth, texWidth_file)); ix++) { //MIN: If stored texture is larger, don't read extra cells
					int index = iy*f->sh.outgassingMapWidth + ix;
					//Read original values
					llong MC = file->ReadLLong();
					double cellArea = 1.0;
					if (version2 >= 7) cellArea = file->ReadDouble();
					if (cellArea < 1E-10) cellArea = 1.0; //to avoid division by zero
					double flux = file->ReadDouble() / no_scans; //not normalized by cell area
					double power = file->ReadDouble() / no_scans; //not normalized by cell area

					if (f->selected) {
						//Calculate dose
						double dose;
						if (source == 0) dose = (double)MC*time;
						else if (source == 1) dose = flux*time / cellArea;
						else if (source == 2) dose = power*time / cellArea;

						double outgassing;
						if (dose == 0) outgassing = 0; //to avoid division by zero later
						else {
							//Convert to outgassing

							if (mode == 0) {
								if (source == 0) outgassing = (double)MC * 0.100 / 1.38E-23 / f->sh.temperature;
								else if (source == 1) outgassing = flux * 0.100 / 1.38E-23 / f->sh.temperature; //Division by 10 because the user will want to see the same outgassing in mbar*l/s
								else if (source == 2) outgassing = power * 0.100 / 1.38E-23 / f->sh.temperature; //(Outgassing is stored internally in Pa*m3/s, for consistent SI unit calculations)
							}
							else if (mode == 1) {
								double moleculePerPhoton = eta0*pow(MAX(1.0,dose/cutoffdose),alpha);
								outgassing = flux*moleculePerPhoton;
							}
							else if (mode == 2) {
								double moleculePerPhoton = InterpolateY(dose, convDistr, FALSE, TRUE);
								outgassing = flux*moleculePerPhoton;
							}
						}
						//Apply outgassing
						//f->outgassingMap[index] = outgassing *0.100; //0.1: mbar*l/s->Pa*m3/s
						f->outgassingMap[index] = outgassing * 1.38E-23 * f->sh.temperature; //1[Pa*m3/s] = kT [particles/sec]

						//Facet diagnostic info
						f->totalDose += flux*time;
						f->totalFlux += flux;
						f->sh.totalOutgassing += f->outgassingMap[index];
					} //if selected
				}
				for (int ie = 0; ie < texWidth_file - f->sh.outgassingMapWidth; ie++) {//Executed if file texture is bigger than expected texture
					//Read extra cells from file without doing anything
					//Read original values
					file->ReadLLong(); //MC
					if (version2 >= 7) file->ReadDouble(); //area
					file->ReadDouble(); //flux
					file->ReadDouble(); //power
				}
			}
			for (int ie = 0; ie < texHeight_file - f->sh.outgassingMapHeight; ie++) {//Executed if file texture is bigger than expected texture
				//Read extra cells ffrom file without doing anything
				for (int iw = 0; iw < texWidth_file; iw++) {
					//Read original values
					file->ReadLLong(); //MC
					if (version2 >= 7) file->ReadDouble(); //area
					file->ReadDouble(); //flux
					file->ReadDouble(); //power
				}
			}
			file->ReadKeyword("}");


		} //if has texture
	}

	//end
	//UpdateSelection();

}

void Geometry::AnalyzeSYNfile(FileReader *file, GLProgress *progressDlg, int *nbNewFacet,
	int *nbTextured, int *nbDifferent, GLProgress *prg){
	//init
	*nbTextured = 0;
	*nbNewFacet = 0;
	*nbDifferent = 0;

	UnSelectAll();
	//char tmp[512];

	file->ReadKeyword("version"); file->ReadKeyword(":");
	int version2;
	version2 = file->ReadInt();
	if (version2 > SYNVERSION) {
		char errMsg[512];
		sprintf(errMsg, "Unsupported SYN version V%d", version2);
		throw Error(errMsg);
	}

	//now read number of facets
	file->ReadKeyword("totalHit"); file->ReadKeyword(":");
	file->ReadLLong();
	file->ReadKeyword("totalDes"); file->ReadKeyword(":");
	file->ReadLLong();
	if (version2 >= 6) {
		file->ReadKeyword("no_scans"); file->ReadKeyword(":");
		/*no_scans = */file->ReadDouble();
	}
	file->ReadKeyword("totalLeak"); file->ReadKeyword(":");
	file->ReadLLong();
	if (version2 > 2) {
		file->ReadKeyword("totalFlux"); file->ReadKeyword(":");
		file->ReadDouble();
		file->ReadKeyword("totalPower"); file->ReadKeyword(":");
		file->ReadDouble();
	}
	file->ReadKeyword("maxDes"); file->ReadKeyword(":");
	file->ReadLLong();
	file->ReadKeyword("nbVertex"); file->ReadKeyword(":");
	file->ReadInt();
	file->ReadKeyword("nbFacet"); file->ReadKeyword(":");
	*nbNewFacet = file->ReadInt(); //gotcha! :)

	//now go for the facets to get their texture ratio, etc.
	for (int i = 0; i < *nbNewFacet && i < GetNbFacet(); i++) {
		prg->SetProgress((double)i / (double)MIN(*nbNewFacet, GetNbFacet()));
		file->JumpSection("facet");
		// Check idx
		int idx = file->ReadInt();
		if (idx != i + 1) throw Error(file->MakeError("Wrong facet index !"));

		file->JumpSection("mesh"); file->ReadKeyword(":");
		if (file->ReadInt()) { //has mesh
			(*nbTextured)++;
			Select(GetFacet(i));
			/*file->ReadKeyword("texDimX");file->ReadKeyword(":");
			if ((this->GetFacet(i)->sh.texWidthD-file->ReadDouble())>1E-8) {
			(*nbDifferent)++;
			continue;
			}
			file->ReadKeyword("texDimY");file->ReadKeyword(":");
			if ((this->GetFacet(i)->sh.texHeightD-file->ReadDouble())>1E-8) {
			(*nbDifferent)++;
			}*/
		}
	}
	UpdateSelection();
}

void Geometry::SaveXML_geometry(pugi::xml_node saveDoc, Worker *work, GLProgress *prg, BOOL saveSelected){
	//TiXmlDeclaration* decl = new TiXmlDeclaration("1.0")="")="");
	//saveDoc->LinkEndChild(decl);

	xml_node geomNode = saveDoc.append_child("Geometry");

	prg->SetMessage("Writing vertices...");
	geomNode.append_child("Vertices").append_attribute("nb") = sh.nbVertex; //creates Vertices node, adds nb attribute and sets its value to sh.nbVertex
	for (int i = 0; i < sh.nbVertex; i++) {
		prg->SetProgress(0.166*((double)i / (double)sh.nbVertex));
		xml_node v = geomNode.child("Vertices").append_child("Vertex");
		v.append_attribute("id") = i;
		v.append_attribute("x") = vertices3[i].x;
		v.append_attribute("y") = vertices3[i].y;
		v.append_attribute("z") = vertices3[i].z;
	}

	prg->SetMessage("Writing facets...");
	geomNode.append_child("Facets");
	geomNode.child("Facets").append_attribute("nb") = sh.nbFacet;
	for (int i = 0, k = 0; i < sh.nbFacet; i++) {
		prg->SetProgress(0.166 + ((double)i / (double)sh.nbFacet) *0.166);
		if (!saveSelected || facets[i]->selected) {
			xml_node f = geomNode.child("Facets").append_child("Facet");
			f.append_attribute("id") = i;
			facets[i]->SaveXML_geom(f);
		}
	}

	prg->SetMessage("Writing model details...");
	geomNode.append_child("Structures").append_attribute("nb") = sh.nbSuper;
	for (int i = 0, k = 0; i < sh.nbSuper; i++) {
		xml_node s = geomNode.child("Structures").append_child("Structure");
		s.append_attribute("id") = i;
		s.append_attribute("name") = strName[i];
	}

	xml_node interfNode = saveDoc.append_child("Interface");

	xml_node selNode = interfNode.append_child("Selections");
	selNode.append_attribute("nb") = (!saveSelected)*(mApp->nbSelection);
	for (int i = 0; (i < mApp->nbSelection) && !saveSelected; i++) { //don't save selections when exporting part of the geometry (saveSelected)
		xml_node newSel = selNode.append_child("Selection");
		newSel.append_attribute("id") = i;
		newSel.append_attribute("name") = mApp->selections[i].name;
		newSel.append_attribute("nb") = mApp->selections[i].nbSel;
		for (int j = 0; j < mApp->selections[i].nbSel; j++) {
			xml_node newItem = newSel.append_child("selItem");
			newItem.append_attribute("id") = j;
			newItem.append_attribute("facet") = mApp->selections[i].selection[j];
		}
	}

	xml_node viewNode = interfNode.append_child("Views");
	viewNode.append_attribute("nb") = (!saveSelected)*(mApp->nbView);
	for (int i = 0; (i < mApp->nbView) && !saveSelected; i++) { //don't save views when exporting part of the geometry (saveSelected)
		xml_node newView = viewNode.append_child("View");
		newView.append_attribute("id") = i;
		newView.append_attribute("name") = mApp->views[i].name;
		newView.append_attribute("projMode") = mApp->views[i].projMode;
		newView.append_attribute("camAngleOx") = mApp->views[i].camAngleOx;
		newView.append_attribute("camAngleOy") = mApp->views[i].camAngleOy;
		newView.append_attribute("camDist") = mApp->views[i].camDist;
		newView.append_attribute("camOffset.x") = mApp->views[i].camOffset.x;
		newView.append_attribute("camOffset.y") = mApp->views[i].camOffset.y;
		newView.append_attribute("camOffset.z") = mApp->views[i].camOffset.z;
		newView.append_attribute("performXY") = mApp->views[i].performXY;
		newView.append_attribute("vLeft") = mApp->views[i].vLeft;
		newView.append_attribute("vRight") = mApp->views[i].vRight;
		newView.append_attribute("vTop") = mApp->views[i].vTop;
		newView.append_attribute("vBottom") = mApp->views[i].vBottom;
	}

	xml_node formulaNode = interfNode.append_child("Formulas");
	formulaNode.append_attribute("nb") = (!saveSelected)*(mApp->nbFormula);
	for (int i = 0; (i < mApp->nbFormula) && !saveSelected; i++) { //don't save formulas when exporting part of the geometry (saveSelected)
		xml_node newFormula = formulaNode.append_child("Formula");
		newFormula.append_attribute("id") = i;
		newFormula.append_attribute("name") = mApp->formulas[i].parser->GetName();
		newFormula.append_attribute("expression") = mApp->formulas[i].parser->GetExpression();
	}

	if (mApp->profilePlotter) {
		std::vector<int> ppViews = mApp->profilePlotter->GetViews();
		xml_node profilePlotterNode = interfNode.append_child("ProfilePlotter");
		profilePlotterNode.append_child("Parameters").append_attribute("logScale") = mApp->profilePlotter->IsLogScaled();
		xml_node viewsNode = profilePlotterNode.append_child("Views");
		for (int v : ppViews) {
			xml_node view = viewsNode.append_child("View");
			view.append_attribute("facetId") = v;
		}
	}

	xml_node simuParamNode = saveDoc.append_child("MolflowSimuSettings");

	simuParamNode.append_child("Gas").append_attribute("mass") = work->gasMass;
	simuParamNode.child("Gas").append_attribute("halfLife") = work->halfLife;

	xml_node timeSettingsNode = simuParamNode.append_child("TimeSettings");

	xml_node userMomentsNode = timeSettingsNode.append_child("UserMoments");
	userMomentsNode.append_attribute("nb") = work->userMoments.size();
	for (size_t i = 0; i < work->userMoments.size(); i++) {
		xml_node newUserEntry = userMomentsNode.append_child("UserEntry");
		newUserEntry.append_attribute("id") = i;
		newUserEntry.append_attribute("content") = work->userMoments[i].c_str();
	}

	timeSettingsNode.append_attribute("timeWindow") = work->timeWindowSize;
	timeSettingsNode.append_attribute("useMaxwellDistr") = work->useMaxwellDistribution;
	timeSettingsNode.append_attribute("calcConstFlow") = work->calcConstantFlow;

	xml_node motionNode = simuParamNode.append_child("Motion");
	motionNode.append_attribute("type") = work->motionType;
	if (work->motionType == 1) { //fixed motion
		xml_node v = motionNode.append_child("VelocityVector");
		v.append_attribute("vx") = work->motionVector2.x;
		v.append_attribute("vy") = work->motionVector2.y;
		v.append_attribute("vz") = work->motionVector2.z;
	}
	else if (work->motionType == 2) { //rotation
		xml_node v = motionNode.append_child("AxisBasePoint");
		v.append_attribute("x") = work->motionVector1.x;
		v.append_attribute("y") = work->motionVector1.y;
		v.append_attribute("z") = work->motionVector1.z;
		xml_node v2 = motionNode.append_child("RotationVector");
		v2.append_attribute("x") = work->motionVector2.x;
		v2.append_attribute("y") = work->motionVector2.y;
		v2.append_attribute("z") = work->motionVector2.z;
	}

	xml_node paramNode = simuParamNode.append_child("Parameters");
	paramNode.append_attribute("nb") = work->parameters.size();
	for (size_t i = 0; i < work->parameters.size(); i++) {
		xml_node newParameter = paramNode.append_child("Parameter");
		newParameter.append_attribute("id") = i;
		newParameter.append_attribute("name") = work->parameters[i].name.c_str();
		newParameter.append_attribute("nbMoments") = (int)work->parameters[i].values.size();
		for (size_t m = 0; m < work->parameters[i].values.size(); m++) {
			xml_node newMoment = newParameter.append_child("Moment");
			newMoment.append_attribute("id") = m;
			newMoment.append_attribute("t") = work->parameters[i].values[m].first;
			newMoment.append_attribute("value") = work->parameters[i].values[m].second;
		}
	}
}

BOOL Geometry::SaveXML_simustate(xml_node saveDoc, Worker *work, BYTE *buffer, SHGHITS *gHits, int nbLeakSave, int nbHHitSave,
	LEAK *pLeak, HIT *pHits, GLProgress *prg, BOOL saveSelected){
	xml_node resultNode = saveDoc.append_child("MolflowResults");
	prg->SetMessage("Writing simulation results...");
	xml_node momentsNode = resultNode.append_child("Moments");
	momentsNode.append_attribute("nb") = work->moments.size() + 1;
	for (size_t m = 0; m <= mApp->worker.moments.size(); m++){
		prg->SetProgress(0.5 + (double)m / (1.0 + (double)mApp->worker.moments.size()));
		xml_node newMoment = momentsNode.append_child("Moment");
		newMoment.append_attribute("id") = m;
		if (m == 0)
			newMoment.append_attribute("time") = "Constant flow";
		else
			newMoment.append_attribute("time") = work->moments[m - 1];

		if (m == 0) { //Write global results. Later these results will probably be time-dependent as well.
			xml_node globalNode = newMoment.append_child("Global");

			xml_node hitsNode = globalNode.append_child("Hits");
			hitsNode.append_attribute("totalHit") = gHits->total.hit.nbHit;
			hitsNode.append_attribute("totalDes") = gHits->total.hit.nbDesorbed;
			hitsNode.append_attribute("totalAbs") = gHits->total.hit.nbAbsorbed;
			hitsNode.append_attribute("totalDist_total") = gHits->distTraveledTotal_total;
			hitsNode.append_attribute("totalDist_fullHitsOnly") = gHits->distTraveledTotal_fullHitsOnly;
			hitsNode.append_attribute("totalLeak") = gHits->nbLeakTotal;
			hitsNode.append_attribute("maxDesorption") = work->maxDesorption;

			xml_node hitCacheNode = globalNode.append_child("Hit_Cache");
			hitCacheNode.append_attribute("nb") = nbHHitSave;
			for (int i = 0; i < nbHHitSave; i++) {
				xml_node newHit = hitCacheNode.append_child("Hit");
				newHit.append_attribute("id") = i;
				newHit.append_attribute("posX") = pHits[i].pos.x;
				newHit.append_attribute("posY") = pHits[i].pos.y;
				newHit.append_attribute("posZ") = pHits[i].pos.z;
				newHit.append_attribute("type") = pHits[i].type;
			}

			xml_node leakCacheNode = globalNode.append_child("Leak_Cache");
			leakCacheNode.append_attribute("nb") = nbLeakSave;
			for (int i = 0; i < nbLeakSave; i++) {
				xml_node newLeak = leakCacheNode.append_child("Leak");
				newLeak.append_attribute("id") = i;
				newLeak.append_attribute("posX") = pLeak[i].pos.x;
				newLeak.append_attribute("posY") = pLeak[i].pos.y;
				newLeak.append_attribute("posZ") = pLeak[i].pos.z;
				newLeak.append_attribute("dirX") = pLeak[i].dir.x;
				newLeak.append_attribute("dirY") = pLeak[i].dir.y;
				newLeak.append_attribute("dirZ") = pLeak[i].dir.z;
			}
		} //end global node

		xml_node facetResultsNode = newMoment.append_child("FacetResults");

		for (int i = 0; i < sh.nbFacet; i++) {
			Facet *f = GetFacet(i);
			xml_node newFacetResult = facetResultsNode.append_child("Facet");
			newFacetResult.append_attribute("id") = i;
			if (m == 0) { //Now it's a global value, will soon become time-dependent
				xml_node facetHitNode = newFacetResult.append_child("Hits");
				facetHitNode.append_attribute("nbHit") = f->sh.counter.hit.nbHit;
				facetHitNode.append_attribute("nbDes") = f->sh.counter.hit.nbDesorbed;
				facetHitNode.append_attribute("nbAbs") = f->sh.counter.hit.nbAbsorbed;
				facetHitNode.append_attribute("sum_v_ort") = f->sh.counter.hit.sum_v_ort;
				facetHitNode.append_attribute("sum_1_per_v") = f->sh.counter.hit.sum_1_per_ort_velocity;
			}

			if (f->sh.isProfile){
				xml_node profileNode = newFacetResult.append_child("Profile");
				profileNode.append_attribute("size") = PROFILE_SIZE;
				APROFILE *profilePtr = (APROFILE *)(buffer + f->sh.hitOffset + sizeof(SHHITS) + m*sizeof(APROFILE)*PROFILE_SIZE);
				for (int p = 0; p < PROFILE_SIZE; p++){
					xml_node slice = profileNode.append_child("Slice");
					slice.append_attribute("id") = p;
					slice.append_attribute("count") = profilePtr[p].count;
					slice.append_attribute("sum_1_per_v") = profilePtr[p].sum_1_per_ort_velocity;
					slice.append_attribute("sum_v_ort") = profilePtr[p].sum_v_ort;
				}
			}

			int profSize = (f->sh.isProfile) ? (PROFILE_SIZE*sizeof(APROFILE)*(1 + (int)mApp->worker.moments.size())) : 0;
			int h = (f->sh.texHeight);
			int w = (f->sh.texWidth);

			if (f->hasMesh){
				xml_node textureNode = newFacetResult.append_child("Texture");
				textureNode.append_attribute("width") = f->sh.texWidth;
				textureNode.append_attribute("height") = f->sh.texHeight;

				AHIT *hits = (AHIT *)((BYTE *)gHits + (f->sh.hitOffset + sizeof(SHHITS) + profSize + m*w*h*sizeof(AHIT)));
				std::stringstream countText, sum1perText, sumvortText;
				countText << '\n'; //better readability in file
				sum1perText << '\n';
				sumvortText << '\n';
				for (int iy = 0; iy < h; iy++) {
					for (int ix = 0; ix < w; ix++) {
						countText << hits[iy*f->sh.texWidth + ix].count << '\t';
						sum1perText << hits[iy*f->sh.texWidth + ix].sum_1_per_ort_velocity << '\t';
						sumvortText << hits[iy*f->sh.texWidth + ix].sum_v_ort_per_area << '\t';
					}
					countText << '\n';
					sum1perText << '\n';
					sumvortText << '\n';
				}
				textureNode.append_child("count").append_child(node_cdata).set_value(countText.str().c_str());
				textureNode.append_child("sum_1_per_v").append_child(node_cdata).set_value(sum1perText.str().c_str());
				textureNode.append_child("sum_v_ort").append_child(node_cdata).set_value(sumvortText.str().c_str());

			} //end texture

			if (f->sh.countDirection && f->dirCache) {
				xml_node dirNode = newFacetResult.append_child("Directions");
				dirNode.append_attribute("width") = f->sh.texWidth;
				dirNode.append_attribute("height") = f->sh.texHeight;

				VHIT *dirs = (VHIT *)((BYTE *)gHits + f->sh.hitOffset + sizeof(SHHITS) + profSize + (1 + (int)work->moments.size())*w*h*sizeof(AHIT) + m*w*h*sizeof(VHIT));

				std::stringstream dirText, dirCountText;
				dirText << '\n'; //better readability in file
				dirCountText << '\n';

				for (int iy = 0; iy < h; iy++) {
					for (int ix = 0; ix < w; ix++) {
						dirText << dirs[iy*f->sh.texWidth + ix].sumDir.x << ",";
						dirText << dirs[iy*f->sh.texWidth + ix].sumDir.y << ",";
						dirText << dirs[iy*f->sh.texWidth + ix].sumDir.z << "\t";
						dirCountText << dirs[iy*f->sh.texWidth + ix].count << "\t";
					}
					dirText << "\n";
					dirCountText << "\n";
				}
				dirNode.append_child("vel.vectors").append_child(node_cdata).set_value(dirText.str().c_str());
				dirNode.append_child("count").append_child(node_cdata).set_value(dirCountText.str().c_str());
			} //end directions

		}

	}

	//Texture Min/Max
	xml_node minMaxNode = resultNode.append_child("TextureMinMax");
	minMaxNode.append_child("With_constant_flow").append_child("Pressure").append_attribute("min") = gHits->texture_limits[0].min.all;
	minMaxNode.child("With_constant_flow").child("Pressure").append_attribute("max") = gHits->texture_limits[0].max.all;
	minMaxNode.child("With_constant_flow").append_child("Density").append_attribute("min") = gHits->texture_limits[1].min.all;
	minMaxNode.child("With_constant_flow").child("Density").append_attribute("max") = gHits->texture_limits[1].max.all;
	minMaxNode.child("With_constant_flow").append_child("Imp.rate").append_attribute("min") = gHits->texture_limits[2].min.all;
	minMaxNode.child("With_constant_flow").child("Imp.rate").append_attribute("max") = gHits->texture_limits[2].max.all;

	minMaxNode.append_child("Moments_only").append_child("Pressure").append_attribute("min") = gHits->texture_limits[0].min.moments_only;
	minMaxNode.child("Moments_only").child("Pressure").append_attribute("max") = gHits->texture_limits[0].max.moments_only;
	minMaxNode.child("Moments_only").append_child("Density").append_attribute("min") = gHits->texture_limits[1].min.moments_only;
	minMaxNode.child("Moments_only").child("Density").append_attribute("max") = gHits->texture_limits[1].max.moments_only;
	minMaxNode.child("Moments_only").append_child("Imp.rate").append_attribute("min") = gHits->texture_limits[2].min.moments_only;
	minMaxNode.child("Moments_only").child("Imp.rate").append_attribute("max") = gHits->texture_limits[2].max.moments_only;

	return TRUE;
}

void Geometry::LoadXML_geom(pugi::xml_node loadXML, Worker *work, GLProgress *progressDlg){
	//mApp->ClearAllSelections();
	//mApp->ClearAllViews();
	//mApp->ClearFormula();
	Clear();

	xml_node geomNode = loadXML.child("Geometry");

	//Vertices
	sh.nbVertex = geomNode.child("Vertices").select_nodes("Vertex").size();
	vertices3 = (VERTEX3D *)malloc(sh.nbVertex * sizeof(VERTEX3D));
	int idx = 0;
	for (xml_node vertex : geomNode.child("Vertices").children("Vertex")) {
		vertices3[idx].x = vertex.attribute("x").as_double();
		vertices3[idx].y = vertex.attribute("y").as_double();
		vertices3[idx].z = vertex.attribute("z").as_double();
		vertices3[idx].selected = FALSE;
		idx++;
	}

	//Structures
	sh.nbSuper = geomNode.child("Structures").select_nodes("Structure").size();
	idx = 0;
	for (xml_node structure : geomNode.child("Structures").children("Structure")) {
		strName[idx] = _strdup(structure.attribute("name").value());
		// For backward compatibilty with STR
		char tmp[256];
		sprintf(tmp, "%s.txt", strName[idx]);
		strFileName[idx] = _strdup(tmp);
		idx++;
	}

	//Parameters (needs to precede facets)
	xml_node simuParamNode = loadXML.child("MolflowSimuSettings");
	xml_node paramNode = simuParamNode.child("Parameters");
	for (xml_node newParameter : paramNode.children("Parameter")){
		Parameter newPar;
		newPar.name = newParameter.attribute("name").as_string();
		for (xml_node newMoment : newParameter.children("Moment")) {
			newPar.AddValue(std::make_pair(newMoment.attribute("t").as_double(),
				newMoment.attribute("value").as_double()));
		}
		work->parameters.push_back(newPar);
	}

	//Facets
	sh.nbFacet = geomNode.child("Facets").select_nodes("Facet").size();
	facets = (Facet **)malloc(sh.nbFacet * sizeof(Facet *));
	memset(facets, 0, sh.nbFacet * sizeof(Facet *));
	idx = 0;
	for (xml_node facetNode : geomNode.child("Facets").children("Facet")) {
		int nbIndex = facetNode.child("Indices").select_nodes("Indice").size();
		if (nbIndex < 3) {
			char errMsg[128];
			sprintf(errMsg, "Facet %d has only %d vertices. ", idx + 1, nbIndex);
			throw Error(errMsg);
		}

		facets[idx] = new Facet(nbIndex);
		facets[idx]->LoadXML(facetNode, sh.nbVertex);

		//Set param names for interface
		if (facets[idx]->sh.sticking_paramId > -1) facets[idx]->userSticking = work->parameters[facets[idx]->sh.sticking_paramId].name;
		if (facets[idx]->sh.opacity_paramId > -1) facets[idx]->userOpacity = work->parameters[facets[idx]->sh.opacity_paramId].name;
		if (facets[idx]->sh.outgassing_paramId > -1) facets[idx]->userOutgassing = work->parameters[facets[idx]->sh.outgassing_paramId].name;

		idx++;
	}

	xml_node interfNode = loadXML.child("Interface");

	xml_node selNode = interfNode.child("Selections");
	//int nbS = selNode.select_nodes("Selection").size();

	for (xml_node sNode : selNode.children("Selection")) {
		ASELECTION s;
		s.name = _strdup(sNode.attribute("name").as_string());
		s.nbSel = sNode.select_nodes("selItem").size();
		s.selection = (int *)malloc((s.nbSel)*sizeof(int));
		idx = 0;
		for (xml_node iNode : sNode.children("selItem"))
			s.selection[idx++] = iNode.attribute("facet").as_int();
		mApp->AddSelection(s.name, s);
	}


	xml_node viewNode = interfNode.child("Views");
	for (xml_node newView : viewNode.children("View")) {
		AVIEW v;
		v.name = _strdup(newView.attribute("name").as_string());
		v.projMode = newView.attribute("projMode").as_int();
		v.camAngleOx = newView.attribute("camAngleOx").as_double();
		v.camAngleOy = newView.attribute("camAngleOy").as_double();
		v.camDist = newView.attribute("camDist").as_double();
		v.camOffset.x = newView.attribute("camOffset.x").as_double();
		v.camOffset.y = newView.attribute("camOffset.y").as_double();
		v.camOffset.z = newView.attribute("camOffset.z").as_double();
		v.performXY = newView.attribute("performXY").as_int();
		v.vLeft = newView.attribute("vLeft").as_double();
		v.vRight = newView.attribute("vRight").as_double();
		v.vTop = newView.attribute("vTop").as_double();
		v.vBottom = newView.attribute("vBottom").as_double();

		mApp->AddView(v.name, v);
	}

	xml_node formulaNode = interfNode.child("Formulas");
	for (xml_node newFormula : formulaNode.children("Formula")) {
		mApp->AddFormula(newFormula.attribute("name").as_string(),
			newFormula.attribute("expression").as_string());
	}

	xml_node ppNode = interfNode.child("ProfilePlotter");
	if (ppNode) {
		if (!mApp->profilePlotter) mApp->profilePlotter = new ProfilePlotter(&mApp->worker);
		xml_node paramsNode = ppNode.child("Parameters");
		if (paramsNode && paramsNode.attribute("logScale"))
			mApp->profilePlotter->SetLogScaled(paramsNode.attribute("logScale").as_bool());
		xml_node viewsNode = ppNode.child("Views");
		if (viewsNode) {
			std::vector<int> views;
			for (xml_node view : viewsNode.children("View"))
				views.push_back(view.attribute("facetId").as_int());
			mApp->profilePlotter->SetViews(views);
		}
	}

	work->gasMass = simuParamNode.child("Gas").attribute("mass").as_double();
	work->halfLife = simuParamNode.child("Gas").attribute("halfLife").as_double();

	xml_node timeSettingsNode = simuParamNode.child("TimeSettings");

	xml_node userMomentsNode = timeSettingsNode.child("UserMoments");
	for (xml_node newUserEntry : userMomentsNode.children("UserEntry")) {
		char tmpExpr[512];
		strcpy(tmpExpr, newUserEntry.attribute("content").as_string());
		work->userMoments.push_back(tmpExpr);
		work->AddMoment(mApp->worker.ParseMoment(tmpExpr));
	}
	work->timeWindowSize = timeSettingsNode.attribute("timeWindow").as_double();
	work->useMaxwellDistribution = timeSettingsNode.attribute("useMaxwellDistr").as_int();
	work->calcConstantFlow = timeSettingsNode.attribute("calcConstFlow").as_int();

	xml_node motionNode = simuParamNode.child("Motion");
	work->motionType = motionNode.attribute("type").as_int();
	if (work->motionType == 1) { //fixed motion
		xml_node v = motionNode.child("VelocityVector");
		work->motionVector2.x = v.attribute("vx").as_double();
		work->motionVector2.y = v.attribute("vy").as_double();
		work->motionVector2.z = v.attribute("vz").as_double();
	}
	else if (work->motionType == 2) { //rotation
		xml_node v = motionNode.child("AxisBasePoint");
		work->motionVector1.x = v.attribute("x").as_double();
		work->motionVector1.y = v.attribute("y").as_double();
		work->motionVector1.z = v.attribute("z").as_double();
		xml_node v2 = motionNode.child("RotationVector");
		work->motionVector2.x = v2.attribute("x").as_double();
		work->motionVector2.y = v2.attribute("y").as_double();
		work->motionVector2.z = v2.attribute("z").as_double();
	}

	InitializeGeometry();
	//AdjustProfile();
	isLoaded = TRUE;

	// Update mesh
	progressDlg->SetMessage("Building mesh...");
	for (int i = 0; i < sh.nbFacet; i++) {
		double p = (double)i / (double)sh.nbFacet;
		progressDlg->SetProgress(p);
		Facet *f = facets[i];
		if (!f->SetTexture(f->sh.texWidthD, f->sh.texHeightD, f->hasMesh)) {
			char errMsg[512];
			sprintf(errMsg, "Not enough memory to build mesh on Facet %d. ", i + 1);
			throw Error(errMsg);
		}
		BuildFacetList(f);
		double nU = Norme(&(f->sh.U));
		f->tRatio = f->sh.texWidthD / nU;
	}
}

void Geometry::InsertXML(pugi::xml_node loadXML, Worker *work, GLProgress *progressDlg, BOOL newStr){
	//mApp->ClearAllSelections();
	//mApp->ClearAllViews();
	//mApp->ClearFormula();
	//Clear();
	int structId = viewStruct;
	if (structId == -1) structId = 0;
	UnSelectAll();

	xml_node geomNode = loadXML.child("Geometry");

	//Vertices
	
	int nbNewVertex=geomNode.child("Vertices").select_nodes("Vertex").size();
	int nbNewFacets = geomNode.child("Facets").select_nodes("Facet").size();
	
	// reallocate memory
	facets = (Facet **)realloc(facets, (nbNewFacets + sh.nbFacet) * sizeof(Facet **));
	memset(facets + sh.nbFacet, 0, nbNewFacets * sizeof(Facet *));
	VERTEX3D *tmp_vertices3 = (VERTEX3D *)malloc((nbNewVertex + sh.nbVertex) * sizeof(VERTEX3D));
	memmove(tmp_vertices3, vertices3, (sh.nbVertex)*sizeof(VERTEX3D));
	memset(tmp_vertices3 + sh.nbVertex, 0, nbNewVertex * sizeof(VERTEX3D));
	SAFE_FREE(vertices3);
	vertices3 = tmp_vertices3;

	// Read geometry vertices
	int idx = sh.nbVertex;
	for (xml_node vertex : geomNode.child("Vertices").children("Vertex")) {
		vertices3[idx].x = vertex.attribute("x").as_double();
		vertices3[idx].y = vertex.attribute("y").as_double();
		vertices3[idx].z = vertex.attribute("z").as_double();
		vertices3[idx].selected = FALSE;
		idx++;
	}

	//Structures
	sh.nbSuper = geomNode.child("Structures").select_nodes("Structure").size();
	idx = 0;
	for (xml_node structure : geomNode.child("Structures").children("Structure")) {
		strName[idx] = _strdup(structure.attribute("name").value());
		// For backward compatibilty with STR
		char tmp[256];
		sprintf(tmp, "%s.txt", strName[idx]);
		strFileName[idx] = _strdup(tmp);
		idx++;
	}

	//Parameters (needs to precede facets)
	xml_node simuParamNode = loadXML.child("MolflowSimuSettings");
	xml_node paramNode = simuParamNode.child("Parameters");
	for (xml_node newParameter : paramNode.children("Parameter")){
		Parameter newPar;
		newPar.name = newParameter.attribute("name").as_string();
		for (xml_node newMoment : newParameter.children("Moment")) {
			newPar.AddValue(std::make_pair(newMoment.attribute("t").as_double(),
				newMoment.attribute("value").as_double()));
		}
		work->parameters.push_back(newPar);
	}

	//Facets
	idx = sh.nbFacet;
	for (xml_node facetNode : geomNode.child("Facets").children("Facet")) {
		int nbIndex = facetNode.child("Indices").select_nodes("Indice").size();
		if (nbIndex < 3) {
			char errMsg[128];
			sprintf(errMsg, "Facet %d has only %d vertices. ", idx + 1, nbIndex);
			throw Error(errMsg);
		}

		facets[idx] = new Facet(nbIndex);
		facets[idx]->LoadXML(facetNode, sh.nbVertex+nbNewVertex, sh.nbVertex);
		facets[idx]->selected = TRUE;

		//Set param names for interface
		if (facets[idx]->sh.sticking_paramId > -1) facets[idx]->userSticking = work->parameters[facets[idx]->sh.sticking_paramId].name;
		if (facets[idx]->sh.opacity_paramId > -1) facets[idx]->userOpacity = work->parameters[facets[idx]->sh.opacity_paramId].name;
		if (facets[idx]->sh.outgassing_paramId > -1) facets[idx]->userOutgassing = work->parameters[facets[idx]->sh.outgassing_paramId].name;

		idx++;
	}

	xml_node interfNode = loadXML.child("Interface");

	xml_node selNode = interfNode.child("Selections");
	//int nbS = selNode.select_nodes("Selection").size();

	for (xml_node sNode : selNode.children("Selection")) {
		ASELECTION s;
		s.name = _strdup(sNode.attribute("name").as_string());
		s.nbSel = sNode.select_nodes("selItem").size();
		s.selection = (int *)malloc((s.nbSel)*sizeof(int));
		idx = 0;
		for (xml_node iNode : sNode.children("selItem"))
			s.selection[idx++] = iNode.attribute("facet").as_int()+sh.nbFacet; //offset selection numbers
		mApp->AddSelection(s.name, s);
	}


	xml_node viewNode = interfNode.child("Views");
	for (xml_node newView : selNode.children("View")) {
		AVIEW v;
		v.name = _strdup(newView.attribute("name").as_string());
		v.projMode = newView.attribute("projMode").as_int();
		v.camAngleOx = newView.attribute("camAngleOx").as_double();
		v.camAngleOy = newView.attribute("camAngleOy").as_double();
		v.camDist = newView.attribute("camDist").as_double();
		v.camOffset.x = newView.attribute("camOffset.x").as_double();
		v.camOffset.y = newView.attribute("camOffset.y").as_double();
		v.camOffset.z = newView.attribute("camOffset.z").as_double();
		v.performXY = newView.attribute("performXY").as_int();
		v.vLeft = newView.attribute("vLeft").as_double();
		v.vRight = newView.attribute("vRight").as_double();
		v.vTop = newView.attribute("vTop").as_double();
		v.vBottom = newView.attribute("vBottom").as_double();

		mApp->AddView(v.name, v);
	}

	sh.nbVertex += nbNewVertex;
	sh.nbFacet += nbNewFacets; //formulas can refer to newly inserted facets

	xml_node formulaNode = interfNode.child("Formulas");
	for (xml_node newFormula : formulaNode.children("Formula")) {
		char tmpExpr[512];
		strcpy(tmpExpr,newFormula.attribute("expression").as_string());
		mApp->OffsetFormula(tmpExpr, sh.nbFacet);
		mApp->AddFormula(newFormula.attribute("name").as_string(),
			tmpExpr);
	}

	/*work->gasMass = simuParamNode.child("Gas").attribute("mass").as_double();
	work->halfLife = simuParamNode.child("Gas").attribute("halfLife").as_double();*/

	/*
	xml_node timeSettingsNode = simuParamNode.child("TimeSettings");

	xml_node userMomentsNode = timeSettingsNode.child("UserMoments");
	for (xml_node newUserEntry : userMomentsNode.children("UserEntry")) {
		char tmpExpr[512];
		strcpy(tmpExpr, newUserEntry.attribute("content").as_string());
		work->userMoments.push_back(tmpExpr);
		work->AddMoment(mApp->worker.ParseMoment(tmpExpr));
	}
	work->timeWindowSize = timeSettingsNode.attribute("timeWindow").as_double();
	work->useMaxwellDistribution = timeSettingsNode.attribute("useMaxwellDistr").as_int();
	work->calcConstantFlow = timeSettingsNode.attribute("calcConstFlow").as_int();
	*/
	
	if (newStr) AddStruct("Inserted ZIP/XML file");

	InitializeGeometry();
	//AdjustProfile();
	isLoaded = TRUE;

	/*
	// Update mesh
	progressDlg->SetMessage("Building mesh...");
	for (int i = 0; i < sh.nbFacet; i++) {
	double p = (double)i / (double)sh.nbFacet;
	progressDlg->SetProgress(p);
	Facet *f = facets[i];
	if (!f->SetTexture(f->sh.texWidthD, f->sh.texHeightD, f->hasMesh)) {
	char errMsg[512];
	sprintf(errMsg, "Not enough memory to build mesh on Facet %d. ", i + 1);
	throw Error(errMsg);
	}
	BuildFacetList(f);
	double nU = Norme(&(f->sh.U));
	f->tRatio = f->sh.texWidthD / nU;
	}
	*/
}


BOOL Geometry::LoadXML_simustate(pugi::xml_node loadXML, Dataport *dpHit, Worker *work, GLProgress *progressDlg){
	if (!loadXML.child("MolflowResults")) return FALSE; //simu state not saved with file
	AccessDataport(dpHit);
	BYTE* buffer = (BYTE*)dpHit->buff;
	SHGHITS *gHits = (SHGHITS *)buffer;
	xml_node resultNode = loadXML.child("MolflowResults");
	xml_node momentsNode = resultNode.child("Moments");
	double nbMoments = (double)momentsNode.select_nodes("Moment").size();
	size_t m = 0;
	for (xml_node newMoment : momentsNode.children("Moment")) {
		progressDlg->SetProgress((double)m / nbMoments);
		if (m == 0) { //read global results
			xml_node globalNode = newMoment.child("Global");

			xml_node hitsNode = globalNode.child("Hits");
			work->nbHit = hitsNode.attribute("totalHit").as_llong();
			work->nbDesorption = hitsNode.attribute("totalDes").as_llong();
			work->nbAbsorption = hitsNode.attribute("totalAbs").as_llong();
			if (hitsNode.attribute("totalDist_total")) { //if it's in the new format where total/partial are separated
				work->distTraveledTotal_total = hitsNode.attribute("totalDist_total").as_double();
				work->distTraveledTotal_fullHitsOnly = hitsNode.attribute("totalDist_fullHitsOnly").as_double();
			} else
			work->distTraveledTotal_total = work->distTraveledTotal_fullHitsOnly = hitsNode.attribute("totalDist").as_double();
			work->nbLeakTotal = hitsNode.attribute("totalLeak").as_llong();
			//work->maxDesorption=hitsNode.attribute("maxDesorption").as_llong();

			HIT pHits[NBHHIT]; //hits temp storage for loading
			work->nbHHit = 0;
			xml_node hitCacheNode = globalNode.child("Hit_Cache");
			for (xml_node newHit : hitCacheNode.children("Hit")) {
				pHits[work->nbHHit].pos.x = newHit.attribute("posX").as_double();
				pHits[work->nbHHit].pos.y = newHit.attribute("posY").as_double();
				pHits[work->nbHHit].pos.z = newHit.attribute("posZ").as_double();
				pHits[work->nbHHit].type = newHit.attribute("type").as_int();
				work->nbHHit++;
			}
			work->SetHHit(pHits, &work->nbHHit, gHits);

			LEAK pLeak[NBHLEAK]; //leak temp storage for loading
			work->nbLastLeaks = 0;
			xml_node leakCacheNode = globalNode.child("Leak_Cache");
			for (xml_node newLeak : leakCacheNode.children("Leak")) {
				pLeak[work->nbLastLeaks].pos.x = newLeak.attribute("posX").as_double();
				pLeak[work->nbLastLeaks].pos.y = newLeak.attribute("posY").as_double();
				pLeak[work->nbLastLeaks].pos.z = newLeak.attribute("posZ").as_double();
				pLeak[work->nbLastLeaks].dir.x = newLeak.attribute("dirX").as_double();
				pLeak[work->nbLastLeaks].dir.y = newLeak.attribute("dirY").as_double();
				pLeak[work->nbLastLeaks].dir.z = newLeak.attribute("dirZ").as_double();
				work->nbLastLeaks++;
			}
			work->SetLeak(pLeak, &work->nbLastLeaks, gHits);
		} //end global node

		xml_node facetResultsNode = newMoment.child("FacetResults");
		for (xml_node newFacetResult : facetResultsNode.children("Facet")) {
			int facetId = newFacetResult.attribute("id").as_int();
			Facet* f = GetFacet(facetId);
			if (m == 0) { //Now it's a global value, will soon become time-dependent
				xml_node facetHitNode = newFacetResult.child("Hits");
				f->sh.counter.hit.nbHit = facetHitNode.attribute("nbHit").as_llong();
				f->sh.counter.hit.nbDesorbed = facetHitNode.attribute("nbDes").as_llong();
				f->sh.counter.hit.nbAbsorbed = facetHitNode.attribute("nbAbs").as_llong();
				f->sh.counter.hit.sum_v_ort = facetHitNode.attribute("sum_v_ort").as_double();
				f->sh.counter.hit.sum_1_per_ort_velocity = facetHitNode.attribute("sum_1_per_v").as_double();
			}

			//Profiles
			if (f->sh.isProfile){
				xml_node profileNode = newFacetResult.child("Profile");
				APROFILE *profilePtr = (APROFILE *)(buffer + f->sh.hitOffset + sizeof(SHHITS) + m*sizeof(APROFILE)*PROFILE_SIZE);
				size_t id = 0;
				for (xml_node slice : profileNode.children("Slice")){
					profilePtr[id].count = slice.attribute("count").as_llong();
					profilePtr[id].sum_1_per_ort_velocity = slice.attribute("sum_1_per_v").as_double();
					profilePtr[id].sum_v_ort = slice.attribute("sum_v_ort").as_double();
					id++;
				}
			}

			//Textures
			int ix, iy;
			int profSize = (f->sh.isProfile) ? (PROFILE_SIZE*sizeof(APROFILE)*(1 + (int)mApp->worker.moments.size())) : 0;
			/*int h = (f->sh.texHeight);
			int w = (f->sh.texWidth);*/

			if (f->hasMesh){
				xml_node textureNode = newFacetResult.child("Texture");
				int texWidth_file = textureNode.attribute("width").as_int();
				int texHeight_file = textureNode.attribute("height").as_int();
				
				/*if (textureNode.attribute("width").as_int() != f->sh.texWidth ||
					textureNode.attribute("height").as_int() != f->sh.texHeight) {
					std::stringstream msg;
					msg << "Texture size mismatch on facet " << facetId + 1 << ".\nExpected: " << f->sh.texWidth << "x" << f->sh.texHeight << "\n"
						<< "In file: " << textureNode.attribute("width").as_int() << "x" << textureNode.attribute("height").as_int();
					throw Error(msg.str().c_str());
				}*/ //We'll treat texture size mismatch, see below

				AHIT *hits = (AHIT *)((BYTE *)gHits + (f->sh.hitOffset + sizeof(SHHITS) + profSize + m*f->sh.texWidth*f->sh.texHeight*sizeof(AHIT)));
				std::stringstream countText, sum1perText, sumvortText;
				countText << textureNode.child_value("count");
				sum1perText << textureNode.child_value("sum_1_per_v");
				sumvortText << textureNode.child_value("sum_v_ort");

				for (iy = 0; iy<(MIN(f->sh.texHeight, texHeight_file)); iy++) { //MIN: If stored texture is larger, don't read extra cells
					for (ix = 0; ix<(MIN(f->sh.texWidth, texWidth_file)); ix++) { //MIN: If stored texture is larger, don't read extra cells
						countText >> hits[iy*f->sh.texWidth + ix].count;
						sum1perText >> hits[iy*f->sh.texWidth + ix].sum_1_per_ort_velocity;
						sumvortText >> hits[iy*f->sh.texWidth + ix].sum_v_ort_per_area;
					}
					for (int ie = 0; ie < texWidth_file - f->sh.texWidth; ie++) {//Executed if file texture is bigger than expected texture
						//Read extra cells from file without doing anything
						llong dummy_ll;
						double dummy_d;
						countText >> dummy_ll;
						sum1perText >> dummy_d;
						sumvortText >> dummy_d;
					}
				}
				for (int ie = 0; ie < texHeight_file - f->sh.texHeight; ie++) {//Executed if file texture is bigger than expected texture
					//Read extra cells ffrom file without doing anything
					for (int iw = 0; iw < texWidth_file; iw++) {
						llong dummy_ll;
						double dummy_d;
						countText >> dummy_ll;
						sum1perText >> dummy_d;
						sumvortText >> dummy_d;
					}
				}
			} //end texture

			if (f->sh.countDirection && f->dirCache) {
				xml_node dirNode = newFacetResult.child("Directions");
				if (dirNode.attribute("width").as_int() != f->sh.texWidth ||
					dirNode.attribute("height").as_int() != f->sh.texHeight) {
					std::stringstream msg;
					msg << "Direction texture size mismatch on facet " << facetId + 1 << ".\nExpected: " << f->sh.texWidth << "x" << f->sh.texHeight << "\n"
						<< "In file: " << dirNode.attribute("width").as_int() << "x" << dirNode.attribute("height").as_int();
					throw Error(msg.str().c_str());
				}

				VHIT *dirs = (VHIT *)((BYTE *)gHits + f->sh.hitOffset + sizeof(SHHITS) 
					+ profSize + (1 + (int)work->moments.size())*f->sh.texWidth*f->sh.texHeight*sizeof(AHIT) 
					+ m*f->sh.texWidth*f->sh.texHeight*sizeof(VHIT));

				std::stringstream dirText, dirCountText;
				dirText << dirNode.child_value("vel.vectors");
				dirCountText << dirNode.child_value("count");

				for (int iy = 0; iy < f->sh.texHeight; iy++) {
					for (int ix = 0; ix < f->sh.texWidth; ix++) {
						std::string component;
						std::getline(dirText, component, ',');
						dirs[iy*f->sh.texWidth + ix].sumDir.x = std::stod(component);
						std::getline(dirText, component, ',');
						dirs[iy*f->sh.texWidth + ix].sumDir.y = std::stod(component);
						dirText >> dirs[iy*f->sh.texWidth + ix].sumDir.z;
						dirCountText >> dirs[iy*f->sh.texWidth + ix].count;
					}
				}
			} //end directions

		} //end facetResult
		m++;
	} //end moment

	xml_node minMaxNode = resultNode.child("TextureMinMax");
	gHits->texture_limits[0].min.all = minMaxNode.child("With_constant_flow").child("Pressure").attribute("min").as_double();
	gHits->texture_limits[0].max.all = minMaxNode.child("With_constant_flow").child("Pressure").attribute("max").as_double();
	gHits->texture_limits[1].min.all = minMaxNode.child("With_constant_flow").child("Density").attribute("min").as_double();
	gHits->texture_limits[1].max.all = minMaxNode.child("With_constant_flow").child("Density").attribute("max").as_double();
	gHits->texture_limits[2].min.all = minMaxNode.child("With_constant_flow").child("Imp.rate").attribute("min").as_double();
	gHits->texture_limits[2].max.all = minMaxNode.child("With_constant_flow").child("Imp.rate").attribute("max").as_double();

	gHits->texture_limits[0].min.moments_only = minMaxNode.child("Moments_only").child("Pressure").attribute("min").as_double();
	gHits->texture_limits[0].max.moments_only = minMaxNode.child("Moments_only").child("Pressure").attribute("max").as_double();
	gHits->texture_limits[1].min.moments_only = minMaxNode.child("Moments_only").child("Density").attribute("min").as_double();
	gHits->texture_limits[1].max.moments_only = minMaxNode.child("Moments_only").child("Density").attribute("max").as_double();
	gHits->texture_limits[2].min.moments_only = minMaxNode.child("Moments_only").child("Imp.rate").attribute("min").as_double();
	gHits->texture_limits[2].max.moments_only = minMaxNode.child("Moments_only").child("Imp.rate").attribute("max").as_double();

	ReleaseDataport(dpHit);
	return true;
}
