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

/*
//Leak detection
#ifdef _DEBUG
#define _CRTDBG_MAP_ALLOC
#define DEBUG_NEW new(_NORMAL_BLOCK, __FILE__, __LINE__)
#define new DEBUG_NEW
#endif
*/

extern GLApplication *theApp;
extern double totalOutgassing;
extern double gasMass;
extern int nonIsothermal;
//extern HANDLE molflowHandle;
extern int needsReload;

// -----------------------------------------------------------

Geometry::Geometry() {

	facets = NULL;
	vertices3 = NULL;
	polyList = 0;
	selectList = 0;
	selectList2 = 0;
	selectList3 = 0;
	arrowList = 0;
	sphereList = 0;
	texMin=0;
	texMax=1;
	autoNorme = TRUE;
	centerNorme = TRUE;
	normeRatio = 1.0f;
	texAutoScale=TRUE;
	texAutoScaleIncludeConstantFlow=TRUE;
	texLogScale=FALSE;

	sh.nbSuper=0;
	viewStruct = -1;
	strcpy(strPath,"");
	Clear();

	//Debug memory check
	 //_ASSERTE (!_CrtDumpMemoryLeaks());;
	 _ASSERTE(_CrtCheckMemory());
}

// -----------------------------------------------------------

Geometry::~Geometry() {
	Clear();
}

// -----------------------------------------------------------

void Geometry::Clear() {

	// Free memory
	if(facets) {
		for(int i=0;i<sh.nbFacet;i++)
			SAFE_DELETE(facets[i]);
		free(facets);
	}
	if(vertices3) free(vertices3);
	for(int i=0;i<sh.nbSuper;i++) {
		SAFE_FREE(strName[i]);
		SAFE_FREE(strFileName[i]);
	}
	memset(strName,0,MAX_SUPERSTR * sizeof(char *));
	memset(strFileName,0,MAX_SUPERSTR * sizeof(char *));
	DeleteGLLists(TRUE,TRUE);

	// Init default
	facets = NULL;         // Facets array
	vertices3 = NULL;      // Facets vertices in (x,y,z) space
	sh.nbFacet = 0;        // Number of facets
	sh.nbVertex = 0;       // Number of vertex
	isLoaded = FALSE;      // isLoaded flag
	texMinAutoscale=0;             // Current minimum
	texMaxAutoscale=1;             // Current maximum
	sh.nbSuper=0;          // Structure number
	nbSelected=0;          // Number of selected facet

	memset(lineList,0,sizeof(lineList));
	memset(strName,0,sizeof(strName));
	memset(strFileName,0,sizeof(strFileName));

	// Init OpenGL material
	memset (&whiteMaterial, 0, sizeof (GLMATERIAL));
	whiteMaterial.Diffuse.r = 0.9f;
	whiteMaterial.Diffuse.g = 0.9f;
	whiteMaterial.Diffuse.b = 0.9f;
	whiteMaterial.Ambient.r = 0.9f;
	whiteMaterial.Ambient.g = 0.9f;
	whiteMaterial.Ambient.b = 0.9f;

	memset (&fillMaterial, 0, sizeof (GLMATERIAL));
	fillMaterial.Diffuse.r = 0.6f;
	fillMaterial.Diffuse.g = 0.65f;
	fillMaterial.Diffuse.b = 0.65f;
	fillMaterial.Ambient.r = 0.45f;
	fillMaterial.Ambient.g = 0.41f;
	fillMaterial.Ambient.b = 0.41f;

	memset (&arrowMaterial, 0, sizeof (GLMATERIAL));
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

// -----------------------------------------------------------
void Geometry::CalcTotalOutGassing() {
	totalOutgassing=0.0;
	for (int i=0;i<sh.nbFacet;i++) {
		if (facets[i]->sh.desorbType>0)
			totalOutgassing+=facets[i]->sh.flow;
		if (facets[i]->sh.useOutgassingFile) 
			for (int l=0;l<(facets[i]->sh.outgassingMapWidth*facets[i]->sh.outgassingMapHeight);l++)
				totalOutgassing+=facets[i]->outgassingMap[l];
	}
}

// -----------------------------------------------------------

int Geometry::GetNbVertex() {
	return sh.nbVertex;
}

// -----------------------------------------------------------

VERTEX3D Geometry::GetFacetCenter(int facet) {

	return facets[facet]->sh.center;

}

// -----------------------------------------------------------

int Geometry::GetNbStructure() {
	return sh.nbSuper;
}

// -----------------------------------------------------------

char *Geometry::GetStructureName(int idx) {
	return strName[idx];
}

// -----------------------------------------------------------

void Geometry::CalculateFacetParam(int facet) {

	Facet *f = facets[facet];

	// Calculate facet normal
	VERTEX3D p0 = vertices3[facets[facet]->indices[0]];
	VERTEX3D v1;
	VERTEX3D v2;
	BOOL consecutive = TRUE;
	int ind=2;

	// TODO: Handle possible collinear consequtive vectors
	int i0 = facets[facet]->indices[0];
	int i1 = facets[facet]->indices[1];
	while (ind<f->sh.nbIndex && consecutive) {
		int i2 = facets[facet]->indices[ind++];

		Sub(&v1 , vertices3+i1 , vertices3+i0); // v1 = P0P1
		Sub(&v2 , vertices3+i2 , vertices3+i1); // v2 = P1P2
		Cross(&(f->sh.N),&v1,&v2);              // Cross product
		consecutive=(Norme(&(f->sh.N))<1e-11);
	}
	f->collinear=consecutive; //mark for later that this facet was on a line
	Normalize(&(f->sh.N));                  // Normalize

	// Calculate Axis Aligned Bounding Box
	f->sh.bb.min.x =  1e100;
	f->sh.bb.min.y =  1e100;
	f->sh.bb.min.z =  1e100;
	f->sh.bb.max.x = -1e100;
	f->sh.bb.max.y = -1e100;
	f->sh.bb.max.z = -1e100;

	for(int i=0;i<f->sh.nbIndex;i++) {
		VERTEX3D p = vertices3[f->indices[i]];
		if( p.x < f->sh.bb.min.x ) f->sh.bb.min.x = p.x;
		if( p.y < f->sh.bb.min.y ) f->sh.bb.min.y = p.y;
		if( p.z < f->sh.bb.min.z ) f->sh.bb.min.z = p.z;
		if( p.x > f->sh.bb.max.x ) f->sh.bb.max.x = p.x;
		if( p.y > f->sh.bb.max.y ) f->sh.bb.max.y = p.y;
		if( p.z > f->sh.bb.max.z ) f->sh.bb.max.z = p.z;
	}

	// Facet center (AABB center)
	f->sh.center.x = (f->sh.bb.max.x + f->sh.bb.min.x)/2.0;
	f->sh.center.y = (f->sh.bb.max.y + f->sh.bb.min.y)/2.0;
	f->sh.center.z = (f->sh.bb.max.z + f->sh.bb.min.z)/2.0;

	// Plane equation
	double A = f->sh.N.x;
	double B = f->sh.N.y;
	double C = f->sh.N.z;
	double D = -Dot(&(f->sh.N),&p0);

	// Facet planeity
	int nb = f->sh.nbIndex;
	double max = 0.0;
	for(int i=1;i<nb;i++) {
		VERTEX3D p = vertices3[f->indices[i]];
		double d = A * p.x + B * p.y + C * p.z + D;
		if(fabs(d)>fabs(max)) max = d;
	}

	// Plane coef
	f->a = A;
	f->b = B;
	f->c = C;
	f->d = D;
	f->err = max;

	f->sh.maxSpeed = 4.0 * sqrt(2.0*8.31*f->sh.temperature/0.001/gasMass);
}

// -----------------------------------------------------------
void Geometry::CheckCollinear() {
	char tmp[256];
	// Check collinear polygon
	int nbCollinear = 0;
	for(int i=0;i<GetNbFacet();i++) {
		if (GetFacet(i)->collinear) nbCollinear++;
	}
	BOOL ok=FALSE;
	if( nbCollinear ) {
		sprintf(tmp,"%d null polygon(s) found !\nThese polygons have all vertices on a single line, thus they do nothing.\nDelete them?",nbCollinear);
		ok = GLMessageBox::Display(tmp,"Info",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONINFO)==GLDLG_OK;
	}
	if (ok) RemoveCollinear();
}

// -----------------------------------------------------------
void Geometry::CheckNonSimple() {
	char tmp[256];
	// Check non simple polygon
	int *nonSimpleList = (int *)malloc(GetNbFacet()*sizeof(int));
	int nbNonSimple = 0;
	for(int i=0;i<GetNbFacet();i++) {
		if( GetFacet(i)->sh.sign==0.0 )
			nonSimpleList[nbNonSimple++]=i;
	}
	BOOL ok=FALSE;
	if( nbNonSimple ) {
		sprintf(tmp,"%d non simple (or null) polygon(s) found !\nSome tasks may not work properly\n"
			"Should I try to correct them (vertex shifting)?",nbNonSimple);
		ok = GLMessageBox::Display(tmp,"Warning",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONWARNING)==GLDLG_OK;
	}

	if (ok) CorrectNonSimple(nonSimpleList,nbNonSimple);
	SAFE_FREE(nonSimpleList);
}

// -----------------------------------------------------------
void Geometry::CheckIsolatedVertex() {
	int nbI = HasIsolatedVertices();
	if(nbI) {
		char tmp[256];
		sprintf(tmp,"Remove %d isolated vertices ?",nbI);
		if( GLMessageBox::Display(tmp,"Question",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONINFO)==GLDLG_OK ) {
			DeleteIsolatedVertices(FALSE);
		}
	}
}

// -----------------------------------------------------------
void Geometry::CorrectNonSimple(int *nonSimpleList,int nbNonSimple) {
	changedSinceSave=TRUE;
	Facet *f;    
	for(int i=0;i<nbNonSimple;i++) {
		f=GetFacet(nonSimpleList[i]);
		if( f->sh.sign==0.0 ) {
			int j=0;
			while ((j<f->sh.nbIndex) && (f->sh.sign==0.0)) {
				f->ShiftVertex();
				InitializeGeometry(nonSimpleList[i],TRUE);
				//f->DetectOrientation();
				j++;
			}
		}
	}
	//BuildGLList();
}

// -----------------------------------------------------------

void Geometry::CreatePolyFromVertices_Convex() {
	//creates facet from selected vertices

	MolFlow *mApp = (MolFlow *)theApp;
	changedSinceSave=TRUE;
	nbSelectedVertex = 0;

	int *vIdx = (int *)malloc(sh.nbVertex*sizeof(int));
	memset(vIdx,0xFF,sh.nbVertex*sizeof(int));
	for(int i=0;i<sh.nbVertex;i++ ) {
		//VERTEX3D *v = GetVertex(i);
		if( vertices3[i].selected ) {
			vIdx[nbSelectedVertex] = i;
			nbSelectedVertex++;
		}
	}

	if (nbSelectedVertex<3) {
		char errMsg[512];
		sprintf(errMsg,"Select at least 3 vertices.");
		throw Error(errMsg);
		return;
	}//at least three vertices

	VERTEX3D U,V,N;
	U.x = vertices3[vIdx[0]].x - vertices3[vIdx[1]].x;
	U.y = vertices3[vIdx[0]].y - vertices3[vIdx[1]].y;
	U.z = vertices3[vIdx[0]].z - vertices3[vIdx[1]].z;
	double nU = Norme(&U);
	ScalarMult(&U,1.0/nU); // Normalize U

	int i2=2;
	double nV;
	do {
		V.x = vertices3[vIdx[0]].x - vertices3[vIdx[i2]].x;
		V.y = vertices3[vIdx[0]].y - vertices3[vIdx[i2]].y;
		V.z = vertices3[vIdx[0]].z - vertices3[vIdx[i2]].z;
		nV = Norme(&V);
		ScalarMult(&V,1.0/nV); // Normalize V
		i2++;
	} while (Dot(&U,&V)>0.99 && i2<nbSelectedVertex); //if U and V are almost the same, the projection would be inaccurate


	//Now we have the U,V plane, let's define it by computing the normal vector:
	Cross(&N,&V,&U); //We have a normal vector
	double nN = Norme(&N);
	ScalarMult(&N,1.0/nN); // Normalize N

	Cross(&V,&N,&U); //Make V perpendicular to U and N (and still in the U,V plane)
	nV = Norme(&V);
	ScalarMult(&V,1.0/nV); // Normalize V

	VERTEX2D *projected = (VERTEX2D *)malloc(nbSelectedVertex * sizeof(VERTEX2D));
	VERTEX2D *debug = (VERTEX2D *)malloc(nbSelectedVertex * sizeof(VERTEX2D));

	//Get coordinates in the U,V system
	for (int i=0;i<nbSelectedVertex;i++) {
		ProjectVertex(&(vertices3[vIdx[i]]),&(projected[i]),U,V,vertices3[vIdx[0]]);
	}

	//Graham scan here on the projected[] array
	int *returnList = (int *)malloc(nbSelectedVertex*sizeof(int));
	grahamMain(projected,nbSelectedVertex,returnList);
	int ii,loopLength;
	for (ii=0;ii<nbSelectedVertex;ii++) {
		if (returnList[ii]==returnList[0] && ii>0) break;
	}
	loopLength=ii;
	//End graham scan


	//a new facet
	sh.nbFacet++;
	facets = (Facet **)realloc(facets, sh.nbFacet * sizeof(Facet *));
	facets[sh.nbFacet-1] = new Facet(loopLength);
	facets[sh.nbFacet-1]->sh.sticking = 0.0;
	facets[sh.nbFacet-1]->sh.sticking = DES_NONE;
	//set selection
	UnSelectAll();
	facets[sh.nbFacet-1]->selected = TRUE;
	nbSelected=1;
	for (int i=0;i<loopLength;i++) {
		facets[sh.nbFacet-1]->indices[i]=vIdx[returnList[i]];
	}
	SAFE_FREE(vIdx);
	SAFE_FREE(projected);
	SAFE_FREE(debug);
	SAFE_FREE(returnList);

	InitializeGeometry();
	mApp->UpdateFacetParams(TRUE);
	UpdateSelection();
	mApp->facetList->SetSelectedRow(sh.nbFacet-1);
	mApp->facetList->ScrollToVisible(sh.nbFacet-1,1,FALSE);
}

void Geometry::CreatePolyFromVertices_Order() {
	//creates facet from selected vertices

	MolFlow *mApp = (MolFlow *)theApp;
	changedSinceSave=TRUE;
	/*nbSelectedVertex = 0;

	int *vIdx = (int *)malloc(sh.nbVertex*sizeof(int));
	memset(vIdx,0xFF,sh.nbVertex*sizeof(int));
	for(int i=0;i<sh.nbVertex;i++ ) {
		//VERTEX3D *v = GetVertex(i);
		if( vertices3[i].selected ) {
			vIdx[nbSelectedVertex] = i;
			nbSelectedVertex++;
		}
	}
	
	if (selectedVertexList.size()<3) {
		char errMsg[512];
		sprintf(errMsg,"Select at least 3 vertices.");
		throw Error(errMsg);
		return;
	}//at least three vertices

	VERTEX3D U,V,N;
	U.x = vertices3[selectedVertexList[0]].x - vertices3[selectedVertexList[1]].x;
	U.y = vertices3[selectedVertexList[0]].y - vertices3[selectedVertexList[1]].y;
	U.z = vertices3[selectedVertexList[0]].z - vertices3[selectedVertexList[1]].z;
	double nU = Norme(&U);
	ScalarMult(&U,1.0/nU); // Normalize U

	int i2=2;
	double nV;
	do {
		V.x = vertices3[selectedVertexList[0]].x - vertices3[selectedVertexList[i2]].x;
		V.y = vertices3[selectedVertexList[0]].y - vertices3[selectedVertexList[i2]].y;
		V.z = vertices3[selectedVertexList[0]].z - vertices3[selectedVertexList[i2]].z;
		nV = Norme(&V);
		ScalarMult(&V,1.0/nV); // Normalize V
		i2++;
	} while (Dot(&U,&V)>0.99 && i2<nbSelectedVertex); //if U and V are almost the same, the projection would be inaccurate


	//Now we have the U,V plane, let's define it by computing the normal vector:
	Cross(&N,&V,&U); //We have a normal vector
	double nN = Norme(&N);
	ScalarMult(&N,1.0/nN); // Normalize N

	Cross(&V,&N,&U); //Make V perpendicular to U and N (and still in the U,V plane)
	nV = Norme(&V);
	ScalarMult(&V,1.0/nV); // Normalize V

	VERTEX2D *projected = (VERTEX2D *)malloc(nbSelectedVertex * sizeof(VERTEX2D));
	VERTEX2D *debug = (VERTEX2D *)malloc(nbSelectedVertex * sizeof(VERTEX2D));

	//Get coordinates in the U,V system
	for (int i=0;i<nbSelectedVertex;i++) {
		ProjectVertex(&(vertices3[vIdx[i]]),&(projected[i]),U,V,vertices3[vIdx[0]]);
	}

	//Graham scan here on the projected[] array
	int *returnList = (int *)malloc(nbSelectedVertex*sizeof(int));
	grahamMain(projected,nbSelectedVertex,returnList);
	int ii,loopLength;
	for (ii=0;ii<nbSelectedVertex;ii++) {
		if (returnList[ii]==returnList[0] && ii>0) break;
	}
	loopLength=ii;
	//End graham scan
	*/

	//a new facet
	sh.nbFacet++;
	facets = (Facet **)realloc(facets, sh.nbFacet * sizeof(Facet *));
	facets[sh.nbFacet-1] = new Facet((int)selectedVertexList.size());
	facets[sh.nbFacet-1]->sh.sticking = 0.0;
	facets[sh.nbFacet-1]->sh.sticking = DES_NONE;
	//set selection
	UnSelectAll();
	facets[sh.nbFacet-1]->selected = TRUE;
	nbSelected=1;
	for (size_t i=0;i<selectedVertexList.size();i++) {
		facets[sh.nbFacet-1]->indices[i]=selectedVertexList[i];
	}
	/*SAFE_FREE(vIdx);
	SAFE_FREE(projected);
	SAFE_FREE(debug);
	SAFE_FREE(returnList);*/

	InitializeGeometry();
	mApp->UpdateFacetParams(TRUE);
	UpdateSelection();
	mApp->facetList->SetSelectedRow(sh.nbFacet-1);
	mApp->facetList->ScrollToVisible(sh.nbFacet-1,1,FALSE);
}


void Geometry::CreateDifference() {
	//creates facet from selected vertices

	MolFlow *mApp = (MolFlow *)theApp;
	changedSinceSave=TRUE;
	nbSelectedVertex = 0;

	if (nbSelected!=2) {
		char errMsg[512];
		sprintf(errMsg,"Select exactly 2 facets.");
		throw Error(errMsg);
		return;
	}//at least three vertices

	int firstFacet=-1;
	int secondFacet=-1;;
	for (int i=0;i<sh.nbFacet && secondFacet<0;i++)
	{
		if (facets[i]->selected) {
			if (firstFacet<0) firstFacet=i;
			else (secondFacet=i);
		}
	}

	//TO DO:
	//swap if normals not collinear
	//shift vertex to make nice cut

	//a new facet
	sh.nbFacet++;
	facets = (Facet **)realloc(facets, sh.nbFacet * sizeof(Facet *));
	facets[sh.nbFacet-1] = new Facet(facets[firstFacet]->sh.nbIndex+facets[secondFacet]->sh.nbIndex+2);
	facets[sh.nbFacet-1]->sh.sticking = 0.0;
	facets[sh.nbFacet-1]->sh.sticking = DES_NONE;
	//set selection
	UnSelectAll();
	facets[sh.nbFacet-1]->selected = TRUE;
	nbSelected=1;
	//one circle on first facet
	int counter=0;
	for (int i=0;i<facets[firstFacet]->sh.nbIndex;i++)
		facets[sh.nbFacet-1]->indices[counter++]=facets[firstFacet]->indices[i];
	//close circle by adding the first vertex again
	facets[sh.nbFacet-1]->indices[counter++]=facets[firstFacet]->indices[0];
	//reverse circle on second facet
	for (int i=facets[secondFacet]->sh.nbIndex-1;i>=0;i--)
		facets[sh.nbFacet-1]->indices[counter++]=facets[secondFacet]->GetIndex(i+1);
	//close circle by adding the first vertex again
	facets[sh.nbFacet-1]->indices[counter++]=facets[secondFacet]->indices[0];

	InitializeGeometry();
	mApp->UpdateFacetParams(TRUE);
	UpdateSelection();
	mApp->facetList->SetSelectedRow(sh.nbFacet-1);
	mApp->facetList->ScrollToVisible(sh.nbFacet-1,1,FALSE);
}


// -----------------------------------------------------------
void Geometry::SelectCoplanar(int width,int height,double tolerance) {


	nbSelectedVertex = 0;

	int *vIdx = (int *)malloc(sh.nbVertex*sizeof(int));
	memset(vIdx,0xFF,sh.nbVertex*sizeof(int));
	for(int i=0;i<sh.nbVertex;i++ ) {
		//VERTEX3D *v = GetVertex(i);
		if( vertices3[i].selected ) {
			vIdx[nbSelectedVertex] = i;
			nbSelectedVertex++;
		}
	}

	VERTEX3D U,V,N;
	U.x = vertices3[vIdx[0]].x - vertices3[vIdx[1]].x;
	U.y = vertices3[vIdx[0]].y - vertices3[vIdx[1]].y;
	U.z = vertices3[vIdx[0]].z - vertices3[vIdx[1]].z;
	double nU = Norme(&U);
	ScalarMult(&U,1.0/nU); // Normalize U

	V.x = vertices3[vIdx[0]].x - vertices3[vIdx[2]].x;
	V.y = vertices3[vIdx[0]].y - vertices3[vIdx[2]].y;
	V.z = vertices3[vIdx[0]].z - vertices3[vIdx[2]].z;
	double nV = Norme(&V);
	ScalarMult(&V,1.0/nV); // Normalize V

	Cross(&N,&V,&U); //We have a normal vector
	double nN = Norme(&N);
	if (nN<1e-8) {
		GLMessageBox::Display("Sorry, the 3 selected vertices are on a line.","Can't define plane",GLDLG_OK,GLDLG_ICONERROR);
		return;
	}
	ScalarMult(&N,1.0/nN); // Normalize N

	// Plane equation
	double A = N.x;
	double B = N.y;
	double C = N.z;
	VERTEX3D p0=vertices3[vIdx[0]];
	double D = -Dot(&N,&p0);

	//double denominator=sqrt(pow(A,2)+pow(B,2)+pow(C,2));
	double distance;

	int outX,outY;

	for(int i=0;i<sh.nbVertex;i++ ) {
		VERTEX3D *v = GetVertex(i);
		BOOL onScreen = GLToolkit::Get2DScreenCoord((float)v->x,(float)v->y,(float)v->z,&outX,&outY); //To improve
		onScreen=(onScreen && outX>=0 && outY>=0 && outX<=(width) && (outY<=height));
		if (onScreen) {
			distance=abs(A*v->x+B*v->y+C*v->z+D);
			if (distance<tolerance) { //vertex is on the plane

				vertices3[i].selected = TRUE;

			}
			else {
				vertices3[i].selected = FALSE;
			}
		}
		else {
			vertices3[i].selected = FALSE;
		}
	}
	SAFE_FREE(vIdx);
}

// -----------------------------------------------------------

VERTEX3D *Geometry::GetVertex(int idx) {
	return vertices3+idx;
}

// -----------------------------------------------------------

Facet *Geometry::GetFacet(int facet) {
	if (facet>=sh.nbFacet || facet<0) {	
		MolFlow *mApp = (MolFlow *)theApp;
		char errMsg[512];
		sprintf(errMsg,"Geometry::GetFacet()\nA process tried to access facet #%d that doesn't exist.\nAutoSaving and probably crashing...",facet+1);
		GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
		mApp->AutoSave();
		throw Error(errMsg);
	}
	return facets[facet];
}

// -----------------------------------------------------------

int Geometry::GetNbFacet() {
	return sh.nbFacet;
}

// -----------------------------------------------------------

AABB Geometry::GetBB() {

	if( viewStruct<0 ) {

		return bb;

	} else {

		// BB of selected struture
		AABB sbb;

		sbb.min.x =  1e100;
		sbb.min.y =  1e100;
		sbb.min.z =  1e100;
		sbb.max.x = -1e100;
		sbb.max.y = -1e100;
		sbb.max.z = -1e100;

		// Axis Aligned Bounding Box
		for(int i=0;i<sh.nbFacet;i++) {
			Facet *f = facets[i];
			if( f->sh.superIdx==viewStruct ) {
				for(int j=0;j<f->sh.nbIndex;j++) {
					VERTEX3D p = vertices3[f->indices[j]];
					if( p.x < sbb.min.x ) sbb.min.x = p.x;
					if( p.y < sbb.min.y ) sbb.min.y = p.y;
					if( p.z < sbb.min.z ) sbb.min.z = p.z;
					if( p.x > sbb.max.x ) sbb.max.x = p.x;
					if( p.y > sbb.max.y ) sbb.max.y = p.y;
					if( p.z > sbb.max.z ) sbb.max.z = p.z;
				}
			}
		}

		return sbb;
	}

}

// -----------------------------------------------------------

VERTEX3D Geometry::GetCenter() {

	if( viewStruct<0 ) {

		return center;

	} else {

		VERTEX3D r;
		AABB sbb = GetBB();
		r.x = (sbb.max.x + sbb.min.x)/2.0;
		r.y = (sbb.max.y + sbb.min.y)/2.0;
		r.z = (sbb.max.z + sbb.min.z)/2.0;

		return r;

	}
}

// -----------------------------------------------------------

void Geometry::InitializeGeometry(int facet_number,BOOL noVertexShift) {

	MolFlow *mApp = (MolFlow *)theApp;
	// Perform various precalculation here for a faster simulation.

	//GLProgress *initGeoPrg = new GLProgress("Initializing geometry...","Please wait");
	//initGeoPrg->SetProgress(0.0);
	//initGeoPrg->SetVisible(TRUE);
	if (facet_number==-1) { //bounding box for all vertices
		bb.min.x =  1e100;
		bb.min.y =  1e100;
		bb.min.z =  1e100;
		bb.max.x = -1e100;
		bb.max.y = -1e100;
		bb.max.z = -1e100;

		// Axis Aligned Bounding Box
		for(int i=0;i<sh.nbVertex;i++) {
			VERTEX3D p = vertices3[i];
			if(!(vertices3[i].selected==FALSE || vertices3[i].selected==TRUE)) vertices3[i].selected=FALSE; //initialize selection
			if( p.x < bb.min.x ) bb.min.x = p.x;
			if( p.y < bb.min.y ) bb.min.y = p.y;
			if( p.z < bb.min.z ) bb.min.z = p.z;
			if( p.x > bb.max.x ) bb.max.x = p.x;
			if( p.y > bb.max.y ) bb.max.y = p.y;
			if( p.z > bb.max.z ) bb.max.z = p.z;
		}
	} else { //bounding box only for the changed facet
		for(int i=0;i<facets[facet_number]->sh.nbIndex;i++) {
			VERTEX3D p = vertices3[facets[facet_number]->indices[i]];
			//if(!(vertices3[i].selected==FALSE || vertices3[i].selected==TRUE)) vertices3[i].selected=FALSE; //initialize selection
			if( p.x < bb.min.x ) bb.min.x = p.x;
			if( p.y < bb.min.y ) bb.min.y = p.y;
			if( p.z < bb.min.z ) bb.min.z = p.z;
			if( p.x > bb.max.x ) bb.max.x = p.x;
			if( p.y > bb.max.y ) bb.max.y = p.y;
			if( p.z > bb.max.z ) bb.max.z = p.z;
		}
	}

	center.x = (bb.max.x + bb.min.x)/2.0;
	center.y = (bb.max.y + bb.min.y)/2.0;
	center.z = (bb.max.z + bb.min.z)/2.0;

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
	for(int i=0;i<sh.nbFacet;i++) {
		//initGeoPrg->SetProgress((double)i/(double)sh.nbFacet);
		if((facet_number==-1)||(i==facet_number)) { //permits to initialize only one facet
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
			VERTEX3D U,V;

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
			ScalarMult(&U,1.0/nU); // Normalize U

			// Construct a normal vector V
			Cross(&V,&(f->sh.N),&U); // |U|=1 and |N|=1 => |V|=1

			// u,v vertices (we start with p0 at 0,0)
			f->vertices2[0].u = 0.0;
			f->vertices2[0].v = 0.0;
			VERTEX2D min; min.u=0.0; min.v=0.0;
			VERTEX2D max; max.u=0.0; max.v=0.0;

			for(int j=1;j<f->sh.nbIndex;j++) {

				VERTEX3D p = vertices3[f->indices[j]];
				VERTEX3D v;
				Sub(&v,&p,&p0); // v = p0p
				f->vertices2[j].u = Dot(&U,&v);  // Project p on U along the V direction
				f->vertices2[j].v = Dot(&V,&v);  // Project p on V along the U direction

				// Bounds
				if( f->vertices2[j].u > max.u ) max.u = f->vertices2[j].u;
				if( f->vertices2[j].v > max.v ) max.v = f->vertices2[j].v;
				if( f->vertices2[j].u < min.u ) min.u = f->vertices2[j].u;
				if( f->vertices2[j].v < min.v ) min.v = f->vertices2[j].v;

			}

			// Calculate facet area (Meister/Gauss formula)
			double A=0.0;
			for(int j=0;j<f->sh.nbIndex;j++) {
				int j1 = IDX(j+1,f->sh.nbIndex);
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
			ScalarMult(&U,uD);
			f->sh.U = U;

			f->sh.nV = V;
			ScalarMult(&V,vD);
			f->sh.V = V;

			Cross(&(f->sh.Nuv),&U,&V);

			// Rescale u,v coordinates
			for(int j=0;j<f->sh.nbIndex;j++) {

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


	if (facet_number==-1) {
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

// -----------------------------------------------------------
void Geometry::RebuildLists() {
	BuildGLList();
}

// -----------------------------------------------------------

DWORD Geometry::GetGeometrySize(std::vector<double> *moments) {

	// Compute number of bytes allocated
	DWORD memoryUsage = 0;
	memoryUsage += sizeof(SHGEOM);
	memoryUsage += sh.nbVertex * sizeof(VERTEX3D);
	for(int i=0;i<sh.nbFacet;i++) 
		memoryUsage += facets[i]->GetGeometrySize();
	memoryUsage+=sizeof(double)*(int)(*moments).size(); //time values
	memoryUsage+=3*sizeof(double)+2*sizeof(BOOL)+sizeof(double); //gas pulse parameters
	return memoryUsage;
}

// -----------------------------------------------------------

void Geometry::CopyGeometryBuffer(BYTE *buffer,std::vector<double> *moments) {
	MolFlow *mApp = (MolFlow *)theApp;
	// Build shared buffer for geometry (see Shared.h)
	int fOffset = sizeof(SHGHITS);
	SHGEOM *shGeom = (SHGEOM *)buffer;
	this->sh.nbMoments=(int)(*moments).size();
	this->sh.gasMass=gasMass;
	this->sh.nonIsothermal=nonIsothermal;
	//this->sh.molflowHandle=molflowHandle;
	//this->sh.totalOutgassing=totalOutgassing;
	memcpy(shGeom,&(this->sh),sizeof(SHGEOM));
	buffer += sizeof(SHGEOM);
	memcpy(buffer,vertices3,sizeof(VERTEX3D)*sh.nbVertex);
	buffer += sizeof(VERTEX3D)*sh.nbVertex;
	for(int i=0;i<sh.nbFacet;i++) {
		Facet *f = facets[i];
		f->sh.hitOffset = fOffset;
		fOffset += f->GetHitsSize((int)mApp->worker.moments.size());
		memcpy(buffer,&(f->sh),sizeof(SHFACET));
		buffer += sizeof(SHFACET);
		memcpy(buffer,f->indices,sizeof(int)*f->sh.nbIndex);
		buffer += sizeof(int)*f->sh.nbIndex;
		memcpy(buffer,f->vertices2,sizeof(VERTEX2D)*f->sh.nbIndex);
		buffer += sizeof(VERTEX2D)*f->sh.nbIndex;
		if (f->sh.useOutgassingFile) {
			memcpy(buffer,f->outgassingMap,sizeof(double)*f->sh.outgassingMapWidth*f->sh.outgassingMapHeight);
			buffer += sizeof(double)*f->sh.outgassingMapWidth*f->sh.outgassingMapHeight;
		}
	}

	// Add surface elements area (reciprocal)
	for(int k=0;k<sh.nbFacet;k++) {
		Facet *f = facets[k];
		DWORD add = 0;
		if(f->sh.isTextured) {

			if( f->mesh ) {

				for(int j=0;j<f->sh.texHeight;j++) {
					for(int i=0;i<f->sh.texWidth;i++) {
						float area = f->mesh[add].area;
						if( area>0.0f ) {
							// Use the sign bit to store isFull flag
							if( f->mesh[add].full )
								*((AHIT *)buffer) = -1.0f / area;
							else
								*((AHIT *)buffer) = 1.0f / area;
						} else {
							*((AHIT *)buffer) = 0.0f;
						}
						add++;
						buffer+=sizeof(AHIT);
					}
				}

			} else {

				double rw = Norme(&(f->sh.U)) / (double)(f->sh.texWidthD);
				double rh = Norme(&(f->sh.V)) / (double)(f->sh.texHeightD);
				float area = (float)(rw*rh);

				for(int j=0;j<f->sh.texHeight;j++) {
					for(int i=0;i<f->sh.texWidth;i++) {
						if( area>0.0f ) {
							*((AHIT *)buffer) = 1.0f / area;
						} else {
							*((AHIT *)buffer) = 0.0f;
						}
						buffer+=sizeof(AHIT);
					}
				}

			}
		}
	}

	//Time moments
	for (int i=0;i<(int)(*moments).size();i++) {
		*((double *)buffer)=(*moments)[i];
		buffer+=sizeof(double);
	}
	//Gas pulse properties
	*((double *)buffer)=mApp->worker.desorptionStartTime; buffer+=sizeof(double);
	*((double *)buffer)=mApp->worker.desorptionStopTime;   buffer+=sizeof(double);
	*((double *)buffer)=mApp->worker.timeWindowSize;      buffer+=sizeof(double);
	*((BOOL *)buffer)=mApp->worker.useMaxwellDistribution;   buffer+=sizeof(BOOL);
	*((BOOL *)buffer)=mApp->worker.calcConstantFlow;   buffer+=sizeof(BOOL);
	*((double *)buffer)=mApp->worker.valveOpenMoment;      buffer+=sizeof(double);
}

// -----------------------------------------------------------

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

// -----------------------------------------------------------

DWORD Geometry::GetHitsSize(std::vector<double> *moments) {
	MolFlow *mApp = (MolFlow *)theApp;
	// Compute number of bytes allocated
	DWORD memoryUsage = 0;
	memoryUsage += sizeof(SHGHITS);
	for(int i=0;i<sh.nbFacet;i++) {
		memoryUsage += facets[i]->GetHitsSize((int)mApp->worker.moments.size());
	}

	return memoryUsage;
}

// -----------------------------------------------------------

DWORD Geometry::GetMaxElemNumber() {

	int nbElem = 0;
	for(int i=0;i<sh.nbFacet;i++) {
		Facet *f = facets[i];
		if( f->mesh ) nbElem += f->sh.texWidth*f->sh.texHeight;
		else          return 0;
	}
	return nbElem;

}

// -----------------------------------------------------------

void Geometry::CopyElemBuffer(BYTE *buffer) {

	int idx = 0;
	for(int i=0;i<sh.nbFacet;i++) {
		Facet *f = facets[i];
		int sz = f->sh.texWidth * f->sh.texHeight * sizeof(SHELEM);
		memcpy(buffer+idx,f->mesh,sz);
		idx+=sz;
	}

}

// -----------------------------------------------------------

void Geometry::BuildShapeList() {

	// Shapes used for direction field rendering

	// 3D arrow (direction field)
	int nbDiv = 10;
	double alpha = 2.0*PI/(double)nbDiv;

	arrowList = glGenLists(1);
	glNewList(arrowList,GL_COMPILE);

	glEnable(GL_CULL_FACE);
	glCullFace(GL_FRONT);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glEnable(GL_LIGHTING);
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_BLEND);

	glBegin(GL_TRIANGLES);

	// Arrow
	for(int i=0;i<nbDiv;i++) {

		double y1 = sin(alpha*(double)i);
		double z1 = cos(alpha*(double)i);
		double y2 = sin(alpha*(double)((i+1)%nbDiv));
		double z2 = cos(alpha*(double)((i+1)%nbDiv));

		glNormal3d( 0.0 , y1 , z1 );
		glVertex3d(-0.5 , 0.5*y1 , 0.5*z1 );
		glNormal3d( 1.0 , 0.0 , 0.0 );
		glVertex3d( 0.5 , 0.0 , 0.0 );
		glNormal3d( 0.0 , y2 , z2 );
		glVertex3d(-0.5 , 0.5*y2 , 0.5*z2 );

	}

	// Cap facets
	for(int i=0;i<nbDiv;i++) {

		double y1 = sin(alpha*(double)i);
		double z1 = cos(alpha*(double)i);
		double y2 = sin(alpha*(double)((i+1)%nbDiv));
		double z2 = cos(alpha*(double)((i+1)%nbDiv));

		glNormal3d( -1.0 , 0.0 , 0.0 );
		glVertex3d( -0.5 , 0.5*y1 , 0.5*z1 );
		glNormal3d( -1.0 , 0.0 , 0.0 );
		glVertex3d( -0.5 , 0.5*y2 , 0.5*z2 );
		glNormal3d( -1.0 , 0.0 , 0.0 );
		glVertex3d( -0.5 , 0.0 , 0.0 );

	}

	glEnd();
	glEndList();

	// Shpere list (isotropic case)
	int nbPhi   = 16;
	int nbTetha = 7;
	double dphi   = 2.0*PI / (double)(nbPhi);
	double dtetha = PI / (double)(nbTetha+1);

	sphereList = glGenLists(1);
	glNewList(sphereList,GL_COMPILE);

	glEnable(GL_CULL_FACE);
	glCullFace(GL_FRONT);
	glDisable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glEnable(GL_LIGHTING);
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_BLEND);

	glBegin(GL_TRIANGLES);

	for(int i=0;i<=nbTetha;i++) {
		for(int j=0;j<nbPhi;j++) {

			VERTEX3D v1,v2,v3,v4;

			v1.x = sin(dtetha*(double)i)*cos(dphi*(double)j);
			v1.y = sin(dtetha*(double)i)*sin(dphi*(double)j);
			v1.z = cos(dtetha*(double)i);

			v2.x = sin(dtetha*(double)(i+1))*cos(dphi*(double)j);
			v2.y = sin(dtetha*(double)(i+1))*sin(dphi*(double)j);
			v2.z = cos(dtetha*(double)(i+1));

			v3.x = sin(dtetha*(double)(i+1))*cos(dphi*(double)(j+1));
			v3.y = sin(dtetha*(double)(i+1))*sin(dphi*(double)(j+1));
			v3.z = cos(dtetha*(double)(i+1));

			v4.x = sin(dtetha*(double)i)*cos(dphi*(double)(j+1));
			v4.y = sin(dtetha*(double)i)*sin(dphi*(double)(j+1));
			v4.z = cos(dtetha*(double)i);

			if( i<nbTetha ) {
				glNormal3d(v1.x , v1.y , v1.z);
				glVertex3d(v1.x , v1.y , v1.z);
				glNormal3d(v2.x , v2.y , v2.z);
				glVertex3d(v2.x , v2.y , v2.z);
				glNormal3d(v3.x , v3.y , v3.z);
				glVertex3d(v3.x , v3.y , v3.z);
			}

			if( i>0 ) {
				glNormal3d(v1.x , v1.y , v1.z);
				glVertex3d(v1.x , v1.y , v1.z);
				glNormal3d(v3.x , v3.y , v3.z);
				glVertex3d(v3.x , v3.y , v3.z);
				glNormal3d(v4.x , v4.y , v4.z);
				glVertex3d(v4.x , v4.y , v4.z);
			}

		}
	}

	glEnd();
	glEndList();


}

// -----------------------------------------------------------

void Geometry::BuildSelectList() {

	nbSelected = 0;

	selectList = glGenLists(1);
	glNewList(selectList,GL_COMPILE);
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
	glNewList(selectList2,GL_COMPILE);
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
	glNewList(selectList3,GL_COMPILE);

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
		if( f->selected ) {
			//DrawFacet(f,FALSE,TRUE,TRUE);
			DrawFacet(f,1,1,1);
			nbSelected++;
		}
	}  
	glLineWidth(1.0f);
	if (antiAliasing) {
		glDisable(GL_BLEND);
		glDisable(GL_LINE_SMOOTH);
	}
	glEndList();


}



// -----------------------------------------------------------

void Geometry::UpdateSelection() {

	DeleteGLLists();
	BuildSelectList();

}

// -----------------------------------------------------------

void Geometry::BuildGLList() {

	// Compile geometry for OpenGL
	for(int j=0;j<sh.nbSuper;j++) {
		lineList[j] = glGenLists(1);
		glNewList(lineList[j],GL_COMPILE);
		for(int i=0;i<sh.nbFacet;i++) {
			if(facets[i]->sh.superIdx==j)
				DrawFacet(facets[i]);
		}
		glEndList();
	}

	polyList = glGenLists(1);
	glNewList(polyList,GL_COMPILE);
	DrawPolys();
	glEndList();

	BuildSelectList();

}


// -----------------------------------------------------------
// Collapse stuff
// -----------------------------------------------------------

int Geometry::AddRefVertex(VERTEX3D *p,VERTEX3D *refs,int *nbRef) {

	BOOL found = FALSE;
	int i=0;
	//VERTEX3D n;
	double v2 = vThreshold*vThreshold;

	while(i<*nbRef && !found) {
		//Sub(&n,p,refs + i);
		double dx=abs(p->x-(refs+i)->x);
		if (dx<vThreshold) {
			double dy=abs(p->y-(refs+i)->y);
			if (dy<vThreshold) {
				double dz=abs(p->z-(refs+i)->z);
				if (dz<vThreshold) {
					found= (dx*dx+dy*dy+dz*dz<v2);
				}
			}
		}
		if(!found) i++;
	}

	if(!found) {
		// Add a new reference vertex
		refs[*nbRef] = *p;
		*nbRef = *nbRef + 1;
	}

	return i;

}


// -----------------------------------------------------------

void Geometry::CollapseVertex(GLProgress *prg,double totalWork) {
	changedSinceSave=TRUE;
	if(!isLoaded) return;
	// Collapse neighbor vertices
	VERTEX3D *refs = (VERTEX3D *)malloc(sh.nbVertex*sizeof(VERTEX3D));
	int      *idx  = (int *)malloc(sh.nbVertex*sizeof(int));
	int       nbRef=0;

	// Collapse
	prg->SetMessage("Collapsing vertices...");
	for(int i=0;i<sh.nbVertex;i++) {
		prg->SetProgress(((double)i/(double)sh.nbVertex) /totalWork);
		idx[i] = AddRefVertex(vertices3 + i,refs,&nbRef);
	}

	// Create the new vertex array
	SAFE_FREE(vertices3);
	vertices3 = (VERTEX3D *)malloc(nbRef * sizeof(VERTEX3D));
	//UnselectAllVertex();

	memcpy(vertices3,refs,nbRef * sizeof(VERTEX3D));
	sh.nbVertex = nbRef;

	// Update facets indices
	for(int i=0;i<sh.nbFacet;i++) {
		Facet *f = facets[i];
		prg->SetProgress(((double)i/(double)sh.nbFacet) * 0.05 + 0.45);
		for(int j=0;j<f->sh.nbIndex;j++)
			f->indices[j] = idx[f->indices[j]];
	}

	free(idx);
	free(refs);

}


// -----------------------------------------------------------

BOOL Geometry::GetCommonEdges(Facet *f1,Facet *f2,int *c1,int *c2,int *chainLength) {

	// Detect common edge between facet
	int p11,p12,p21,p22,lgth,si,sj;
	int maxLength = 0;

	for(int i=0;i<f1->sh.nbIndex;i++) {

		p11 = f1->GetIndex(i);
		p12 = f1->GetIndex(i+1);

		for(int j=0;j<f2->sh.nbIndex;j++) {

			p21 = f2->GetIndex(j);
			p22 = f2->GetIndex(j+1);

			if( p11==p22 && p12==p21 ) {

				// Common edge found
				si = i;
				sj = j;
				lgth = 1;

				// Loop until the end of the common edge chain
				i+=2;
				j-=1;
				p12 = f1->GetIndex(i);
				p21 = f2->GetIndex(j);
				BOOL ok = (p12==p21);
				while( lgth<f1->sh.nbIndex && lgth<f2->sh.nbIndex && ok ) {
					p12 = f1->GetIndex(i);
					p21 = f2->GetIndex(j);
					ok = (p12==p21);
					if( ok ) {
						i++;j--;
						lgth++;
					}
				}

				if( lgth>maxLength ) {
					*c1 = si;
					*c2 = sj;
					maxLength = lgth;
				}

			}

		}

	}

	if( maxLength>0 ) {
		*chainLength = maxLength;
		return TRUE;
	}

	return FALSE;

}

// -----------------------------------------------------------

Facet *Geometry::MergeFacet(Facet *f1,Facet *f2) {
	changedSinceSave = TRUE;
	// Merge 2 facets into 1 when possible and create a new facet
	// otherwise return NULL.
	int  c1;
	int  c2;
	int  l;
	BOOL end = FALSE;
	Facet *nF = NULL;

	if( GetCommonEdges(f1,f2,&c1,&c2,&l) ) {
		int commonNo=f1->sh.nbIndex+f2->sh.nbIndex-2*l;
		if (commonNo==0) { //two identical facets, so return a copy of f1
			nF = new Facet(f1->sh.nbIndex);
			nF->Copy(f1);
			for(int i=0;i<f1->sh.nbIndex;i++)
				nF->indices[i] = f1->GetIndex(i);
			return nF;
		}
		int nbI=0;
		nF = new Facet(commonNo);
		// Copy params from f1
		//nF->Copy(f1);
		nF->Copy(f1);

		if( l==f1->sh.nbIndex ) {

			// f1 absorbed, copy indices from f2
			for(int i=0;i<f2->sh.nbIndex - l;i++)
				nF->indices[nbI++] = f2->GetIndex(c2+2+i);

		} else if ( l==f2->sh.nbIndex ) {

			// f2 absorbed, copy indices from f1
			for(int i=0;i<f1->sh.nbIndex - l;i++)
				nF->indices[nbI++] = f1->GetIndex(c1+l+i);

		} else {    

			// Copy indices from f1
			for(int i=0;i<f1->sh.nbIndex - (l-1);i++)
				nF->indices[nbI++] = f1->GetIndex(c1+l+i);
			// Copy indices from f2
			for(int i=0;i<f2->sh.nbIndex - (l+1);i++)
				nF->indices[nbI++] = f2->GetIndex(c2+2+i);

		}

	}

	return nF;

}

// -----------------------------------------------------------
void Geometry::Collapse(double vT,double fT,double lT,BOOL doSelectedOnly,GLProgress *prg) {	
	changedSinceSave = TRUE;
	Facet *fi,*fj;
	Facet *merged;
	MolFlow *mApp = (MolFlow *)theApp;

	vThreshold = vT;
	double totalWork=(1.0 + (double)(fT>0.0) + (double)(lT>0.0)); //for progress indicator
	// Collapse vertex
	if (vT>0.0) {
		CollapseVertex(prg,totalWork);
		RemoveCollinear();
		RemoveNullFacet();
		InitializeGeometry();
	}


	if( fT>0.0 ) {

		// Collapse facets
		int i = 0;
		prg->SetMessage("Collapsing facets...");
		while(i<sh.nbFacet) {
			prg->SetProgress((1.0+((double)i/(double)sh.nbFacet)) / totalWork);
			fi = facets[i];
			// Search a coplanar facet
			int j = i+1;
			while((!doSelectedOnly || fi->selected) && j<sh.nbFacet) {
				fj = facets[j];
				merged = NULL;
				if( (!doSelectedOnly || fj->selected) && fi->IsCoplanar(fj,fT) ) {
					// Collapse
					merged = MergeFacet(fi,fj);
					if( merged ) {
						// Replace the old 2 facets by the new one
						SAFE_DELETE(fi);
						SAFE_DELETE(fj);
						for(int k=j;k<sh.nbFacet-1;k++)
							facets[k] = facets[k+1];
						sh.nbFacet--;
						facets[i] = merged;
						//InitializeGeometry(i);
						//SetFacetTexture(i,facets[i]->tRatio,facets[i]->hasMesh);  //rebuild mesh
						fi = facets[i];
						mApp->RenumberSelections(j);
						mApp->RenumberFormulas(j);
						j=i+1;
					}
				}
				if(!merged) j++;
			}
			i++;
		}
	}
	//Collapse collinear sides. Takes some time, so only if threshold>0
	prg->SetMessage("Collapsing collinear sides...");
	if(lT>0.0) {
		for(int i=0;i<sh.nbFacet;i++) {
			prg->SetProgress((1.0+(double)(fT>0.0)+((double)i/(double)sh.nbFacet)) / totalWork);
			if (!doSelectedOnly || facets[i]->selected)
				MergecollinearSides(facets[i],lT);
		}
	}
	prg->SetMessage("Rebuilding geometry...");
	for(int i=0;i<sh.nbFacet;i++) {

		Facet *f = facets[i];

		// Revert orientation if normal has been swapped
		// This happens when the second vertex is no longer convex
		VERTEX3D n,v1,v2;
		double   d;
		int i0 = facets[i]->indices[0];
		int i1 = facets[i]->indices[1];
		int i2 = facets[i]->indices[2];

		Sub(&v1 , vertices3+i1 , vertices3+i0); // v1 = P0P1
		Sub(&v2 , vertices3+i2 , vertices3+i1); // v2 = P1P2
		Cross(&n,&v1,&v2);                      // Cross product
		d = Dot(&n,&(f->sh.N));
		if( d<0.0 ) f->SwapNormal();

	}



	// Delete old resources
	for(int i=0;i<sh.nbSuper;i++)
		DeleteGLLists(TRUE,TRUE);

	// Reinitialise geom
	InitializeGeometry();

}

void Geometry::MergecollinearSides(Facet *f,double lT) {
	changedSinceSave = TRUE;
	BOOL collinear;
	double linTreshold=cos(lT*PI/180);
	// Merge collinear sides
	for(int k=0;(k<f->sh.nbIndex&&f->sh.nbIndex>3);k++){
		k=k;
		do {
			collinear=FALSE;
			int p0=f->indices[k];
			int p1=f->indices[(k+1)%f->sh.nbIndex];
			int p2=f->indices[(k+2)%f->sh.nbIndex]; //to compare last side with first too
			VERTEX3D p0p1;
			VERTEX3D p0p2;
			Sub(&p0p1,&vertices3[p1],&vertices3[p0]);
			Sub(&p0p2,&vertices3[p2],&vertices3[p1]);
			Normalize(&p0p1);
			Normalize(&p0p2);
			collinear=(Dot(&p0p1,&p0p2)>=linTreshold);
			if (collinear&&f->sh.nbIndex>3) { //collinear
				for (int l=(k+1)%f->sh.nbIndex;l<f->sh.nbIndex-1;l++){
					f->indices[l]=f->indices[l+1];
				}
				f->sh.nbIndex--;
			}
		} while (collinear&&f->sh.nbIndex>3);
	}
}

// -----------------------------------------------------------

void Geometry::Rebuild() {

	// Rebuild internal structure on geometry change

	// Remove texture (improvement TODO)
	for(int i=0;i<sh.nbFacet;i++)
		if (facets[i]->sh.isTextured)
			facets[i]->SetTexture(0.0,0.0,FALSE);

	// Delete old resources
	DeleteGLLists(TRUE,TRUE);

	// Reinitialise geom
	InitializeGeometry();

}

// -----------------------------------------------------------

void Geometry::MoveVertexTo(int idx,double x,double y,double z) {
	vertices3[idx].x = x;
	vertices3[idx].y = y;
	vertices3[idx].z = z;
}

// -----------------------------------------------------------

void Geometry::SwapNormal() {

	if(!IsLoaded()) {
		GLMessageBox::Display("No geometry loaded.","No geometry",GLDLG_OK,GLDLG_ICONERROR);
		return;
	}
	if(GetNbSelected()<=0) return;
	changedSinceSave=TRUE;
	for(int i=0;i<sh.nbFacet;i++) {
		Facet *f = facets[i];
		if( f->selected ) f->SwapNormal();
	}

	DeleteGLLists(TRUE,TRUE);
	BuildGLList();

}

// -----------------------------------------------------------

void Geometry::ShiftVertex() {

	if(!IsLoaded()) {
		GLMessageBox::Display("No geometry loaded.","No geometry",GLDLG_OK,GLDLG_ICONERROR);
		return;
	}
	if(GetNbSelected()<=0) return;
	changedSinceSave = TRUE;
	DeleteGLLists(TRUE,TRUE);
	for(int i=0;i<sh.nbFacet;i++) {
		Facet *f = facets[i];
		if( f->selected ) {
			f->ShiftVertex();
			InitializeGeometry(i,TRUE);// Reinitialise geom
		}
	}
	// Delete old resource
	BuildGLList();
}

// -----------------------------------------------------------

int Geometry::InvalidateDeviceObjects() {

	DeleteGLLists(TRUE,TRUE);
	DELETE_LIST(arrowList);
	DELETE_LIST(sphereList);
	for(int i=0;i<sh.nbFacet;i++)
		facets[i]->InvalidateDeviceObjects();

	return GL_OK;

}

// -----------------------------------------------------------

int Geometry::RestoreDeviceObjects() {

	if(!IsLoaded()) return GL_OK;

	for(int i=0;i<sh.nbFacet;i++) {
		Facet *f = facets[i];
		f->RestoreDeviceObjects();
		BuildFacetList(f);
	}

	BuildGLList();

	return GL_OK;

}

// -----------------------------------------------------------

void Geometry::BuildFacetList(Facet *f) {

	// Rebuild OpenGL geomtetry with texture
	if( f->sh.isTextured && f->textureVisible ) {

		// Facet geometry
		glNewList(f->glList,GL_COMPILE);
		if( f->sh.nbIndex == 3 ) {
			glBegin(GL_TRIANGLES);
			FillFacet(f,TRUE);
			glEnd();
		} else if( f->sh.nbIndex == 4 ) {
			glBegin(GL_QUADS);
			FillFacet(f,TRUE);
			glEnd();
		} else {
			glBegin(GL_TRIANGLES);
			Triangulate(f,TRUE);
			glEnd();
		}
		glEndList();

	}

}

// -----------------------------------------------------------

void Geometry::SetFacetTexture(int facet,double ratio,BOOL mesh) {

	Facet *f = facets[facet];
	double nU = Norme(&(f->sh.U));
	double nV = Norme(&(f->sh.V));
	if (!f->SetTexture(nU*ratio,nV*ratio,mesh)) {
		char errMsg[512];
		sprintf(errMsg,"Not enough memory to build mesh on Facet %d. ",facet+1);
		throw Error(errMsg);
	}
	f->tRatio = ratio;
	BuildFacetList(f);

}

// -----------------------------------------------------------
// Testing purpose function, construct a PIPE
// -----------------------------------------------------------
void  Geometry::BuildPipe(double L,double R,double s,int step) {
	needsReload=TRUE;
	Clear();

	MolFlow *mApp = (MolFlow *)theApp;
	mApp->ClearAllSelections();
	mApp->ClearAllViews();
	sprintf(sh.name,"PIPE%g",L/R);

	int nbDecade = 0;
	int nbTF = 9*nbDecade;
	int nbTV = 4*nbTF;

	sh.nbVertex  = 2*step + nbTV;
	if(!(vertices3 = (VERTEX3D *)malloc(sh.nbVertex * sizeof(VERTEX3D))))
		throw Error("Couldn't allocate memory for vertices");

	sh.nbFacet   = step + 2 + nbTF;
	sh.nbSuper = 1;
	if (!(facets    = (Facet **)malloc(sh.nbFacet * sizeof(Facet *))))
		throw Error("Couldn't allocate memory for facets");
	memset(facets,0,sh.nbFacet * sizeof(Facet *));
	// Vertices
	for(int i=0;i<step;i++) {
		double angle = (double)i/(double)step * 2*PI;
		vertices3[2*i+nbTV].x = R*cos(angle);
		vertices3[2*i+nbTV].y = R*sin(angle);
		vertices3[2*i+nbTV].z = 0.0;
		vertices3[2*i+1+nbTV].x = R*cos(angle);
		vertices3[2*i+1+nbTV].y = R*sin(angle);
		vertices3[2*i+1+nbTV].z = L;
	}

	try {
		// Cap facet
	facets[0+nbTF] = new Facet(step);
	facets[0+nbTF]->sh.sticking = 1.0;
	facets[0+nbTF]->sh.desorbType = DES_COSINE;
	facets[0+nbTF]->sh.flow = 1.0;
	for(int i=0;i<step;i++)
		facets[0+nbTF]->indices[i] = 2*i+nbTV;

	facets[1+nbTF] = new Facet(step);
	facets[1+nbTF]->sh.sticking = 1.0;
	facets[1+nbTF]->sh.desorbType = DES_NONE;
	for(int i=0;i<step;i++)
		facets[1+nbTF]->indices[step-i-1] = 2*i+1+nbTV;

	// Wall facet
	for(int i=0;i<step;i++) {
		facets[i+2+nbTF] = new Facet(4);
		facets[i+2+nbTF]->sh.reflectType = REF_DIFFUSE;
		facets[i+2+nbTF]->sh.sticking = s;
		facets[i+2+nbTF]->indices[0]  = 2*i+nbTV;
		facets[i+2+nbTF]->indices[1]  = 2*i+1+nbTV;
		if( i<step-1 ) {
			facets[i+2+nbTF]->indices[2] = 2*(i+1)+1+nbTV;
			facets[i+2+nbTF]->indices[3] = 2*(i+1)+nbTV;
		} else {
			facets[i+2+nbTF]->indices[2] = 1+nbTV;
			facets[i+2+nbTF]->indices[3] = 0+nbTV;
		}
	}

	// Volatile facet
	for(int d=0;d<nbDecade;d++) {
		for(int i=0;i<9;i++) {

			double z = (double)(i+1) * pow(10,(double)d);
			int idx = d*36+i*4;

			vertices3[idx+0].x = -R;
			vertices3[idx+0].y =  R;
			vertices3[idx+0].z =  z;
			vertices3[idx+1].x =  R;
			vertices3[idx+1].y =  R;
			vertices3[idx+1].z =  z;
			vertices3[idx+2].x =  R;
			vertices3[idx+2].y = -R;
			vertices3[idx+2].z =  z;
			vertices3[idx+3].x = -R;
			vertices3[idx+3].y = -R;
			vertices3[idx+3].z =  z;

			facets[9*d+i] = new Facet(4);
			facets[9*d+i]->sh.sticking = 0.0;
			facets[9*d+i]->sh.opacity = 0.0;
			facets[9*d+i]->sh.isVolatile = TRUE;
			facets[9*d+i]->indices[0] = idx+0;
			facets[9*d+i]->indices[1] = idx+1;
			facets[9*d+i]->indices[2] = idx+2;
			facets[9*d+i]->indices[3] = idx+3;

		}
	}
	} catch (std::bad_alloc) {
		Clear();
		throw Error("Couldn't reserve memory for the facets");
	} catch (...) {
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

void Geometry::Merge(int nbV,int nbF,VERTEX3D *nV,Facet **nF) {
	changedSinceSave = TRUE;
	// Merge the current geometry with the specified one
	if( !nbV || !nbF ) return;

	// Reallocate mem
	Facet   **nFacets    = (Facet **)malloc((sh.nbFacet+nbF) * sizeof(Facet *));
	VERTEX3D *nVertices3 = (VERTEX3D *)malloc((sh.nbVertex+nbV) * sizeof(VERTEX3D));


	if(sh.nbFacet) memcpy(nFacets,facets,sizeof(Facet *) * sh.nbFacet);
	memcpy(nFacets+sh.nbFacet,nF,sizeof(Facet *) * nbF);

	if(sh.nbVertex) memcpy(nVertices3,vertices3,sizeof(VERTEX3D) * sh.nbVertex);
	memcpy(nVertices3+sh.nbVertex,nV,sizeof(VERTEX3D) * nbV);

	SAFE_FREE(facets);
	SAFE_FREE(vertices3);
	facets = nFacets;
	vertices3 = nVertices3;
	//UnselectAllVertex();

	// Shift indices
	for(int i=0;i<nbF;i++) {
		Facet *f = facets[sh.nbFacet + i];
		for(int j=0;j<f->sh.nbIndex;j++)
			f->indices[j] += sh.nbVertex;
	}

	sh.nbVertex += nbV;
	sh.nbFacet += nbF;

}

// -----------------------------------------------------------

void Geometry::UpdateName(FileReader *file) {

	char *p = strrchr(file->GetName(),'\\');
	if(!p) p = strrchr(file->GetName(),'/');

	if(!p)
		strcpy(sh.name,file->GetName());
	else
		strcpy(sh.name,p+1);

}

// -----------------------------------------------------------

void Geometry::AdjustProfile() {

	// Backward compatibily with TXT profile (To be improved)
	for(int i=0;i<sh.nbFacet;i++) {
		Facet *f = facets[i];
		if( f->sh.profileType==REC_PRESSUREU ) {
			VERTEX3D v0;
			Sub(&v0 , vertices3+(f->indices[1]) , vertices3+(f->indices[0])); // v0 = P0P1
			double n0 = Norme(&v0);
			double nU = Norme(&(f->sh.U));
			if(IS_ZERO(n0-nU)) f->sh.profileType = REC_PRESSUREU; // Select U
			else               f->sh.profileType = REC_PRESSUREV; // Select V
		}
	}

}

// -----------------------------------------------------------

void Geometry::LoadASE(FileReader *file,GLProgress *prg) {

	Clear();

	MolFlow *mApp = (MolFlow *)theApp;
	mApp->ClearAllSelections();
	mApp->ClearAllViews();
	ASELoader ase(file);
	ase.Load();

	// Compute total of facet
	sh.nbFacet = 0;
	for(int i=0;i<ase.nbObj;i++) sh.nbFacet += ase.OBJ[i].nb_face;

	// Allocate mem
	sh.nbVertex = 3*sh.nbFacet;
	facets = (Facet **)malloc(sh.nbFacet * sizeof(Facet *));
	memset(facets,0,sh.nbFacet * sizeof(Facet *));
	vertices3 = (VERTEX3D *)malloc(sh.nbVertex * sizeof(VERTEX3D));
	memset(vertices3,0,sh.nbVertex * sizeof(VERTEX3D));

	// Fill 
	int nb = 0;
	for(int i=0;i<ase.nbObj;i++) {

		for(int j=0;j<ase.OBJ[i].nb_face;j++ ) {
			vertices3[3*nb+0] = ase.OBJ[i].pts[ase.OBJ[i].face[j].v1];
			vertices3[3*nb+1] = ase.OBJ[i].pts[ase.OBJ[i].face[j].v2];
			vertices3[3*nb+2] = ase.OBJ[i].pts[ase.OBJ[i].face[j].v3];
			facets[nb] = new Facet(3);
			facets[nb]->indices[0]=3*nb+0;
			facets[nb]->indices[1]=3*nb+1;
			facets[nb]->indices[2]=3*nb+2;
			nb++;
		}

	}

	UpdateName(file);
	sh.nbSuper = 1;
	strName[0] = _strdup(sh.name);
	strFileName[0] = _strdup(file->GetName());
	char *e = strrchr(strName[0],'.');
	if(e) *e=0;
	InitializeGeometry();
	isLoaded = TRUE;

}

// -----------------------------------------------------------

void Geometry::LoadSTR(FileReader *file,GLProgress *prg) {

	char nPath[512];
	char fPath[512];
	char fName[512];
	char sName[512];
	int nF,nV;
	Facet **F;
	VERTEX3D *V;
	FileReader *fr;

	Clear();

	MolFlow *mApp = (MolFlow *)theApp;
	mApp->ClearAllSelections();
	mApp->ClearAllViews();
	// Load multiple structure file
	sh.nbSuper = file->ReadInt();

	strcpy(fPath,file->ReadLine());
	strcpy(nPath,FileUtils::GetPath(file->GetName()));

	for(int n=0;n<sh.nbSuper;n++) {

		int i1 = file->ReadInt();
		int i2 = file->ReadInt();
		fr = NULL;
		strcpy(sName,file->ReadWord());
		strName[n] = _strdup(sName);
		char *e = strrchr(strName[n],'.');
		if(e) *e=0;

		sprintf(fName,"%s\\%s",nPath,sName);
		if( FileUtils::Exist(fName) ) {
			fr = new FileReader(fName);
			strcpy(strPath,nPath);
		} else {
			sprintf(fName,"%s\\%s",fPath,sName);
			if( FileUtils::Exist(fName) ) {
				fr = new FileReader(fName);
				strcpy(strPath,fPath);
			}
		}

		if( !fr ) {
			char errMsg[512];
			sprintf(errMsg,"Cannot find %s",sName);
			throw Error(errMsg);
		}

		strFileName[n] = _strdup(sName);
		LoadTXTGeom(fr,&nV,&nF,&V,&F,n);
		Merge(nV,nF,V,F);
		SAFE_FREE(V);
		SAFE_FREE(F);
		delete fr;

	}

	UpdateName(file);
	InitializeGeometry();
	AdjustProfile();
	isLoaded = TRUE;

}

void Geometry::LoadSTL(FileReader *file,GLProgress *prg,double scaleFactor) {

	MolFlow *mApp = (MolFlow *)theApp;
	mApp->ClearAllSelections();
	mApp->ClearAllViews();
	char *w;

	prg->SetMessage("Clearing current geometry...");
	Clear();

	// First pass
	prg->SetMessage("Counting facets in STL file...");
	file->ReadKeyword("solid");
	file->ReadLine(); // solid name
	w = file->ReadWord();
	while(strcmp(w,"facet")==0) {
		sh.nbFacet++;
		file->JumpSection("endfacet");
		w = file->ReadWord();
	}
	if( strcmp(w,"endsolid")!=0 ) throw Error("Unexpected or not supported STL keyword, 'endsolid' required");

	// Allocate mem
	sh.nbVertex = 3*sh.nbFacet;
	facets = (Facet **)malloc(sh.nbFacet * sizeof(Facet *));
	memset(facets,0,sh.nbFacet * sizeof(Facet *));
	vertices3 = (VERTEX3D *)malloc(sh.nbVertex * sizeof(VERTEX3D));
	memset(vertices3,0,sh.nbVertex * sizeof(VERTEX3D));

	// Second pass
	prg->SetMessage("Reading facets...");
	file->SeekStart();
	file->ReadKeyword("solid");
	file->ReadLine();
	for(int i=0;i<sh.nbFacet;i++) {

		double p = (double)i/(double)(sh.nbFacet);
		prg->SetProgress(p);

		file->ReadKeyword("facet");
		file->ReadKeyword("normal");
		file->ReadDouble();
		file->ReadDouble();
		file->ReadDouble();
		file->ReadKeyword("outer");
		file->ReadKeyword("loop");

		file->ReadKeyword("vertex");
		vertices3[3*i+0].x = file->ReadDouble()*scaleFactor;
		vertices3[3*i+0].y = file->ReadDouble()*scaleFactor;
		vertices3[3*i+0].z = file->ReadDouble()*scaleFactor;

		file->ReadKeyword("vertex");
		vertices3[3*i+1].x = file->ReadDouble()*scaleFactor;
		vertices3[3*i+1].y = file->ReadDouble()*scaleFactor;
		vertices3[3*i+1].z = file->ReadDouble()*scaleFactor;

		file->ReadKeyword("vertex");
		vertices3[3*i+2].x = file->ReadDouble()*scaleFactor;
		vertices3[3*i+2].y = file->ReadDouble()*scaleFactor;
		vertices3[3*i+2].z = file->ReadDouble()*scaleFactor;

		file->ReadKeyword("endloop");
		file->ReadKeyword("endfacet");

		facets[i] = new Facet(3);
		facets[i]->indices[0]=3*i+0;
		facets[i]->indices[1]=3*i+2;
		facets[i]->indices[2]=3*i+1;

	}

	sh.nbSuper = 1;
	UpdateName(file);
	strName[0] = _strdup(sh.name);
	strFileName[0] = _strdup(file->GetName());
	char *e = strrchr(strName[0],'.');
	if(e) *e=0;
	prg->SetMessage("Initializing geometry...");
	InitializeGeometry();
	isLoaded = TRUE;

}

// -----------------------------------------------------------

void Geometry::LoadTXT(FileReader *file,GLProgress *prg) {

	MolFlow *mApp = (MolFlow *)theApp;
	mApp->ClearAllSelections();
	mApp->ClearAllViews();
	Clear();
	LoadTXTGeom(file,&(sh.nbVertex),&(sh.nbFacet),&vertices3,&facets);
	UpdateName(file);
	sh.nbSuper = 1;
	strName[0] = _strdup(sh.name);
	strFileName[0] = _strdup(file->GetName());
	char *e = strrchr(strName[0],'.');
	if(e) *e=0;
	InitializeGeometry();
	AdjustProfile();
	isLoaded = TRUE;

}

// -----------------------------------------------------------

void Geometry::InsertTXT(FileReader *file,GLProgress *prg,BOOL newStr) {

	//Clear();
	int structId=viewStruct;
	if (structId==-1) structId=0;
	InsertTXTGeom(file,&(sh.nbVertex),&(sh.nbFacet),&vertices3,&facets,structId,newStr);
	//UpdateName(file);
	//sh.nbSuper = 1;
	//strName[0] = _strdup(sh.name);
	//strFileName[0] = _strdup(file->GetName());
	char *e = strrchr(strName[0],'.');
	if(e) *e=0;
	InitializeGeometry();
	AdjustProfile();
	isLoaded = TRUE;

}

// -----------------------------------------------------------

void Geometry::InsertSTL(FileReader *file,GLProgress *prg,double scaleFactor,BOOL newStr) {

	//Clear();
	int structId=viewStruct;
	if (structId==-1) structId=0;
	InsertSTLGeom(file,&(sh.nbVertex),&(sh.nbFacet),&vertices3,&facets,structId,scaleFactor,newStr);
	//UpdateName(file);
	//sh.nbSuper = 1;
	//strName[0] = _strdup(sh.name);
	//strFileName[0] = _strdup(file->GetName());
	char *e = strrchr(strName[0],'.');
	if(e) *e=0;
	InitializeGeometry();
	AdjustProfile();
	isLoaded = TRUE;

}

// -----------------------------------------------------------

void Geometry::InsertGEO(FileReader *file,GLProgress *prg,BOOL newStr) {

	//Clear();
	int structId=viewStruct;
	if (structId==-1) structId=0;
	InsertGEOGeom(file,&(sh.nbVertex),&(sh.nbFacet),&vertices3,&facets,structId,newStr);
	//UpdateName(file);
	//sh.nbSuper = 1;
	//strName[0] = _strdup(sh.name);
	//strFileName[0] = _strdup(file->GetName());
	char *e = strrchr(strName[0],'.');
	if(e) *e=0;
	InitializeGeometry();
	AdjustProfile();
	isLoaded = TRUE;

}

void Geometry::InsertSYN(FileReader *file,GLProgress *prg,BOOL newStr) {

	int structId=viewStruct;
	if (structId==-1) structId=0;
	InsertSYNGeom(file,&(sh.nbVertex),&(sh.nbFacet),&vertices3,&facets,structId,newStr);
	char *e = strrchr(strName[0],'.');
	if(e) *e=0;
	InitializeGeometry();
	AdjustProfile();
}

// -----------------------------------------------------------

void Geometry::LoadTXTGeom(FileReader *file,int *nbV,int *nbF,VERTEX3D **V,Facet ***F,int strIdx) {

	file->ReadInt(); // Unused
	tNbHit = file->ReadLLong();
	tNbLeak = file->ReadLLong();
	tNbDesorption = file->ReadLLong();
	tNbDesorptionMax = file->ReadLLong(); 

	int nV = file->ReadInt();
	int nF = file->ReadInt();

	// Allocate memory
	Facet   **f = (Facet **)malloc(nF * sizeof(Facet *));
	memset(f,0,nF * sizeof(Facet *));
	VERTEX3D *v = (VERTEX3D *)malloc(nV * sizeof(VERTEX3D));

	// Read geometry vertices
	for(int i=0;i<nV;i++) {
		v[i].x = file->ReadDouble();
		v[i].y = file->ReadDouble();
		v[i].z = file->ReadDouble();
	}

	// Read geometry facets (indexed from 1)
	for(int i=0;i<nF;i++) {
		int nb = file->ReadInt();
		f[i] = new Facet(nb);
		for(int j=0;j<nb;j++)
			f[i]->indices[j]=file->ReadInt()-1;
	}

	// Read facets params
	for(int i=0;i<nF;i++) {
		f[i]->LoadTXT(file);
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

void Geometry::InsertTXTGeom(FileReader *file,int *nbVertex,int *nbFacet,VERTEX3D **vertices3,Facet ***facets,int strIdx,BOOL newStruct) {

	UnSelectAll();

	//tNbHit = file->ReadLLong();
	//tNbLeak = file->ReadInt();
	//tNbDesorption = file->ReadLLong();
	//tNbDesorptionMax = file->ReadLLong(); 
	for (int i=0;i<5;i++) file->ReadInt(); //leading lines

	int nbNewVertex = file->ReadInt();
	int nbNewFacets = file->ReadInt();

	// Allocate memory
	*facets = (Facet **)realloc(*facets,(nbNewFacets+*nbFacet) * sizeof(Facet **));
	memset(*facets+*nbFacet,0,nbNewFacets * sizeof(Facet *));
	//*vertices3 = (VERTEX3D*)realloc(*vertices3,(nbNewVertex+*nbVertex) * sizeof(VERTEX3D));
	VERTEX3D *tmp_vertices3 = (VERTEX3D *)malloc((nbNewVertex+*nbVertex) * sizeof(VERTEX3D));
	memmove(tmp_vertices3,*vertices3,(*nbVertex)*sizeof(VERTEX3D));
	memset(tmp_vertices3+*nbVertex,0,nbNewVertex * sizeof(VERTEX3D));
	SAFE_FREE(*vertices3);
	*vertices3=tmp_vertices3;

	// Read geometry vertices
	for(int i=*nbVertex;i<(*nbVertex+nbNewVertex);i++) {
		(*vertices3+i)->x = file->ReadDouble();
		(*vertices3+i)->y = file->ReadDouble();
		(*vertices3+i)->z = file->ReadDouble();
		(*vertices3+i)->selected = FALSE;
	}

	// Read geometry facets (indexed from 1)
	for(int i=*nbFacet;i<(*nbFacet+nbNewFacets);i++) {
		int nb = file->ReadInt();
		*(*facets+i) = new Facet(nb);
		(*facets)[i]->selected = TRUE;
		for(int j=0;j<nb;j++)
			(*facets)[i]->indices[j]=file->ReadInt()-1+*nbVertex;
	}

	// Read facets params
	for(int i=*nbFacet;i<(*nbFacet+nbNewFacets);i++) {
		(*facets)[i]->LoadTXT(file);
		if (newStruct) {
			(*facets)[i]->sh.superIdx = sh.nbSuper;
		} else {
			(*facets)[i]->sh.superIdx = strIdx;
		}
	}

	*nbVertex += nbNewVertex;
	*nbFacet += nbNewFacets;
	if (newStruct) sh.nbSuper++;

}

// -----------------------------------------------------------

void Geometry::InsertGEOGeom(FileReader *file,int *nbVertex,int *nbFacet,VERTEX3D **vertices3,Facet ***facets,int strIdx,BOOL newStruct) {

	MolFlow *mApp = (MolFlow *)theApp;
	UnSelectAll();
	char tmp[512];

	file->ReadKeyword("version");file->ReadKeyword(":");
	int version2;
	version2 = file->ReadInt();

	file->ReadKeyword("totalHit");file->ReadKeyword(":");
	file->ReadLLong();
	file->ReadKeyword("totalDes");file->ReadKeyword(":");
	file->ReadLLong();
	file->ReadKeyword("totalLeak");file->ReadKeyword(":");
	file->ReadLLong();
	if (version2>=12) {
		file->ReadKeyword("totalAbs");file->ReadKeyword(":");
		file->ReadLLong();
		file->ReadKeyword("totalDist");file->ReadKeyword(":");
		file->ReadDouble();
	}
	file->ReadKeyword("maxDes");file->ReadKeyword(":");
	file->ReadLLong(); 
	file->ReadKeyword("nbVertex");file->ReadKeyword(":");
	int nbNewVertex = file->ReadInt();
	file->ReadKeyword("nbFacet");file->ReadKeyword(":");
	int nbNewFacets = file->ReadInt();
	file->ReadKeyword("nbSuper");file->ReadKeyword(":");
	int nbNewSuper = file->ReadInt();
	int nbF=0;
	int nbV=0;
	if (version2>=2) {
		file->ReadKeyword("nbFormula");file->ReadKeyword(":");
		nbF = file->ReadInt();
		file->ReadKeyword("nbView");file->ReadKeyword(":");
		nbV = file->ReadInt();
	}
	int nbS=0;
	if (version2>=8) {
		file->ReadKeyword("nbSelection");file->ReadKeyword(":");
		nbS = file->ReadInt();
	}
	if (version2>=7) {
		file->ReadKeyword("gasMass");file->ReadKeyword(":");
		gasMass = file->ReadDouble();
	}
	if (version2>=10) { //time-dependent version
		file->ReadKeyword("userMoments");file->ReadKeyword("{");
		file->ReadKeyword("nb");file->ReadKeyword(":");
		int nb=file->ReadInt();

		for(int i=0;i<nb;i++) 
			file->ReadString();
		file->ReadKeyword("}");
	}
	if (version2>=11) { //gas pulse parameters
		file->ReadKeyword("desorptionStart");file->ReadKeyword(":");
		file->ReadDouble();
		file->ReadKeyword("desorptionStop");file->ReadKeyword(":");
		file->ReadDouble();
		file->ReadKeyword("timeWindow");file->ReadKeyword(":");
		file->ReadDouble();
		file->ReadKeyword("useMaxwellian");file->ReadKeyword(":");
		file->ReadInt();
	}

	if (version2>=12) { //2013.aug.22
		file->ReadKeyword("calcConstantFlow");file->ReadKeyword(":");
		file->ReadInt();
	}

	if (version2>=2) {
		file->ReadKeyword("formulas");file->ReadKeyword("{");
		for(int i=0;i<nbF;i++) {
			char tmpName[256];
			char tmpExpr[512];
			strcpy(tmpName,file->ReadString());
			strcpy(tmpExpr,file->ReadString());
			mApp->OffsetFormula(tmpExpr,sh.nbFacet);
			mApp->AddFormula(tmpName,tmpExpr);
		}
		file->ReadKeyword("}");

		file->ReadKeyword("views");file->ReadKeyword("{");
		for(int i=0;i<nbV;i++) {
			char tmpName[256];
			AVIEW v;
			strcpy(tmpName,file->ReadString());
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
			mApp->AddView(tmpName,v);
		}
		file->ReadKeyword("}");
	}

	if (version2>=8) {
		file->ReadKeyword("selections");file->ReadKeyword("{");
		for(int i=0;i<nbS;i++) {
			ASELECTION s;
			char tmpName[256];
			strcpy(tmpName,file->ReadString());
			s.name=_strdup(tmpName);
			s.nbSel = file->ReadInt();
			s.selection = (int *)malloc((s.nbSel)*sizeof(int));

			for (int j=0;j<s.nbSel;j++) {
				s.selection[j]=file->ReadInt()+sh.nbFacet; //offset facet number by current number of facets
			}
			mApp->AddSelection(s.name,s);
		}
		file->ReadKeyword("}");
	}

	file->ReadKeyword("structures");file->ReadKeyword("{");
	for(int i=0;i<nbNewSuper;i++) {
		strName[i] = _strdup(file->ReadString());
		// For backward compatibilty with STR
		sprintf(tmp,"%s.txt",strName[i]);
		strFileName[i] = _strdup(tmp);
	}
	file->ReadKeyword("}");

	// Reallocate memory
	*facets = (Facet **)realloc(*facets,(nbNewFacets+*nbFacet) * sizeof(Facet **));
	memset(*facets+*nbFacet,0,nbNewFacets * sizeof(Facet *));
	//*vertices3 = (VERTEX3D*)realloc(*vertices3,(nbNewVertex+*nbVertex) * sizeof(VERTEX3D));
	VERTEX3D *tmp_vertices3 = (VERTEX3D *)malloc((nbNewVertex+*nbVertex) * sizeof(VERTEX3D));
	memmove(tmp_vertices3,*vertices3,(*nbVertex)*sizeof(VERTEX3D));
	memset(tmp_vertices3+*nbVertex,0,nbNewVertex * sizeof(VERTEX3D));
	SAFE_FREE(*vertices3);
	*vertices3=tmp_vertices3;

	// Read geometry vertices
	file->ReadKeyword("vertices");file->ReadKeyword("{");
	for(int i=*nbVertex;i<(*nbVertex+nbNewVertex);i++) {
		// Check idx
		int idx = file->ReadInt();
		if( idx != i-*nbVertex+1 ) throw Error(file->MakeError("Wrong vertex index !"));
		(*vertices3+i)->x = file->ReadDouble();
		(*vertices3+i)->y = file->ReadDouble();
		(*vertices3+i)->z = file->ReadDouble();
		(*vertices3+i)->selected = FALSE;
	}
	file->ReadKeyword("}");

	if (version2>=6) {
		// Read leaks
		file->ReadKeyword("leaks");file->ReadKeyword("{");
		file->ReadKeyword("nbLeak");file->ReadKeyword(":");
		int nbleak2 = file->ReadInt();
		for (int i=0;i<nbleak2;i++) {
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
		file->ReadKeyword("hits");file->ReadKeyword("{");
		file->ReadKeyword("nbHHit");file->ReadKeyword(":");
		int nbHHit2 = file->ReadInt();
		for (int i=0;i<nbHHit2;i++) {
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
	for(int i=*nbFacet;i<(*nbFacet+nbNewFacets);i++) {
		file->ReadKeyword("facet");  
		// Check idx
		int idx = file->ReadInt();
		if( idx != i+1-*nbFacet ) throw Error(file->MakeError("Wrong facet index !"));
		file->ReadKeyword("{");
		file->ReadKeyword("nbIndex");
		file->ReadKeyword(":");
		int nb = file->ReadInt();

		if( nb<3 ) {
			char errMsg[512];
			sprintf(errMsg,"Facet %d has only %d vertices. ",i,nb);
			throw Error(errMsg);
		}

		*(*facets+i) = new Facet(nb);
		(*facets)[i]->LoadGEO(file,version2,nbNewVertex);
		(*facets)[i]->selected = TRUE;
		for(int j=0;j<nb;j++)
			(*facets)[i]->indices[j]+=*nbVertex;
		file->ReadKeyword("}");
		if (newStruct) {
			(*facets)[i]->sh.superIdx += sh.nbSuper;
			if ((*facets)[i]->sh.superDest>0) (*facets)[i]->sh.superDest+= sh.nbSuper;
		} else {
			(*facets)[i]->sh.superIdx += strIdx;
			if ((*facets)[i]->sh.superDest>0) (*facets)[i]->sh.superDest+= strIdx;
		}
	}


	*nbVertex += nbNewVertex;
	*nbFacet += nbNewFacets;
	if (newStruct) sh.nbSuper+=nbNewSuper;
	else if (sh.nbSuper<strIdx+nbNewSuper) sh.nbSuper=strIdx+nbNewSuper;

}


void Geometry::InsertSYNGeom(FileReader *file,int *nbVertex,int *nbFacet,VERTEX3D **vertices3,Facet ***facets,int strIdx,BOOL newStruct) {

	MolFlow *mApp = (MolFlow *)theApp;
	UnSelectAll();
	char tmp[512];

	file->ReadKeyword("version");file->ReadKeyword(":");
	int version2;
	version2 = file->ReadInt();
	if( version2>SYNVERSION ) {
		char errMsg[512];
		sprintf(errMsg,"Unsupported SYN version V%d",version2);
		throw Error(errMsg);
	}

	file->ReadKeyword("totalHit");file->ReadKeyword(":");
	file->ReadLLong();
	file->ReadKeyword("totalDes");file->ReadKeyword(":");
	file->ReadLLong();
	file->ReadKeyword("totalLeak");file->ReadKeyword(":");
	file->ReadLLong();
	if (version2>2) {
		file->ReadKeyword("totalFlux");file->ReadKeyword(":");
		file->ReadDouble();
		file->ReadKeyword("totalPower");file->ReadKeyword(":");
		file->ReadDouble();
	}
	file->ReadKeyword("maxDes");file->ReadKeyword(":");
	file->ReadLLong(); 
	file->ReadKeyword("nbVertex");file->ReadKeyword(":");
	int nbNewVertex = file->ReadInt();
	file->ReadKeyword("nbFacet");file->ReadKeyword(":");
	int nbNewFacets = file->ReadInt();
	file->ReadKeyword("nbSuper");file->ReadKeyword(":");
	int nbNewSuper = file->ReadInt();
	int nbF=0;
	int nbV=0;

	file->ReadKeyword("nbFormula");file->ReadKeyword(":");
	nbF = file->ReadInt();
	file->ReadKeyword("nbView");file->ReadKeyword(":");
	nbV = file->ReadInt();
	int nbS=0;

	file->ReadKeyword("nbSelection");file->ReadKeyword(":");
	nbS = file->ReadInt();
	if (version2>1) {
		file->ReadKeyword("nbRegions");file->ReadKeyword(":");
		int nbR = file->ReadInt();
		//result=PARfileList(nbR);

		file->ReadKeyword("PARfiles");file->ReadKeyword("{");
		for (int i=0;i<nbR;i++) {
			/*char tmp[512];
			strcpy(tmp,file->ReadString());
			result.fileNames[i]=_strdup(tmp);*/
			file->ReadString();
		}
		file->ReadKeyword("}");
	}

	file->ReadKeyword("formulas");file->ReadKeyword("{");
	for(int i=0;i<nbF;i++) {
		char tmpName[256];
		char tmpExpr[512];
		strcpy(tmpName,file->ReadString());
		strcpy(tmpExpr,file->ReadString());
		//mApp->AddFormula(tmpName,tmpExpr); //we don't add SynRad formulas to MolFlow

	}
	file->ReadKeyword("}");

	file->ReadKeyword("views");file->ReadKeyword("{");
	for(int i=0;i<nbV;i++) {
		char tmpName[256];
		AVIEW v;
		strcpy(tmpName,file->ReadString());
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
		mApp->AddView(tmpName,v);
	}
	file->ReadKeyword("}");

	file->ReadKeyword("selections");file->ReadKeyword("{");
	for(int i=0;i<nbS;i++) {
		ASELECTION s;
		char tmpName[256];
		strcpy(tmpName,file->ReadString());
		s.name=_strdup(tmpName);
		s.nbSel = file->ReadInt();
		s.selection = (int *)malloc((s.nbSel)*sizeof(int));

		for (int j=0;j<s.nbSel;j++) {
			s.selection[j]=file->ReadInt()+sh.nbFacet;
		}
		mApp->AddSelection(s.name,s);
	}
	file->ReadKeyword("}");

	file->ReadKeyword("structures");file->ReadKeyword("{");
	for(int i=0;i<nbNewSuper;i++) {
		strName[i] = _strdup(file->ReadString());
		// For backward compatibilty with STR
		sprintf(tmp,"%s.txt",strName[i]);
		strFileName[i] = _strdup(tmp);
	}
	file->ReadKeyword("}");

	// Reallocate memory
	*facets = (Facet **)realloc(*facets,(nbNewFacets+*nbFacet) * sizeof(Facet **));
	memset(*facets+*nbFacet,0,nbNewFacets * sizeof(Facet *));
	//*vertices3 = (VERTEX3D*)realloc(*vertices3,(nbNewVertex+*nbVertex) * sizeof(VERTEX3D));
	VERTEX3D *tmp_vertices3 = (VERTEX3D *)malloc((nbNewVertex+*nbVertex) * sizeof(VERTEX3D));
	memmove(tmp_vertices3,*vertices3,(*nbVertex)*sizeof(VERTEX3D));
	memset(tmp_vertices3+*nbVertex,0,nbNewVertex * sizeof(VERTEX3D));
	SAFE_FREE(*vertices3);
	*vertices3=tmp_vertices3;

	// Read geometry vertices
	file->ReadKeyword("vertices");file->ReadKeyword("{");
	for(int i=*nbVertex;i<(*nbVertex+nbNewVertex);i++) {
		// Check idx
		int idx = file->ReadInt();
		if( idx != i-*nbVertex+1 ) throw Error(file->MakeError("Wrong vertex index !"));
		(*vertices3+i)->x = file->ReadDouble();
		(*vertices3+i)->y = file->ReadDouble();
		(*vertices3+i)->z = file->ReadDouble();
		(*vertices3+i)->selected = FALSE;
	}
	file->ReadKeyword("}");


	// Read leaks
	file->ReadKeyword("leaks");file->ReadKeyword("{");
	file->ReadKeyword("nbLeak");file->ReadKeyword(":");
	int nbleak_local = file->ReadInt();
	for (int i=0;i<nbleak_local;i++) {
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
	file->ReadKeyword("hits");file->ReadKeyword("{");
	file->ReadKeyword("nbHHit");file->ReadKeyword(":");
	int nbHHit_local = file->ReadInt();
	for (int i=0;i<nbHHit_local;i++) {
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
	for(int i=*nbFacet;i<(*nbFacet+nbNewFacets);i++) {
		file->ReadKeyword("facet");  
		// Check idx
		int idx = file->ReadInt();
		if( idx != i+1-*nbFacet ) throw Error(file->MakeError("Wrong facet index !"));
		file->ReadKeyword("{");
		file->ReadKeyword("nbIndex");
		file->ReadKeyword(":");
		int nb = file->ReadInt();

		if( nb<3 ) {
			char errMsg[512];
			sprintf(errMsg,"Facet %d has only %d vertices. ",i,nb);
			throw Error(errMsg);
		}

		*(*facets+i) = new Facet(nb);
		(*facets)[i]->LoadSYN(file,version2,nbNewVertex);
		(*facets)[i]->selected = TRUE;
		for(int j=0;j<nb;j++)
			(*facets)[i]->indices[j]+=*nbVertex;
		file->ReadKeyword("}");
		if (newStruct) {
			(*facets)[i]->sh.superIdx += sh.nbSuper;
		} else {
			(*facets)[i]->sh.superIdx = strIdx;
		}
	}


	*nbVertex += nbNewVertex;
	*nbFacet += nbNewFacets;
	if (newStruct) sh.nbSuper+=nbNewSuper;
	//return result;
}

// -----------------------------------------------------------

void Geometry::InsertSTLGeom(FileReader *file,int *nbVertex,int *nbFacet,VERTEX3D **vertices3,Facet ***facets,int strIdx,double scaleFactor,BOOL newStruct) {

	UnSelectAll();
	char *w;

	int nbNewFacets=0;
	// First pass
	file->ReadKeyword("solid");
	file->ReadLine(); // solid name
	w = file->ReadWord();
	while(strcmp(w,"facet")==0) {
		nbNewFacets++;
		file->JumpSection("endfacet");
		w = file->ReadWord();
	}
	if( strcmp(w,"endsolid")!=0 ) throw Error("Unexpected or not supported STL keyword, 'endsolid' required");

	// Allocate memory
	int nbNewVertex=3*nbNewFacets;
	*facets = (Facet **)realloc(*facets,(nbNewFacets+*nbFacet) * sizeof(Facet **));
	memset(*facets+*nbFacet,0,nbNewFacets * sizeof(Facet *));
	//*vertices3 = (VERTEX3D*)realloc(*vertices3,(nbNewVertex+*nbVertex) * sizeof(VERTEX3D));
	VERTEX3D *tmp_vertices3 = (VERTEX3D *)malloc((nbNewVertex+*nbVertex) * sizeof(VERTEX3D));
	memmove(tmp_vertices3,*vertices3,(*nbVertex)*sizeof(VERTEX3D));
	memset(tmp_vertices3+*nbVertex,0,nbNewVertex * sizeof(VERTEX3D));
	SAFE_FREE(*vertices3);
	*vertices3=tmp_vertices3;

	// Second pass
	file->SeekStart();
	file->ReadKeyword("solid");
	file->ReadLine();
	for(int i=0;i<nbNewFacets;i++) {

		file->ReadKeyword("facet");
		file->ReadKeyword("normal");
		file->ReadDouble();
		file->ReadDouble();
		file->ReadDouble();
		file->ReadKeyword("outer");
		file->ReadKeyword("loop");

		file->ReadKeyword("vertex");
		(*vertices3)[*nbVertex+3*i+0].x = file->ReadDouble()*scaleFactor;
		(*vertices3)[*nbVertex+3*i+0].y = file->ReadDouble()*scaleFactor;
		(*vertices3)[*nbVertex+3*i+0].z = file->ReadDouble()*scaleFactor;

		file->ReadKeyword("vertex");
		(*vertices3)[*nbVertex+3*i+1].x = file->ReadDouble()*scaleFactor;
		(*vertices3)[*nbVertex+3*i+1].y = file->ReadDouble()*scaleFactor;
		(*vertices3)[*nbVertex+3*i+1].z = file->ReadDouble()*scaleFactor;

		file->ReadKeyword("vertex");
		(*vertices3)[*nbVertex+3*i+2].x = file->ReadDouble()*scaleFactor;
		(*vertices3)[*nbVertex+3*i+2].y = file->ReadDouble()*scaleFactor;
		(*vertices3)[*nbVertex+3*i+2].z = file->ReadDouble()*scaleFactor;

		file->ReadKeyword("endloop");
		file->ReadKeyword("endfacet");

		*(*facets+i+*nbFacet) = new Facet(3);
		(*facets)[i+*nbFacet]->selected = TRUE;
		(*facets)[i+*nbFacet]->indices[0] = *nbVertex+3*i+0;
		(*facets)[i+*nbFacet]->indices[1] = *nbVertex+3*i+1;
		(*facets)[i+*nbFacet]->indices[2] = *nbVertex+3*i+2;

		if (newStruct) {
			(*facets)[i+*nbFacet]->sh.superIdx = sh.nbSuper;
		} else {
			(*facets)[i+*nbFacet]->sh.superIdx = strIdx;
		}
	}

	*nbVertex += nbNewVertex;
	*nbFacet += nbNewFacets;
	if (newStruct) sh.nbSuper++;

}

// -----------------------------------------------------------

void Geometry::SaveProfile(FileWriter *file,Dataport *dpHit,int super,BOOL saveSelected,BOOL txtFormat) {
	MolFlow *mApp = (MolFlow *)theApp;
	BYTE *buffer = (BYTE *)dpHit->buff;
	if (!txtFormat) file->Write("profiles {\n");
	// Profiles
	int nbProfile=0;
	int *profileFacet = (int *)malloc((sh.nbFacet)*sizeof(int));
	for(int i=0;i<sh.nbFacet;i++) 
		if (((!saveSelected)||(facets[i]->selected)) && facets[i]->sh.isProfile)
			profileFacet[nbProfile++] = i;
	if (!txtFormat) {
		file->Write(" number: ");file->WriteInt(nbProfile,"\n");
		file->Write(" facets: ");
		for (int i=0;i<nbProfile;i++) 
			file->WriteInt(profileFacet[i],"\t");
	}
	file->Write("\n");
	for (size_t m=0;(m<=mApp->worker.moments.size()) || (m==0 && !txtFormat);m++){
		char tmp[128];
		sprintf(tmp," moment %d {\n",m);
		file->Write(tmp);
		for(int j=0;j<PROFILE_SIZE;j++) {
			for(int i=0;i<nbProfile;i++) {
				Facet *f = GetFacet(profileFacet[i]);
				llong *profilePtr = (llong *)(buffer + f->sh.hitOffset + sizeof(SHHITS)+m*sizeof(llong)*PROFILE_SIZE);
				char tmp2[128];
				sprintf(tmp2,"%9I64d",profilePtr[j]);
				file->Write(tmp2);
				if (!txtFormat) file->Write("\t"); else file->Write(" ");
			}
			if (nbProfile>0) file->Write("\n");
		}		
		file->Write(" }\n");
	}

	if (!txtFormat) file->Write("}\n");
	SAFE_FREE(profileFacet);
}

// -----------------------------------------------------------

void Geometry::LoadProfile(FileReader *file,Dataport *dpHit,int version) {
	MolFlow *mApp = (MolFlow *)theApp;
	AccessDataport(dpHit);
	BYTE *buffer = (BYTE *)dpHit->buff;
	file->ReadKeyword("profiles");file->ReadKeyword("{");
	// Profiles
	int nbProfile;
	file->ReadKeyword("number");file->ReadKeyword(":");nbProfile=file->ReadInt();
	int *profileFacet = (int *)malloc((nbProfile)*sizeof(int));
	file->ReadKeyword("facets");file->ReadKeyword(":");
	for (int i=0;i<nbProfile;i++) 
		profileFacet[i]=file->ReadInt();
	for (size_t m=0;m<=mApp->worker.moments.size() || (version<10 && m==0);m++){
		if (version>=10) {
			file->ReadKeyword("moment");
			if(m!=file->ReadInt()) {
				char errMsg[512];
				sprintf(errMsg,"Unexpected profile moment");
				throw Error(errMsg);
				break;
			}
			file->ReadKeyword("{");
		}
		
		for(int j=0;j<PROFILE_SIZE;j++) {
			for(int i=0;i<nbProfile;i++) {
				Facet *f = GetFacet(profileFacet[i]);
				llong *profilePtr = (llong *)(buffer + f->sh.hitOffset + sizeof(SHHITS)+m*PROFILE_SIZE*sizeof(llong));
				profilePtr[j]=file->ReadLLong();
			}
		}
		if (version>=10) file->ReadKeyword("}");
	}
	ReleaseDataport(dpHit);
	SAFE_FREE(profileFacet);
}


// -----------------------------------------------------------
extern GLApplication *theApp;

void Geometry::LoadGEO(FileReader *file,GLProgress *prg,LEAK *pleak,int *nbleak,HIT *pHits,int *nbHHit,int *version,Worker *worker) {

	MolFlow *mApp = (MolFlow *)theApp;
	mApp->ClearAllSelections();
	mApp->ClearAllViews();
	prg->SetMessage("Clearing current geometry...");
	Clear();
	mApp->ClearFormula();


	// Globals
	char tmp[512];
	prg->SetMessage("Reading GEO file header...");
	file->ReadKeyword("version");file->ReadKeyword(":");
	*version = file->ReadInt();
	if( *version>GEOVERSION ) {
		char errMsg[512];
		sprintf(errMsg,"Unsupported GEO version V%d",*version);
		throw Error(errMsg);
	}

	file->ReadKeyword("totalHit");file->ReadKeyword(":");
	tNbHit = file->ReadLLong();
	file->ReadKeyword("totalDes");file->ReadKeyword(":");
	tNbDesorption = file->ReadLLong();
	file->ReadKeyword("totalLeak");file->ReadKeyword(":");
	tNbLeak = file->ReadLLong();
	if (*version>=12) {
		file->ReadKeyword("totalAbs");file->ReadKeyword(":");
		tNbAbsorption = file->ReadLLong();
		file->ReadKeyword("totalDist");file->ReadKeyword(":");
		distTraveledTotal = file->ReadDouble();
	} else {
		tNbAbsorption=0;
		distTraveledTotal=0.0;
	}
	file->ReadKeyword("maxDes");file->ReadKeyword(":");
	tNbDesorptionMax = file->ReadLLong(); 
	file->ReadKeyword("nbVertex");file->ReadKeyword(":");
	sh.nbVertex = file->ReadInt();
	file->ReadKeyword("nbFacet");file->ReadKeyword(":");
	sh.nbFacet = file->ReadInt();
	file->ReadKeyword("nbSuper");file->ReadKeyword(":");
	sh.nbSuper = file->ReadInt();
	int nbF=0;
	int nbV=0;
	if (*version>=2) {
		file->ReadKeyword("nbFormula");file->ReadKeyword(":");
		nbF = file->ReadInt();
		file->ReadKeyword("nbView");file->ReadKeyword(":");
		nbV = file->ReadInt();
	}
	int nbS=0;
	if (*version>=8) {
		file->ReadKeyword("nbSelection");file->ReadKeyword(":");
		nbS = file->ReadInt();
	}
	if (*version>=7) {
		file->ReadKeyword("gasMass");file->ReadKeyword(":");
		gasMass = file->ReadDouble();
	}
	if (*version>=10) { //time-dependent version
		file->ReadKeyword("userMoments");file->ReadKeyword("{");
		file->ReadKeyword("nb");file->ReadKeyword(":");
		int nb=file->ReadInt();

		for(int i=0;i<nb;i++) {
			char tmpExpr[512];
			strcpy(tmpExpr,file->ReadString());
			mApp->worker.userMoments.push_back(tmpExpr);
			mApp->worker.AddMoment(mApp->worker.ParseMoment(tmpExpr));
		}
		file->ReadKeyword("}");
	}
	if (*version>=11) { //pulse version
		file->ReadKeyword("desorptionStart");file->ReadKeyword(":");
		worker->desorptionStartTime=file->ReadDouble();
		file->ReadKeyword("desorptionStop");file->ReadKeyword(":");
		worker->desorptionStopTime=file->ReadDouble();
		file->ReadKeyword("timeWindow");file->ReadKeyword(":");
		worker->timeWindowSize=file->ReadDouble();
		file->ReadKeyword("useMaxwellian");file->ReadKeyword(":");
		worker->useMaxwellDistribution=file->ReadInt();
	}
	if (*version>=12) { //2013.aug.22
		file->ReadKeyword("calcConstantFlow");file->ReadKeyword(":");
		worker->calcConstantFlow=file->ReadInt();
	}
	if (*version>=2) {
		file->ReadKeyword("formulas");file->ReadKeyword("{");
		for(int i=0;i<nbF;i++) {
			char tmpName[256];
			char tmpExpr[512];
			strcpy(tmpName,file->ReadString());
			strcpy(tmpExpr,file->ReadString());
			mApp->AddFormula(tmpName,tmpExpr);
		}
		file->ReadKeyword("}");

		file->ReadKeyword("views");file->ReadKeyword("{");
		for(int i=0;i<nbV;i++) {
			char tmpName[256];
			AVIEW v;
			strcpy(tmpName,file->ReadString());
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
			mApp->AddView(tmpName,v);
		}
		file->ReadKeyword("}");
	}

	if (*version>=8) {
		file->ReadKeyword("selections");file->ReadKeyword("{");
		for(int i=0;i<nbS;i++) {
			ASELECTION s;
			char tmpName[256];
			strcpy(tmpName,file->ReadString());
			s.name=_strdup(tmpName);
			s.nbSel = file->ReadInt();
			s.selection = (int *)malloc((s.nbSel)*sizeof(int));

			for (int j=0;j<s.nbSel;j++) {
				s.selection[j]=file->ReadInt();
			}
			mApp->AddSelection(s.name,s);
		}
		file->ReadKeyword("}");
	}

	file->ReadKeyword("structures");file->ReadKeyword("{");
	for(int i=0;i<sh.nbSuper;i++) {
		strName[i] = _strdup(file->ReadString());
		// For backward compatibilty with STR
		sprintf(tmp,"%s.txt",strName[i]);
		strFileName[i] = _strdup(tmp);
	}
	file->ReadKeyword("}");

	// Allocate memory
	facets = (Facet **)malloc(sh.nbFacet * sizeof(Facet *));
	memset(facets,0,sh.nbFacet * sizeof(Facet *));
	vertices3 = (VERTEX3D *)malloc(sh.nbVertex * sizeof(VERTEX3D));

	// Read vertices
	prg->SetMessage("Reading vertices...");
	file->ReadKeyword("vertices");file->ReadKeyword("{");
	for(int i=0;i<sh.nbVertex;i++) {
		// Check idx
		int idx = file->ReadInt();
		if( idx != i+1 ) throw Error(file->MakeError("Wrong vertex index !"));
		vertices3[i].x = file->ReadDouble();
		vertices3[i].y = file->ReadDouble();
		vertices3[i].z = file->ReadDouble();
	}
	file->ReadKeyword("}");

	if (*version>=6) {
		prg->SetMessage("Reading leaks and hits...");
		// Read leaks
		file->ReadKeyword("leaks");file->ReadKeyword("{");
		file->ReadKeyword("nbLeak");file->ReadKeyword(":");
		*nbleak = file->ReadInt();
		for (int i=0;i<*nbleak;i++) {
			int idx = file->ReadInt();
			if( idx != i ) throw Error(file->MakeError("Wrong leak index !"));
			(pleak+i)->pos.x = file->ReadDouble();
			(pleak+i)->pos.y = file->ReadDouble();
			(pleak+i)->pos.z = file->ReadDouble();

			(pleak+i)->dir.x = file->ReadDouble();
			(pleak+i)->dir.y = file->ReadDouble();
			(pleak+i)->dir.z = file->ReadDouble();   
		}
		file->ReadKeyword("}");

		// Read hit cache
		file->ReadKeyword("hits");file->ReadKeyword("{");
		file->ReadKeyword("nbHHit");file->ReadKeyword(":");
		*nbHHit = file->ReadInt();
		for (int i=0;i<*nbHHit;i++) {
			int idx = file->ReadInt();
			if( idx != i ) throw Error(file->MakeError("Wrong hit cache index !"));
			(pHits+i)->pos.x = file->ReadDouble();
			(pHits+i)->pos.y = file->ReadDouble();
			(pHits+i)->pos.z = file->ReadDouble();

			(pHits+i)->type = file->ReadInt();
		}
		file->ReadKeyword("}");
	}

	// Read facets
	prg->SetMessage("Reading facets...");
	for(int i=0;i<sh.nbFacet;i++) {
		file->ReadKeyword("facet");
		// Check idx
		int idx = file->ReadInt();
		if( idx != i+1 ) throw Error(file->MakeError("Wrong facet index !"));
		file->ReadKeyword("{");
		file->ReadKeyword("nbIndex");
		file->ReadKeyword(":");
		int nbI = file->ReadInt();
		if( nbI<3 ) {
			char errMsg[512];
			sprintf(errMsg,"Facet %d has only %d vertices. ",i+1,nbI);
			throw Error(errMsg);
		}
		prg->SetProgress((float)i/sh.nbFacet);
		facets[i] = new Facet(nbI);
		facets[i]->LoadGEO(file,*version,sh.nbVertex);
		file->ReadKeyword("}");
	}

	InitializeGeometry();
	AdjustProfile();
	isLoaded = TRUE;
	UpdateName(file);

	// Update mesh
	prg->SetMessage("Building mesh...");
	for(int i=0;i<sh.nbFacet;i++) {
		double p = (double)i/(double)sh.nbFacet;
		prg->SetProgress(p);
		Facet *f = facets[i];
		if (!f->SetTexture(f->sh.texWidthD,f->sh.texHeightD,f->hasMesh)) {
			char errMsg[512];
			sprintf(errMsg,"Not enough memory to build mesh on Facet %d. ",i+1);
			throw Error(errMsg);
		}
		BuildFacetList(f);
		double nU = Norme(&(f->sh.U));
		f->tRatio = f->sh.texWidthD / nU;
	}

}


void Geometry::LoadSYN(FileReader *file,GLProgress *prg,LEAK *pleak,int *nbleak,HIT *pHits,int *nbHHit,int *version) {


	MolFlow *mApp = (MolFlow *)theApp;
	mApp->ClearAllSelections();
	mApp->ClearAllViews();
	prg->SetMessage("Clearing current geometry...");
	Clear();
	mApp->ClearFormula();

	// Globals
	char tmp[512];
	prg->SetMessage("Reading SYN file header...");
	file->ReadKeyword("version");file->ReadKeyword(":");
	*version = file->ReadInt();
	if( *version>SYNVERSION ) {
		char errMsg[512];
		sprintf(errMsg,"Unsupported SYN version V%d",*version);
		throw Error(errMsg);
	}

	file->ReadKeyword("totalHit");file->ReadKeyword(":");
	tNbHit = 0;file->ReadLLong();
	file->ReadKeyword("totalDes");file->ReadKeyword(":");
	tNbDesorption = 0;file->ReadLLong();
	file->ReadKeyword("totalLeak");file->ReadKeyword(":");
	tNbLeak = 0;file->ReadLLong();
	if (*version>2) {
		file->ReadKeyword("totalFlux");file->ReadKeyword(":");
		file->ReadDouble();
		file->ReadKeyword("totalPower");file->ReadKeyword(":");
		file->ReadDouble();
	}
	file->ReadKeyword("maxDes");file->ReadKeyword(":");
	tNbDesorptionMax = 0;file->ReadLLong(); 
	file->ReadKeyword("nbVertex");file->ReadKeyword(":");
	sh.nbVertex = file->ReadInt();
	file->ReadKeyword("nbFacet");file->ReadKeyword(":");
	sh.nbFacet = file->ReadInt();
	file->ReadKeyword("nbSuper");file->ReadKeyword(":");
	sh.nbSuper = file->ReadInt();
	int nbF=0;
	int nbV=0;

	file->ReadKeyword("nbFormula");file->ReadKeyword(":");
	nbF = file->ReadInt();
	file->ReadKeyword("nbView");file->ReadKeyword(":");
	nbV = file->ReadInt();
	int nbS=0;
	file->ReadKeyword("nbSelection");file->ReadKeyword(":");
	nbS = file->ReadInt();

	if (*version>1) {
		file->ReadKeyword("nbRegions");file->ReadKeyword(":");
		int nbR = file->ReadInt();
		//result=PARfileList(nbR);

		file->ReadKeyword("PARfiles");file->ReadKeyword("{");
		for (int i=0;i<nbR;i++) {
			/*char tmp[512];
			strcpy(tmp,file->ReadString());
			result.fileNames[i]=_strdup(tmp);*/
			file->ReadString();
		}
		file->ReadKeyword("}");
	}

	file->ReadKeyword("formulas");file->ReadKeyword("{");
	for(int i=0;i<nbF;i++) {
		char tmpName[256];
		char tmpExpr[512];
		strcpy(tmpName,file->ReadString());
		strcpy(tmpExpr,file->ReadString());
		mApp->AddFormula(tmpName,tmpExpr);
	}
	file->ReadKeyword("}");

	file->ReadKeyword("views");file->ReadKeyword("{");
	for(int i=0;i<nbV;i++) {
		char tmpName[256];
		AVIEW v;
		strcpy(tmpName,file->ReadString());
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
		mApp->AddView(tmpName,v);
	}
	file->ReadKeyword("}");

	file->ReadKeyword("selections");file->ReadKeyword("{");
	for(int i=0;i<nbS;i++) {
		ASELECTION s;
		char tmpName[256];
		strcpy(tmpName,file->ReadString());
		s.name=_strdup(tmpName);
		s.nbSel = file->ReadInt();
		s.selection = (int *)malloc((s.nbSel)*sizeof(int));

		for (int j=0;j<s.nbSel;j++) {
			s.selection[j]=file->ReadInt();
		}
		mApp->AddSelection(s.name,s);
	}
	file->ReadKeyword("}");

	file->ReadKeyword("structures");file->ReadKeyword("{");
	for(int i=0;i<sh.nbSuper;i++) {
		strName[i] = _strdup(file->ReadString());
		// For backward compatibilty with STR
		sprintf(tmp,"%s.txt",strName[i]);
		strFileName[i] = _strdup(tmp);
	}
	file->ReadKeyword("}");

	// Allocate memory
	facets = (Facet **)malloc(sh.nbFacet * sizeof(Facet *));
	memset(facets,0,sh.nbFacet * sizeof(Facet *));
	vertices3 = (VERTEX3D *)malloc(sh.nbVertex * sizeof(VERTEX3D));

	// Read vertices
	prg->SetMessage("Reading vertices...");
	file->ReadKeyword("vertices");file->ReadKeyword("{");
	for(int i=0;i<sh.nbVertex;i++) {
		// Check idx
		int idx = file->ReadInt();
		if( idx != i+1 ) throw Error(file->MakeError("Wrong vertex index !"));
		vertices3[i].x = file->ReadDouble();
		vertices3[i].y = file->ReadDouble();
		vertices3[i].z = file->ReadDouble();
	}
	file->ReadKeyword("}");
	prg->SetMessage("Reading leaks and hits...");
	// Read leaks
	file->ReadKeyword("leaks");file->ReadKeyword("{");
	file->ReadKeyword("nbLeak");file->ReadKeyword(":");
	*nbleak = 0;
	int nbleak_local=file->ReadInt();
	for (int i=0;i<nbleak_local;i++) {
		int idx = file->ReadInt();
		if( idx != i ) throw Error(file->MakeError("Wrong leak index !"));
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
	file->ReadKeyword("hits");file->ReadKeyword("{");
	file->ReadKeyword("nbHHit");file->ReadKeyword(":");
	*nbHHit = 0;
	int nbHHit_local=file->ReadInt();
	for (int i=0;i<nbHHit_local;i++) {
		int idx = file->ReadInt();
		if( idx != i ) throw Error(file->MakeError("Wrong hit cache index !"));
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
	for(int i=0;i<sh.nbFacet;i++) {
		file->ReadKeyword("facet");
		// Check idx
		int idx = file->ReadInt();
		if( idx != i+1 ) throw Error(file->MakeError("Wrong facet index !"));
		file->ReadKeyword("{");
		file->ReadKeyword("nbIndex");
		file->ReadKeyword(":");
		int nbI = file->ReadInt();
		if( nbI<3 ) {
			char errMsg[512];
			sprintf(errMsg,"Facet %d has only %d vertices. ",i,nbI);
			throw Error(errMsg);
		}
		prg->SetProgress((float)i/sh.nbFacet);
		facets[i] = new Facet(nbI);
		facets[i]->LoadSYN(file,*version,sh.nbVertex);
		file->ReadKeyword("}");
	}

	prg->SetMessage("Initalizing geometry and building mesh...");
	InitializeGeometry();
	AdjustProfile();
	isLoaded = TRUE;
	UpdateName(file);

	// Update mesh
	prg->SetMessage("Drawing textures...");
	for(int i=0;i<sh.nbFacet;i++) {
		double p = (double)i/(double)sh.nbFacet;
		prg->SetProgress(p);
		Facet *f = facets[i];
		//f->SetTexture(f->sh.texWidthD,f->sh.texHeightD,f->hasMesh);
		BuildFacetList(f);
		//double nU = Norme(&(f->sh.U));
		//f->tRatio = f->sh.texWidthD / nU;
	}
	//return result;
}


// -----------------------------------------------------------
bool Geometry::LoadTextures(FileReader *file,GLProgress *prg,Dataport *dpHit,int version) {
	MolFlow *mApp = (MolFlow *)theApp;
	if(file->SeekFor("{textures}")) {
		char tmp[256];
		//versions 3+
		// Block dpHit during the whole disc reading

		AccessDataport(dpHit);
		AHIT readVal;

		// Globals
		BYTE *buffer = (BYTE *)dpHit->buff;
		SHGHITS *gHits = (SHGHITS *)buffer;

		gHits->total.hit.nbHit=tNbHit;
		gHits->total.hit.nbDesorbed=tNbDesorption;
		gHits->total.hit.nbAbsorbed=tNbAbsorption;
		gHits->nbLeakTotal = tNbLeak;
		gHits->distTraveledTotal = distTraveledTotal;

		// Read facets
		file->ReadKeyword("minHit");file->ReadKeyword(":");
		gHits->minHit = (AHIT)file->ReadDouble();
		file->ReadKeyword("maxHit");file->ReadKeyword(":");
		gHits->maxHit = (AHIT)file->ReadDouble();

		for (size_t m=0;m<=mApp->worker.moments.size() || (m==0 && version<10);m++){
			if (version>=10) {
				file->ReadKeyword("moment");
				if(m!=file->ReadInt()) {
					throw Error("Unexpected profile moment");
					break;
				}
			file->ReadKeyword("{");
			}

			for(int i=0;i<sh.nbFacet;i++) {
				Facet *f = facets[i];
				if (f->hasMesh || version<8) {
					prg->SetProgress((double)(i+m*sh.nbFacet)/(double)(mApp->worker.moments.size()*sh.nbFacet)*0.33+0.66);
					file->ReadKeyword("texture_facet");
					// Check idx
					int idx = file->ReadInt();

					if( idx != i+1 ) {
						sprintf(tmp,"Wrong facet index. Expected %d, read %d.",i+1,idx);
						throw Error(file->MakeError(tmp));
					}
					file->ReadKeyword("{");

					int ix,iy;

					///Load textures, for GEO file version 3+

					int profSize = (f->sh.isProfile)?((1+(int)mApp->worker.moments.size())*(PROFILE_SIZE*sizeof(llong))):0;
					int h=(f->sh.texHeight);
					int w=(f->sh.texWidth);
					AHIT *hits = (AHIT *)((BYTE *)gHits + (f->sh.hitOffset + sizeof(SHHITS) + profSize + m*w*h*sizeof(AHIT)));

					for (iy=0;iy<h;iy++) {
						for (ix=0;ix<w;ix++) {
							readVal=(AHIT)file->ReadDouble();
							*(hits+(iy*(f->sh.texWidth)+ix))=readVal;
						}
					}
					file->ReadKeyword("}");
				}
			}
			if (version>=10) file->ReadKeyword("}");
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

// -----------------------------------------------------------

void Geometry::SaveGEO(FileWriter *file,GLProgress *prg,Dataport *dpHit,std::vector<std::string> userMoments,Worker *worker,
					   BOOL saveSelected,LEAK *pleak,int *nbleakSave,HIT *pHits,int *nbHHitSave,BOOL crashSave) {

	MolFlow *mApp = (MolFlow *)theApp;
	prg->SetMessage("Counting hits...");
	if(!IsLoaded()) throw Error("Nothing to save !");


	// Block dpHit during the whole disc writing
	AccessDataport(dpHit);

	// Globals
	BYTE *buffer = (BYTE *)dpHit->buff;
	SHGHITS *gHits = (SHGHITS *)buffer;


	float dCoef = 1.0f;
	int ix,iy;

	switch(gHits->mode) {

	case MC_MODE:
		if( gHits->total.hit.nbDesorbed>0 ) {
			dCoef = (float)totalOutgassing / (float)gHits->total.hit.nbDesorbed;
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

	}

	prg->SetMessage("Writing geometry details...");
	file->Write("version:");file->WriteInt(GEOVERSION,"\n");
	file->Write("totalHit:");file->WriteLLong(gHits->total.hit.nbHit,"\n");
	file->Write("totalDes:");file->WriteLLong(gHits->total.hit.nbDesorbed,"\n");
	file->Write("totalLeak:");file->WriteLLong(gHits->nbLeakTotal,"\n");
	file->Write("totalAbs:");file->WriteLLong(gHits->total.hit.nbAbsorbed,"\n");
	file->Write("totalDist:");file->WriteDouble(gHits->distTraveledTotal,"\n");
	file->Write("maxDes:");file->WriteLLong(tNbDesorptionMax,"\n");
	file->Write("nbVertex:");file->WriteInt(sh.nbVertex,"\n");
	file->Write("nbFacet:");file->WriteInt(saveSelected?nbSelected:sh.nbFacet,"\n");
	file->Write("nbSuper:");file->WriteInt(sh.nbSuper,"\n");
	file->Write("nbFormula:");file->WriteInt(mApp->nbFormula,"\n");
	file->Write("nbView:");file->WriteInt(mApp->nbView,"\n");
	file->Write("nbSelection:");file->WriteInt(mApp->nbSelection,"\n");
	file->Write("gasMass:");file->WriteDouble(gasMass,"\n");  

	file->Write("userMoments {\n");
	file->Write(" nb:");file->WriteInt((int)userMoments.size());
	for(size_t u=0;u<userMoments.size();u++) {
		file->Write("\n \"");
		file->Write(userMoments[u].c_str());
		file->Write("\"");
	}
	file->Write("\n}\n");
	
	file->Write("desorptionStart:");file->WriteDouble(worker->desorptionStartTime,"\n");
	file->Write("desorptionStop:");file->WriteDouble(worker->desorptionStopTime,"\n");
	file->Write("timeWindow:");file->WriteDouble(worker->timeWindowSize,"\n");
	file->Write("useMaxwellian:");file->WriteInt(worker->useMaxwellDistribution,"\n");
	file->Write("calcConstantFlow:");file->WriteInt(worker->calcConstantFlow,"\n");

	file->Write("formulas {\n");
	for(int i=0;i<mApp->nbFormula;i++) {
		file->Write("  \"");
		file->Write(mApp->formulas[i].parser->GetName());
		file->Write("\" \"");
		file->Write(mApp->formulas[i].parser->GetExpression());
		file->Write("\"\n");
	}
	file->Write("}\n");

	file->Write("views {\n");
	for(int i=0;i<mApp->nbView;i++) {
		file->Write("  \"");
		file->Write(mApp->views[i].name);
		file->Write("\"\n");
		file->WriteInt(mApp->views[i].projMode," ");
		file->WriteDouble(mApp->views[i].camAngleOx," ");
		file->WriteDouble(mApp->views[i].camAngleOy," ");
		file->WriteDouble(mApp->views[i].camDist," ");
		file->WriteDouble(mApp->views[i].camOffset.x," ");
		file->WriteDouble(mApp->views[i].camOffset.y," ");
		file->WriteDouble(mApp->views[i].camOffset.z," ");
		file->WriteInt(mApp->views[i].performXY," ");
		file->WriteDouble(mApp->views[i].vLeft," ");
		file->WriteDouble(mApp->views[i].vRight," ");
		file->WriteDouble(mApp->views[i].vTop," ");
		file->WriteDouble(mApp->views[i].vBottom,"\n");
	}
	file->Write("}\n");

	file->Write("selections {\n");
	for(int i=0;i<mApp->nbSelection;i++) {
		file->Write("  \"");
		file->Write(mApp->selections[i].name);
		file->Write("\"\n ");
		file->WriteInt(mApp->selections[i].nbSel,"\n");
		for (int j=0;j<mApp->selections[i].nbSel;j++) {
			file->Write("  ");
			file->WriteInt(mApp->selections[i].selection[j],"\n");
		}
		//file->Write("\n");
	}
	file->Write("}\n");

	file->Write("structures {\n");
	for(int i=0;i<sh.nbSuper;i++) {
		file->Write("  \"");
		file->Write(strName[i]);
		file->Write("\"\n");
	}
	file->Write("}\n");
	//vertices
	prg->SetMessage("Writing vertices...");
	file->Write("vertices {\n");
	for(int i=0;i<sh.nbVertex;i++) {
		prg->SetProgress(0.33*((double)i/(double)sh.nbVertex));
		file->Write("  ");
		file->WriteInt(i+1," ");
		file->WriteDouble(vertices3[i].x," ");
		file->WriteDouble(vertices3[i].y," ");
		file->WriteDouble(vertices3[i].z,"\n");
	}
	file->Write("}\n");

	//leaks
	prg->SetMessage("Writing leaks...");
	file->Write("leaks {\n");
	file->Write("  nbLeak:");file->WriteInt(*nbleakSave,"\n");
	for (int i=0;i<*nbleakSave;i++) {

		file->Write("  ");
		file->WriteInt(i," ");
		file->WriteDouble((pleak+i)->pos.x," ");
		file->WriteDouble((pleak+i)->pos.y," ");
		file->WriteDouble((pleak+i)->pos.z," ");

		file->WriteDouble((pleak+i)->dir.x," ");
		file->WriteDouble((pleak+i)->dir.y," ");
		file->WriteDouble((pleak+i)->dir.z,"\n");    
	}
	file->Write("}\n");

	//hit cache (lines and dots)
	prg->SetMessage("Writing hit cache...");
	file->Write("hits {\n");
	file->Write("  nbHHit:");file->WriteInt(*nbHHitSave,"\n");
	for (int i=0;i<*nbHHitSave;i++) {

		file->Write("  ");
		file->WriteInt(i," ");
		file->WriteDouble((pHits+i)->pos.x," ");
		file->WriteDouble((pHits+i)->pos.y," ");
		file->WriteDouble((pHits+i)->pos.z," ");

		file->WriteInt((pHits+i)->type,"\n");
	}
	file->Write("}\n");

	//facets

	prg->SetMessage("Writing facets..."); 

	for(int i=0,k=0;i<sh.nbFacet;i++) {
		prg->SetProgress(0.33+((double)i/(double)sh.nbFacet) *0.33);
		if( saveSelected ) {
			if( facets[i]->selected ) { facets[i]->SaveGEO(file,k);k++; }
		} else {
			facets[i]->SaveGEO(file,i);
		}
	}

	prg->SetMessage("Writing profiles...");
	SaveProfile(file,dpHit,-1,saveSelected);

	///Save textures, for GEO file version 3
	char tmp[256];
	sprintf(tmp,"{textures}\n");
	file->Write(tmp);
	file->Write("minHit:");file->WriteDouble(gHits->minHit,"\n");
	file->Write("maxHit:");file->WriteDouble(gHits->maxHit,"\n");

	//Selections
	//SaveSelections();

	prg->SetMessage("Writing textures..."); 
	for (size_t m=0;m<=mApp->worker.moments.size();m++){
		sprintf(tmp,"moment %d {\n",m);
		file->Write(tmp);
		for(int i=0;i<sh.nbFacet;i++) {
			prg->SetProgress((double)(i+m*sh.nbFacet)/(double)(mApp->worker.moments.size()*sh.nbFacet)*0.33+0.66);
			Facet *f = facets[i];
			if (f->hasMesh) {
				int h=(f->sh.texHeight);
				int w=(f->sh.texWidth);
				int profSize = (f->sh.isProfile)?(PROFILE_SIZE*sizeof(llong)*(1+(int)mApp->worker.moments.size())):0;
				AHIT *hits = (AHIT *)((BYTE *)gHits + (f->sh.hitOffset + sizeof(SHHITS) + profSize + m*w*h*sizeof(AHIT)));

				//char tmp[256];
				sprintf(tmp," texture_facet %d {\n",i+1);
				file->Write(tmp);

				for (iy=0;iy<h;iy++) {
					for (ix=0;ix<w;ix++) {
						file->WriteDouble(*(hits+(iy*(f->sh.texWidth)+ix)),"\t");
					}
					file->Write("\n");
				}
				file->Write(" }\n");
			}
		}
		file->Write("}\n");
	}

	ReleaseDataport(dpHit);

	//Debug memory check
	 //_ASSERTE (!_CrtDumpMemoryLeaks());;
	 _ASSERTE(_CrtCheckMemory());
}

// -----------------------------------------------------------

void Geometry::SaveTXT(FileWriter *file,Dataport *dpHit,BOOL saveSelected) {

	if(!IsLoaded()) throw Error("Nothing to save !");

	// Unused
	file->WriteInt(0,"\n");

	// Block dpHit during the whole disc writing
	AccessDataport(dpHit);

	// Globals
	BYTE *buffer = (BYTE *)dpHit->buff;
	SHGHITS *gHits = (SHGHITS *)buffer;

	// Unused
	file->WriteLLong(gHits->total.hit.nbHit,"\n");
	file->WriteLLong(gHits->nbLeakTotal,"\n");
	file->WriteLLong(gHits->total.hit.nbDesorbed,"\n");
	file->WriteLLong(tNbDesorptionMax,"\n");

	file->WriteInt(sh.nbVertex,"\n");
	file->WriteInt(saveSelected?nbSelected:sh.nbFacet,"\n");

	// Read geometry vertices
	for(int i=0;i<sh.nbVertex;i++) {
		file->WriteDouble(vertices3[i].x," ");
		file->WriteDouble(vertices3[i].y," ");
		file->WriteDouble(vertices3[i].z,"\n");
	}

	// Facets
	for(int i=0;i<sh.nbFacet;i++) {
		Facet *f = facets[i];
		int j;
		if( saveSelected ) {
			if( f->selected ) {
				file->WriteInt(f->sh.nbIndex," ");
				for(j=0;j<f->sh.nbIndex-1;j++)
					file->WriteInt(f->indices[j]+1," ");
				file->WriteInt(f->indices[j]+1,"\n");
			}
		} else {
			file->WriteInt(f->sh.nbIndex," ");
			for(j=0;j<f->sh.nbIndex-1;j++)
				file->WriteInt(f->indices[j]+1," ");
			file->WriteInt(f->indices[j]+1,"\n");
		}
	}

	// Params
	for(int i=0;i<sh.nbFacet;i++) {

		// Update facet hits from shared mem
		Facet *f = facets[i];
		SHHITS *shF = (SHHITS *)(buffer + f->sh.hitOffset);
		memcpy(&(f->sh.counter),shF,sizeof(SHHITS));
		if( saveSelected ) {
			if(f->selected) f->SaveTXT(file);
		} else {
			f->SaveTXT(file);
		}

	}

	SaveProfile(file,dpHit,-1,saveSelected,TRUE);

	ReleaseDataport(dpHit);

}

// -----------------------------------------------------------------------
void Geometry::ExportTextures(FILE *file,int mode,Dataport *dpHit,BOOL saveSelected) {
	MolFlow *mApp = (MolFlow *)theApp;
	//if(!IsLoaded()) throw Error("Nothing to save !");

	// Block dpHit during the whole disc writing
	BYTE *buffer = NULL;
	if( dpHit )
		if( AccessDataport(dpHit) ) 
			buffer = (BYTE *)dpHit->buff;

	// Globals
	//BYTE *buffer = (BYTE *)dpHit->buff;
	//SHGHITS *gHits = (SHGHITS *)buffer;

	for (size_t m=0;m<=mApp->worker.moments.size();m++){
		if (m==0) fprintf(file," moment 0 (Constant Flow){\n");
		else fprintf(file," moment %d (%g s){\n",m,mApp->worker.moments[m-1]);
		// Facets

		for(int i=0;i<sh.nbFacet;i++) {
			Facet *f = facets[i];
			

			if( f->selected ) {
				AHIT *hits=NULL;
				VHIT *dirs=NULL;
				fprintf(file,"FACET%d\n",i+1);
				
				if( f->mesh || f->sh.countDirection ) {
					char tmp[256];
					float dCoef = 1.0f;
					if(!buffer) return;
					SHGHITS *shGHit = (SHGHITS *)buffer;
					int nbMoments = (int)mApp->worker.moments.size();
					int profSize = (f->sh.isProfile)?(PROFILE_SIZE*sizeof(llong)*(1+nbMoments)):0;
					int w = f->sh.texWidth;
					int h = f->sh.texHeight;
					int tSize=w*h*sizeof(AHIT);
					int dSize=w*h*sizeof(VHIT);
					if (f->mesh) hits = (AHIT *)((BYTE *)buffer + (f->sh.hitOffset + sizeof(SHHITS) + profSize + m*tSize));
					if (f->sh.countDirection) dirs = (VHIT *)((BYTE *)buffer + (f->sh.hitOffset + sizeof(SHHITS) + profSize*(1+nbMoments) + tSize*(1+nbMoments) + m*dSize));
				

					switch(mode) {

					case 0: // Element area
						for(int j=0;j<h;j++) {
							for(int i=0;i<w;i++) {						
								sprintf(tmp,"%g",f->mesh[i+j*w].area);
								if( tmp ) fprintf(file,"%s",tmp);
								if( j<w-1 ) 
									fprintf(file,"\t");
							}
							fprintf(file,"\n");
						}
						break;


					case 1:  // Pressure

						// Lock during update					
						dCoef = 1.0f;
						if( shGHit->mode == MC_MODE ) dCoef = (float)totalOutgassing *((m==0)?1.0f:(float)((mApp->worker.desorptionStopTime-mApp->worker.desorptionStartTime)
							/mApp->worker.timeWindowSize))/ (float)shGHit->total.hit.nbDesorbed;
						for(int j=0;j<h;j++) {
							for(int i=0;i<w;i++) {						
								sprintf(tmp,"%g",hits[i+j*w]*dCoef);
								fprintf(file,"%s",tmp);
								fprintf(file,"\t");
							}
							fprintf(file,"\n");						
						}
						break;


					case 2: // Average velocity
						for(int j=0;j<h;j++) {
							for(int i=0;i<w;i++) {						
								if( f->sh.countDirection ) {
									double n = dirs[i+j*w].sumSpeed/(double)dirs[i+j*w].count;
									sprintf(tmp,"%g",n);
								} else {
									sprintf(tmp,"None");
								}
								fprintf(file,"%s",tmp);
								fprintf(file,"\t");
							}
							fprintf(file,"\n");
						}
						break;

					case 3: // Velocity
						for(int j=0;j<h;j++) {
							for(int i=0;i<w;i++) {						
								if( f->sh.countDirection ) {
									sprintf(tmp,"%g,%g,%g",
										dirs[i+j*w].dir.x,
										dirs[i+j*w].dir.y,
										dirs[i+j*w].dir.z);
								} else {
									sprintf(tmp,"None");
								}
								fprintf(file,"%s",tmp);
								fprintf(file,"\t");
							}
							fprintf(file,"\n");
						}
						break;

					case 4: // Velocity vector Count
						for(int j=0;j<h;j++) {
							for(int i=0;i<w;i++) {
								if( f->sh.countDirection ) {
									sprintf(tmp,"%I64d",dirs[i+j*w].count);
								} else {
									sprintf(tmp,"None");
								}
								fprintf(file,"%s",tmp);
								fprintf(file,"\t");
							}
							fprintf(file,"\n");
						}
						break;

					case 5:  // Hits

						// Lock during update					
						//if( shGHit->mode == MC_MODE ) dCoef = 1.0f / (float)shGHit->total.hit.nbDesorbed;
						for(int j=0;j<h;j++) {
							for(int i=0;i<w;i++) {					
								float scale = (float)(40.0/(sqrt(8*8.31*f->sh.temperature/(PI*(gasMass*1e-3))))); //was divided by 11.7 //0.274974=1/sqrt(8*R/PI/1000)/40
								if(f->sh.opacity>0.0 && f->sh.opacity!=1.0) scale = scale/(float)f->sh.opacity;
								if(f->sh.is2sided) scale = scale / 2.0f;
								float val=(float)((hits[i+j*w]/(float)scale*f->mesh[i+j*w].area));
								sprintf(tmp,"%g",val);
								fprintf(file,"%s",tmp);
								fprintf(file,"\t");
							}
							fprintf(file,"\n");						
						}
						break;
					}
				}
				else {
					fprintf(file,"No mesh.\n");
				}
				fprintf(file,"\n"); //Current facet exported. 
			}
		}
		fprintf(file," }\n");
	}
	ReleaseDataport(dpHit);
}

void Geometry::ImportDesorption(FileReader *file,Dataport *dpHit) {

	//if(!IsLoaded()) throw Error("Nothing to save !");

	// Block dpHit during the whole disc writing

	for (int i=0;i<sh.nbFacet;i++) { //clear previous desorption maps
		facets[i]->hasOutgassingMap = FALSE;
		facets[i]->sh.useOutgassingFile = FALSE;
		facets[i]->sh.desorbType = DES_NONE; //clear previously set desorptions
		facets[i]->selected = FALSE;
		facets[i]->UnselectElem();
	}
	// Facets
	while (!file->IsEof()) {
		file->ReadKeyword("facet");
		int facetId=file->ReadInt()-1;
		file->ReadKeyword("{");
		if (!(facetId>=0 && facetId<sh.nbFacet)) {
			file->MakeError("Invalid facet Id (loaded desorption file for a different geometry?)");
			return;
		}

		Facet *f = facets[facetId];
		f->hasOutgassingMap=TRUE;
		f->sh.useOutgassingFile=TRUE; //turn on file usage by default
		f->sh.desorbType=DES_COSINE; //auto-set to cosine
		Select(f);
		file->ReadKeyword("cell_size_cm");file->ReadKeyword(":");
		double ratio=f->sh.outgassingFileRatio=file->ReadDouble();
		if (f->sh.outgassingFileRatio!=0.0) {
			f->sh.outgassingFileRatio=1.0/f->sh.outgassingFileRatio; //cell size -> samples per cm
			ratio=f->sh.outgassingFileRatio;
		}
		double nU = Norme(&(f->sh.U));
		double nV = Norme(&(f->sh.V));
		int w = f->sh.outgassingMapWidth = (int)ceil(nU*ratio*0.999999); //double precision written to file
		int h = f->sh.outgassingMapHeight= (int)ceil(nV*ratio*0.999999); //double precision written to file
		f->outgassingMap=(double*)malloc(w*h*sizeof(double));
		for(int i=0;i<w;i++) {
			for(int j=0;j<h;j++) {
				f->outgassingMap[i+j*w]=file->ReadDouble();
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

// -----------------------------------------------------------

void Geometry::SaveSTR(Dataport *dpHit,BOOL saveSelected) {

	if(!IsLoaded()) throw Error("Nothing to save !");
	if(sh.nbSuper<1) throw Error("Cannot save single structure in STR format");

	// Block dpHit during the whole disc writting
	AccessDataport(dpHit);
	for(int i=0;i<sh.nbSuper;i++)
		SaveSuper(dpHit,i);
	ReleaseDataport(dpHit);

}


// -----------------------------------------------------------

void Geometry::SaveSuper(Dataport *dpHit,int s) {

	char fName[512];
	sprintf(fName,"%s/%s",strPath,strFileName[s]);
	FileWriter *file = new FileWriter(fName);

	// Unused
	file->WriteInt(0,"\n");

	// Globals
	BYTE *buffer = (BYTE *)dpHit->buff;
	SHGHITS *gHits = (SHGHITS *)buffer;

	//Extract data of the specified super structure
	llong totHit  = 0;
	llong totAbs  = 0;
	llong totDes  = 0;
	int *refIdx = (int *)malloc(sh.nbVertex*sizeof(int));
	memset(refIdx,0xFF,sh.nbVertex*sizeof(int));
	int nbV=0;
	int nbF=0;

	for(int i=0;i<sh.nbFacet;i++) {
		Facet *f = facets[i];
		if( f->sh.superIdx==s ) {
			totHit += f->sh.counter.hit.nbHit;
			totAbs += f->sh.counter.hit.nbAbsorbed;
			totDes += f->sh.counter.hit.nbDesorbed;
			for(int j=0;j<f->sh.nbIndex;j++)
				refIdx[f->indices[j]] = 1;
			nbF++;
		}
	}

	for(int i=0;i<sh.nbVertex;i++) {
		if(refIdx[i]>=0) {
			refIdx[i] = nbV;
			nbV++;
		}
	}

	file->WriteLLong(totHit,"\n");
	file->WriteLLong(gHits->nbLeakTotal,"\n");
	file->WriteLLong(totDes,"\n");
	file->WriteLLong(tNbDesorptionMax,"\n");

	file->WriteInt(nbV,"\n");
	file->WriteInt(nbF,"\n");

	// Read geometry vertices
	for(int i=0;i<sh.nbVertex;i++) {
		if( refIdx[i]>=0 ) {
			file->WriteDouble(vertices3[i].x," ");
			file->WriteDouble(vertices3[i].y," ");
			file->WriteDouble(vertices3[i].z,"\n");
		}
	}

	// Facets
	for(int i=0;i<sh.nbFacet;i++) {
		Facet *f = facets[i];
		int j;
		if( f->sh.superIdx==s ) {
			file->WriteInt(f->sh.nbIndex," ");
			for(j=0;j<f->sh.nbIndex-1;j++)
				file->WriteInt(refIdx[f->indices[j]]+1," ");
			file->WriteInt(refIdx[f->indices[j]]+1,"\n");
		}
	}

	// Params
	for(int i=0;i<sh.nbFacet;i++) {

		// Update facet hits from shared mem
		Facet *f = facets[i];
		if( f->sh.superIdx==s ) {
			SHHITS *shF = (SHHITS *)(buffer + f->sh.hitOffset);
			memcpy(&(f->sh.counter),shF,sizeof(SHHITS));
			f->SaveTXT(file);
		}

	}

	SaveProfile(file,dpHit,s);

	SAFE_DELETE(file);
	free(refIdx);

}


// -----------------------------------------------------------

void Geometry::RemoveLinkFacet() { //unused

	// Remove facet used as link for superstructure (not needed)
	int nb = 0;
	for(int i=0;i<sh.nbFacet;i++)
		if(facets[i]->IsLinkFacet()) nb++;

	if(nb==0) return;

	Facet **f = (Facet **)malloc((sh.nbFacet-nb) * sizeof(Facet *));

	nb=0;
	for(int i=0;i<sh.nbFacet;i++) {
		if(facets[i]->IsLinkFacet()) {
			delete facets[i];
		} else {
			f[nb++] = facets[i];
		}
	}

	SAFE_FREE(facets);
	facets = f;
	sh.nbFacet = nb;

}

// -----------------------------------------------------------

int Geometry::HasIsolatedVertices() {

	// Check if there are unused vertices
	int *check = (int *)malloc(sh.nbVertex*sizeof(int));
	memset(check,0,sh.nbVertex*sizeof(int));

	for(int i=0;i<sh.nbFacet;i++) {
		Facet *f = facets[i];
		for(int j=0;j<f->sh.nbIndex;j++) {
			check[f->indices[j]]++;
		}
	}

	int nbUnused = 0;
	for(int i=0;i<sh.nbVertex;i++) {
		if(!check[i]) nbUnused++;
	}

	SAFE_FREE(check);
	return nbUnused;

}

// -----------------------------------------------------------

void  Geometry::DeleteIsolatedVertices(BOOL selectedOnly) {
	changedSinceSave = TRUE;
	// Remove unused vertices
	int *check = (int *)malloc(sh.nbVertex*sizeof(int));
	memset(check,0,sh.nbVertex*sizeof(int));

	for(int i=0;i<sh.nbFacet;i++) {
		Facet *f = facets[i];
		for(int j=0;j<f->sh.nbIndex;j++) {
			check[f->indices[j]]++;
		}
	}

	int nbUnused = 0;
	for(int i=0;i<sh.nbVertex;i++) {
		if(!check[i]&&!(selectedOnly&&!vertices3[i].selected)) nbUnused++;
	}

	int nbVert = sh.nbVertex - nbUnused;


	if(nbVert==0) {
		// Remove all
		SAFE_FREE(check);
		Clear();
		return;
	}


	// Update facet indices
	int *newId = (int *)malloc(sh.nbVertex*sizeof(int));
	memset(newId,0,sh.nbVertex*sizeof(int));
	for(int i=0,n=0;i<sh.nbVertex;i++) {
		if(check[i] || (selectedOnly && !vertices3[i].selected)) {
			newId[i]=n;
			n++;
		}
	}
	for(int i=0;i<sh.nbFacet;i++) {
		Facet *f = facets[i];
		for(int j=0;j<f->sh.nbIndex;j++) {
			f->indices[j] = newId[f->indices[j]];
		}
	}



	VERTEX3D *nVert = (VERTEX3D *)malloc(nbVert*sizeof(VERTEX3D));
	_ASSERTE(nVert);
	for(int i=0,n=0;i<sh.nbVertex;i++) {
		if(check[i] || (selectedOnly && !vertices3[i].selected)) {
			nVert[n]=vertices3[i];
			n++;
		}
	}

	SAFE_FREE(vertices3);
	vertices3 = nVert;
	sh.nbVertex = nbVert;

	SAFE_FREE(check);
	SAFE_FREE(newId);

	// Delete old resources
	//DeleteGLLists(TRUE,TRUE);

	//InitializeGeometry();

}

// -----------------------------------------------------------

void  Geometry::SelectIsolatedVertices() {

	UnselectAllVertex();	
	// Select unused vertices
	int *check = (int *)malloc(sh.nbVertex*sizeof(int));
	memset(check,0,sh.nbVertex*sizeof(int));

	for(int i=0;i<sh.nbFacet;i++) {
		Facet *f = facets[i];
		for(int j=0;j<f->sh.nbIndex;j++) {
			check[f->indices[j]]++;
		}
	}

	for(int i=0;i<sh.nbVertex;i++) {
		if(!check[i]) vertices3[i].selected=TRUE;
	}

	SAFE_FREE(check);
}


// -----------------------------------------------------------

void Geometry::RemoveCollinear() {
	MolFlow *mApp = (MolFlow *)theApp;
	changedSinceSave=TRUE;
	int nb = 0;
	for(int i=0;i<sh.nbFacet;i++)
		if(facets[i]->collinear) nb++;

	if(nb==0) return;
	/*
	if(sh.nbFacet-nb==0) {
	// Remove all
	Clear();
	return;
	}
	*/

	Facet   **f = (Facet **)malloc((sh.nbFacet-nb) * sizeof(Facet *));

	nb=0;
	for(int i=0;i<sh.nbFacet;i++) {
		if(facets[i]->collinear) {
			delete facets[i];
			mApp->RenumberSelections(i);
			mApp->RenumberFormulas(i);
		} else {
			f[nb++] = facets[i];
		}
	}

	SAFE_FREE(facets);
	facets = f;
	sh.nbFacet = nb;

	// Delete old resources
	DeleteGLLists(TRUE,TRUE);

	BuildGLList();
	mApp->UpdateModelParams();

}

// -----------------------------------------------------------

void Geometry::RemoveFromStruct(int numToDel) {
	MolFlow *mApp = (MolFlow *)theApp;
	changedSinceSave=TRUE;
	int nb = 0;
	for(int i=0;i<sh.nbFacet;i++)
		if(facets[i]->sh.superIdx==numToDel) nb++;

	if(nb==0) return;
	/*
	if(sh.nbFacet-nb==0) {
	// Remove all
	Clear();
	return;
	}
	*/

	Facet   **f = (Facet **)malloc((sh.nbFacet-nb) * sizeof(Facet *));

	nb=0;
	for(int i=0;i<sh.nbFacet;i++) {
		if(facets[i]->sh.superIdx==numToDel) {
			delete facets[i];
			mApp->RenumberSelections(i);
			mApp->RenumberFormulas(i);
		} else {
			f[nb++] = facets[i];
		}
	}

	SAFE_FREE(facets);
	facets = f;
	sh.nbFacet = nb;
	CalcTotalOutGassing();

}

// -----------------------------------------------------------

void Geometry::RemoveSelectedVertex() {
	MolFlow *mApp = (MolFlow *)theApp;
	changedSinceSave=TRUE;
	
	//Analyze facets
	std::vector<int> facetsToRemove,facetsToChange;
	for (int f=0;f<sh.nbFacet;f++) {
		int nbSelVertex=0;
		for (int i=0;i<facets[f]->sh.nbIndex;i++)
			if (vertices3[facets[f]->indices[i]].selected)
				nbSelVertex++;
		if (nbSelVertex) {
			facetsToChange.push_back(f);
			if ((facets[f]->sh.nbIndex-nbSelVertex)<=2)
				facetsToRemove.push_back(f);
		}
	}

	for (size_t c=0;c<facetsToChange.size();c++) {
		Facet* f=facets[facetsToChange[c]];
		int nbRemove=0;
		for (size_t i=0;(int)i<f->sh.nbIndex;i++) //count how many to remove			
			if (vertices3[f->indices[i]].selected)
				nbRemove++;
		int *newIndices = (int *)malloc((f->sh.nbIndex-nbRemove)*sizeof(int));
		int nb=0;
		for (size_t i=0;(int)i<f->sh.nbIndex;i++)
			if (!vertices3[f->indices[i]].selected) newIndices[nb++]=f->indices[i];

		SAFE_FREE(f->indices);f->indices=newIndices;
		SAFE_FREE(f->vertices2);
		SAFE_FREE(f->visible);
		f->sh.nbIndex -= nbRemove;
		f->vertices2 = (VERTEX2D *)malloc(f->sh.nbIndex*sizeof(VERTEX2D));
		memset(f->vertices2,0,f->sh.nbIndex * sizeof(VERTEX2D));
		f->visible = (BOOL *)malloc(f->sh.nbIndex*sizeof(BOOL));
		_ASSERTE(f->visible);
		memset(f->visible,0xFF,f->sh.nbIndex*sizeof(BOOL));
	}
	
	if (facetsToRemove.size()) {
		Facet   **newFacets = (Facet **)malloc((sh.nbFacet-facetsToRemove.size()) * sizeof(Facet *));
		size_t nextToRemove=0;
		size_t nextToAdd=0;
		for (size_t f=0;(int)f<sh.nbFacet;f++) {
			if (nextToRemove<facetsToRemove.size() && f==facetsToRemove[nextToRemove]) {
				delete facets[f];
				mApp->RenumberSelections(f);
				mApp->RenumberFormulas(f);
				nextToRemove++;
			} else {
				newFacets[nextToAdd++]=facets[f];
			}
		}
		SAFE_DELETE(facets);
		facets=newFacets;
		sh.nbFacet-=facetsToRemove.size();
	}

	DeleteIsolatedVertices(TRUE);

	DeleteGLLists(TRUE,TRUE);

	BuildGLList();
	CalcTotalOutGassing();

	//Debug memory check
	 //_ASSERTE (!_CrtDumpMemoryLeaks());;
	 _ASSERTE(_CrtCheckMemory());
}

void Geometry::RemoveSelected() {
	MolFlow *mApp = (MolFlow *)theApp;
	changedSinceSave=TRUE;
	int nb = 0;
	for(int i=0;i<sh.nbFacet;i++)
		if(facets[i]->selected) nb++;

	if(nb==0) return;
	/*
	if(sh.nbFacet-nb==0) {
	// Remove all
	Clear();
	return;
	}
	*/

	Facet   **f = (Facet **)malloc((sh.nbFacet-nb) * sizeof(Facet *));

	nb=0;
	for(int i=0;i<sh.nbFacet;i++) {
		if(facets[i]->selected) {
			delete facets[i];
			mApp->RenumberSelections(i);
			mApp->RenumberFormulas(i);
		} else {
			f[nb++] = facets[i];
		}
	}

	SAFE_FREE(facets);
	facets = f;
	sh.nbFacet = nb;

	// Delete old resources
	DeleteGLLists(TRUE,TRUE);

	BuildGLList();
	CalcTotalOutGassing();

	//Debug memory check
	 //_ASSERTE (!_CrtDumpMemoryLeaks());;
	 _ASSERTE(_CrtCheckMemory());
}

// -----------------------------------------------------------
/*BOOL AskToReset_Geom(Worker *work) {
MolFlow *mApp = (MolFlow *)theApp;
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

// -----------------------------------------------------------

int Geometry::ExplodeSelected(BOOL toMap,int desType,double exponent,double *values) {

	MolFlow *mApp = (MolFlow *)theApp;
	changedSinceSave = TRUE;
	if( nbSelected==0 ) return -1;

	// Check that all facet has a mesh
	BOOL ok = TRUE;
	int idx=0;
	while(ok && idx<sh.nbFacet) {
		if(facets[idx]->selected)
			ok = facets[idx]->hasMesh;
		idx++;
	}
	if( !ok ) return -2;

	int nb = 0;
	int FtoAdd = 0;
	int VtoAdd = 0;
	Facet::FACETGROUP *blocks = (Facet::FACETGROUP *)malloc(nbSelected * sizeof(Facet::FACETGROUP));

	for(int i=0;i<sh.nbFacet;i++) {
		if(facets[i]->selected) {
			facets[i]->Explode(blocks+nb);
			FtoAdd += blocks[nb].nbF;
			VtoAdd += blocks[nb].nbV;
			nb++;
		}
	}

	// Update vertex array
	VERTEX3D *ptrVert;
	int       vIdx;
	VERTEX3D *nVert = (VERTEX3D *) malloc( (sh.nbVertex + VtoAdd)*sizeof(VERTEX3D) );
	memcpy(nVert,vertices3,sh.nbVertex*sizeof(VERTEX3D));

	ptrVert = nVert + sh.nbVertex;
	vIdx = sh.nbVertex;
	nb=0;
	for(int i=0;i<sh.nbFacet;i++) {
		if(facets[i]->selected) {
			facets[i]->FillVertexArray(ptrVert);
			for(int j=0;j<blocks[nb].nbF;j++) {
				for(int k=0;k<blocks[nb].facets[j]->sh.nbIndex;k++) {
					blocks[nb].facets[j]->indices[k] = vIdx + k;
				}
				vIdx+=blocks[nb].facets[j]->sh.nbIndex;
			}
			ptrVert += blocks[nb].nbV;
			nb++;
		}
	}
	SAFE_FREE(vertices3);
	vertices3=nVert;
	for (int i=sh.nbVertex;i<sh.nbVertex+VtoAdd;i++)
		vertices3[i].selected=FALSE;
	sh.nbVertex += VtoAdd;

	// Update facet
	Facet   **f = (Facet **)malloc((sh.nbFacet+FtoAdd-nbSelected) * sizeof(Facet *));

	// Delete selected
	nb=0;
	for(int i=0;i<sh.nbFacet;i++) {
		if(facets[i]->selected) {
			delete facets[i];
			mApp->RenumberSelections(i);
			mApp->RenumberFormulas(i);
		} else {
			f[nb++] = facets[i];
		}
	}

	// Add new facets
	int count=0;
	for(int i=0;i<nbSelected;i++) {
		for(int j=0;j<blocks[i].nbF;j++) {
			f[nb++] = blocks[i].facets[j];
			if (toMap) { //set outgassing values
				f[nb-1]->sh.flow=*(values+count++);
				if (f[nb-1]->sh.flow>0.0) {
					f[nb-1]->sh.desorbType=desType+1;
					f[nb-1]->selected=TRUE;
					if (f[nb-1]->sh.desorbType==DES_COSINE_N) f[nb-1]->sh.desorbTypeN=exponent;
				} else {
					f[nb-1]->sh.desorbType=DES_NONE;
					f[nb-1]->selected=FALSE;
				}
			}
		}
	}

	// Free allocated memory
	for(int i=0;i<nbSelected;i++) {
		SAFE_FREE(blocks[i].facets);
	}
	SAFE_FREE(blocks);

	SAFE_FREE(facets);
	facets = f;
	sh.nbFacet = nb;

	// Delete old resources
	DeleteGLLists(TRUE,TRUE);

	InitializeGeometry();

	return 0;

}

// -----------------------------------------------------------

void Geometry::RemoveNullFacet() {
	MolFlow *mApp = (MolFlow *)theApp;
	// Remove degenerated facet (area~0.0)
	int nb = 0;
	double areaMin = vThreshold*vThreshold;
	for(int i=0;i<sh.nbFacet;i++)
		if(facets[i]->sh.area<areaMin) nb++;

	if(nb==0) return;

	Facet   **f = (Facet **)malloc((sh.nbFacet-nb) * sizeof(Facet *));

	nb=0;
	for(int i=0;i<sh.nbFacet;i++) {
		if(facets[i]->sh.area<areaMin) {
			delete facets[i];
			mApp->RenumberSelections(i);
			mApp->RenumberFormulas(i);
		} else {
			f[nb++] = facets[i];
		}
	}

	SAFE_FREE(facets);
	facets = f;
	sh.nbFacet = nb;

	// Delete old resources
	DeleteGLLists(TRUE,TRUE);

	BuildGLList();

}

// -----------------------------------------------------------

BOOL Geometry::IsLoaded() {
	return isLoaded;
}

// -----------------------------------------------------------
void Geometry::AlignFacets(int* selection,int nbSelected,int Facet_source,int Facet_dest,int Anchor_source,int Anchor_dest,
						   int Aligner_source,int Aligner_dest,BOOL invertNormal,BOOL invertDir1,BOOL invertDir2,BOOL copy,Worker *worker) {
							   MolFlow *mApp = (MolFlow *)theApp;
							   double counter=0.0;	
							   double selected=(double)GetNbSelected();
							   if (selected==0.0) return;
							   GLProgress *prgAlign = new GLProgress("Aligning facets...","Please wait");
							   prgAlign->SetProgress(0.0);
							   prgAlign->SetVisible(TRUE);
							   if (!mApp->AskToReset(worker)) return;
							   if (copy) CloneSelectedFacets(); //move
							   BOOL *alreadyMoved=(BOOL*)malloc(sh.nbVertex*sizeof(BOOL*));
							   memset(alreadyMoved,FALSE,sh.nbVertex*sizeof(BOOL*));

							   //Translating facets to align anchors
							   int temp;
							   if (invertDir1) { //change anchor and direction on source
								   temp=Aligner_source;
								   Aligner_source=Anchor_source;
								   Anchor_source=temp;
							   }
							   if (invertDir2) { //change anchor and direction on destination
								   temp=Aligner_dest;
								   Aligner_dest=Anchor_dest;
								   Anchor_dest=temp;
							   }
							   VERTEX3D Translation;
							   Sub(&Translation,GetVertex(Anchor_dest),GetVertex(Anchor_source));

							   int nb = 0;
							   for(int i=0;i<nbSelected;i++) {
								   counter+=0.333;
								   prgAlign->SetProgress(counter/selected);
								   for(int j=0;j<facets[selection[i]]->sh.nbIndex;j++) {
									   if (!alreadyMoved[facets[selection[i]]->indices[j]]) {
										   vertices3[facets[selection[i]]->indices[j]].x+=Translation.x;
										   vertices3[facets[selection[i]]->indices[j]].y+=Translation.y;
										   vertices3[facets[selection[i]]->indices[j]].z+=Translation.z;
										   alreadyMoved[facets[selection[i]]->indices[j]] = TRUE;
									   }
								   }
							   }

							   SAFE_FREE(alreadyMoved);

							   //Rotating to match normal vectors
							   VERTEX3D Axis;
							   VERTEX3D Normal;
							   double angle;
							   Normal=facets[Facet_dest]->sh.N;
							   if (invertNormal) ScalarMult(&Normal,-1.0);
							   Cross(&Axis,&(facets[Facet_source]->sh.N),&Normal);
							   if (Norme(&Axis)<1e-5) { //The two normals are either collinear or the opposite
								   if ((Dot(&(facets[Facet_dest]->sh.N),&(facets[Facet_source]->sh.N))>0.99999 && (!invertNormal)) || 
									   (Dot(&(facets[Facet_dest]->sh.N),&(facets[Facet_source]->sh.N))<0.00001 && invertNormal)) { //no rotation needed
										   Axis.x=1.0;
										   Axis.y=0.0;
										   Axis.z=0.0;
										   angle=0.0;
								   } else { //180deg rotation needed
									   Cross(&Axis,&(facets[Facet_source]->sh.U),&(facets[Facet_source]->sh.N));
									   Normalize(&Axis);
									   angle=180.0;
								   }
							   } else {
								   Normalize(&Axis);
								   angle=Dot(&(facets[Facet_source]->sh.N),&(facets[Facet_dest]->sh.N))/Norme(&(facets[Facet_source]->sh.N))/Norme(&(facets[Facet_dest]->sh.N));
								   //BOOL opposite=(angle<0.0);
								   angle=acos(angle);
								   angle=angle/PI*180;
								   if (invertNormal) angle=180.0-angle;
							   }


							   BOOL *alreadyRotated=(BOOL*)malloc(sh.nbVertex*sizeof(BOOL*));
							   memset(alreadyRotated,FALSE,sh.nbVertex*sizeof(BOOL*));


							   nb = 0;
							   for(int i=0;i<nbSelected;i++) {
								   counter+=0.333;
								   prgAlign->SetProgress(counter/selected);
								   for(int j=0;j<facets[selection[i]]->sh.nbIndex;j++) {
									   if (!alreadyRotated[facets[selection[i]]->indices[j]]) {
										   //rotation comes here
										   Rotate(&(vertices3[facets[selection[i]]->indices[j]]),*(GetVertex(Anchor_dest)),Axis,angle);
										   alreadyRotated[facets[selection[i]]->indices[j]] = TRUE;
									   }
								   }
							   }

							   SAFE_FREE(alreadyRotated);

							   //Rotating to match direction points

							   VERTEX3D Dir1,Dir2;
							   Sub(&Dir1,GetVertex(Aligner_dest),GetVertex(Anchor_dest));
							   Sub(&Dir2,GetVertex(Aligner_source),GetVertex(Anchor_source));
							   Cross(&Axis,&Dir2,&Dir1);
							   if (Norme(&Axis)<1e-5) { //The two directions are either collinear or the opposite
								   if (Dot(&Dir1,&Dir2)>0.99999) { //no rotation needed
									   Axis.x=1.0;
									   Axis.y=0.0;
									   Axis.z=0.0;
									   angle=0.0;
								   } else { //180deg rotation needed
									   //construct a vector perpendicular to the normal
									   Axis=facets[Facet_source]->sh.N;
									   Normalize(&Axis);
									   angle=180.0;
								   }
							   } else {
								   Normalize(&Axis);
								   angle=Dot(&Dir1,&Dir2)/Norme(&Dir1)/Norme(&Dir2);
								   //BOOL opposite=(angle<0.0);
								   angle=acos(angle);
								   angle=angle/PI*180;
								   //if (invertNormal) angle=180.0-angle;
							   }

							   BOOL *alreadyRotated2=(BOOL*)malloc(sh.nbVertex*sizeof(BOOL*));
							   memset(alreadyRotated2,FALSE,sh.nbVertex*sizeof(BOOL*));


							   nb = 0;
							   for(int i=0;i<nbSelected;i++) {
								   counter+=0.333;
								   prgAlign->SetProgress(counter/selected);
								   for(int j=0;j<facets[selection[i]]->sh.nbIndex;j++) {
									   if (!alreadyRotated2[facets[selection[i]]->indices[j]]) {
										   //rotation comes here
										   Rotate(&(vertices3[facets[selection[i]]->indices[j]]),*(GetVertex(Anchor_dest)),Axis,angle);
										   alreadyRotated2[facets[selection[i]]->indices[j]] = TRUE;
									   }
								   }
							   }

							   SAFE_FREE(alreadyRotated2);

							   InitializeGeometry();
							   //update textures
							   try {
								   for(int i=0;i<nbSelected;i++)
								   SetFacetTexture(selection[i],facets[selection[i]]->tRatio,facets[selection[i]]->hasMesh);	   
							   } catch (Error &e) {
								   GLMessageBox::Display((char *)e.GetMsg(),"Error",GLDLG_OK,GLDLG_ICONERROR);
								   return;
							   }
							   prgAlign->SetVisible(FALSE);
							   SAFE_DELETE(prgAlign);
}


// -------------------------------------------------------------
void Geometry::MoveSelectedFacets(double dX,double dY,double dZ,BOOL copy,Worker *worker) {
	MolFlow *mApp = (MolFlow *)theApp;
	GLProgress *prgMove = new GLProgress("Moving selected facets...","Please wait");
	prgMove->SetProgress(0.0);
	prgMove->SetVisible(TRUE);
	if (!(dX==0.0&&dY==0.0&&dZ==0.0)) {
		if (!mApp->AskToReset(worker)) return;
		int nbSelFacet=0;
		if (copy) CloneSelectedFacets(); //move
		double counter=1.0;
		double selected=(double)GetNbSelected();
		if (selected==0.0) return;

		BOOL *alreadyMoved=(BOOL*)malloc(sh.nbVertex*sizeof(BOOL*));
		memset(alreadyMoved,FALSE,sh.nbVertex*sizeof(BOOL*));


		int nb = 0;
		for(int i=0;i<sh.nbFacet;i++) {
			if(facets[i]->selected) {
				counter+=1.0;
				prgMove->SetProgress(counter/selected);
				for(int j=0;j<facets[i]->sh.nbIndex;j++) {
					if (!alreadyMoved[facets[i]->indices[j]]) {
						vertices3[facets[i]->indices[j]].x+=dX;
						vertices3[facets[i]->indices[j]].y+=dY;
						vertices3[facets[i]->indices[j]].z+=dZ;
						alreadyMoved[facets[i]->indices[j]] = TRUE;
					}
				}
			}
		}

		SAFE_FREE(alreadyMoved);

		InitializeGeometry();
		//update textures
		try {
			for(int i=0;i<sh.nbFacet;i++) if(facets[i]->selected) SetFacetTexture(i,facets[i]->tRatio,facets[i]->hasMesh);	   
		} catch (Error &e) {
			GLMessageBox::Display((char *)e.GetMsg(),"Error",GLDLG_OK,GLDLG_ICONERROR);
			return;
		}
	}
	prgMove->SetVisible(FALSE);
	SAFE_DELETE(prgMove);
}

// -------------------------------------------------------------
void Geometry::MirrorSelectedFacets(VERTEX3D P0,VERTEX3D N,BOOL copy,Worker *worker) {
	MolFlow *mApp = (MolFlow *)theApp;
	double selected=(double)GetNbSelected();
	double counter=0.0;
	if (selected==0.0) return;
	GLProgress *prgMirror = new GLProgress("Mirroring selected facets...","Please wait");
	prgMirror->SetProgress(0.0);
	prgMirror->SetVisible(TRUE);

	if (!mApp->AskToReset(worker)) return;
	int nbSelFacet=0;
	if (copy) CloneSelectedFacets(); //move
	BOOL *alreadyMirrored=(BOOL*)malloc(sh.nbVertex*sizeof(BOOL*));
	memset(alreadyMirrored,FALSE,sh.nbVertex*sizeof(BOOL*));


	int nb = 0;
	for(int i=0;i<sh.nbFacet;i++) {
		if(facets[i]->selected) {
			counter+=1.0;
			prgMirror->SetProgress(counter/selected);
			nbSelFacet++;
			for(int j=0;j<facets[i]->sh.nbIndex;j++) {
				if (!alreadyMirrored[facets[i]->indices[j]]) {
					//mirroring comes here
					Mirror(&(vertices3[facets[i]->indices[j]]),P0,N);
					alreadyMirrored[facets[i]->indices[j]] = TRUE;
				}
			}
		}
	}

	SAFE_FREE(alreadyMirrored);
	if (nbSelFacet==0) return;
	SwapNormal();
	InitializeGeometry();
	//update textures
	try {
		for(int i=0;i<sh.nbFacet;i++) if(facets[i]->selected) SetFacetTexture(i,facets[i]->tRatio,facets[i]->hasMesh);		   
	} catch (Error &e) {
		GLMessageBox::Display((char *)e.GetMsg(),"Error",GLDLG_OK,GLDLG_ICONERROR);
		return;
	}

	prgMirror->SetVisible(FALSE);
	SAFE_DELETE(prgMirror);
}

// -------------------------------------------------------------
void Geometry::RotateSelectedFacets(VERTEX3D AXIS_P0,VERTEX3D AXIS_DIR,double theta,BOOL copy,Worker *worker) {
	MolFlow *mApp = (MolFlow *)theApp;
	double selected=(double)GetNbSelected();
	double counter=0.0;
	if (selected==0.0) return;
	GLProgress *prgRotate = new GLProgress("Rotating selected facets...","Please wait");
	prgRotate->SetProgress(0.0);
	prgRotate->SetVisible(TRUE);
	if (theta!=0.0) {
		if (!mApp->AskToReset(worker)) return;
		if (copy) CloneSelectedFacets(); //move
		BOOL *alreadyRotated=(BOOL*)malloc(sh.nbVertex*sizeof(BOOL*));
		memset(alreadyRotated,FALSE,sh.nbVertex*sizeof(BOOL*));


		int nb = 0;
		for(int i=0;i<sh.nbFacet;i++) {
			if(facets[i]->selected) {
				counter+=1.0;
				prgRotate->SetProgress(counter/selected);
				for(int j=0;j<facets[i]->sh.nbIndex;j++) {
					if (!alreadyRotated[facets[i]->indices[j]]) {
						//rotation comes here
						Rotate(&(vertices3[facets[i]->indices[j]]),AXIS_P0,AXIS_DIR,theta);
						alreadyRotated[facets[i]->indices[j]] = TRUE;
					}
				}
			}
		}

		SAFE_FREE(alreadyRotated);
		InitializeGeometry();
		//update textures
		try {
			for(int i=0;i<sh.nbFacet;i++) if(facets[i]->selected) SetFacetTexture(i,facets[i]->tRatio,facets[i]->hasMesh);		   
		} catch (Error &e) {
			GLMessageBox::Display((char *)e.GetMsg(),"Error",GLDLG_OK,GLDLG_ICONERROR);
			return;
		}
	}
	prgRotate->SetVisible(FALSE);
	SAFE_DELETE(prgRotate);
}


// -------------------------------------------------------------
void Geometry::CloneSelectedFacets() { //create clone of selected facets
	double counter=0.0;
	double selected=(double)GetNbSelected();
	if (selected==0.0) return;
	int *copyId=(int*)malloc(sh.nbVertex*sizeof(int*));
	memset(copyId,-1,sh.nbVertex*sizeof(int*));

	//count how many new vertices to create

	int nb = sh.nbVertex-1;
	for(int i=0;i<sh.nbFacet;i++) {
		if(facets[i]->selected) {
			counter+=1.0;
			for(int j=0;j<facets[i]->sh.nbIndex;j++) {
				if (copyId[facets[i]->indices[j]]==-1) {
					nb++;
					copyId[facets[i]->indices[j]]=nb;
				}
			}
		}
	}
	/*
	VERTEX3D *newVertices=(VERTEX3D*)malloc((nb+1)*sizeof(VERTEX3D*));
	memcpy(newVertices,vertices3,sh.nbVertex*sizeof(VERTEX3D));
	SAFE_FREE(vertices3);
	*/
	//vertices3=(VERTEX3D *)realloc(vertices3,(nb+1)*sizeof(VERTEX3D));
	VERTEX3D *tmp_vertices3 = (VERTEX3D *)malloc((nb+1) * sizeof(VERTEX3D));
	memmove(tmp_vertices3,vertices3,(sh.nbVertex)*sizeof(VERTEX3D));
	memset(tmp_vertices3+sh.nbVertex,0,(nb+1-sh.nbVertex) * sizeof(VERTEX3D));
	SAFE_FREE(vertices3);
	vertices3=tmp_vertices3;
	for (int i=0;i<sh.nbVertex;i++) {
		if (copyId[i]!=-1) {
			vertices3[copyId[i]].x=vertices3[i].x;
			vertices3[copyId[i]].y=vertices3[i].y;
			vertices3[copyId[i]].z=vertices3[i].z;
			vertices3[copyId[i]].selected=vertices3[i].selected;
		}
	}
	sh.nbVertex=nb+1;
	sh.nbFacet+=(int)selected;
	facets = (Facet **)realloc(facets, sh.nbFacet * sizeof(Facet *));
	int nb2=sh.nbFacet-(int)selected-1; //copy new facets
	for(int i=0;i<sh.nbFacet-(int)selected;i++) {
		if(facets[i]->selected) {
			nb2++;
			facets[nb2] = new Facet(facets[i]->sh.nbIndex);
			facets[nb2]->Copy(facets[i],TRUE);
			//copy indices
			for (int j=0;j<facets[i]->sh.nbIndex;j++) {
				facets[nb2]->indices[j]=facets[i]->indices[j];
			} 
			facets[i]->selected = FALSE;
		}
	}
	for(int i=(sh.nbFacet-(int)selected);i<sh.nbFacet;i++) {
		for (int j=0;j<facets[i]->sh.nbIndex;j++) {
			facets[i]->indices[j]=copyId[facets[i]->indices[j]];
		}
	}
	CalcTotalOutGassing();

	//Debug memory check
	 //_ASSERTE (!_CrtDumpMemoryLeaks());;
	 _ASSERTE(_CrtCheckMemory());
}

// -------------------------------------------------------------
void Geometry::MoveSelectedVertex(double dX,double dY,double dZ,BOOL copy,Worker *worker) {

	MolFlow *mApp = (MolFlow *)theApp;
	if (!(dX==0.0&&dY==0.0&&dZ==0.0)) {
		if (!mApp->AskToReset(worker)) return;
		changedSinceSave = TRUE;
		if (!copy) { //move
			for(int i=0;i<sh.nbVertex;i++) {
				if (vertices3[i].selected) {
					vertices3[i].x+=dX;
					vertices3[i].y+=dY;
					vertices3[i].z+=dZ;
				}
			}
			InitializeGeometry();
		} else { //copy
			int nbVertexOri=sh.nbVertex;
			for(int i=0;i<nbVertexOri;i++) {
				if (vertices3[i].selected) {
					AddVertex(vertices3[i].x+dX,vertices3[i].y+dY,vertices3[i].z+dZ);
				}
			}
		}

	}
}

// -------------------------------------------------------------
void Geometry::AddVertex(double X,double Y,double Z) {

	changedSinceSave = TRUE;

	//a new vertex
	sh.nbVertex++;
	VERTEX3D *verticesNew = (VERTEX3D *)malloc(sh.nbVertex * sizeof(VERTEX3D));
	memcpy(verticesNew,vertices3,(sh.nbVertex-1) * sizeof(VERTEX3D)); //copy old vertices
	SAFE_FREE(vertices3);
	verticesNew[sh.nbVertex-1].x=X;
	verticesNew[sh.nbVertex-1].y=Y;
	verticesNew[sh.nbVertex-1].z=Z;
	verticesNew[sh.nbVertex-1].selected=TRUE;
	vertices3=verticesNew;

	//InitializeGeometry();

}

void Geometry::GetSelection(int **selection,int *nbSel){
	int sel=0;
	*nbSel=GetNbSelected();
	int *selected = (int *)malloc((*nbSel)*sizeof(int));
	for (int i=0;i<sh.nbFacet;i++)
		if (facets[i]->selected) selected[sel++]=i;
	*selection=selected;
}

void Geometry::SetSelection(int **selection,int *nbSel){
	MolFlow *mApp = (MolFlow *)theApp;
	UnSelectAll();
	for (int i=0;i<*nbSel;i++) {
		int toSelect=(*selection)[i];
		if (toSelect<sh.nbFacet) facets[toSelect]->selected = TRUE;
	}
	UpdateSelection();
	mApp->facetList->ScrollToVisible((*selection)[*nbSel-1],0,TRUE); //in facet list, select the last facet of selection group
	mApp->UpdateFacetParams(TRUE);
}

void Geometry::AddStruct(char *name) {
	strName[sh.nbSuper++]=_strdup(name);
	BuildGLList();
}

void Geometry::DelStruct(int numToDel) {
	MolFlow *mApp = (MolFlow *)theApp;
	RemoveFromStruct(numToDel);
	CheckIsolatedVertex();
	mApp->UpdateModelParams();

	for (int i=0;i<sh.nbFacet;i++) {
		if (facets[i]->sh.superIdx>numToDel) facets[i]->sh.superIdx--;
		if (facets[i]->sh.superDest>numToDel) facets[i]->sh.superDest--;
	}
	for (int j=numToDel;j<(sh.nbSuper-1);j++)
	{
		strName[j]=_strdup(strName[j+1]);
	}
	sh.nbSuper--;
	BuildGLList();

	//Debug memory check
	 //_ASSERTE (!_CrtDumpMemoryLeaks());;
	 _ASSERTE(_CrtCheckMemory());
}


void Geometry::ScaleSelectedVertices(VERTEX3D invariant,double factor,BOOL copy,Worker *worker) {

	MolFlow *mApp = (MolFlow *)theApp;

	if (!mApp->AskToReset(worker)) return;
	changedSinceSave = TRUE;
	if (!copy) { //scale
		for(int i=0;i<sh.nbVertex;i++) {
			if (vertices3[i].selected) {
				vertices3[i].x=invariant.x+factor*(vertices3[i].x-invariant.x);
				vertices3[i].y=invariant.y+factor*(vertices3[i].y-invariant.y);
				vertices3[i].z=invariant.z+factor*(vertices3[i].z-invariant.z);
			}
		}
	} else { //scale and copy
		int nbVertexOri=sh.nbVertex;
		for(int i=0;i<nbVertexOri;i++) {
			if (vertices3[i].selected) {
				AddVertex(invariant.x+factor*(vertices3[i].x-invariant.x),
					invariant.y+factor*(vertices3[i].y-invariant.y),
					invariant.z+factor*(vertices3[i].z-invariant.z));
			}
		}
	}
	InitializeGeometry();
}

void Geometry::ScaleSelectedFacets(VERTEX3D invariant,double factorX,double factorY,double factorZ,BOOL copy,Worker *worker) {

	MolFlow *mApp = (MolFlow *)theApp;
	GLProgress *prgMove = new GLProgress("Scaling selected facets...","Please wait");
	prgMove->SetProgress(0.0);
	prgMove->SetVisible(TRUE);

	if (!mApp->AskToReset(worker)) return;
	int nbSelFacet=0;
	if (copy) CloneSelectedFacets(); //move
	double counter=1.0;
	double selected=(double)GetNbSelected();
	if (selected==0.0) return;

	BOOL *alreadyMoved=(BOOL*)malloc(sh.nbVertex*sizeof(BOOL*));
	memset(alreadyMoved,FALSE,sh.nbVertex*sizeof(BOOL*));


	int nb = 0;
	for(int i=0;i<sh.nbFacet;i++) {
		if(facets[i]->selected) {
			counter+=1.0;
			prgMove->SetProgress(counter/selected);
			for(int j=0;j<facets[i]->sh.nbIndex;j++) {
				if (!alreadyMoved[facets[i]->indices[j]]) {
					vertices3[facets[i]->indices[j]].x=invariant.x+factorX*(vertices3[facets[i]->indices[j]].x-invariant.x);
					vertices3[facets[i]->indices[j]].y=invariant.y+factorY*(vertices3[facets[i]->indices[j]].y-invariant.y);
					vertices3[facets[i]->indices[j]].z=invariant.z+factorZ*(vertices3[facets[i]->indices[j]].z-invariant.z);
					alreadyMoved[facets[i]->indices[j]] = TRUE;
				}
			}
		}
	}

	SAFE_FREE(alreadyMoved);

	InitializeGeometry();
	//update textures
	//for(int i=0;i<sh.nbFacet;i++) if(facets[i]->selected) SetFacetTexture(i,facets[i]->tRatio,facets[i]->hasMesh);	   

	prgMove->SetVisible(FALSE);
	SAFE_DELETE(prgMove);
}

  void Geometry::EmptySelectedVertexList(){
	  selectedVertexList=std::vector<int>();
  }

  void Geometry::RemoveFromSelectedVertexList(int vertexId){
	  for (size_t j=0;j<selectedVertexList.size();j++)
			  if (selectedVertexList[j]==vertexId)
				  selectedVertexList.erase(selectedVertexList.begin()+j);
  }

  void Geometry::AddToSelectedVertexList(int vertexId){
	  selectedVertexList.push_back(vertexId);
  }