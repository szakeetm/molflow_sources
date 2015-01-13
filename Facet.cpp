/*
  File:        Facet.cpp
  Description: Facet class (memory management)
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

#include "Facet.h"
#include "Utils.h"
#include <malloc.h>
#include <crtdbg.h>
#include <string.h>
#include <math.h>
#include "GLToolkit.h"


using namespace pugi;
/*
//Leak detection
#ifdef _DEBUG
#define _CRTDBG_MAP_ALLOC
#define DEBUG_NEW new(_NORMAL_BLOCK, __FILE__, __LINE__)
#define new DEBUG_NEW
#endif
*/

#define MAX(x,y) (((x)<(y))?(y):(x))
#define MIN(x,y) (((x)<(y))?(x):(y))

// Colormap stuff
extern COLORREF rainbowCol[];

static int colorMap[65536];
static BOOL colorInited = FALSE;

// -----------------------------------------------------------

Facet::Facet(int nbIndex) {

	indices = (int *)malloc(nbIndex*sizeof(int));                    // Ref to Geometry VERTEX3D
	vertices2 = (VERTEX2D *)malloc(nbIndex * sizeof(VERTEX2D));      // Local U,V coordinates
	memset(vertices2, 0, nbIndex * sizeof(VERTEX2D));

	sh.nbIndex = nbIndex;
	memset(&sh.counter, 0, sizeof(sh.counter));
	/*sh.counter.hit.nbDesorbed = 0;
	sh.counter.hit.nbAbsorbed = 0;
	sh.counter.hit.nbHit = 0;*/

	sh.sticking = 0.0;
	sh.opacity = 1.0;
	sh.temperature = 293.15; // 20degC
	sh.flow = 0.0;           // 1 unit*l/s //will be outgasssing
	sh.mass = 28.0;          // Nitrogen
	sh.desorbType = DES_NONE;
	sh.desorbTypeN = 0.0;
	sh.reflectType = REF_DIFFUSE;
	sh.profileType = REC_NONE;
	sh.texWidth = 0;
	sh.texHeight = 0;
	sh.texWidthD = 0.0;
	sh.texHeightD = 0.0;
	sh.center.x = 0.0;
	sh.center.y = 0.0;
	sh.center.z = 0.0;
	sh.is2sided = FALSE;
	sh.isProfile = FALSE;
	//sh.isOpaque = TRUE;
	sh.isTextured = FALSE;
	sh.sign = 0.0;
	sh.countDes = FALSE;
	sh.countAbs = FALSE;
	sh.countRefl = FALSE;
	sh.countTrans = FALSE;
	sh.countACD = FALSE;
	sh.countDirection = FALSE;
	sh.superIdx = 0;
	sh.superDest = 0;
	sh.teleportDest = 0;
	sh.isVolatile = FALSE;
	sh.useOutgassingFile = FALSE;
	sh.accomodationFactor = 1.0;

	sh.outgassing_paramId=-1;
	sh.opacity_paramId=-1;
	sh.sticking_paramId=-1;

	sh.isMoving = FALSE;

	hasOutgassingMap = FALSE;
	outgassingMap = NULL;

	textureVisible = TRUE;
	volumeVisible = TRUE;

	texDimW = 0;
	texDimH = 0;
	tRatio = 0.0;

	mesh = NULL;
	meshPts = NULL;
	hasMesh = FALSE;
	nbElem = 0;
	selectedElem.u = 0;
	selectedElem.v = 0;
	selectedElem.width = 0;
	selectedElem.height = 0;
	dirCache = NULL;
	textureError = FALSE;
	hasOutgassingMap = FALSE;

	userOutgassing = "";
	userOpacity = "";
	userSticking = "";

	// Init the colormap at the first facet construction
	for (int i = 0; i < 65536 && !colorInited; i++) {

		double r1, g1, b1;
		double r2, g2, b2;
		int colId = i / 8192;
		//int colId = i/10923;

		r1 = (double)((rainbowCol[colId] >> 16) & 0xFF);
		g1 = (double)((rainbowCol[colId] >> 8) & 0xFF);
		b1 = (double)((rainbowCol[colId] >> 0) & 0xFF);

		r2 = (double)((rainbowCol[colId + 1] >> 16) & 0xFF);
		g2 = (double)((rainbowCol[colId + 1] >> 8) & 0xFF);
		b2 = (double)((rainbowCol[colId + 1] >> 0) & 0xFF);

		double rr = (double)(i - colId * 8192) / 8192.0;
		//double rr = (double)(i-colId*10923) / 10923;
		SATURATE(rr, 0.0, 1.0);
		colorMap[i] = (COLORREF)((int)(r1 + (r2 - r1)*rr) +
			(int)(g1 + (g2 - g1)*rr) * 256 +
			(int)(b1 + (b2 - b1)*rr) * 65536);

	}
	colorMap[65535] = 0xFFFFFF; // Saturation color
	colorInited = TRUE;

	glTex = 0;
	glList = 0;
	glElem = 0;
	glSelElem = 0;
	selected = FALSE;
	visible = (BOOL *)malloc(nbIndex*sizeof(BOOL));
	memset(visible, 0xFF, nbIndex*sizeof(BOOL));
	//visible[5]=1; //Troll statement to corrupt heap (APPVERIF debug test)
}

// -----------------------------------------------------------

Facet::~Facet() {
	SAFE_FREE(indices);
	SAFE_FREE(vertices2);
	SAFE_FREE(mesh);
	SAFE_FREE(dirCache);
	DELETE_TEX(glTex);
	DELETE_LIST(glList);
	DELETE_LIST(glElem);
	DELETE_LIST(glSelElem);
	SAFE_FREE(visible);
	for (int i = 0; i < nbElem; i++)
		SAFE_FREE(meshPts[i].pts);
	SAFE_FREE(meshPts);
	SAFE_FREE(outgassingMap);
}

// -----------------------------------------------------------

BOOL Facet::IsLinkFacet() {
	return ((sh.opacity == 0.0) && (sh.sticking >= 1.0));
}

// -----------------------------------------------------------

void Facet::LoadGEO(FileReader *file, int version, int nbVertex) {

	file->ReadKeyword("indices"); file->ReadKeyword(":");
	for (int i = 0; i < sh.nbIndex; i++) {
		indices[i] = file->ReadInt() - 1;
		if (indices[i] >= nbVertex)
			throw Error(file->MakeError("Facet index out of bounds"));
	}

	file->ReadKeyword("sticking"); file->ReadKeyword(":");
	sh.sticking = file->ReadDouble();
	file->ReadKeyword("opacity"); file->ReadKeyword(":");
	sh.opacity = file->ReadDouble();
	file->ReadKeyword("desorbType"); file->ReadKeyword(":");
	sh.desorbType = file->ReadInt();
	if (version >= 9) {
		file->ReadKeyword("desorbTypeN"); file->ReadKeyword(":");
		sh.desorbTypeN = file->ReadDouble();
	}
	else {
		ConvertOldDesorbType();
	}
	file->ReadKeyword("reflectType"); file->ReadKeyword(":");
	sh.reflectType = file->ReadInt();
	file->ReadKeyword("profileType"); file->ReadKeyword(":");
	sh.profileType = file->ReadInt();
	file->ReadKeyword("superDest"); file->ReadKeyword(":");
	sh.superDest = file->ReadInt();
	file->ReadKeyword("superIdx"); file->ReadKeyword(":");
	sh.superIdx = file->ReadInt();
	file->ReadKeyword("is2sided"); file->ReadKeyword(":");
	sh.is2sided = file->ReadInt();
	if (version < 8) {
		file->ReadKeyword("area"); file->ReadKeyword(":");
		sh.area = file->ReadDouble();
	}
	file->ReadKeyword("mesh"); file->ReadKeyword(":");
	hasMesh = file->ReadInt();
	if (version >= 7) {
		file->ReadKeyword("outgassing"); file->ReadKeyword(":");
		sh.flow = file->ReadDouble()*0.100; //mbar*l/s -> Pa*m3/s
	}
	file->ReadKeyword("texDimX"); file->ReadKeyword(":");
	sh.texWidthD = file->ReadDouble();
	file->ReadKeyword("texDimY"); file->ReadKeyword(":");
	sh.texHeightD = file->ReadDouble();
	file->ReadKeyword("countDes"); file->ReadKeyword(":");
	sh.countDes = file->ReadInt();
	file->ReadKeyword("countAbs"); file->ReadKeyword(":");
	sh.countAbs = file->ReadInt();
	file->ReadKeyword("countRefl"); file->ReadKeyword(":");
	sh.countRefl = file->ReadInt();
	file->ReadKeyword("countTrans"); file->ReadKeyword(":");
	sh.countTrans = file->ReadInt();
	file->ReadKeyword("acMode"); file->ReadKeyword(":");
	sh.countACD = file->ReadInt();
	file->ReadKeyword("nbAbs"); file->ReadKeyword(":");
	sh.counter.hit.nbAbsorbed = file->ReadLLong();
	file->ReadKeyword("nbDes"); file->ReadKeyword(":");
	sh.counter.hit.nbDesorbed = file->ReadLLong();
	file->ReadKeyword("nbHit"); file->ReadKeyword(":");
	sh.counter.hit.nbHit = file->ReadLLong();
	if (version >= 2) {
		// Added in GEO version 2
		file->ReadKeyword("temperature"); file->ReadKeyword(":");
		sh.temperature = file->ReadDouble();
		file->ReadKeyword("countDirection"); file->ReadKeyword(":");
		sh.countDirection = file->ReadInt();
	}
	if (version >= 4) {
		// Added in GEO version 4
		file->ReadKeyword("textureVisible"); file->ReadKeyword(":");
		textureVisible = file->ReadInt();
		file->ReadKeyword("volumeVisible"); file->ReadKeyword(":");
		volumeVisible = file->ReadInt();
	}

	if (version >= 5) {
		// Added in GEO version 5
		file->ReadKeyword("teleportDest"); file->ReadKeyword(":");
		sh.teleportDest = file->ReadInt();
	}

	if (version >= 13) {
		// Added in GEO version 13
		file->ReadKeyword("accomodationFactor"); file->ReadKeyword(":");
		sh.accomodationFactor = file->ReadDouble();
	}

	UpdateFlags();

}

void Facet::LoadXML(xml_node f,int nbVertex) {
	int idx = 0;
	for (xml_node indice : f.child("Indices").children("Indice")) {
		indices[idx] = indice.attribute("vertex").as_int();
		if (indices[idx] >= nbVertex) {
			char err[128];
			sprintf(err, "Facet %d refers to vertex %d which doesn't exist", f.attribute("id").as_int() + 1, idx + 1);
			throw Error(err);
		}
		idx++;
	}
	sh.sticking = f.child("Sticking").attribute("constValue").as_double();
	sh.sticking_paramId = f.child("Sticking").attribute("parameterId").as_int();
	sh.opacity = f.child("Opacity").attribute("constValue").as_double();
	sh.opacity_paramId = f.child("Opacity").attribute("parameterId").as_int();
	sh.is2sided = f.child("Opacity").attribute("is2sided").as_int();
	sh.flow = f.child("Outgassing").attribute("constValue").as_double();
	sh.desorbType = f.child("Outgassing").attribute("desType").as_int();
	sh.desorbTypeN = f.child("Outgassing").attribute("desExponent").as_double();
	sh.outgassing_paramId = f.child("Outgassing").attribute("parameterId").as_int();
	//TO DO: store dynamic outgassing
	sh.temperature = f.child("Temperature").attribute("value").as_double();
	sh.accomodationFactor = f.child("Temperature").attribute("accFactor").as_double();
	sh.reflectType = f.child("Reflection").attribute("type").as_int();
	sh.superIdx = f.child("Structure").attribute("inStructure").as_int();
	sh.superDest = f.child("Structure").attribute("linksTo").as_int();
	sh.teleportDest = f.child("Teleport").attribute("target").as_int();
	sh.isMoving = f.child("Motion").attribute("isMoving").as_bool();
	xml_node recNode = f.child("Recordings");
	sh.profileType = recNode.child("Profile").attribute("type").as_int();
	xml_node texNode = recNode.child("Texture");
	hasMesh = texNode.attribute("hasMesh").as_bool();
	sh.texWidthD = texNode.attribute("texDimX").as_double();
	sh.texHeightD = texNode.attribute("texDimY").as_double();
	sh.countDes = texNode.attribute("countDes").as_int();
	sh.countAbs = texNode.attribute("countAbs").as_int();
	sh.countRefl = texNode.attribute("countRefl").as_int();
	sh.countTrans = texNode.attribute("countTrans").as_int();
	sh.countDirection = texNode.attribute("countDir").as_int();
	sh.countACD = texNode.attribute("countAC").as_int();

	textureVisible = f.child("ViewSettings").attribute("textureVisible").as_int();
	volumeVisible = f.child("ViewSettings").attribute("volumeVisible").as_int();
	
	UpdateFlags();
}


void Facet::LoadSYN(FileReader *file, int version, int nbVertex) {

	file->ReadKeyword("indices"); file->ReadKeyword(":");
	for (int i = 0; i < sh.nbIndex; i++) {
		indices[i] = file->ReadInt() - 1;
		if (indices[i] >= nbVertex)
			throw Error(file->MakeError("Facet index out of bounds"));
	}

	file->ReadKeyword("sticking"); file->ReadKeyword(":");
	sh.sticking = file->ReadDouble();
	if (version >= 4) {
		file->ReadKeyword("roughness"); file->ReadKeyword(":");
		file->ReadDouble();
	}
	file->ReadKeyword("opacity"); file->ReadKeyword(":");
	sh.opacity = file->ReadDouble();
	file->ReadKeyword("reflectType"); file->ReadKeyword(":");
	sh.reflectType = file->ReadInt();
	if (sh.reflectType > REF_MIRROR) sh.reflectType = REF_DIFFUSE; //treat material reflection
	file->ReadKeyword("profileType"); file->ReadKeyword(":");
	sh.profileType = 0; file->ReadInt();
	file->ReadKeyword("hasSpectrum"); file->ReadKeyword(":");
	file->ReadInt();
	file->ReadKeyword("superDest"); file->ReadKeyword(":");
	sh.superDest = file->ReadInt();
	file->ReadKeyword("superIdx"); file->ReadKeyword(":");
	sh.superIdx = file->ReadInt();
	file->ReadKeyword("is2sided"); file->ReadKeyword(":");
	sh.is2sided = file->ReadInt();
	file->ReadKeyword("mesh"); file->ReadKeyword(":");
	hasMesh = FALSE; file->ReadInt();
	file->ReadKeyword("texDimX"); file->ReadKeyword(":");
	sh.texWidthD = 0.0; file->ReadDouble();
	file->ReadKeyword("texDimY"); file->ReadKeyword(":");
	sh.texHeightD = 0.0; file->ReadDouble();
	if (version < 3) {
		file->ReadKeyword("countDes"); file->ReadKeyword(":");
		file->ReadInt();
	}
	file->ReadKeyword("countAbs"); file->ReadKeyword(":");
	sh.countAbs = FALSE; file->ReadInt();
	file->ReadKeyword("countRefl"); file->ReadKeyword(":");
	sh.countRefl = FALSE; file->ReadInt();
	file->ReadKeyword("countTrans"); file->ReadKeyword(":");
	sh.countTrans = FALSE; file->ReadInt();
	file->ReadKeyword("nbAbs"); file->ReadKeyword(":");
	sh.counter.hit.nbAbsorbed = 0; file->ReadLLong();
	if (version < 3) {
		file->ReadKeyword("nbDes"); file->ReadKeyword(":");
		sh.counter.hit.nbDesorbed = 0;
		file->ReadLLong();
	}
	file->ReadKeyword("nbHit"); file->ReadKeyword(":");
	sh.counter.hit.nbHit = 0; file->ReadLLong();
	if (version >= 3) {
		file->ReadKeyword("fluxAbs"); file->ReadKeyword(":");
		file->ReadDouble();
		file->ReadKeyword("powerAbs"); file->ReadKeyword(":");
		file->ReadDouble();
	}
	file->ReadKeyword("countDirection"); file->ReadKeyword(":");
	sh.countDirection = FALSE; file->ReadInt();
	file->ReadKeyword("textureVisible"); file->ReadKeyword(":");
	textureVisible = file->ReadInt();
	file->ReadKeyword("volumeVisible"); file->ReadKeyword(":");
	volumeVisible = file->ReadInt();
	file->ReadKeyword("teleportDest"); file->ReadKeyword(":");
	sh.teleportDest = file->ReadInt();

	UpdateFlags();

}


// -----------------------------------------------------------

void Facet::LoadTXT(FileReader *file) {

	// Opacity parameters descripton (TXT format)
	// -4    => Pressure profile (1 sided)
	// -3    => Desorption distribution
	// -2    => Angular profile
	// -1    => Pressure profile (2 sided)
	// [0,1] => Partial opacity (1 sided)
	// [1,2] => Partial opacity (2 sided)

	// Read facet parameters from TXT format
	sh.sticking = file->ReadDouble();
	double o = file->ReadDouble();
	sh.area = file->ReadDouble();
	sh.counter.hit.nbDesorbed = (llong)(file->ReadDouble() + 0.5);
	sh.counter.hit.nbHit = (llong)(file->ReadDouble() + 0.5);
	sh.counter.hit.nbAbsorbed = (llong)(file->ReadDouble() + 0.5);
	sh.desorbType = (int)(file->ReadDouble() + 0.5);

	// Convert opacity
	sh.profileType = REC_NONE;
	if (o < 0.0) {

		sh.opacity = 0.0;
		if (IS_ZERO(o + 1.0)) {
			sh.profileType = REC_PRESSUREU;
			sh.is2sided = TRUE;
		}
		if (IS_ZERO(o + 2.0))
			sh.profileType = REC_ANGULAR;
		if (IS_ZERO(o + 4.0)) {
			sh.profileType = REC_PRESSUREU;
			sh.is2sided = FALSE;
		}

	}
	else {
		if (o >= 1.0000001) {
			sh.opacity = o - 1.0;
			sh.is2sided = TRUE;
		}
		else
			sh.opacity = o;
	}

	// Convert desorbType
	switch (sh.desorbType) {
	case 0:
		sh.desorbType = DES_COSINE;
		break;
	case 1:
		sh.desorbType = DES_UNIFORM;
		break;
	case 2:
	case 3:
	case 4:
		sh.desorbType = sh.desorbType + 1; // cos^n
		break;
	}
	ConvertOldDesorbType();
	sh.reflectType = (int)(file->ReadDouble() + 0.5);

	// Convert reflectType
	switch (sh.reflectType) {
	case 0:
		sh.reflectType = REF_DIFFUSE;
		break;
	case 1:
		sh.reflectType = REF_MIRROR;
		break;
	default:
		sh.reflectType = REF_DIFFUSE;
		break;
	}

	file->ReadDouble(); // Unused

	if (sh.counter.hit.nbDesorbed == 0)
		sh.desorbType = DES_NONE;

	if (IsLinkFacet())
		sh.superDest = (int)(sh.sticking + 0.5);

	UpdateFlags();

}

void Facet::SaveTXT(FileWriter *file) {

	if (!sh.superDest)
		file->WriteDouble(sh.sticking, "\n");
	else {
		file->WriteDouble((double)sh.superDest, "\n");
		sh.opacity = 0.0;
	}

	if (sh.is2sided)
		file->WriteDouble(sh.opacity + 1.0, "\n");
	else
		file->WriteDouble(sh.opacity, "\n");

	file->WriteDouble(sh.area, "\n");

	if (sh.desorbType != DES_NONE)
		file->WriteDouble(1.0, "\n");
	else
		file->WriteDouble(0.0, "\n");
	file->WriteDouble(0.0, "\n"); //nbHit
	file->WriteDouble(0.0, "\n"); //nbAbsorbed

	file->WriteDouble(0.0, "\n"); //no desorption

	switch (sh.reflectType) {
	case REF_DIFFUSE:
		file->WriteDouble(0.0, "\n");
		break;
	case REF_MIRROR:
		file->WriteDouble(1.0, "\n");
		break;
	case REF_UNIFORM:
		file->WriteDouble(2.0, "\n");
	default:
		file->WriteDouble((double)(sh.reflectType), "\n");
		break;
	}

	file->WriteDouble(0.0, "\n"); // Unused
}

void Facet::SaveGEO(FileWriter *file, int idx) {

	char tmp[256];

	sprintf(tmp, "facet %d {\n", idx + 1);
	file->Write(tmp);
	file->Write("  nbIndex:"); file->WriteInt(sh.nbIndex, "\n");
	file->Write("  indices:\n");
	for (int i = 0; i < sh.nbIndex; i++) {
		file->Write("    ");
		file->WriteInt(indices[i] + 1, "\n");
	}
	//file->Write("\n");
	file->Write("  sticking:"); file->WriteDouble(sh.sticking, "\n");
	file->Write("  opacity:"); file->WriteDouble(sh.opacity, "\n");
	file->Write("  desorbType:"); file->WriteInt(sh.desorbType, "\n");
	file->Write("  desorbTypeN:"); file->WriteDouble(sh.desorbTypeN, "\n");
	file->Write("  reflectType:"); file->WriteInt(sh.reflectType, "\n");
	file->Write("  profileType:"); file->WriteInt(sh.profileType, "\n");
	file->Write("  superDest:"); file->WriteInt(sh.superDest, "\n");
	file->Write("  superIdx:"); file->WriteInt(sh.superIdx, "\n");
	file->Write("  is2sided:"); file->WriteInt(sh.is2sided, "\n");
	file->Write("  mesh:"); file->WriteInt((mesh != NULL), "\n");
	file->Write("  outgassing:"); file->WriteDouble(sh.flow*10.00, "\n"); //Pa*m3/s -> mbar*l/s for compatibility with old versions
	file->Write("  texDimX:"); file->WriteDouble(sh.texWidthD, "\n");
	file->Write("  texDimY:"); file->WriteDouble(sh.texHeightD, "\n");
	file->Write("  countDes:"); file->WriteInt(sh.countDes, "\n");
	file->Write("  countAbs:"); file->WriteInt(sh.countAbs, "\n");
	file->Write("  countRefl:"); file->WriteInt(sh.countRefl, "\n");
	file->Write("  countTrans:"); file->WriteInt(sh.countTrans, "\n");
	file->Write("  acMode:"); file->WriteInt(sh.countACD, "\n");
	file->Write("  nbAbs:"); file->WriteLLong(sh.counter.hit.nbAbsorbed, "\n");
	file->Write("  nbDes:"); file->WriteLLong(sh.counter.hit.nbDesorbed, "\n");
	file->Write("  nbHit:"); file->WriteLLong(sh.counter.hit.nbHit, "\n");

	// Version 2
	file->Write("  temperature:"); file->WriteDouble(sh.temperature, "\n");
	file->Write("  countDirection:"); file->WriteInt(sh.countDirection, "\n");

	// Version 4
	file->Write("  textureVisible:"); file->WriteInt(textureVisible, "\n");
	file->Write("  volumeVisible:"); file->WriteInt(volumeVisible, "\n");

	// Version 5
	file->Write("  teleportDest:"); file->WriteInt(sh.teleportDest, "\n");

	// Version 13
	file->Write("  accomodationFactor:"); file->WriteDouble(sh.accomodationFactor, "\n");

	file->Write("}\n");
}

// -----------------------------------------------------------

void Facet::DetectOrientation() {

	// Detect polygon orientation (clockwise or counter clockwise)
	// p= 1.0 => The second vertex is convex and vertex are counter clockwise.
	// p=-1.0 => The second vertex is concave and vertex are clockwise.
	// p= 0.0 => The polygon is not a simple one and orientation cannot be detected.

	POLYGON p;
	p.nbPts = sh.nbIndex;
	p.pts = vertices2;
	p.sign = 1.0;

	BOOL convexFound = FALSE;
	int i = 0;
	while (i < p.nbPts && !convexFound) {
		VERTEX2D c;
		BOOL empty = EmptyTriangle(&p, i - 1, i, i + 1, &c);
		if (empty || sh.nbIndex == 3) {
			int _i1 = IDX(i - 1, p.nbPts);
			int _i2 = IDX(i, p.nbPts);
			int _i3 = IDX(i + 1, p.nbPts);
			if (IsInPoly(c.u, c.v, p.pts, p.nbPts)) {
				convexFound = TRUE;
				// Orientation
				if (IsConvex(&p, i)) p.sign = 1.0;
				else                 p.sign = -1.0;
			}
		}
		i++;
	}

	if (!convexFound) {
		// Not a simple polygon
		sh.sign = 0.0;
	}
	else {
		sh.sign = p.sign;
	}

}

// -----------------------------------------------------------

void Facet::UpdateFlags() {

	sh.isProfile = (sh.profileType != REC_NONE);
	//sh.isOpaque = (sh.opacity != 0.0);
	sh.isTextured = ((texDimW*texDimH) > 0);

}

// -----------------------------------------------------------

int Facet::RestoreDeviceObjects() {

	// Initialize scene objects (OpenGL)
	if (sh.isTextured > 0) {
		glGenTextures(1, &glTex);
		glList = glGenLists(1);
	}

	BuildMeshList();
	BuildSelElemList();

	return GL_OK;

}

// -----------------------------------------------------------

int Facet::InvalidateDeviceObjects() {

	// Free all alocated resource (OpenGL)
	DELETE_TEX(glTex);
	DELETE_LIST(glList);
	DELETE_LIST(glElem);
	DELETE_LIST(glSelElem);
	return GL_OK;

}

// -----------------------------------------------------------

BOOL Facet::SetTexture(double width, double height, BOOL useMesh) {

	BOOL dimOK = (width*height > 0.0000001);

	if (dimOK) {
		sh.texWidthD = width;
		sh.texHeightD = height;
		sh.texWidth = (int)ceil(width - (double)1e-9); //double precision written to file
		sh.texHeight = (int)ceil(height - (double)1e-9); //double precision written to file
		dimOK = (sh.texWidth > 0 && sh.texHeight > 0);
	}
	else {
		sh.texWidth = 0;
		sh.texHeight = 0;
		sh.texWidthD = 0.0;
		sh.texHeightD = 0.0;
	}

	texDimW = 0;
	texDimH = 0;
	hasMesh = FALSE;
	SAFE_FREE(mesh);
	SAFE_FREE(dirCache);
	DELETE_TEX(glTex);
	DELETE_LIST(glList);
	DELETE_LIST(glElem);
	for (int i = 0; i < nbElem; i++)
		SAFE_FREE(meshPts[i].pts);
	SAFE_FREE(meshPts);
	nbElem = 0;
	UnselectElem();

	if (dimOK) {

		// Add a 1 texel border for bilinear filtering (rendering purpose)
		texDimW = GetPower2(sh.texWidth + 2);
		texDimH = GetPower2(sh.texHeight + 2);
		if (texDimW < 4) texDimW = 4;
		if (texDimH < 4) texDimH = 4;
		glGenTextures(1, &glTex);
		glList = glGenLists(1);
		if (useMesh)
			if (!BuildMesh()) return FALSE;
		if (sh.countDirection) {
			int dirSize = sh.texWidth*sh.texHeight*sizeof(VHIT);
			dirCache = (VHIT *)calloc(sh.texWidth*sh.texHeight, sizeof(VHIT));
			if (!dirCache) return FALSE;
			//memset(dirCache,0,dirSize); //already done by calloc
		}

	}

	UpdateFlags();
	return TRUE;
}

// -----------------------------------------------------------

void Facet::glVertex2u(double u, double v) {

	glVertex3d(sh.O.x + sh.U.x*u + sh.V.x*v,
		sh.O.y + sh.U.y*u + sh.V.y*v,
		sh.O.z + sh.U.z*u + sh.V.z*v);

}

// -----------------------------------------------------------

BOOL Facet::BuildMesh() {

	mesh = (SHELEM *)malloc(sh.texWidth * sh.texHeight * sizeof(SHELEM));
	if (!mesh) {
		//Couldn't allocate memory
		return FALSE;
		//throw Error("malloc failed on Facet::BuildMesh()");
	}
	meshPts = (MESH *)malloc(sh.texWidth * sh.texHeight * sizeof(MESH));
	if (!meshPts) {
		return FALSE;
	}
	hasMesh = TRUE;
	memset(mesh, 0, sh.texWidth * sh.texHeight * sizeof(SHELEM));
	memset(meshPts, 0, sh.texWidth * sh.texHeight * sizeof(MESH));

	POLYGON P1, P2;
	double sx, sy, A, tA;
	double iw = 1.0 / (double)sh.texWidthD;
	double ih = 1.0 / (double)sh.texHeightD;
	double rw = Norme(&sh.U) * iw;
	double rh = Norme(&sh.V) * ih;
	double *vList;
	double fA = iw*ih;
	int    nbv;

	P1.pts = (VERTEX2D *)malloc(4 * sizeof(VERTEX2D));
	_ASSERTE(P1.pts);
	P1.nbPts = 4;
	P1.sign = 1.0;
	P2.nbPts = sh.nbIndex;
	P2.pts = vertices2;
	P2.sign = -sh.sign;
	tA = 0.0;
	nbElem = 0;

	for (int j = 0; j < sh.texHeight; j++) {
		sy = (double)j;
		for (int i = 0; i < sh.texWidth; i++) {
			sx = (double)i;

			BOOL allInside = FALSE;
			double u0 = sx * iw;
			double v0 = sy * ih;
			double u1 = (sx + 1.0) * iw;
			double v1 = (sy + 1.0) * ih;
			float  uC, vC;
			mesh[i + j*sh.texWidth].elemId = -1;

			if (sh.nbIndex <= 4) {

				// Optimization for quad and triangle
				allInside = IsInPoly(u0, v0, vertices2, sh.nbIndex);
				allInside = allInside && IsInPoly(u0, v1, vertices2, sh.nbIndex);
				allInside = allInside && IsInPoly(u1, v0, vertices2, sh.nbIndex);
				allInside = allInside && IsInPoly(u1, v1, vertices2, sh.nbIndex);

			}

			if (!allInside) {

				// Intersect element with the facet (facet boundaries)
				P1.pts[0].u = u0;
				P1.pts[0].v = v0;
				P1.pts[1].u = u1;
				P1.pts[1].v = v0;
				P1.pts[2].u = u1;
				P1.pts[2].v = v1;
				P1.pts[3].u = u0;
				P1.pts[3].v = v1;
				A = GetInterArea(&P1, &P2, visible, &uC, &vC, &nbv, &vList);
				if (!IS_ZERO(A)) {

					if (A > (fA + 1e-10)) {

						// Polyon intersection error !
						// Switch back to brute force
						A = GetInterAreaBF(&P2, u0, v0, u1, v1, &uC, &vC);
						mesh[i + j*sh.texWidth].area = (float)(A*(rw*rh) / (iw*ih));
						mesh[i + j*sh.texWidth].uCenter = uC;
						mesh[i + j*sh.texWidth].vCenter = vC;
						mesh[i + j*sh.texWidth].full = IS_ZERO(fA - A);

					}
					else {

						// !! P1 and P2 are in u,v coordinates !!
						mesh[i + j*sh.texWidth].area = (float)(A*(rw*rh) / (iw*ih));
						mesh[i + j*sh.texWidth].uCenter = uC;
						mesh[i + j*sh.texWidth].vCenter = vC;
						mesh[i + j*sh.texWidth].full = IS_ZERO(fA - A);
						mesh[i + j*sh.texWidth].elemId = nbElem;

						// Mesh coordinates
						meshPts[nbElem].nbPts = nbv;
						meshPts[nbElem].pts = (VERTEX2D *)malloc(nbv * sizeof(VERTEX2D));
						_ASSERTE(meshPts[nbElem].pts);
						for (int n = 0; n < nbv; n++) {
							meshPts[nbElem].pts[n].u = vList[2 * n];
							meshPts[nbElem].pts[n].v = vList[2 * n + 1];
						}
						nbElem++;

					}

				}
				SAFE_FREE(vList);

			}
			else {

				mesh[i + j*sh.texWidth].area = (float)(rw*rh);
				mesh[i + j*sh.texWidth].uCenter = (float)(u0 + u1) / 2.0f;
				mesh[i + j*sh.texWidth].vCenter = (float)(v0 + v1) / 2.0f;
				mesh[i + j*sh.texWidth].full = TRUE;
				mesh[i + j*sh.texWidth].elemId = nbElem;

				// Mesh coordinates
				meshPts[nbElem].nbPts = 4;
				meshPts[nbElem].pts = (VERTEX2D *)malloc(4 * sizeof(VERTEX2D));
				_ASSERTE(meshPts[nbElem].pts);
				meshPts[nbElem].pts[0].u = u0;
				meshPts[nbElem].pts[0].v = v0;
				meshPts[nbElem].pts[1].u = u1;
				meshPts[nbElem].pts[1].v = v0;
				meshPts[nbElem].pts[2].u = u1;
				meshPts[nbElem].pts[2].v = v1;
				meshPts[nbElem].pts[3].u = u0;
				meshPts[nbElem].pts[3].v = v1;
				nbElem++;

			}

			tA += mesh[i + j*sh.texWidth].area;

		}
	}

	// Check meshing accuracy (TODO)
	/*
	int p = (int)(ceil(log10(sh.area)));
	double delta = pow(10.0,(double)(p-5));
	if( fabs(sh.area - tA)>delta ) {
	}
	*/

	free(P1.pts);
	BuildMeshList();
	return TRUE;
}

// -----------------------------------------------------------

void Facet::BuildMeshList() {

	if (!meshPts)
		return;

	DELETE_LIST(glElem);

	// Build OpenGL geometry for meshing
	glElem = glGenLists(1);
	glNewList(glElem, GL_COMPILE);


	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	for (int i = 0; i < nbElem; i++) {
		glBegin(GL_POLYGON);
		for (int n = 0; n < meshPts[i].nbPts; n++) {
			glEdgeFlag(TRUE);
			glVertex2u(meshPts[i].pts[n].u, meshPts[i].pts[n].v);
		}
		glEnd();

	}

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glEndList();

}

// -----------------------------------------------------------

void Facet::BuildSelElemList() {

	DELETE_LIST(glSelElem);
	int nbSel = 0;

	if (mesh && selectedElem.width != 0 && selectedElem.height != 0) {

		glSelElem = glGenLists(1);
		glNewList(glSelElem, GL_COMPILE);
		glColor3f(1.0f, 1.0f, 1.0f);
		glLineWidth(1.0f);
		glEnable(GL_LINE_SMOOTH);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glEnable(GL_POLYGON_OFFSET_LINE);
		glPolygonOffset(-1.0f, -1.0f);
		for (int i = 0; i < selectedElem.width; i++) {
			for (int j = 0; j < selectedElem.height; j++) {

				int add = (selectedElem.u + i) + (selectedElem.v + j)*sh.texWidth;
				int elId = mesh[add].elemId;
				if (elId >= 0) {
					glBegin(GL_POLYGON);
					for (int n = 0; n < meshPts[elId].nbPts; n++) {
						glEdgeFlag(TRUE);
						glVertex2u(meshPts[elId].pts[n].u, meshPts[elId].pts[n].v);
					}
					glEnd();
					nbSel++;
				}

			}
		}


		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glDisable(GL_POLYGON_OFFSET_LINE);
		glDisable(GL_LINE_SMOOTH);
		glEndList();

		// Empty selection
		if (nbSel == 0) UnselectElem();

	}

}

// -----------------------------------------------------------

void Facet::UnselectElem() {

	DELETE_LIST(glSelElem);
	selectedElem.width = 0;
	selectedElem.height = 0;

}

// -----------------------------------------------------------

void Facet::SelectElem(int u, int v, int width, int height) {

	UnselectElem();

	if (mesh && u >= 0 && u < sh.texWidth && v >= 0 && v < sh.texHeight) {

		int maxW = sh.texWidth - u;
		int maxH = sh.texHeight - v;
		selectedElem.u = u;
		selectedElem.v = v;
		selectedElem.width = MIN(maxW, width);
		selectedElem.height = MIN(maxH, height);
		BuildSelElemList();

	}

}

// -----------------------------------------------------------

void Facet::RenderSelectedElem() {
	if (glSelElem) glCallList(glSelElem);
}

// -----------------------------------------------------------

void Facet::Explode(FACETGROUP *group) {

	int nb = 0;
	if (!(group->facets = (Facet **)malloc(nbElem*sizeof(Facet *))))
		throw Error("Not enough memory to create new facets");
	for (int i = 0; i < nbElem; i++) {
		try {
			Facet *f = new Facet(meshPts[i].nbPts);
			f->Copy(this);
			group->facets[i] = f;
		}
		catch (...) {
			for (int d = 0; d < i; d++)
				SAFE_DELETE(group->facets[d]);
			throw Error("Cannot reserve memory for new facet(s)");
		}

		nb += meshPts[i].nbPts;
	}

	group->nbF = nbElem;
	group->nbV = nb;

}

// -----------------------------------------------------------

void Facet::FillVertexArray(VERTEX3D *v) {

	int nb = 0;
	for (int i = 0; i < nbElem; i++) {
		for (int j = 0; j < meshPts[i].nbPts; j++) {
			v[nb].x = sh.O.x + sh.U.x*meshPts[i].pts[j].u + sh.V.x*meshPts[i].pts[j].v;
			v[nb].y = sh.O.y + sh.U.y*meshPts[i].pts[j].u + sh.V.y*meshPts[i].pts[j].v;
			v[nb].z = sh.O.z + sh.U.z*meshPts[i].pts[j].u + sh.V.z*meshPts[i].pts[j].v;
			nb++;
		}
	}

}

// -----------------------------------------------------------

DWORD Facet::GetGeometrySize() {

	DWORD s = sizeof(SHFACET)
		+(sh.nbIndex * sizeof(int))
		+ (sh.nbIndex * sizeof(VERTEX2D));

	// Size of the 'element area' array passed to the geomotry buffer
	if (sh.isTextured) s += sizeof(AHIT)*sh.texWidth*sh.texHeight;
	if (sh.useOutgassingFile == 1) s += sizeof(double)*sh.outgassingMapHeight*sh.outgassingMapWidth;
	return s;
}

// -----------------------------------------------------------

DWORD Facet::GetHitsSize(int nbMoments) {

	return   sizeof(SHHITS)
		+(sh.texWidth*sh.texHeight*sizeof(AHIT))*(1 + nbMoments)
		+ (sh.isProfile ? (PROFILE_SIZE*sizeof(APROFILE)*(1 + nbMoments)) : 0)
		+ (sh.countDirection ? (sh.texWidth*sh.texHeight*sizeof(VHIT)*(1 + nbMoments)) : 0);

}

// -----------------------------------------------------------

DWORD Facet::GetTexSwapSize(BOOL useColormap) {

	DWORD tSize = texDimW*texDimH;
	if (useColormap) tSize = tSize * 4;
	return tSize;

}

// -----------------------------------------------------------

DWORD Facet::GetTexSwapSizeForRatio(double ratio, BOOL useColor) {

	double nU = Norme(&(sh.U));
	double nV = Norme(&(sh.V));
	double width = nU*ratio;
	double height = nV*ratio;

	BOOL dimOK = (width*height > 0.0000001);

	if (dimOK) {

		int iwidth = (int)ceil(width);
		int iheight = (int)ceil(height);
		int m = MAX(iwidth, iheight);
		int tDim = GetPower2(m);
		if (tDim < 16) tDim = 16;
		DWORD tSize = tDim*tDim;
		if (useColor) tSize = tSize * 4;
		return tSize;

	}
	else {

		return 0;

	}

}

// --------------------------------------------------------------------

DWORD Facet::GetNbCell() {
	return sh.texHeight * sh.texWidth;
}

// --------------------------------------------------------------------

DWORD Facet::GetNbCellForRatio(double ratio) {

	double nU = Norme(&(sh.U));
	double nV = Norme(&(sh.V));
	double width = nU*ratio;
	double height = nV*ratio;

	BOOL dimOK = (width*height > 0.0000001);

	if (dimOK) {
		int iwidth = (int)ceil(width);
		int iheight = (int)ceil(height);
		return iwidth*iheight;
	}
	else {
		return 0;
	}

}

// -----------------------------------------------------------

DWORD Facet::GetTexRamSize(size_t nbMoments) {

	size_t size = sizeof(AHIT)*nbMoments;
	if (mesh) size += sizeof(SHELEM);
	if (sh.countDirection) size += sizeof(VHIT)*nbMoments;
	return (sh.texWidth*sh.texHeight*size);

}

// -----------------------------------------------------------

DWORD Facet::GetTexRamSizeForRatio(double ratio, BOOL useMesh, BOOL countDir, size_t nbMoments) {

	double nU = Norme(&(sh.U));
	double nV = Norme(&(sh.V));
	double width = nU*ratio;
	double height = nV*ratio;

	BOOL dimOK = (width*height > 0.0000001);

	if (dimOK) {

		int iwidth = (int)ceil(width);
		int iheight = (int)ceil(height);
		int size = sizeof(AHIT)*nbMoments;
		if (useMesh) size += sizeof(SHELEM);
		if (countDir) size += sizeof(VHIT)*nbMoments;
		return iwidth * iheight * size;

	}
	else {

		return 0;

	}

}

// -----------------------------------------------------------

#define SUM_NEIGHBOR(i,j,we)                      \
	if( (i)>=0 && (i)<=w && (j)>=0 && (j)<=h ) {    \
	add = (i)+(j)*sh.texWidth;                    \
	if( mesh[add].area>0.0 ) {                   \
	if (textureMode==0) sum += we*(texBuffer[add].count*scaleF);          \
	  else if (textureMode==1) sum += we*(texBuffer[add].sum_1_per_speed*scaleF);          \
	  else if (textureMode==2) sum += we*(texBuffer[add].sum_v_ort_per_area*scaleF);          \
	  W=W+we;                                     \
	}                                             \
	}

double Facet::GetSmooth(int i, int j, AHIT *texBuffer, int textureMode, double scaleF) {

	float W = 0.0f;
	float sum = 0.0;
	int w = sh.texWidth - 1;
	int h = sh.texHeight - 1;
	int add;

	SUM_NEIGHBOR(i - 1, j - 1, 1.0f);
	SUM_NEIGHBOR(i - 1, j + 1, 1.0f);
	SUM_NEIGHBOR(i + 1, j - 1, 1.0f);
	SUM_NEIGHBOR(i + 1, j + 1, 1.0f);
	SUM_NEIGHBOR(i, j - 1, 2.0f);
	SUM_NEIGHBOR(i, j + 1, 2.0f);
	SUM_NEIGHBOR(i - 1, j, 2.0f);
	SUM_NEIGHBOR(i + 1, j, 2.0f);

	if (W == 0.0f)
		return 0.0f;
	else
		return sum / W;


}

// -----------------------------------------------------------
#define LOG10(x) log10f((float)x)

void Facet::BuildTexture(AHIT *texBuffer, int textureMode, double min, double max, BOOL useColorMap,
	double dCoeff1, double dCoeff2, double dCoeff3, BOOL doLog) {

	int size = sh.texWidth*sh.texHeight;
	int tSize = texDimW*texDimH;
	if (size == 0 || tSize == 0) return;

	double scaleFactor = 1.0;
	int val;

	glBindTexture(GL_TEXTURE_2D, glTex);
	if (useColorMap) {

		// -------------------------------------------------------
		// 16 Bit rainbow colormap
		// -------------------------------------------------------

		// Scale
		if (min < max) {
			if (doLog) {
				if (min < 1e-20) min = 1e-20;
				scaleFactor = 65534.0 / (log10(max) - log10(min)); // -1 for saturation color
			}
			else {
				scaleFactor = 65534.0 / (max - min); // -1 for saturation color
			}
		}
		else {
			doLog = FALSE;
			min = 0;
		}

		int *buff32 = (int *)malloc(tSize * 4);
		if (!buff32) throw Error("Out of memory in Facet::BuildTexture()");
		memset(buff32, 0, tSize * 4);
		for (int j = 0; j < sh.texHeight; j++) {
			for (int i = 0; i < sh.texWidth; i++) {
				int idx = i + j*sh.texWidth;
				double physicalValue;
				switch (textureMode) {
				case 0: //pressure
					physicalValue = texBuffer[idx].sum_v_ort_per_area*dCoeff1;
					break;
				case 1: //impingement rate
					physicalValue = (double)texBuffer[idx].count / (this->mesh[idx].area*(sh.is2sided ? 2.0 : 1.0))*dCoeff2;
					break;
				case 2: //particle density
					physicalValue = 2.0 * texBuffer[idx].sum_1_per_speed / (this->mesh[idx].area*(sh.is2sided ? 2.0 : 1.0))*dCoeff3;
					break;
				}
				if (doLog) {
					val = (int)((log10(physicalValue) - log10(min))*scaleFactor + 0.5);
				}
				else {
					val = (int)((physicalValue - min)*scaleFactor + 0.5);
				}
				SATURATE(val, 0, 65535);
				buff32[(i + 1) + (j + 1)*texDimW] = colorMap[val];
				if (texBuffer[idx].count == 0.0) buff32[(i + 1) + (j + 1)*texDimW] = (COLORREF)(65535 + 256 + 1); //show unset value as white
			}
		}

		/*
		// Perform edge smoothing (only with mesh)
		if( mesh ) {
		for(int j=-1;j<=sh.texHeight;j++) {
		for(int i=-1;i<=sh.texWidth;i++) {
		BOOL doSmooth = (i<0) || (i>=sh.texWidth) ||
		(j<0) || (j>=sh.texHeight) ||
		mesh[i+j*sh.texWidth].area==0.0f;
		if( doSmooth ) {
		if( doLog ) {
		val = (int)((log10(GetSmooth(i,j,texBuffer,dCoeff))-log10(min))*scaleFactor+0.5f);
		} else {
		val = (int)((GetSmooth(i,j,texBuffer,dCoeff)-min)*scaleFactor+0.5f);
		}
		SATURATE(val,0,65535);
		buff32[(i+1) + (j+1)*texDimW] = colorMap[val];
		}
		}
		}
		}
		*/



		glTexImage2D(
			GL_TEXTURE_2D,       // Type
			0,                   // No Mipmap
			GL_RGBA,             // Format RGBA
			texDimW,             // Width
			texDimH,             // Height
			0,                   // Border
			GL_RGBA,             // Format RGBA
			GL_UNSIGNED_BYTE,    // 8 Bit/pixel
			buff32               // Data
			);

		GLToolkit::CheckGLErrors("Facet::BuildTexture()");

		free(buff32);
	}
	else {

		// -------------------------------------------------------
		// 8 bit Luminance
		// -------------------------------------------------------
		if (min < max) {
			if (doLog) {
				if (min < 1e-20) min = 1e-20;
				scaleFactor = 255.0 / (log10(max) - log10(min)); // -1 for saturation color
			}
			else {
				scaleFactor = 255.0 / (max - min); // -1 for saturation color
			}
		}
		else {
			doLog = FALSE;
			min = 0;
		}

		unsigned char *buff8 = (unsigned char *)malloc(tSize*sizeof(unsigned char));
		if (!buff8) throw Error("Out of memory in Facet::BuildTexture()");
		memset(buff8, 0, tSize*sizeof(unsigned char));
		float fmin = (float)min;

		for (int j = 0; j < sh.texHeight; j++) {
			for (int i = 0; i < sh.texWidth; i++) {
				int idx = i + j*sh.texWidth;
				double physicalValue;
				switch (textureMode) {
				case 0: //pressure
					physicalValue = texBuffer[idx].sum_v_ort_per_area*dCoeff1;
					break;
				case 1: //impingement rate
					physicalValue = (double)texBuffer[idx].count / (this->mesh[idx].area*(sh.is2sided ? 2.0 : 1.0))*dCoeff2;
					break;
				case 2: //particle density
					physicalValue = 2.0 * texBuffer[idx].sum_1_per_speed / (this->mesh[idx].area*(sh.is2sided ? 2.0 : 1.0))*dCoeff2;
					break;
				}
				if (doLog) {
					val = (int)((log10(physicalValue) - log10(min))*scaleFactor + 0.5f);
				}
				else {
					val = (int)((physicalValue - min)*scaleFactor + 0.5f);
				}
				SATURATE(val, 0, 255);
				buff8[(i + 1) + (j + 1)*texDimW] = val;
			}
		}

		/*
		// Perform edge smoothing (only with mesh)
		if( mesh ) {
		for(int j=-1;j<=sh.texHeight;j++) {
		for(int i=-1;i<=sh.texWidth;i++) {
		BOOL doSmooth = (i<0) || (i>=sh.texWidth) ||
		(j<0) || (j>=sh.texHeight) ||
		mesh[i+j*sh.texWidth].area==0.0;
		if( doSmooth ) {
		if( doLog ) {
		val = (int)((LOG10(GetSmooth(i,j,texBuffer,dCoeff))-LOG10(min))*scaleFactor+0.5f);
		} else {
		val = (int)((GetSmooth(i,j,texBuffer,dCoeff)-min)*scaleFactor+0.5f);
		}
		SATURATE(val,0,255);
		buff8[(i+1) + (j+1)*texDimW] = val;
		}
		}
		}
		}*/

		glTexImage2D(
			GL_TEXTURE_2D,       // Type
			0,                   // No Mipmap
			GL_LUMINANCE,        // Format luminance
			texDimW,             // Width
			texDimH,             // Height
			0,                   // Border
			GL_LUMINANCE,        // Format luminance
			GL_UNSIGNED_BYTE,    // 8 Bit/pixel
			buff8                // Data
			);

		free(buff8);

	}
	GLToolkit::CheckGLErrors("Facet::BuildTexture()");
}

// -----------------------------------------------------------

void Facet::SwapNormal() {

	// Revert vertex order (around the second point)

	int *tmp = (int *)malloc(sh.nbIndex*sizeof(int));
	for (int i = sh.nbIndex, j = 0; i > 0; i--, j++)
		tmp[(i + 1) % sh.nbIndex] = GetIndex(j + 1);
	free(indices);
	indices = tmp;

	/* normal recalculated at reinitialize
	// Invert normal
	sh.N.x = -sh.N.x;
	sh.N.y = -sh.N.y;
	sh.N.z = -sh.N.z;*/

}

// -----------------------------------------------------------

void Facet::ShiftVertex() {

	// Shift vertex

	int *tmp = (int *)malloc(sh.nbIndex*sizeof(int));
	for (int i = 0; i < sh.nbIndex; i++)
		tmp[i] = GetIndex(i + 1);
	free(indices);
	indices = tmp;

}

// -----------------------------------------------------------

void Facet::InitVisibleEdge() {

	// Detect non visible edge (for polygon which contains holes)
	memset(visible, 0xFF, sh.nbIndex*sizeof(BOOL));

	for (int i = 0; i < sh.nbIndex; i++) {

		int p11 = GetIndex(i);
		int p12 = GetIndex(i + 1);

		for (int j = i + 1; j < sh.nbIndex; j++) {

			int p21 = GetIndex(j);
			int p22 = GetIndex(j + 1);

			if ((p11 == p22 && p12 == p21) || (p11 == p21 && p12 == p22)) {
				// Invisible edge found
				visible[i] = FALSE;
				visible[j] = FALSE;
			}

		}

	}

}

// -----------------------------------------------------------

BOOL Facet::IsCoplanar(Facet *f, double threshold) {

	// Detect if 2 facets are in the same plane (orientation preserving)
	// and have same parameters (used by collapse)

	return (fabs(a - f->a) < threshold) &&
		(fabs(b - f->b) < threshold) &&
		(fabs(c - f->c) < threshold) &&
		(fabs(d - f->d) < threshold) &&
		(sh.desorbType == f->sh.desorbType) &&
		(sh.sticking == f->sh.sticking) &&
		(sh.flow == f->sh.flow) &&
		(sh.opacity == f->sh.opacity) &&
		(sh.is2sided == f->sh.is2sided) &&
		(sh.reflectType == f->sh.reflectType) &&
		(sh.temperature == f->sh.temperature);
	//TODO: Add other properties!

}

// -----------------------------------------------------------

void Facet::Copy(Facet *f, BOOL copyMesh) {

	sh.sticking = f->sh.sticking;
	sh.opacity = f->sh.opacity;
	sh.area = f->sh.area;
	sh.desorbType = f->sh.desorbType;
	sh.desorbTypeN = f->sh.desorbTypeN;
	sh.reflectType = f->sh.reflectType;
	if (copyMesh) {
		sh.profileType = f->sh.profileType;
	}
	else {
		sh.profileType = REC_NONE;
	}
	sh.is2sided = f->sh.is2sided;
	//sh.bb = f->sh.bb;
	//sh.center = f->sh.center;

	//sh.counter = f->sh.counter;
	sh.flow = f->sh.flow;

	sh.mass = f->sh.mass;
	//sh.nbIndex = f->sh.nbIndex;
	//sh.nU = f->sh.nU;
	//sh.Nuv = f->sh.Nuv;
	//sh.nV = f->sh.nV;
	sh.superIdx = f->sh.superIdx;
	sh.superDest = f->sh.superDest;
	sh.teleportDest = f->sh.teleportDest;
	sh.temperature = f->sh.temperature;
	//sh.texHeight = f->sh.texHeight;
	//sh.texHeightD = f->sh.texHeightD;
	//sh.texWidth = f->sh.texWidth;
	//sh.texWidthD = f->sh.texWidthD;
	//sh.U = f->sh.U;
	//sh.V = f->sh.V;
	//dirCache = f->dirCache;
	if (copyMesh) {
		sh.countAbs = f->sh.countAbs;
		sh.countRefl = f->sh.countRefl;
		sh.countTrans = f->sh.countTrans;
		sh.countDes = f->sh.countDes;
		sh.countACD = f->sh.countACD;
		sh.countDirection = f->sh.countDirection;
		sh.isTextured = f->sh.isTextured;
		hasMesh = f->hasMesh;
		tRatio = f->tRatio;
	}
	this->UpdateFlags();
	//nbElem = f->nbElem;
	//texDimH = f->texDimH;
	//texDimW = f->texDimW;
	textureVisible = f->textureVisible;

	//visible = f->visible; //Dragons ahead!
	volumeVisible = f->volumeVisible;
	a = f->a;
	b = f->b;
	c = f->c;
	d = f->d;
	err = f->err;
	sh.N = f->sh.N;
	selected = f->selected;
}

// -----------------------------------------------------------

int Facet::GetIndex(int idx) {

	if (idx < 0) {
		return indices[(sh.nbIndex + idx) % sh.nbIndex];
	}
	else {
		return indices[idx % sh.nbIndex];
	}

}

void Facet::ConvertOldDesorbType() {
	if (sh.desorbType >= 3 && sh.desorbType <= 5) {
		sh.desorbTypeN = (double)(sh.desorbType - 1);
		sh.desorbType = DES_COSINE_N;
	}
}

void  Facet::SaveXML_geom(pugi::xml_node f){
	xml_node e=f.append_child("Sticking");
	e.append_attribute("constValue")= sh.sticking;
	e.append_attribute("parameterId")= sh.sticking_paramId;

	e = f.append_child("Opacity");
	e.append_attribute("constValue")= sh.opacity;
	e.append_attribute("parameterId")= sh.opacity_paramId;
	e.append_attribute("is2sided")= sh.is2sided;

	e = f.append_child("Outgassing");
	e.append_attribute("constValue") = sh.flow;
	e.append_attribute("parameterId") = sh.outgassing_paramId;
	e.append_attribute("desType") = sh.desorbType;
	e.append_attribute("desExponent") = sh.desorbTypeN;
	e.append_attribute("dynamicFromFile") = sh.useOutgassingFile;

	e = f.append_child("Temperature");
	e.append_attribute("value") = sh.temperature;
	e.append_attribute("accFactor") = sh.accomodationFactor;

	e = f.append_child("Reflection");
	e.append_attribute("type") = sh.reflectType;

	e = f.append_child("Structure");
	e.append_attribute("inStructure") = sh.superIdx;
	e.append_attribute("linksTo") = sh.superDest;

	e = f.append_child("Teleport");
	e.append_attribute("target") = sh.teleportDest;

	e = f.append_child("Motion");
	e.append_attribute("isMoving") = sh.isMoving;

	e = f.append_child("Recordings");
	xml_node t = e.append_child("Profile");
	t.append_attribute("type") = sh.profileType;
	switch (sh.profileType) {
	case 0:
		t.append_attribute("name") = "none";
		break;
	case 1:
		t.append_attribute("name") = "pressure u";
		break;
	case 2:
		t.append_attribute("name") = "pressure v";
		break;
	case 3:
		t.append_attribute("name") = "angular";
		break;
	case 4:
		t.append_attribute("name") = "speed";
		break;
	case 5:
		t.append_attribute("name") = "ortho.v";
		break;
	}
	t = e.append_child("Texture");
	t.append_attribute("hasMesh") = mesh != NULL;
	t.append_attribute("texDimX") = sh.texWidthD;
	t.append_attribute("texDimY") = sh.texHeightD;
	t.append_attribute("countDes") = sh.countDes;
	t.append_attribute("countAbs") = sh.countAbs;
	t.append_attribute("countRefl") = sh.countRefl;
	t.append_attribute("countTrans") = sh.countTrans;
	t.append_attribute("countDir") = sh.countDirection;
	t.append_attribute("countAC") = sh.countACD;
	
	e = f.append_child("ViewSettings");
	e.append_attribute("textureVisible") = textureVisible;
	e.append_attribute("volumeVisible") = volumeVisible;

	f.append_child("Indices").append_attribute("nb") = sh.nbIndex;
	for (size_t i = 0; i < sh.nbIndex; i++) {
		xml_node indice = f.child("Indices").append_child("Indice");
		indice.append_attribute("id") = i;
		indice.append_attribute("vertex") = indices[i];
	}
}