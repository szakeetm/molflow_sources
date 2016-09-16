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
#ifdef MOLFLOW
#include "MolFlow.h"
#endif

#ifdef SYNRAD
#include "SynRad.h"
#endif
#include "Facet.h"
#include "Utils.h"
#include <malloc.h>

#include <string.h>
#include <math.h>
#include "GLApp/GLToolkit.h"
#include <sstream>

using namespace pugi;

#define MAX(x,y) (((x)<(y))?(y):(x))
#define MIN(x,y) (((x)<(y))?(x):(y))

// Colormap stuff
extern COLORREF rainbowCol[]; //defined in GLGradient.cpp



#ifdef MOLFLOW
extern MolFlow *mApp;
#endif

#ifdef SYNRAD
extern SynRad*mApp;
#endif
static int colorMap[65536];
static BOOL colorInited = FALSE;

// -----------------------------------------------------------

Facet::Facet(int nbIndex) {

	indices = (int *)malloc(nbIndex*sizeof(int));                    // Ref to Geometry VERTEX3D
	vertices2 = (VERTEX2D *)malloc(nbIndex * sizeof(VERTEX2D));      // Local U,V coordinates
	memset(vertices2, 0, nbIndex * sizeof(VERTEX2D));

	sh.nbIndex = nbIndex;

	//ResizeCounter(nbMoments);

	/*memset(&sh.counter, 0, sizeof(sh.counter));
	sh.counter.hit.nbDesorbed = 0;
	sh.counter.hit.nbAbsorbed = 0;
	sh.counter.hit.nbHit = 0;*/
	memset(&counterCache, 0, sizeof(SHHITS));


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

	sh.enableSojournTime = FALSE;
	sh.sojournFreq = 1E13;
	sh.sojournE = 100;

	sh.outgassing_paramId = -1;
	sh.opacity_paramId = -1;
	sh.sticking_paramId = -1;

	sh.isMoving = FALSE;

	hasOutgassingFile = FALSE;
	outgassingMap = NULL;
	totalFlux = sh.totalOutgassing = totalDose = 0.0;

	textureVisible = TRUE;
	volumeVisible = TRUE;

	texDimW = 0;
	texDimH = 0;
	tRatio = 0.0;

	//mesh = NULL;
	//meshPts = NULL;
	cellPropertiesIds = NULL;
	meshvector = NULL;
	meshvectorsize = 0;
	hasMesh = FALSE;
	//nbElem = 0;
	selectedElem.u = 0;
	selectedElem.v = 0;
	selectedElem.width = 0;
	selectedElem.height = 0;
	dirCache = NULL;
	textureError = FALSE;
	hasOutgassingFile = FALSE;

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



Facet::~Facet() {
	SAFE_FREE(indices);
	SAFE_FREE(vertices2);
	SAFE_FREE(cellPropertiesIds);
	SAFE_FREE(dirCache);
	DELETE_TEX(glTex);
	DELETE_LIST(glList);
	DELETE_LIST(glElem);
	DELETE_LIST(glSelElem);
	SAFE_FREE(visible);
 	for (size_t i = 0; i < meshvectorsize; i++)
 		SAFE_FREE(meshvector[i].points);
 	SAFE_FREE(meshvector);
	SAFE_FREE(outgassingMap);
}


/*
void Facet::ResetCounter() {
	SHHITS zeroes;memset(&zeroes, 0, sizeof(zeroes)); //A new zero-value vector
	std::fill(counter.begin(), counter.end(), zeroes); //Initialize each moment with 0 values
}

void Facet::ResizeCounter(size_t nbMoments) {
	SHHITS zeroes;memset(&zeroes, 0, sizeof(zeroes)); //A new zero-value vector
	counter = std::vector<SHHITS>(nbMoments + 1); //Reserve a counter for each moment, plus an extra for const. flow
	std::fill(counter.begin(), counter.end(), zeroes); //Initialize each moment with 0 values
}
*/


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
	counterCache.hit.nbAbsorbed = file->ReadLLong();

	file->ReadKeyword("nbDes"); file->ReadKeyword(":");
	counterCache.hit.nbDesorbed = file->ReadLLong();

	file->ReadKeyword("nbHit"); file->ReadKeyword(":");

	counterCache.hit.nbHit = file->ReadLLong();
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

void Facet::LoadXML(xml_node f, int nbVertex, BOOL isMolflowFile, int vertexOffset) {
	int idx = 0;
	for (xml_node indice : f.child("Indices").children("Indice")) {
		indices[idx] = indice.attribute("vertex").as_int() + vertexOffset;
		if (indices[idx] >= nbVertex) {
			char err[128];
			sprintf(err, "Facet %d refers to vertex %d which doesn't exist", f.attribute("id").as_int() + 1, idx + 1);
			throw Error(err);
		}
		idx++;
	}
	sh.opacity = f.child("Opacity").attribute("constValue").as_double();
	sh.is2sided = f.child("Opacity").attribute("is2sided").as_int();
	sh.superIdx = f.child("Structure").attribute("inStructure").as_int();
	sh.superDest = f.child("Structure").attribute("linksTo").as_int();
	sh.teleportDest = f.child("Teleport").attribute("target").as_int();

	if (isMolflowFile) {
		sh.sticking = f.child("Sticking").attribute("constValue").as_double();
		sh.sticking_paramId = f.child("Sticking").attribute("parameterId").as_int();
		sh.opacity_paramId = f.child("Opacity").attribute("parameterId").as_int();
		sh.flow = f.child("Outgassing").attribute("constValue").as_double();
		sh.desorbType = f.child("Outgassing").attribute("desType").as_int();
		sh.desorbTypeN = f.child("Outgassing").attribute("desExponent").as_double();
		sh.outgassing_paramId = f.child("Outgassing").attribute("parameterId").as_int();
		hasOutgassingFile = f.child("Outgassing").attribute("hasOutgassingFile").as_int();
		sh.useOutgassingFile = f.child("Outgassing").attribute("useOutgassingFile").as_int();
		sh.temperature = f.child("Temperature").attribute("value").as_double();
		sh.accomodationFactor = f.child("Temperature").attribute("accFactor").as_double();
		xml_node reflNode = f.child("Reflection");
		sh.reflectType = reflNode.attribute("type").as_int();
		if (reflNode.attribute("enableSojournTime")) {
			sh.enableSojournTime = reflNode.attribute("enableSojournTime").as_bool();
			if (!reflNode.attribute("sojournFreq")) {//Backward compatibility with ver. before 2.6.25
				sh.sojournFreq = 1.0 / reflNode.attribute("sojournTheta0").as_double();
				sh.sojournE = 8.31 * reflNode.attribute("sojournE").as_double();
			}
			else {
				sh.sojournFreq = reflNode.attribute("sojournFreq").as_double();
				sh.sojournE = reflNode.attribute("sojournE").as_double();
			}
		}
		else {
			//Already set to default when calling Molflow::LoadFile()
		}
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

		xml_node outgNode = f.child("DynamicOutgassing");
		if ((hasOutgassingFile) && outgNode && outgNode.child("map")) {
			sh.outgassingMapWidth = outgNode.attribute("width").as_int();
			sh.outgassingMapHeight = outgNode.attribute("height").as_int();
			sh.outgassingFileRatio = outgNode.attribute("ratio").as_double();
			totalDose = outgNode.attribute("totalDose").as_double();
			sh.totalOutgassing = outgNode.attribute("totalOutgassing").as_double();
			totalFlux = outgNode.attribute("totalFlux").as_double();

			std::stringstream outgText;
			outgText << outgNode.child_value("map");
			outgassingMap = (double*)malloc(sh.outgassingMapWidth*sh.outgassingMapHeight*sizeof(double));

			for (int iy = 0; iy < sh.outgassingMapHeight; iy++) {
				for (int ix = 0; ix < sh.outgassingMapWidth; ix++) {
					outgText >> outgassingMap[iy*sh.outgassingMapWidth + ix];
				}
			}

		}
		else hasOutgassingFile = sh.useOutgassingFile = 0; //if outgassing map was incorrect, don't use it
	} //else use default values at Facet() constructor

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

	if (version >= 9) { //new reflection model
		file->ReadKeyword("reflectType"); file->ReadKeyword(":");
		sh.reflectType = REF_DIFFUSE; int reflType = file->ReadInt(); //Discard Synrad diffuse
		file->ReadKeyword("sticking"); file->ReadKeyword(":");
		sh.sticking = 0; file->ReadDouble(); //Discard Synrad sticking

		if (reflType >= 2) { //Material reflection: update index from the material's name
			file->ReadKeyword("materialName"); file->ReadKeyword(":"); file->ReadWord();
		}
		file->ReadKeyword("doScattering"); file->ReadKeyword(":");
		file->ReadInt();
		file->ReadKeyword("rmsRoughness"); file->ReadKeyword(":");
		file->ReadDouble();
		file->ReadKeyword("autoCorrLength"); file->ReadKeyword(":");
		file->ReadDouble();
		file->ReadKeyword("opacity"); file->ReadKeyword(":");
		sh.opacity = file->ReadDouble();
	}
	else { //legacy reflection model
		file->ReadKeyword("sticking"); file->ReadKeyword(":");
		sh.sticking = 0; file->ReadDouble(); //Discard Synrad sticking
		if (version >= 4) {
			file->ReadKeyword("roughness"); file->ReadKeyword(":");
			file->ReadDouble(); //roughness
		}
		file->ReadKeyword("opacity"); file->ReadKeyword(":");
		sh.opacity = file->ReadDouble();
		file->ReadKeyword("reflectType"); file->ReadKeyword(":");
		sh.reflectType = REF_DIFFUSE; file->ReadInt(); //Discard Synrad diffuse
	}

	file->ReadKeyword("profileType"); file->ReadKeyword(":");
	sh.profileType = 0; file->ReadInt(); //Discard Synrad profile
	file->ReadKeyword("hasSpectrum"); file->ReadKeyword(":");
	file->ReadInt();
	file->ReadKeyword("superDest"); file->ReadKeyword(":");
	sh.superDest = file->ReadInt();
	file->ReadKeyword("superIdx"); file->ReadKeyword(":");
	sh.superIdx = file->ReadInt();
	file->ReadKeyword("is2sided"); file->ReadKeyword(":");
	sh.is2sided = file->ReadInt();
	file->ReadKeyword("mesh"); file->ReadKeyword(":");
	hasMesh = FALSE; file->ReadInt(); //Discard synrad texture
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
	counterCache.hit.nbAbsorbed = 0; file->ReadLLong();
	if (version < 3) {
		file->ReadKeyword("nbDes"); file->ReadKeyword(":");
		counterCache.hit.nbDesorbed = 0;
		file->ReadLLong();
	}
	file->ReadKeyword("nbHit"); file->ReadKeyword(":");
	counterCache.hit.nbHit = 0; file->ReadLLong();
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
	counterCache.hit.nbDesorbed = (llong)(file->ReadDouble() + 0.5);
	counterCache.hit.nbHit = (llong)(file->ReadDouble() + 0.5);
	counterCache.hit.nbAbsorbed = (llong)(file->ReadDouble() + 0.5);
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

	if (counterCache.hit.nbDesorbed == 0)
		sh.desorbType = DES_NONE;

	if (IsLinkFacet()) {
		sh.superDest = (int)(sh.sticking + 0.5);
		sh.sticking = 0;
	}

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
	file->Write("  mesh:"); file->WriteInt((cellPropertiesIds != NULL), "\n");


	file->Write("  outgassing:"); file->WriteDouble(sh.flow*10.00, "\n"); //Pa*m3/s -> mbar*l/s for compatibility with old versions
	file->Write("  texDimX:"); file->WriteDouble(sh.texWidthD, "\n");
	file->Write("  texDimY:"); file->WriteDouble(sh.texHeightD, "\n");


	file->Write("  countDes:"); file->WriteInt(sh.countDes, "\n");
	file->Write("  countAbs:"); file->WriteInt(sh.countAbs, "\n");
	file->Write("  countRefl:"); file->WriteInt(sh.countRefl, "\n");
	file->Write("  countTrans:"); file->WriteInt(sh.countTrans, "\n");
	file->Write("  acMode:"); file->WriteInt(sh.countACD, "\n");
	file->Write("  nbAbs:"); file->WriteLLong(0/*sh.counter.hit.nbAbsorbed*/, "\n");
	file->Write("  nbDes:"); file->WriteLLong(0/*sh.counter.hit.nbDesorbed*/, "\n");
	file->Write("  nbHit:"); file->WriteLLong(0/*sh.counter.hit.nbHit*/, "\n");

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



// -----------------------------------------------------------

void Facet::UpdateFlags() {

	sh.isProfile = (sh.profileType != REC_NONE);
	//sh.isOpaque = (sh.opacity != 0.0);
	sh.isTextured = ((texDimW*texDimH) > 0);

}





size_t Facet::GetGeometrySize() {

	size_t s = sizeof(SHFACET)
		+ (sh.nbIndex * sizeof(int))
		+ (sh.nbIndex * sizeof(VERTEX2D));

	// Size of the 'element area' array passed to the geometry buffer
	if (sh.isTextured) s += sizeof(AHIT)*sh.texWidth*sh.texHeight;
	if (sh.useOutgassingFile == 1) s += sizeof(double)*sh.outgassingMapHeight*sh.outgassingMapWidth;
	return s;

}

// -----------------------------------------------------------

size_t Facet::GetHitsSize(size_t nbMoments) {

	return   (1+nbMoments)*(
		sizeof(SHHITS)+
		+ (sh.texWidth*sh.texHeight*sizeof(AHIT))
		+ (sh.isProfile ? (PROFILE_SIZE*sizeof(APROFILE)) : 0)
		+ (sh.countDirection ? (sh.texWidth*sh.texHeight*sizeof(VHIT)) : 0)
		);


}





size_t Facet::GetTexRamSize(size_t nbMoments) {

	int sizePerCell = 220; //estimation
	int sizePerFacet = 9000; //estimation
	return (sh.texWidth*sh.texHeight*sizePerCell + sh.isTextured*sizePerFacet);

}

size_t Facet::GetTexRamSizeForRatio(double ratio, BOOL useMesh, BOOL countDir, size_t nbMoments) {
	double nU = Norme(sh.U);
	double nV = Norme(sh.V);
	double width = nU*ratio;
	double height = nV*ratio;

	BOOL dimOK = (width*height > 0.0000001);

	if (dimOK) {
		int iwidth = (int)ceil(width);
		int iheight = (int)ceil(height);
		/*int size = 2*sizeof(double)+sizeof(llong);
		if(useMesh) size += sizeof(SHELEM)+sizeof(MESH);
		if(countDir) size += sizeof(VHIT);*/
		int size = 220; //estimation
		return iwidth * iheight * size + 9000;
	}
	else {
		return 0;
	}
}

#define SUM_NEIGHBOR(i,j,we)                      \
	if( (i)>=0 && (i)<=w && (j)>=0 && (j)<=h ) {    \
	add = (i)+(j)*sh.texWidth;                    \
	if( GetMeshArea(add)>0.0 ) {                   \
	if (textureMode==0) sum += we*(texBuffer[add].count*scaleF);          \
				  		  	  	  else if (textureMode==1) sum += we*(texBuffer[add].sum_1_per_ort_velocity*scaleF);          \
				  		  	  	  else if (textureMode==2) sum += we*(texBuffer[add].sum_v_ort_per_area*scaleF);          \
		  W=W+we;                                     \
				}                                             \
				}


double Facet::GetSmooth(int i, int j, AHIT *texBuffer, int textureMode, double scaleF) {

	double W = 0.0f;
	double sum = 0.0;
	int w = sh.texWidth - 1;
	int h = sh.texHeight - 1;
	int add;

	SUM_NEIGHBOR(i - 1, j - 1, 1.0);
	SUM_NEIGHBOR(i - 1, j + 1, 1.0);
	SUM_NEIGHBOR(i + 1, j - 1, 1.0);
	SUM_NEIGHBOR(i + 1, j + 1, 1.0);
	SUM_NEIGHBOR(i, j - 1, 2.0);
	SUM_NEIGHBOR(i, j + 1, 2.0);
	SUM_NEIGHBOR(i - 1, j, 2.0);
	SUM_NEIGHBOR(i + 1, j, 2.0);

	if (W == 0.0f)
		return 0.0f;
	else
		return sum / W;
}

// -----------------------------------------------------------
#define LOG10(x) log10f((float)x)

void Facet::BuildTexture(AHIT *texBuffer, int textureMode, double min, double max, BOOL useColorMap,
	double dCoeff1, double dCoeff2, double dCoeff3, BOOL doLog, size_t m) {
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
					physicalValue = (double)texBuffer[idx].count / (this->GetMeshArea(idx)*(sh.is2sided ? 2.0 : 1.0))*dCoeff2;
					break;
				case 2: //particle density
					physicalValue = texBuffer[idx].sum_1_per_ort_velocity / (this->GetMeshArea(idx)*(sh.is2sided ? 2.0 : 1.0))*dCoeff3;

					//Correction for double-density effect (measuring density on desorbing/absorbing facets):
					if (counterCache.hit.nbHit>0 || counterCache.hit.nbDesorbed>0)
						if (counterCache.hit.nbAbsorbed > 0 || counterCache.hit.nbDesorbed > 0) //otherwise save calculation time
							physicalValue *= 1.0 - ((double)counterCache.hit.nbAbsorbed + (double)counterCache.hit.nbDesorbed) / ((double)counterCache.hit.nbHit + (double)counterCache.hit.nbDesorbed) / 2.0;

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


		GLint width, height, format;
		glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, &width);
		glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, &height);
		glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_INTERNAL_FORMAT, &format);
		if (format == GL_RGBA && width == texDimW && height == texDimH) {
			//Update texture
			glTexSubImage2D(
				GL_TEXTURE_2D,       // Type


				0,                   // No Mipmap
				0,					// X offset
				0,					// Y offset
				texDimW,             // Width
				texDimH,             // Height

				GL_RGBA,             // Format RGBA
				GL_UNSIGNED_BYTE,    // 8 Bit/pixel
				buff32              // Data
			);
		}
		else {
			//Rebuild texture
			glTexImage2D(
				GL_TEXTURE_2D,       // Type
				0,                   // No Mipmap
				GL_RGBA,             // Format RGBA
				texDimW,             // Width
				texDimH,             // Height
				0,                   // Border
				GL_RGBA,             // Format RGBA
				GL_UNSIGNED_BYTE,    // 8 Bit/pixel
				buff32              // Data
			);
		}
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
					physicalValue = (double)texBuffer[idx].count / (this->GetMeshArea(idx)*(sh.is2sided ? 2.0 : 1.0))*dCoeff2;
					break;
				case 2: //particle density
					physicalValue = texBuffer[idx].sum_1_per_ort_velocity / (this->GetMeshArea(idx)*(sh.is2sided ? 2.0 : 1.0))*dCoeff3;

					//Correction for double-density effect (measuring density on desorbing/absorbing facets):
					if (counterCache.hit.nbHit>0 || counterCache.hit.nbDesorbed>0)
						if (counterCache.hit.nbAbsorbed > 0 || counterCache.hit.nbDesorbed > 0) //otherwise save calculation time
							physicalValue *= 1.0 - ((double)counterCache.hit.nbAbsorbed + (double)counterCache.hit.nbDesorbed) / ((double)counterCache.hit.nbHit + (double)counterCache.hit.nbDesorbed) / 2.0;

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


		GLint width, height, format;
		glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, &width);
		glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, &height);
		glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_INTERNAL_FORMAT, &format);
		if (format == GL_LUMINANCE && width == texDimW && height == texDimH) {
			//Update texture
			glTexSubImage2D(
				GL_TEXTURE_2D,       // Type
				0,                   // No Mipmap
				0,					// X offset
				0,					// Y offset
				texDimW,             // Width
				texDimH,             // Height
				GL_LUMINANCE,         // Format RGBA
				GL_UNSIGNED_BYTE,    // 8 Bit/pixel
				buff8                // Data
			);
		}



		else {
			//Rebuild texture
			glTexImage2D(
				GL_TEXTURE_2D,       // Type

				0,                   // No Mipmap
				GL_LUMINANCE,         // Format RGBA
				texDimW,             // Width
				texDimH,             // Height

				0,                   // Border
				GL_LUMINANCE,         // Format RGBA
				GL_UNSIGNED_BYTE,    // 8 Bit/pixel
				buff8                // Data
			);
		}
		free(buff8);

	}
	GLToolkit::CheckGLErrors("Facet::BuildTexture()");
}

BOOL Facet::IsCoplanarAndEqual(Facet *f, double threshold) {

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



void Facet::ConvertOldDesorbType() {
	if (sh.desorbType >= 3 && sh.desorbType <= 5) {
		sh.desorbTypeN = (double)(sh.desorbType - 1);
		sh.desorbType = DES_COSINE_N;
	}
}

void  Facet::SaveXML_geom(pugi::xml_node f){
	xml_node e = f.append_child("Sticking");
	e.append_attribute("constValue") = sh.sticking;
	e.append_attribute("parameterId") = sh.sticking_paramId;

	e = f.append_child("Opacity");
	e.append_attribute("constValue") = sh.opacity;
	e.append_attribute("parameterId") = sh.opacity_paramId;
	e.append_attribute("is2sided") = sh.is2sided;

	e = f.append_child("Outgassing");
	e.append_attribute("constValue") = sh.flow;
	e.append_attribute("parameterId") = sh.outgassing_paramId;
	e.append_attribute("desType") = sh.desorbType;
	e.append_attribute("desExponent") = sh.desorbTypeN;
	e.append_attribute("hasOutgassingFile") = hasOutgassingFile;
	e.append_attribute("useOutgassingFile") = sh.useOutgassingFile;

	e = f.append_child("Temperature");
	e.append_attribute("value") = sh.temperature;
	e.append_attribute("accFactor") = sh.accomodationFactor;

	e = f.append_child("Reflection");
	e.append_attribute("type") = sh.reflectType;
	e.append_attribute("enableSojournTime") = sh.enableSojournTime;
	e.append_attribute("sojournFreq") = sh.sojournFreq;
	e.append_attribute("sojournE") = sh.sojournE;

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
	t.append_attribute("hasMesh") = cellPropertiesIds != NULL;
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


	if (hasOutgassingFile){
		xml_node textureNode = f.append_child("DynamicOutgassing");
		textureNode.append_attribute("width") = sh.outgassingMapWidth;
		textureNode.append_attribute("height") = sh.outgassingMapHeight;
		textureNode.append_attribute("ratio") = sh.outgassingFileRatio;
		textureNode.append_attribute("totalDose") = totalDose;
		textureNode.append_attribute("totalOutgassing") = sh.totalOutgassing;
		textureNode.append_attribute("totalFlux") = totalFlux;

		std::stringstream outgText;
		outgText << '\n'; //better readability in file
		for (int iy = 0; iy < sh.outgassingMapHeight; iy++) {
			for (int ix = 0; ix < sh.outgassingMapWidth; ix++) {
				outgText << outgassingMap[iy*sh.outgassingMapWidth + ix] << '\t';
			}
			outgText << '\n';
		}
		textureNode.append_child("map").append_child(node_cdata).set_value(outgText.str().c_str());

	} //end texture
}