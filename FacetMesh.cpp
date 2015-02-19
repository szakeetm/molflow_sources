/*
File:        FacetMesh.cpp
Description: Facet mesh configuration dialog
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

#include "FacetMesh.h"
#include "GLApp/GLLabel.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLMessageBox.h"
#include "Utils.h" 
#include "MolFlow.h"
extern MolFlow *mApp;

//-----------------------------------------------------------------------------

FacetMesh::FacetMesh(Worker *w):GLWindow() {

	worker = w;
	geom = w->GetGeometry();

	SetIconfiable(TRUE);

	int wD = 320;
	int hD = 487;
	iPanel = new GLTitledPanel("Facet Dimesion");
	iPanel->SetBounds(3, 1, 309, 45);
	Add(iPanel);
	aPanel = new GLTitledPanel("Texture properties");
	aPanel->SetBounds(3, 52, 309, 119);
	Add(aPanel);
	mPanel = new GLTitledPanel("Texture cell / memory");
	mPanel->SetBounds(3, 177, 309, 44);
	Add(mPanel);
	vPanel = new GLTitledPanel("View settings");
	vPanel->SetBounds(3, 227, 309, 44);
	Add(vPanel);
	desPanel = new GLTitledPanel("Dynamic desorption");
	desPanel->SetBounds(3, 415, 309, 44);
	Add(desPanel);
	paramPanel = new GLTitledPanel("Additional parameters");
	paramPanel->SetBounds(3, 277, 309, 132);
	Add(paramPanel);
	l1 = new GLLabel("\201 length:");
	iPanel->SetCompBounds(l1, 19, 18, 40, 12);
	iPanel->Add(l1);

	l2 = new GLLabel("\202 length:");
	iPanel->SetCompBounds(l2, 153, 18, 41, 12);
	iPanel->Add(l2);

	vLength = new GLTextField(0, "");
	vLength->SetEditable(FALSE);
	iPanel->SetCompBounds(vLength, 196, 17, 72, 18);
	iPanel->Add(vLength);

	uLength = new GLTextField(0, "");
	uLength->SetEditable(FALSE);
	iPanel->SetCompBounds(uLength, 65, 17, 72, 18);
	iPanel->Add(uLength);

	lengthText = new GLTextField(0, "");
	aPanel->SetCompBounds(lengthText, 196, 35, 72, 18);
	aPanel->Add(lengthText);

	perCm = new GLLabel("cells/cm");
	aPanel->SetCompBounds(perCm, 139, 38, 40, 12);
	aPanel->Add(perCm);

	resolutionText = new GLTextField(0, "");
	aPanel->SetCompBounds(resolutionText, 75, 35, 62, 18);
	aPanel->Add(resolutionText);

	l5 = new GLLabel("Resolution:");
	aPanel->SetCompBounds(l5, 13, 38, 52, 12);
	aPanel->Add(l5);

	enableBtn = new GLToggle(0, "Enable texture");
	aPanel->SetCompBounds(enableBtn, 9, 18, 83, 16);
	aPanel->Add(enableBtn);

	recordDesBtn = new GLToggle(0, "Count desorption");
	aPanel->SetCompBounds(recordDesBtn, 8, 62, 94, 16);
	aPanel->Add(recordDesBtn);

	perCell = new GLLabel("cm/cell");
	aPanel->SetCompBounds(perCell, 268, 38, 35, 12);
	aPanel->Add(perCell);

	recordDirBtn = new GLToggle(0, "Record direction vectors");
	aPanel->SetCompBounds(recordDirBtn, 166, 100, 125, 16);
	aPanel->Add(recordDirBtn);

	recordTransBtn = new GLToggle(0, "Count transparent pass");
	aPanel->SetCompBounds(recordTransBtn, 166, 81, 120, 16);
	aPanel->Add(recordTransBtn);

	recordReflBtn = new GLToggle(0, "Count reflection");
	aPanel->SetCompBounds(recordReflBtn, 166, 62, 89, 16);
	aPanel->Add(recordReflBtn);

	recordACBtn = new GLToggle(0, "Angular coefficient");
	aPanel->SetCompBounds(recordACBtn, 8, 100, 101, 16);
	aPanel->Add(recordACBtn);

	recordAbsBtn = new GLToggle(0, "Count absorption");
	aPanel->SetCompBounds(recordAbsBtn, 8, 81, 94, 16);
	aPanel->Add(recordAbsBtn);

	showTexture = new GLToggle(0, "Draw Texture");
	vPanel->SetCompBounds(showTexture, 9, 19, 80, 16);
	vPanel->Add(showTexture);

	showVolume = new GLToggle(0, "Draw Volume");
	vPanel->SetCompBounds(showVolume, 113, 19, 81, 16);
	vPanel->Add(showVolume);

	cellText = new GLTextField(0, "");
	mPanel->SetCompBounds(cellText, 195, 19, 107, 18);
	cellText->SetEditable(FALSE);
	mPanel->Add(cellText);

	l8 = new GLLabel("Cells:");
	mPanel->SetCompBounds(l8, 166, 22, 29, 12);
	mPanel->Add(l8);

	ramText = new GLTextField(0, "");
	ramText->SetEditable(FALSE);
	mPanel->SetCompBounds(ramText, 58, 19, 100, 18);
	mPanel->Add(ramText);

	l7 = new GLLabel("Memory:");
	mPanel->SetCompBounds(l7, 10, 22, 43, 12);
	mPanel->Add(l7);

	quickApply = new GLButton(0, "Quick Apply");
	vPanel->SetCompBounds(quickApply, 203, 16, 99, 20);
	vPanel->Add(quickApply);

	fileDesText = new GLTextField(0, "");
	desPanel->SetCompBounds(fileDesText, 200, 19, 63, 18);
	fileDesText->SetEditable(FALSE);
	desPanel->Add(fileDesText);

	label3 = new GLLabel("mbar.l/s");
	desPanel->SetCompBounds(label3, 264, 22, 39, 12);
	desPanel->Add(label3);

	label1 = new GLLabel("Desorption:");
	desPanel->SetCompBounds(label1, 142, 22, 53, 12);
	desPanel->Add(label1);

	label2 = new GLLabel("Use file:");
	desPanel->SetCompBounds(label2, 10, 22, 39, 12);
	desPanel->Add(label2);

	facetMovingToggle = new GLToggle(0, "Moving part");
	paramPanel->SetCompBounds(facetMovingToggle, 10, 109, 74, 16);
	paramPanel->Add(facetMovingToggle);

	facetSuperDest = new GLTextField(0, "");
	paramPanel->SetCompBounds(facetSuperDest, 182, 88, 81, 18);
	paramPanel->Add(facetSuperDest);

	label8 = new GLLabel("Link to:");
	paramPanel->SetCompBounds(label8, 142, 90, 35, 12);
	paramPanel->Add(label8);

	facetStructure = new GLTextField(0, "");
	paramPanel->SetCompBounds(facetStructure, 60, 88, 76, 18);
	paramPanel->Add(facetStructure);

	label7 = new GLLabel("Structure:");
	paramPanel->SetCompBounds(label7, 10, 90, 46, 12);
	paramPanel->Add(label7);

	facetTeleport = new GLTextField(0, "");
	paramPanel->SetCompBounds(facetTeleport, 141, 66, 122, 18);
	paramPanel->Add(facetTeleport);

	label4 = new GLLabel("Teleport to facet:");
	paramPanel->SetCompBounds(label4, 10, 68, 74, 12);
	paramPanel->Add(label4);

	label5 = new GLLabel("Temp.accom. coefficient:");
	paramPanel->SetCompBounds(label5, 10, 46, 126, 13);
	paramPanel->Add(label5);

	label6 = new GLLabel("Reflection:");
	paramPanel->SetCompBounds(label6, 10, 21, 50, 12);
	paramPanel->Add(label6);

	facetReflType = new GLCombo(0);
	paramPanel->SetCompBounds(facetReflType, 141, 19, 122, 20);
	paramPanel->Add(facetReflType);
	facetReflType->SetSize(3);
	facetReflType->SetValueAt(0, "Diffuse");
	facetReflType->SetValueAt(1, "Mirror");
	facetReflType->SetValueAt(2, "Uniform");

	facetUseDesFile = new GLCombo(0);
	desPanel->SetCompBounds(facetUseDesFile, 55, 18, 82, 20);
	desPanel->Add(facetUseDesFile);
	facetUseDesFile->SetSize(1);
	facetUseDesFile->SetValueAt(0, "No file imported");

	facetAccFactor = new GLTextField(0, "");
	paramPanel->SetCompBounds(facetAccFactor, 141, 44, 122, 18);
	paramPanel->Add(facetAccFactor);

	SetTitle("Advanced facet parameters");

	Refresh(0, NULL);
	// Position dialog next to Facet parameters
	int facetX, facetY, facetW, facetH;
	mApp->facetPanel->GetBounds(&facetX, &facetY, &facetW, &facetH);
	SetBounds(facetX - wD - 10, facetY + 20, wD, hD);

	RestoreDeviceObjects();
}

//-----------------------------------------------------------------------------


void FacetMesh::UpdateSize() {
	
	char tmp[64];

	if( enableBtn->GetState() ) {

		llong ram = 0;
		llong cell = 0;
		int nbFacet = geom->GetNbFacet();

		if( recordACBtn->GetState() ) {

			for(int i=0;i<nbFacet;i++) {
				Facet *f = geom->GetFacet(i);
				if(f->sh.opacity==1.0) {
					cell += (llong)f->GetNbCell();
					ram += (llong)f->GetTexRamSize(1+worker->moments.size());
				}
			}
			ram += (((cell-1)*cell)/2 + 8*cell)*((llong)sizeof(ACFLOAT));

		} else {

			for(int i=0;i<nbFacet;i++) {
				Facet *f = geom->GetFacet(i);
				cell += (llong)f->GetNbCell();
				ram += (llong)f->GetTexRamSize(1+worker->moments.size());
			}

		}
		ramText->SetText(FormatMemoryLL(ram));
		sprintf(tmp,"%d",(int)cell);
		cellText->SetText(tmp);

	} else {

		ramText->SetText("0 bytes");
		cellText->SetText("0");

	}

}

//-----------------------------------------------------------------------------

void FacetMesh::UpdateSizeForRatio() {
	if (!geom->IsLoaded()) return;
	double ratio;
	char tmp[64];
	BOOL boundMap = TRUE;// boundaryBtn->GetState();
	BOOL recordDir = recordDirBtn->GetState();

	if( !enableBtn->GetState() ) {
		ramText->SetText(FormatMemory(0));
		cellText->SetText("0");
		return;
	}

	if( sscanf(resolutionText->GetText(),"%lf",&ratio)==0 ) {
		ramText->SetText("");
		cellText->SetText("");
		return;
	}

	llong ram = 0;
	llong cell = 0;
	int nbFacet = geom->GetNbFacet();
	if( recordACBtn->GetState() ) {

		for(int i=0;i<nbFacet;i++) {
			Facet *f = geom->GetFacet(i);
			//if(f->sh.opacity==1.0) {
				if(f->selected) {
					cell += (llong)f->GetNbCellForRatio(ratio);
					ram += (llong)f->GetTexRamSizeForRatio(ratio,boundMap,FALSE,1+worker->moments.size());
				} else {
					cell += (llong)f->GetNbCell();
					ram += (llong)f->GetTexRamSize(1+worker->moments.size());
				}
			//}
		}
		ram += (((cell-1)*cell)/2 + 8*cell)*((llong)sizeof(ACFLOAT));

	} else {

		for(int i=0;i<nbFacet;i++) {
			Facet *f = geom->GetFacet(i);
			if(f->selected) {
				cell += (llong)f->GetNbCellForRatio(ratio);
				ram += (llong)f->GetTexRamSizeForRatio(ratio,boundMap,recordDir,1+worker->moments.size());
			} else {
				cell += (llong)f->GetNbCell();
				ram += (llong)f->GetTexRamSize(1+worker->moments.size());
			}
		}

	}

	ramText->SetText(FormatMemoryLL(ram));
	sprintf(tmp,"%d",(int)cell);
	cellText->SetText(tmp);

}

//-----------------------------------------------------------------------------

void FacetMesh::Refresh(int nbSel, int* selection) {

	char tmp[128];
	double maxU=0.0;
	double maxV=0.0;
	double minU=1.0e100;
	double minV=1.0e100;	

	BOOL somethingSelected = nbSel>0;
	enableBtn->SetEnabled(somethingSelected);
	recordDesBtn->SetEnabled(somethingSelected);
	recordAbsBtn->SetEnabled(somethingSelected);
	recordReflBtn->SetEnabled(somethingSelected);
	recordTransBtn->SetEnabled(somethingSelected);
	recordACBtn->SetEnabled(somethingSelected);
	recordDirBtn->SetEnabled(somethingSelected);
	showTexture->SetEnabled(somethingSelected);
	showVolume->SetEnabled(somethingSelected);
	resolutionText->SetEditable(somethingSelected);
	lengthText->SetEditable(somethingSelected);
	facetReflType->SetEditable(somethingSelected);
	facetAccFactor->SetEditable(somethingSelected);
	facetTeleport->SetEditable(somethingSelected);
	facetStructure->SetEditable(somethingSelected);
	facetSuperDest->SetEditable(somethingSelected);
	facetUseDesFile->SetEditable(somethingSelected);
	facetMovingToggle->SetEnabled(somethingSelected);

	if (!geom->IsLoaded()) return;

	if (!somethingSelected) { //Empty selection, clear
		enableBtn->SetState(0);
		resolutionText->SetText("");
		lengthText->SetText("");
		recordDesBtn->SetState(0);
		recordAbsBtn->SetState(0);
		recordReflBtn->SetState(0);
		recordTransBtn->SetState(0);
		recordACBtn->SetState(0);
		recordDirBtn->SetState(0);
		showTexture->SetState(0);
		showVolume->SetState(0);
		facetUseDesFile->SetSelectedValue("");
		facetReflType->SetSelectedValue("");
		facetAccFactor->Clear();
		facetSuperDest->Clear();
		facetMovingToggle->SetState(0);
		facetStructure->Clear();
		facetTeleport->Clear();
		return;
	}

	Facet* f0 = geom->GetFacet(selection[0]);

	BOOL isEnabledE = TRUE;
	BOOL isBoundE = TRUE;
	BOOL CountDesE = TRUE;
	BOOL CountAbsE = TRUE;
	BOOL CountReflE = TRUE;
	BOOL CountTransE = TRUE;
	BOOL CountACE = TRUE;
	BOOL CountDirE = TRUE;
	BOOL TexVisibleE = TRUE;
	BOOL VolVisibleE = TRUE;
	BOOL ratioE = TRUE;
	BOOL teleportE = TRUE;
	BOOL accFactorE = TRUE;
	BOOL superDestE = TRUE;
	BOOL superIdxE = TRUE;
	BOOL reflectTypeE = TRUE;
	BOOL hasOutgMapE = TRUE;
	BOOL useOutgMapE = TRUE;
	BOOL isMovingE = TRUE;

	for(int i=1;i<nbSel;i++) {

		Facet *f = geom->GetFacet(selection[i]);
		//if( f->selected ) {
			double nU = Norme(&(f->sh.U));
			double nV = Norme(&(f->sh.V));
			maxU = MAX(maxU,nU);
			maxV = MAX(maxV,nV);
			minU = MIN(minU,nU);
			minV = MIN(minV,nV);
			isEnabledE=isEnabledE && (f0->sh.isTextured == f->sh.isTextured);
			isBoundE = isBoundE && (f0->hasMesh == f->hasMesh);
			CountDesE = CountDesE && f0->sh.countDes==f->sh.countDes;
			CountAbsE = CountAbsE && f0->sh.countAbs == f->sh.countAbs;
			CountReflE = CountReflE && f0->sh.countRefl == f->sh.countRefl;
			CountTransE = CountTransE && f0->sh.countTrans == f->sh.countTrans;
			CountACE = CountACE && f0->sh.countACD == f->sh.countACD;
			CountDirE = CountDirE && f0->sh.countDirection == f->sh.countDirection;
			TexVisibleE = TexVisibleE && f0->textureVisible == f->textureVisible;
			VolVisibleE = VolVisibleE && f0->volumeVisible == f->volumeVisible;
			ratioE = ratioE && abs(f0->tRatio - f->tRatio) < 1E-8;
			teleportE = teleportE && (f0->sh.teleportDest == f->sh.teleportDest);
			accFactorE = accFactorE && (abs(f0->sh.accomodationFactor - f->sh.accomodationFactor)<1e-7);
			superDestE = superDestE && (f0->sh.superDest == f->sh.superDest);
			superIdxE = superIdxE && (f0->sh.superIdx == f->sh.superIdx);
			reflectTypeE = reflectTypeE && (f0->sh.reflectType == f->sh.reflectType);
			hasOutgMapE = hasOutgMapE && (f0->hasOutgassingMap == f->hasOutgassingMap);
			useOutgMapE = useOutgMapE && (f0->sh.useOutgassingFile == f->sh.useOutgassingFile);
			isMovingE = isMovingE && (f0->sh.isMoving == f->sh.isMoving);
	}

	if( nbSel==1 ) {
		Facet *f = f0;
		sprintf(tmp,"Advanced parameters (Facet #%d)",selection[0]+1);
		SetTitle(tmp);
		sprintf(tmp,"%g",maxU);
		uLength->SetText(tmp);
		sprintf(tmp,"%g",maxV);
		vLength->SetText(tmp);
	} else {
		sprintf(tmp,"Advanced parameters (%d selected)",nbSel);
		SetTitle(tmp);
		sprintf(tmp,"%g (MAX)",maxU);
		uLength->SetText(tmp);
		sprintf(tmp,"%g (MAX)",maxV);
		vLength->SetText(tmp);
	}

	enableBtn->AllowMixedState(!isEnabledE); enableBtn->SetState(isEnabledE ? f0->sh.isTextured : 2);
	recordDesBtn->AllowMixedState(!CountDesE); recordACBtn->SetState(CountDesE ? f0->sh.countDes : 2);
	recordAbsBtn->AllowMixedState(!CountAbsE); recordAbsBtn->SetState(CountAbsE ? f0->sh.countAbs : 2);
	recordReflBtn->AllowMixedState(!CountReflE); recordReflBtn->SetState(CountReflE ? f0->sh.countRefl : 2);
	recordTransBtn->AllowMixedState(!CountTransE); recordTransBtn->SetState(CountTransE ? f0->sh.countTrans : 2);
	recordACBtn->AllowMixedState(!CountACE); recordACBtn->SetState(CountACE ? f0->sh.countACD : 2);
	recordDirBtn->AllowMixedState(!CountDirE); recordDirBtn->SetState(CountDirE ? f0->sh.countDirection : 2);
	showTexture->AllowMixedState(!TexVisibleE); showTexture->SetState(TexVisibleE ? f0->textureVisible : 2);
	showVolume->AllowMixedState(!VolVisibleE); showVolume->SetState(VolVisibleE ? f0->volumeVisible : 2);
	facetMovingToggle->AllowMixedState(!isMovingE); facetMovingToggle->SetState(isMovingE ? f0->sh.isMoving : 2);

	if( isEnabledE) {
		if (f0->sh.isTextured) { //All facets have textures
			if (ratioE) { //All facets have textures with same resolution
				resolutionText->SetText(f0->tRatio);
				lengthText->SetText(1.0 / f0->tRatio);
			}
			else { //Mixed resolution
				resolutionText->SetText(isEnabledE ? "..." : "");
				lengthText->SetText(isEnabledE ? "..." : "");
			}
		} else { //None of the facets have textures
				resolutionText->SetText("");
				lengthText->SetText("");
			}
	}
	else { //Mixed state
		resolutionText->SetText("");
		lengthText->SetText("");
	}

	if (teleportE) facetTeleport->SetText(f0->sh.teleportDest); else facetTeleport->SetText("...");
	if (accFactorE) facetAccFactor->SetText(f0->sh.accomodationFactor); else facetAccFactor->SetText("...");
	if (reflectTypeE) facetReflType->SetSelectedIndex(f0->sh.reflectType); else facetReflType->SetSelectedValue("...");
	if (hasOutgMapE) { //all selected equally HAVE or equally DON'T HAVE outgassing maps

		if (!f0->hasOutgassingMap) { //All selected DON'T HAVE outgassing maps
			facetUseDesFile->SetSize(1);
			facetUseDesFile->SetSelectedIndex(0); //no map
			facetUseDesFile->SetSelectedValue("No map loaded");
			facetUseDesFile->SetEditable(FALSE);
		}
		else { //All selected HAVE outgassing maps

			facetUseDesFile->SetSize(2);
			facetUseDesFile->SetValueAt(0, "Use user VALUES");
			facetUseDesFile->SetValueAt(1, "Use desorption FILE");
			facetUseDesFile->SetEditable(TRUE);
			if (useOutgMapE) {
				facetUseDesFile->SetSelectedIndex(f0->sh.useOutgassingFile);
			}
			else {
				facetUseDesFile->SetSelectedValue("...");
			}
		}
	} else { //Mixed: some have it, some don't
		facetUseDesFile->SetSelectedIndex(0);
		facetUseDesFile->SetSize(1);
		facetUseDesFile->SetSelectedValue("...");
		facetUseDesFile->SetEditable(FALSE);
	}

	if (superDestE) {
		if (f0->sh.superDest == 0) {
			facetSuperDest->SetText("no");
		}
		else {
			sprintf(tmp, "%d", f0->sh.superDest);
			facetSuperDest->SetText(tmp);
		}
	}
	else {
		facetSuperDest->SetText("...");
	}
	if (superIdxE) {
		sprintf(tmp, "%d", f0->sh.superIdx + 1);
		facetStructure->SetText(tmp);
	}
	else {
		facetStructure->SetText("...");
	}

	if (isMovingE) {
		facetMovingToggle->SetState(f0->sh.isMoving);
		facetMovingToggle->AllowMixedState(FALSE);
	}
	else {
		facetMovingToggle->SetState(2);
		facetMovingToggle->AllowMixedState(TRUE);
	}

	UpdateSize();
}


//-----------------------------------------------------------------------------

BOOL FacetMesh::Apply() {
	if (!mApp->AskToReset(worker)) return FALSE;
	BOOL boundMap = TRUE; // boundaryBtn->GetState();
	int nbSelected;
	int* selection;
	double ratio;
	geom->GetSelection(&selection, &nbSelected);
	int nbPerformed = 0.0;

	if (enableBtn->GetState() == 1) {

		// Check counting mode
		if (!recordDesBtn->GetState() && !recordAbsBtn->GetState() &&
			!recordReflBtn->GetState() && !recordTransBtn->GetState() &&
			!recordACBtn->GetState() && !recordDirBtn->GetState()) {
			GLMessageBox::Display("Please select counting mode", "Error", GLDLG_OK, GLDLG_ICONINFO);
			Refresh(nbSelected, selection);
			return FALSE;
		}

		// Resolution
		if (!resolutionText->GetNumber(&ratio) || ratio<=0.0) {
			GLMessageBox::Display("Invalid texture resolution", "Error", GLDLG_OK, GLDLG_ICONERROR);
			Refresh(nbSelected, selection);
			return FALSE;
		}

	}

	// Superstructure

	BOOL structChanged = FALSE; //if a facet gets into a new structure, we have to re-render the geometry

	int superStruct;
	BOOL doSuperStruct = FALSE;
	if (sscanf(facetStructure->GetText(), "%d", &superStruct)>0 && superStruct>0 && superStruct <= geom->GetNbStructure()) doSuperStruct = TRUE;
	else {
		if (strcmp(facetStructure->GetText(), "...") == 0) doSuperStruct = FALSE;
		else{
			GLMessageBox::Display("Invalid superstructre number", "Error", GLDLG_OK, GLDLG_ICONERROR);
			Refresh(nbSelected, selection);
			return FALSE;
		}
	}


	// Super structure destination (link)
	int superDest;
	BOOL doSuper = FALSE;
	if (strcmp(facetSuperDest->GetText(), "none") == 0 || strcmp(facetSuperDest->GetText(), "no") == 0 || strcmp(facetSuperDest->GetText(), "0") == 0) {
		doSuper = TRUE;
		superDest = 0;
	}
	else if (sscanf(facetSuperDest->GetText(), "%d", &superDest)>0) {
		if (superDest == superStruct) {
			GLMessageBox::Display("Link and superstructure can't be the same", "Error", GLDLG_OK, GLDLG_ICONERROR);
			Refresh(nbSelected, selection);
			return FALSE;
		}
		else if (superDest>0 && superDest <= geom->GetNbStructure()) doSuper = TRUE;
	}
	else if (strcmp(facetSuperDest->GetText(), "...") == 0) doSuper = FALSE;

	else {

		GLMessageBox::Display("Invalid superstructure destination", "Error", GLDLG_OK, GLDLG_ICONERROR);
		Refresh(nbSelected, selection);
		return FALSE;
	}

	// teleport
	int teleport;
	BOOL doTeleport = FALSE;

	if (facetTeleport->GetNumberInt(&teleport)) {
		if (teleport<0 || teleport>geom->GetNbFacet()) {
			GLMessageBox::Display("Invalid teleport destination\n(If no teleport: set number to 0)", "Error", GLDLG_OK, GLDLG_ICONERROR);
			Refresh(nbSelected, selection);
			return FALSE;
		}
		else if (teleport>0 && geom->GetFacet(teleport - 1)->selected) {
			char tmp[256];
			sprintf(tmp, "The teleport destination of facet #%d can't be itself!", teleport);
			GLMessageBox::Display(tmp, "Error", GLDLG_OK, GLDLG_ICONERROR);
			Refresh(nbSelected, selection);
			return FALSE;
		}
		doTeleport = TRUE;
	}
	else {
		if (strcmp(facetTeleport->GetText(), "...") == 0) doTeleport = FALSE;
		else {
			GLMessageBox::Display("Invalid teleport destination\n(If no teleport: set number to 0)", "Error", GLDLG_OK, GLDLG_ICONERROR);
			Refresh(nbSelected, selection);
			return FALSE;
		}
	}

	// temp.accomodation factor
	double accfactor;
	BOOL doAccfactor = FALSE;
	if (facetAccFactor->GetNumber(&accfactor)) {
		if (accfactor<0.0 || accfactor>1.0) {
			GLMessageBox::Display("Facet accomodation factor must be between 0 and 1", "Error", GLDLG_OK, GLDLG_ICONERROR);
			Refresh(nbSelected, selection);
			return FALSE;
		}
		doAccfactor = TRUE;
	}
	else {
		if (strcmp(facetAccFactor->GetText(), "...") == 0) doAccfactor = FALSE;
		else {
			GLMessageBox::Display("Invalid accomodation factor number", "Error", GLDLG_OK, GLDLG_ICONERROR);
			Refresh(nbSelected, selection);
			return FALSE;
		}
	}

	// Use desorption map
	int useMapA = 0;
	BOOL doUseMapA = FALSE;
	if (strcmp(facetUseDesFile->GetSelectedValue(), "...") == 0) doUseMapA = FALSE;
	else {
		useMapA = (facetUseDesFile->GetSelectedIndex() == 1);
		BOOL missingMap = FALSE;
		int missingMapId;
		if (useMapA) {
			for (int i = 0; i<geom->GetNbFacet(); i++) {
				if (geom->GetFacet(i)->selected && !geom->GetFacet(i)->hasOutgassingMap) {
					missingMap = TRUE;
					missingMapId = i;
				}
			}
		}
		if (missingMap) {
			char tmp[256];
			sprintf(tmp, "%d is selected but doesn't have any outgassing map loaded.", missingMapId + 1);
			GLMessageBox::Display(tmp, "Can't use map on all facets", GLDLG_OK, GLDLG_ICONERROR);
			Refresh(nbSelected, selection);
			return FALSE;
		}
		doUseMapA = TRUE;
	}

	// Reflection type
	int reflType = facetReflType->GetSelectedIndex();







	//Check complete, let's apply
	progressDlg = new GLProgress("Applying mesh settings","Please wait");
	progressDlg->SetVisible(TRUE);
	progressDlg->SetProgress(0.0);
	int count=0;
	for(int i=0;i<nbSelected;i++) {
		Facet *f = geom->GetFacet(selection[i]);
		if (recordDesBtn->GetState()<2) f->sh.countDes = recordDesBtn->GetState()*enableBtn->GetState();
		if (recordAbsBtn->GetState()<2) f->sh.countAbs = recordAbsBtn->GetState()*enableBtn->GetState();
		if (recordReflBtn->GetState()<2) f->sh.countRefl = recordReflBtn->GetState()*enableBtn->GetState();
		if (recordTransBtn->GetState()<2) f->sh.countTrans = recordTransBtn->GetState()*enableBtn->GetState();
		if (recordACBtn->GetState()<2) f->sh.countACD = recordACBtn->GetState()*enableBtn->GetState();
		if (recordDirBtn->GetState()<2) f->sh.countDirection = recordDirBtn->GetState()*enableBtn->GetState();
		BOOL hasAnyTexture = f->sh.countDes || f->sh.countAbs || f->sh.countRefl || f->sh.countTrans || f->sh.countACD || f->sh.countDirection;
		if (doTeleport) f->sh.teleportDest = teleport;
		if (doAccfactor) f->sh.accomodationFactor = accfactor;
		if (reflType >= 0) f->sh.reflectType = reflType;
		if (doSuperStruct) {
			if (f->sh.superIdx != (superStruct - 1)) {
				f->sh.superIdx = superStruct - 1;
				structChanged = TRUE;
			}
		}
		if (doSuper) {
			f->sh.superDest = superDest;
			if (superDest) f->sh.opacity = 1.0; // Force opacity for link facet
		}
		//Moving or not
		if (facetMovingToggle->GetState() < 2) f->sh.isMoving = facetMovingToggle->GetState();

		if (doUseMapA) {
			f->sh.useOutgassingFile = useMapA;
		}

		//set textures
		try {
			if (structChanged) geom->RebuildLists();
			geom->SetFacetTexture(selection[i], hasAnyTexture ? ratio : 0.0, hasAnyTexture ? boundMap : 0.0);
		} catch (Error &e) {
			GLMessageBox::Display((char *)e.GetMsg(),"Error",GLDLG_OK,GLDLG_ICONWARNING);
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			return FALSE;
		} catch (...) {
			GLMessageBox::Display("Unexpected error while setting textures","Error",GLDLG_OK,GLDLG_ICONWARNING);
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			return FALSE;
		}
		if (showTexture->GetState()<2) f->textureVisible = showTexture->GetState();
		if (showVolume->GetState()<2) f->volumeVisible = showVolume->GetState();
		nbPerformed++;
		progressDlg->SetProgress((double)nbPerformed/(double)nbSelected);
	} //main cycle end


	if (progressDlg) progressDlg->SetVisible(FALSE);
	SAFE_DELETE(progressDlg);
	return TRUE;

}

//-----------------------------------------------------------------------------
void FacetMesh::QuickApply() {
	//Apply view settings without stopping the simulation


	double nbSelected = (double)geom->GetNbSelected();
	double nbPerformed = 0.0;

	for(int i=0;i<geom->GetNbFacet();i++) {

		Facet *f = geom->GetFacet(i);
		if( f->selected ) {

			if (showTexture->GetState()<2) f->textureVisible = showTexture->GetState();
			if (showVolume->GetState()<2) f->volumeVisible = showVolume->GetState();

			nbPerformed+=1.0;
			progressDlg->SetProgress(nbPerformed/nbSelected);
		}

	}
	geom->RebuildLists();
}

//-----------------------------------------------------------------------------

void FacetMesh::UpdateToggle(GLComponent *src) {

	/*if (src==boundaryBtn) {
		recordACBtn->SetState(FALSE);
	} else if(src==enableBtn) {
		//boundaryBtn->SetState(enableBtn->GetState());
	} else */

	if(src==recordDesBtn ) {
		enableBtn->SetState(TRUE);
		//boundaryBtn->SetState(TRUE);
		recordACBtn->SetState(FALSE);
	} else if(src==recordAbsBtn ) {
		enableBtn->SetState(TRUE);
		//boundaryBtn->SetState(TRUE);
		recordACBtn->SetState(FALSE);
	} else if(src==recordReflBtn ) {
		enableBtn->SetState(TRUE);
		//boundaryBtn->SetState(TRUE);
		recordACBtn->SetState(FALSE);
	} else if(src==recordTransBtn ) {
		enableBtn->SetState(TRUE);
		//boundaryBtn->SetState(TRUE);
		recordACBtn->SetState(FALSE);
	} else if(src==recordDirBtn ) {
		enableBtn->SetState(TRUE);
		//boundaryBtn->SetState(TRUE);
		recordACBtn->SetState(FALSE);
	} else if(src==recordACBtn) {
		if( recordACBtn->GetState() ) {
			enableBtn->SetState(TRUE);
			//boundaryBtn->SetState(TRUE);
			recordDesBtn->SetState(FALSE);
			recordAbsBtn->SetState(FALSE);
			recordReflBtn->SetState(FALSE);
			recordTransBtn->SetState(FALSE);
			recordDirBtn->SetState(FALSE);
		}
	}

	UpdateSizeForRatio();
}

//-----------------------------------------------------------------------------

void FacetMesh::ProcessMessage(GLComponent *src,int message) {

	switch(message) {

		// -------------------------------------------------------------
	case MSG_BUTTON:
		/*if(src==cancelButton) {

			GLWindow::ProcessMessage(NULL,MSG_CLOSE);

		} else if (src==applyButton) {


			//if (worker->running) worker->Stop_Public();
			if( Apply() )
				GLWindow::ProcessMessage(NULL,MSG_CLOSE);

		} else */if (src==quickApply) {

			progressDlg = new GLProgress("Applying view settings","Please wait");
			progressDlg->SetVisible(TRUE);
			progressDlg->SetProgress(0.5);

			QuickApply();

			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);

		} 
		/*
		else if (src==updateButton) {

		UpdateSizeForRatio();

		}
		*/
		break;

		// -------------------------------------------------------------
	case MSG_TEXT_UPD:
		enableBtn->SetState(TRUE);
		
		if (src == resolutionText) {
			UpdateSizeForRatio();
			mApp->facetApplyBtn->SetEnabled(TRUE);
			double res;
			if (resolutionText->GetNumber(&res) && res != 0.0)
				lengthText->SetText(1.0 / res);
			else
				lengthText->SetText("");
		}
		else if (src == lengthText) {
			double length;
			if (lengthText->GetNumber(&length) && length != 0.0) {
				resolutionText->SetText(1.0 / length);
				UpdateSizeForRatio();
				mApp->facetApplyBtn->SetEnabled(TRUE);
			} else
				resolutionText->SetText("");
		}
		else if (src == facetTeleport) {
			mApp->facetApplyBtn->SetEnabled(TRUE);
		}
		else if (src == facetAccFactor) {
			mApp->facetApplyBtn->SetEnabled(TRUE);
		}
		else if (src == facetSuperDest || src == facetStructure) {
			mApp->facetApplyBtn->SetEnabled(TRUE);
		}
		
		break;

		// -------------------------------------------------------------
	case MSG_TOGGLE:
		UpdateToggle(src);
		mApp->facetApplyBtn->SetEnabled(TRUE);
		break;
	case MSG_TEXT:
		if (src == facetTeleport) {
			mApp->ApplyFacetParams();
		}
		else if (src == facetAccFactor) {
			mApp->ApplyFacetParams();
		}
		else if (src == facetSuperDest || src == facetStructure) {
			mApp->ApplyFacetParams();
		}
		break;
	case MSG_COMBO:
		if (src == facetReflType) {
			mApp->facetApplyBtn->SetEnabled(TRUE);
		}
		else if (src == facetUseDesFile) {
			mApp->facetApplyBtn->SetEnabled(TRUE);
		}
		break;
	}

	GLWindow::ProcessMessage(src,message);
}

