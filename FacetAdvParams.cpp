/*
File:        FacetAdvParams.cpp
Description: Advanced facet settings window (used to be fact mesh settings)
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

#include "FacetAdvParams.h"
#include "Facet.h"
#include "GLApp/GLLabel.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLMessageBox.h"
#include "GLApp/MathTools.h"
#include "GLApp/GLToggle.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLTextField.h"
#include "GLApp/GLTitledPanel.h"
//#include "GLApp/GLProgress.h"
#include "GLApp/GLCombo.h"
//#include "Worker.h"
#include "Geometry.h"

#ifdef MOLFLOW
#include "MolFlow.h"
#endif

#ifdef SYNRAD
#include "SynRad.h"
#endif
#ifdef MOLFLOW
extern MolFlow *mApp;
#endif

#ifdef SYNRAD
extern SynRad*mApp;
#endif

//-----------------------------------------------------------------------------

FacetAdvParams::FacetAdvParams(Worker *w) :GLWindow() {

	worker = w;
	geom = w->GetGeometry();

	SetIconfiable(true);

	int wD = 320;
	int hD = 574;
	aPanel = new GLTitledPanel("Texture properties");
	aPanel->SetBounds(5, 3, 309, 123);
	Add(aPanel);
	mPanel = new GLTitledPanel("Texture cell / memory");
	mPanel->SetBounds(5, 129, 309, 44);
	Add(mPanel);
	vPanel = new GLTitledPanel("View settings");
	vPanel->SetBounds(5, 430, 309, 44);
	Add(vPanel);
	desPanel = new GLTitledPanel("Dynamic desorption");
	desPanel->SetBounds(5, 478, 309, 69);
	Add(desPanel);
	paramPanel = new GLTitledPanel("Additional parameters");
	paramPanel->SetBounds(5, 177, 309, 174);
	Add(paramPanel);
	angleMapPanel = new GLTitledPanel("Incident angle distribution");
	angleMapPanel->SetBounds(5, 357, 309, 67);
	Add(angleMapPanel);
	lengthText = new GLTextField(0, "");
	aPanel->SetCompBounds(lengthText, 180, 36, 72, 18);
	aPanel->Add(lengthText);

	perCm = new GLLabel("cells/cm");
	aPanel->SetCompBounds(perCm, 125, 39, 40, 12);
	aPanel->Add(perCm);

	resolutionText = new GLTextField(0, "");
	aPanel->SetCompBounds(resolutionText, 60, 36, 62, 18);
	aPanel->Add(resolutionText);

	l5 = new GLLabel("Resolution:");
	aPanel->SetCompBounds(l5, 5, 39, 52, 12);
	aPanel->Add(l5);

	enableBtn = new GLToggle(0, "Enable texture");
	aPanel->SetCompBounds(enableBtn, 5, 18, 83, 16);
	aPanel->Add(enableBtn);

	recordDesBtn = new GLToggle(0, "Count desorption");
	aPanel->SetCompBounds(recordDesBtn, 5, 60, 94, 16);
	aPanel->Add(recordDesBtn);

	perCell = new GLLabel("cm/cell");
	aPanel->SetCompBounds(perCell, 260, 39, 35, 12);
	aPanel->Add(perCell);

	recordDirBtn = new GLToggle(0, "Record direction vectors");
	aPanel->SetCompBounds(recordDirBtn, 165, 102, 125, 16);
	aPanel->Add(recordDirBtn);

	recordTransBtn = new GLToggle(0, "Count transparent pass");
	aPanel->SetCompBounds(recordTransBtn, 165, 81, 120, 16);
	aPanel->Add(recordTransBtn);

	recordReflBtn = new GLToggle(0, "Count reflection");
	aPanel->SetCompBounds(recordReflBtn, 165, 60, 89, 16);
	aPanel->Add(recordReflBtn);

	recordACBtn = new GLToggle(0, "Angular coefficient");
	aPanel->SetCompBounds(recordACBtn, 5, 102, 101, 16);
	aPanel->Add(recordACBtn);

	recordAbsBtn = new GLToggle(0, "Count absorption");
	aPanel->SetCompBounds(recordAbsBtn, 5, 81, 94, 16);
	aPanel->Add(recordAbsBtn);

	showTexture = new GLToggle(0, "Draw Texture");
	vPanel->SetCompBounds(showTexture, 10, 18, 80, 16);
	vPanel->Add(showTexture);

	showVolume = new GLToggle(0, "Draw Volume");
	vPanel->SetCompBounds(showVolume, 110, 18, 81, 16);
	vPanel->Add(showVolume);

	cellText = new GLTextField(0, "");
	mPanel->SetCompBounds(cellText, 195, 19, 107, 18);
	mPanel->Add(cellText);

	l8 = new GLLabel("Cells:");
	mPanel->SetCompBounds(l8, 166, 22, 29, 12);
	mPanel->Add(l8);

	ramText = new GLTextField(0, "");
	mPanel->SetCompBounds(ramText, 58, 19, 100, 18);
	mPanel->Add(ramText);

	l7 = new GLLabel("Memory:");
	mPanel->SetCompBounds(l7, 10, 22, 43, 12);
	mPanel->Add(l7);

	quickApply = new GLButton(0, "<- Change draw");
	vPanel->SetCompBounds(quickApply, 200, 15, 99, 20);
	vPanel->Add(quickApply);

	fileYieldText = new GLTextField(0, "");
	desPanel->SetCompBounds(fileYieldText, 205, 18, 58, 18);
	desPanel->Add(fileYieldText);

	label3 = new GLLabel("mol/ph");
	desPanel->SetCompBounds(label3, 265, 21, 33, 12);
	desPanel->Add(label3);

	label1 = new GLLabel("Avg");
	desPanel->SetCompBounds(label1, 155, 21, 48, 12);
	desPanel->Add(label1);

	label2 = new GLLabel("Use file:");
	desPanel->SetCompBounds(label2, 5, 21, 39, 12);
	desPanel->Add(label2);

	facetMovingToggle = new GLToggle(0, "Moving part");
	paramPanel->SetCompBounds(facetMovingToggle, 10, 111, 74, 16);
	paramPanel->Add(facetMovingToggle);

	facetSuperDest = new GLTextField(0, "");
	paramPanel->SetCompBounds(facetSuperDest, 199, 90, 101, 18);
	paramPanel->Add(facetSuperDest);

	label8 = new GLLabel("Link to:");
	paramPanel->SetCompBounds(label8, 163, 96, 35, 12);
	paramPanel->Add(label8);

	facetStructure = new GLTextField(0, "");
	paramPanel->SetCompBounds(facetStructure, 60, 90, 91, 18);
	paramPanel->Add(facetStructure);

	label7 = new GLLabel("Structure:");
	paramPanel->SetCompBounds(label7, 10, 91, 46, 12);
	paramPanel->Add(label7);

	facetTeleport = new GLTextField(0, "");
	paramPanel->SetCompBounds(facetTeleport, 155, 66, 145, 18);
	paramPanel->Add(facetTeleport);

	label4 = new GLLabel("Teleport to facet:");
	paramPanel->SetCompBounds(label4, 10, 69, 74, 12);
	paramPanel->Add(label4);

	label5 = new GLLabel("Accomodation coefficient:");
	paramPanel->SetCompBounds(label5, 10, 47, 126, 13);
	paramPanel->Add(label5);

	label6 = new GLLabel("Reflection:");
	paramPanel->SetCompBounds(label6, 10, 22, 50, 12);
	paramPanel->Add(label6);

	facetReflType = new GLCombo(0);
	paramPanel->SetCompBounds(facetReflType, 155, 18, 147, 20);
	paramPanel->Add(facetReflType);
	facetReflType->SetSize(3);
	facetReflType->SetValueAt(0, "Diffuse");
	facetReflType->SetValueAt(1, "Mirror");
	facetReflType->SetValueAt(2, "Uniform");

	facetUseDesFile = new GLCombo(0);
	desPanel->SetCompBounds(facetUseDesFile, 50, 18, 95, 20);
	desPanel->Add(facetUseDesFile);
	facetUseDesFile->SetSize(1);
	facetUseDesFile->SetValueAt(0, "No file imported");

	facetAccFactor = new GLTextField(0, "");
	paramPanel->SetCompBounds(facetAccFactor, 155, 42, 145, 18);
	paramPanel->Add(facetAccFactor);

	fileDoseText = new GLTextField(0, "");
	desPanel->SetCompBounds(fileDoseText, 205, 42, 58, 18);
	desPanel->Add(fileDoseText);

	fileFluxText = new GLTextField(0, "");
	desPanel->SetCompBounds(fileFluxText, 50, 42, 55, 18);
	desPanel->Add(fileFluxText);

	label11 = new GLLabel("ph/cm2");
	desPanel->SetCompBounds(label11, 265, 45, 36, 12);
	desPanel->Add(label11);

	label9 = new GLLabel("ph/s/cm2");
	desPanel->SetCompBounds(label9, 105, 45, 44, 12);
	desPanel->Add(label9);

	label12 = new GLLabel("Avg");
	desPanel->SetCompBounds(label12, 155, 45, 49, 12);
	desPanel->Add(label12);

	label10 = new GLLabel("Avg");
	desPanel->SetCompBounds(label10, 5, 45, 44, 12);
	desPanel->Add(label10);

	enableSojournTime = new GLToggle(0, "Wall Sojourn time");
	paramPanel->SetCompBounds(enableSojournTime, 10, 132, 95, 16);
	paramPanel->Add(enableSojournTime);

	sojournLabel3 = new GLLabel("J/molecule");
	paramPanel->SetCompBounds(sojournLabel3, 250, 153, 56, 13);
	paramPanel->Add(sojournLabel3);

	sojournE = new GLTextField(0, "");
	paramPanel->SetCompBounds(sojournE, 200, 150, 50, 18);
	paramPanel->Add(sojournE);

	sojournLabel2 = new GLLabel("Hz   Binding E:");
	paramPanel->SetCompBounds(sojournLabel2, 125, 153, 76, 13);
	paramPanel->Add(sojournLabel2);

	sojournLabel1 = new GLLabel("Attempt freq:");
	paramPanel->SetCompBounds(sojournLabel1, 10, 153, 67, 13);
	paramPanel->Add(sojournLabel1);

	sojournFreq = new GLTextField(0, "");
	paramPanel->SetCompBounds(sojournFreq, 75, 150, 50, 18);
	paramPanel->Add(sojournFreq);

	SojournInfoButton = new GLButton(0, "Info");
	paramPanel->SetCompBounds(SojournInfoButton, 229, 126, 69, 20);
	paramPanel->Add(SojournInfoButton);

	label14 = new GLLabel("x");
	angleMapPanel->SetCompBounds(label14, 237, 17, 10, 12);
	angleMapPanel->Add(label14);

	angleMapDistrSizeLabel = new GLLabel("Distribution size:");
	angleMapPanel->SetCompBounds(angleMapDistrSizeLabel, 105, 17, 74, 12);
	angleMapPanel->Add(angleMapDistrSizeLabel);

	angleMapHeightBox = new GLTextField(0, "");
	angleMapPanel->SetCompBounds(angleMapHeightBox, 249, 14, 50, 18);
	angleMapPanel->Add(angleMapHeightBox);

	angleMapWidthBox = new GLTextField(0, "");
	angleMapPanel->SetCompBounds(angleMapWidthBox, 184, 14, 50, 18);
	angleMapPanel->Add(angleMapWidthBox);

	angleMapImportButton = new GLButton(0, "Import CSV");
	angleMapPanel->SetCompBounds(angleMapImportButton, 110, 36, 94, 20);
	angleMapPanel->Add(angleMapImportButton);

	angleMapExportButton = new GLButton(0, "Export to CSV");
	angleMapPanel->SetCompBounds(angleMapExportButton, 10, 36, 94, 20);
	angleMapPanel->Add(angleMapExportButton);

	angleMapRecordCheckbox = new GLToggle(0, "Record");
	angleMapPanel->SetCompBounds(angleMapRecordCheckbox, 10, 16, 54, 16);
	angleMapPanel->Add(angleMapRecordCheckbox);

	angleMapReleaseButton = new GLButton(0, "Release recorded");
	angleMapPanel->SetCompBounds(angleMapReleaseButton, 208, 36, 94, 20);
	angleMapPanel->Add(angleMapReleaseButton);

	remeshButton = new GLButton(0, "Force remesh");
	aPanel->SetCompBounds(remeshButton, 216, 13, 79, 20);
	aPanel->Add(remeshButton);

	SetTitle("Advanced facet parameters");
	// Center dialog
	int wS, hS;
	GLToolkit::GetScreenSize(&wS, &hS);
	int xD = (wS - wD) / 2;
	int yD = (hS - hD) / 2;
	SetBounds(xD, yD, wD, hD);

	cellText->SetEditable(false);
	ramText->SetEditable(false);
	fileDoseText->SetEditable(false);
	fileYieldText->SetEditable(false);
	fileFluxText->SetEditable(false);


	Refresh(std::vector<size_t>());
	Reposition(wD, hD);

	RestoreDeviceObjects();
}

//-----------------------------------------------------------------------------


void FacetAdvParams::UpdateSize() {

	char tmp[64];

	if (enableBtn->GetState()) {

		llong ram = 0;
		llong cell = 0;
		size_t nbFacet = geom->GetNbFacet();

		if (recordACBtn->GetState()) {

			for (size_t i = 0; i < nbFacet; i++) {
				Facet *f = geom->GetFacet(i);
				if (f->sh.opacity == 1.0) {
					cell += (llong)f->GetNbCell();
					ram += (llong)f->GetTexRamSize(1 + worker->moments.size());
				}
			}
			ram += (((cell - 1)*cell) / 2 + 8 * cell)*((llong)sizeof(ACFLOAT));

		}
		else {

			for (size_t i = 0; i < nbFacet; i++) {
				Facet *f = geom->GetFacet(i);
				cell += (llong)f->GetNbCell();
				ram += (llong)f->GetTexRamSize(1 + worker->moments.size());
			}

		}
		ramText->SetText(FormatMemoryLL(ram));
		sprintf(tmp, "%zd", cell);
		cellText->SetText(tmp);

	}
	else {

		ramText->SetText("0 bytes");
		cellText->SetText("0");

	}

}

//-----------------------------------------------------------------------------

void FacetAdvParams::UpdateSizeForRatio() {
	if (!geom->IsLoaded()) return;
	double ratio;
	char tmp[64];
	bool boundMap = true;// boundaryBtn->GetState();
	bool recordDir = recordDirBtn->GetState();

	if (!enableBtn->GetState()) {
		ramText->SetText(FormatMemory(0));
		cellText->SetText("0");
		return;
	}

	if (sscanf(resolutionText->GetText(), "%lf", &ratio) == 0) {
		ramText->SetText("");
		cellText->SetText("");
		return;
	}

	llong ram = 0;
	llong cell = 0;
	size_t nbFacet = geom->GetNbFacet();
	if (recordACBtn->GetState()) {

		for (size_t i = 0; i < nbFacet; i++) {
			Facet *f = geom->GetFacet(i);
			//if(f->sh.opacity==1.0) {
			if (f->selected) {
				cell += (llong)f->GetNbCellForRatio(ratio);
				ram += (llong)f->GetTexRamSizeForRatio(ratio, boundMap, false, 1 + worker->moments.size());
			}
			else {
				cell += (llong)f->GetNbCell();
				ram += (llong)f->GetTexRamSize(1 + worker->moments.size());
			}
			//}
		}
		ram += (((cell - 1)*cell) / 2 + 8 * cell)*((llong)sizeof(ACFLOAT));

	}
	else {

		for (size_t i = 0; i < nbFacet; i++) {
			Facet *f = geom->GetFacet(i);
			if (f->selected) {
				cell += (llong)f->GetNbCellForRatio(ratio);
				ram += (llong)f->GetTexRamSizeForRatio(ratio, boundMap, recordDir, 1 + worker->moments.size());
			}
			else {
				cell += (llong)f->GetNbCell();
				ram += (llong)f->GetTexRamSize(1 + worker->moments.size());
			}
		}

	}

	ramText->SetText(FormatMemoryLL(ram));
	sprintf(tmp, "%zd", cell);
	cellText->SetText(tmp);

}

//-----------------------------------------------------------------------------

void FacetAdvParams::Refresh(std::vector<size_t> selection) {

	sumArea = sumOutgassing = 0;

	bool somethingSelected = selection.size() > 0;
	enableBtn->SetEnabled(somethingSelected);
	recordDesBtn->SetEnabled(somethingSelected);
	recordAbsBtn->SetEnabled(somethingSelected);
	recordReflBtn->SetEnabled(somethingSelected);
	recordTransBtn->SetEnabled(somethingSelected);
	recordACBtn->SetEnabled(somethingSelected);
	recordDirBtn->SetEnabled(somethingSelected);
	angleMapRecordCheckbox->SetEnabled(somethingSelected);
	angleMapWidthBox->SetEditable(somethingSelected);
	angleMapHeightBox->SetEditable(somethingSelected);
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
	enableSojournTime->SetEnabled(somethingSelected);
	sojournFreq->SetEditable(somethingSelected);
	sojournE->SetEditable(somethingSelected);

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
		angleMapRecordCheckbox->SetState(0);
		angleMapWidthBox->SetText("");
		angleMapHeightBox->SetText("");
		showTexture->SetState(0);
		showVolume->SetState(0);
		facetUseDesFile->SetSelectedValue("");
		facetReflType->SetSelectedValue("");
		facetAccFactor->Clear();
		facetSuperDest->Clear();
		facetMovingToggle->SetState(0);
		facetStructure->Clear();
		facetTeleport->Clear();
		enableSojournTime->SetState(0);
		enableSojournTime->SetText("Wall sojourn time");
		sojournFreq->SetText("");
		sojournE->SetText("");
		return;
	}

	Facet* f0 = geom->GetFacet(selection[0]);

	bool isEnabledE = true;
	bool isBoundE = true;
	bool CountDesE = true;
	bool CountAbsE = true;
	bool CountReflE = true;
	bool CountTransE = true;
	bool CountACE = true;
	bool CountDirE = true;
	bool RecordAngleMapE = true;
	bool AngleMapWidthE = true;
	bool AngleMapHeightE = true;
	bool TexVisibleE = true;
	bool VolVisibleE = true;
	bool ratioE = true;
	bool teleportE = true;
	bool accFactorE = true;
	bool superDestE = true;
	bool superIdxE = true;
	bool reflectTypeE = true;
	bool hasOutgMapE = true;
	bool useOutgMapE = true;
	bool yieldEqual = true;
	bool fluxAEqual = true;
	bool doseAEqual = true;
	bool isMovingE = true;
	bool dynOutgEqual = true;
	bool dynOutgAEqual = true;
	bool hasSojournE = true;
	bool sojournFreqE = true;
	bool sojournEE = true;

	double f0Area = f0->GetArea();
	sumArea = f0Area;
	sumOutgassing = f0->sh.totalOutgassing;
	for (size_t i = 1; i < selection.size(); i++) {
		Facet *f = geom->GetFacet(selection[i]);
		double fArea = f->GetArea();
		sumArea += fArea;
		sumOutgassing += f->sh.totalOutgassing;
		isEnabledE = isEnabledE && (f0->sh.isTextured == f->sh.isTextured);
		isBoundE = isBoundE && (f0->hasMesh == f->hasMesh);
		CountDesE = CountDesE && f0->sh.countDes == f->sh.countDes;
		CountAbsE = CountAbsE && f0->sh.countAbs == f->sh.countAbs;
		CountReflE = CountReflE && f0->sh.countRefl == f->sh.countRefl;
		CountTransE = CountTransE && f0->sh.countTrans == f->sh.countTrans;
		CountACE = CountACE && f0->sh.countACD == f->sh.countACD;
		CountDirE = CountDirE && f0->sh.countDirection == f->sh.countDirection;
		RecordAngleMapE = RecordAngleMapE && f0->sh.recordAngleMap == f->sh.recordAngleMap;
		AngleMapWidthE = AngleMapWidthE && (f0->sh.angleMapPhiWidth == f->sh.angleMapPhiWidth);
		AngleMapHeightE = AngleMapHeightE && (f0->sh.angleMapThetaHeight == f->sh.angleMapThetaHeight);
		TexVisibleE = TexVisibleE && f0->textureVisible == f->textureVisible;
		VolVisibleE = VolVisibleE && f0->volumeVisible == f->volumeVisible;
		ratioE = ratioE && abs(f0->tRatio - f->tRatio) < 1E-8;
		teleportE = teleportE && (f0->sh.teleportDest == f->sh.teleportDest);
		accFactorE = accFactorE && IsEqual(f0->sh.accomodationFactor, f->sh.accomodationFactor);
		superDestE = superDestE && (f0->sh.superDest == f->sh.superDest);
		superIdxE = superIdxE && (f0->sh.superIdx == f->sh.superIdx);
		reflectTypeE = reflectTypeE && (f0->sh.reflectType == f->sh.reflectType);
		hasOutgMapE = hasOutgMapE && (f0->hasOutgassingFile == f->hasOutgassingFile);
		useOutgMapE = useOutgMapE && (f0->sh.useOutgassingFile == f->sh.useOutgassingFile);
		dynOutgEqual = dynOutgEqual && IsEqual(f0->sh.totalOutgassing, f->sh.totalOutgassing);
		dynOutgAEqual = dynOutgAEqual && IsEqual(f0->sh.totalOutgassing / f0Area, f->sh.totalOutgassing / fArea);
		yieldEqual = yieldEqual && IsEqual(f0->sh.totalOutgassing / f0->sh.temperature / f0->totalFlux, f->sh.totalOutgassing / f->sh.temperature / f->totalFlux);
		fluxAEqual = fluxAEqual && IsEqual(f0->totalFlux / f0Area, f->totalFlux / fArea);
		doseAEqual = doseAEqual && IsEqual(f0->totalDose / f0Area, f->totalDose / fArea);
		isMovingE = isMovingE && (f0->sh.isMoving == f->sh.isMoving);
		hasSojournE = hasSojournE && (f0->sh.enableSojournTime == f->sh.enableSojournTime);
		sojournFreqE = sojournFreqE && IsEqual(f0->sh.sojournFreq, f->sh.sojournFreq);
		sojournEE = sojournEE && IsEqual(f0->sh.sojournE, f->sh.sojournE);
	}

	/*
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
	}*/

	enableBtn->AllowMixedState(!isEnabledE); enableBtn->SetState(isEnabledE ? f0->sh.isTextured : 2);
	recordDesBtn->AllowMixedState(!CountDesE); recordDesBtn->SetState(CountDesE ? f0->sh.countDes : 2);
	recordAbsBtn->AllowMixedState(!CountAbsE); recordAbsBtn->SetState(CountAbsE ? f0->sh.countAbs : 2);
	recordReflBtn->AllowMixedState(!CountReflE); recordReflBtn->SetState(CountReflE ? f0->sh.countRefl : 2);
	recordTransBtn->AllowMixedState(!CountTransE); recordTransBtn->SetState(CountTransE ? f0->sh.countTrans : 2);
	recordACBtn->AllowMixedState(!CountACE); recordACBtn->SetState(CountACE ? f0->sh.countACD : 2);
	recordDirBtn->AllowMixedState(!CountDirE); recordDirBtn->SetState(CountDirE ? f0->sh.countDirection : 2);
	angleMapRecordCheckbox->AllowMixedState(!RecordAngleMapE); angleMapRecordCheckbox->SetState(RecordAngleMapE ? f0->sh.recordAngleMap : 2);
	showTexture->AllowMixedState(!TexVisibleE); showTexture->SetState(TexVisibleE ? f0->textureVisible : 2);
	showVolume->AllowMixedState(!VolVisibleE); showVolume->SetState(VolVisibleE ? f0->volumeVisible : 2);
	facetMovingToggle->AllowMixedState(!isMovingE); facetMovingToggle->SetState(isMovingE ? f0->sh.isMoving : 2);
	enableSojournTime->AllowMixedState(!hasSojournE); enableSojournTime->SetState(hasSojournE ? f0->sh.enableSojournTime : 2);

	if (isEnabledE) {
		if (f0->sh.isTextured) { //All facets have textures
			if (ratioE) { //All facets have textures with same resolution
				resolutionText->SetText(f0->tRatio);
				lengthText->SetText(1.0 / f0->tRatio);
			}
			else { //Mixed resolution
				resolutionText->SetText(isEnabledE ? "..." : "");
				lengthText->SetText(isEnabledE ? "..." : "");
			}
		}
		else { //None of the facets have textures
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
		//mApp->facetFlow->SetEditable(!f0->hasOutgassingFile);
		//mApp->facetFlowArea->SetEditable(!f0->hasOutgassingFile);
		if (!f0->hasOutgassingFile) { //All selected DON'T HAVE outgassing maps
			facetUseDesFile->SetSize(1);
			facetUseDesFile->SetSelectedIndex(0); //no map
			facetUseDesFile->SetSelectedValue("No map loaded");
			facetUseDesFile->SetEditable(false);
			fileFluxText->SetText("");
			fileDoseText->SetText("");
			fileYieldText->SetText("");
		}
		else { //All selected HAVE outgassing maps

			facetUseDesFile->SetSize(2);
			facetUseDesFile->SetValueAt(0, "Use user values");
			facetUseDesFile->SetValueAt(1, "Use des. file");
			facetUseDesFile->SetEditable(true);
			if (useOutgMapE) {
				facetUseDesFile->SetSelectedIndex(f0->sh.useOutgassingFile);
			}
			else {
				facetUseDesFile->SetSelectedValue("...");
			}
			char tmp[64];
			if (fluxAEqual) sprintf(tmp, "%.2E", f0->totalFlux / f0->sh.area); else sprintf(tmp, "...");
			fileFluxText->SetText(tmp);
			if (doseAEqual) sprintf(tmp, "%.2E", f0->totalDose / f0->sh.area); else sprintf(tmp, "...");
			fileDoseText->SetText(tmp);
			if (yieldEqual) sprintf(tmp, "%.2E", f0->sh.totalOutgassing / (1.38E-23*f0->sh.temperature) / f0->totalFlux); else sprintf(tmp, "...");
			fileYieldText->SetText(tmp);
			if (useOutgMapE) {
				mApp->facetFlow->SetEditable(!f0->sh.useOutgassingFile);
				mApp->facetFlowArea->SetEditable(!f0->sh.useOutgassingFile);
				if (f0->sh.useOutgassingFile) {
					sprintf(tmp, "%.1E", sumOutgassing * 10.00); //10.00: Pa*m3/s -> mbar*l/s
					mApp->facetFlow->SetText(tmp);
					sprintf(tmp, "%.1E", sumOutgassing * 10.00 / sumArea);
					mApp->facetFlowArea->SetText(tmp);
				}
				else {
					//Let the main program handle this
				}
			}
			else { //some use it, some not
				facetUseDesFile->SetSelectedValue("...");
				mApp->facetFlow->SetEditable(false);
				mApp->facetFlowArea->SetEditable(false);
			}
		}
	}
	else { //Mixed: some have it, some don't
		facetUseDesFile->SetSelectedIndex(0);
		facetUseDesFile->SetSize(1);
		facetUseDesFile->SetSelectedValue("...");
		facetUseDesFile->SetEditable(false);
		fileFluxText->SetText("");
		fileDoseText->SetText("");
		fileYieldText->SetText("");
	}

	if (superDestE) {
		if (f0->sh.superDest == 0) {
			facetSuperDest->SetText("no");
		}
		else {
			facetSuperDest->SetText(f0->sh.superDest);
		}
	}
	else {
		facetSuperDest->SetText("...");
	}
	if (superIdxE) {
		facetStructure->SetText(f0->sh.superIdx + 1);
	}
	else {
		facetStructure->SetText("...");
	}

	if (isMovingE) {
		facetMovingToggle->SetState(f0->sh.isMoving);
		facetMovingToggle->AllowMixedState(false);
	}
	else {
		facetMovingToggle->SetState(2);
		facetMovingToggle->AllowMixedState(true);
	}

	if (AngleMapWidthE) {
		angleMapWidthBox->SetText((int)f0->sh.angleMapPhiWidth);
	}
	else {
		angleMapWidthBox->SetText("...");
	}

	if (AngleMapHeightE) {
		angleMapHeightBox->SetText((int)f0->sh.angleMapThetaHeight);
	}
	else {
		angleMapHeightBox->SetText("...");
	}

	if (enableSojournTime->GetState() == 0) {
		enableSojournTime->SetText("Wall sojourn time");
		sojournFreq->SetEditable(false);
		sojournE->SetEditable(false);
	}
	else {
		sojournFreq->SetEditable(true);
		sojournE->SetEditable(true);
	}

	if (sojournFreqE) {
		sojournFreq->SetText(f0->sh.sojournFreq);
	}
	else {
		sojournFreq->SetText("...");
	}
	if (sojournEE) {
		sojournE->SetText(f0->sh.sojournE);
	}
	else {
		sojournE->SetText("...");
	}
	CalcSojournTime();
	UpdateSize();
}

void FacetAdvParams::Reposition(int wD, int hD) {
	if (wD == 0) wD = this->GetWidth();
	if (hD == 0) hD = this->GetHeight();
	// Position dialog next to Facet parameters
	int facetX, facetY, facetW, facetH;
	mApp->facetPanel->GetBounds(&facetX, &facetY, &facetW, &facetH);
	SetBounds(facetX - wD - 10, facetY + 20, wD, hD);
}


bool FacetAdvParams::ApplyTexture(bool force) {
	bool boundMap = true; // boundaryBtn->GetState();
	double ratio = 0.0;
	std::vector<size_t> selectedFacets = geom->GetSelectedFacets();
	int nbPerformed = 0;
	bool doRatio = false;
	if (enableBtn->GetState() == 1) { //check if valid texture settings are to be applied

									  // Check counting mode
		if (!recordDesBtn->GetState() && !recordAbsBtn->GetState() &&
			!recordReflBtn->GetState() && !recordTransBtn->GetState() &&
			!recordACBtn->GetState() && !recordDirBtn->GetState()) {
			GLMessageBox::Display("Please select counting mode", "Error", GLDLG_OK, GLDLG_ICONINFO);
			return false;
		}

		// Resolution
		if (resolutionText->GetNumber(&ratio) && ratio >= 0.0) {
			//Got a valid number
			doRatio = true;
		}
		else if (strcmp(resolutionText->GetText(), "...") != 0) { //Not in mixed "..." state
			GLMessageBox::Display("Invalid texture resolution\nMust be a non-negative number", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return false;
		}
		else {
			//Mixed state: leave doRatio as false
		}
	}

	if (!mApp->AskToReset(worker)) return false;
	progressDlg = new GLProgress("Applying mesh settings", "Please wait");
	progressDlg->SetVisible(true);
	progressDlg->SetProgress(0.0);
	int count = 0;
	for (auto sel : selectedFacets) {
		Facet *f = geom->GetFacet(sel);
		bool hadAnyTexture = f->sh.countDes || f->sh.countAbs || f->sh.countRefl || f->sh.countTrans || f->sh.countACD || f->sh.countDirection;
		bool hadDirCount = f->sh.countDirection;

		if (enableBtn->GetState() == 0 || ratio == 0.0) {
			//Let the user disable textures with the main switch or by typing 0 as resolution
			f->sh.countDes = f->sh.countAbs = f->sh.countRefl = f->sh.countTrans = f->sh.countACD = f->sh.countDirection = false;
		}
		else {
			if (recordDesBtn->GetState() < 2) f->sh.countDes = recordDesBtn->GetState();
			if (recordAbsBtn->GetState() < 2) f->sh.countAbs = recordAbsBtn->GetState();
			if (recordReflBtn->GetState() < 2) f->sh.countRefl = recordReflBtn->GetState();
			if (recordTransBtn->GetState() < 2) f->sh.countTrans = recordTransBtn->GetState();
			if (recordACBtn->GetState() < 2) f->sh.countACD = recordACBtn->GetState();
			if (recordDirBtn->GetState() < 2) f->sh.countDirection = recordDirBtn->GetState();
		}

		bool hasAnyTexture = f->sh.countDes || f->sh.countAbs || f->sh.countRefl || f->sh.countTrans || f->sh.countACD || f->sh.countDirection;

		//set textures
		try {
			bool needsRemeshing = force || (hadAnyTexture != hasAnyTexture) || (hadDirCount != f->sh.countDirection) || (doRatio && (!IS_ZERO(geom->GetFacet(sel)->tRatio - ratio)));
			if (needsRemeshing) geom->SetFacetTexture(sel, hasAnyTexture ? ratio : 0.0, hasAnyTexture ? boundMap : false);
		}
		catch (Error &e) {
			GLMessageBox::Display((char *)e.GetMsg(), "Error", GLDLG_OK, GLDLG_ICONWARNING);
			progressDlg->SetVisible(false);
			SAFE_DELETE(progressDlg);
			return false;
		}
		catch (...) {
			GLMessageBox::Display("Unexpected error while setting textures", "Error", GLDLG_OK, GLDLG_ICONWARNING);
			progressDlg->SetVisible(false);
			SAFE_DELETE(progressDlg);
			return false;
		}
		nbPerformed++;
		progressDlg->SetProgress((double)nbPerformed / (double)selectedFacets.size());
	} //main cycle end

	if (progressDlg) progressDlg->SetVisible(false);
	SAFE_DELETE(progressDlg);
	return true;
}

bool FacetAdvParams::Apply() {
	std::vector<size_t> selectedFacets=geom->GetSelectedFacets();
	int nbPerformed = 0;
	/*
	bool boundMap = true; // boundaryBtn->GetState();
	double ratio = 0.0;
	
	
	bool doRatio = false;
	if (enableBtn->GetState() == 1) { //check if valid texture settings are to be applied

		// Check counting mode
		if (!recordDesBtn->GetState() && !recordAbsBtn->GetState() &&
			!recordReflBtn->GetState() && !recordTransBtn->GetState() &&
			!recordACBtn->GetState() && !recordDirBtn->GetState()) {
			GLMessageBox::Display("Please select counting mode", "Error", GLDLG_OK, GLDLG_ICONINFO);
			return false;
		}

		// Resolution
		if (resolutionText->GetNumber(&ratio) && ratio >= 0.0) {
			//Got a valid number
			doRatio = true;
		}
		else if (strcmp(resolutionText->GetText(), "...") != 0)  { //Not in mixed "..." state
			GLMessageBox::Display("Invalid texture resolution\nMust be a non-negative number", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return false;
		}
		else {
			//Mixed state: leave doRatio as false
		}
	}*/

	// Superstructure

	bool structChanged = false; //if a facet gets into a new structure, we have to re-render the geometry

	int superStruct;
	bool doSuperStruct = false;
	if (sscanf(facetStructure->GetText(), "%d", &superStruct) > 0 && superStruct > 0 && superStruct <= geom->GetNbStructure()) doSuperStruct = true;
	else {
		if (strcmp(facetStructure->GetText(), "...") == 0) doSuperStruct = false;
		else{
			GLMessageBox::Display("Invalid superstructre number", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return false;
		}
	}


	// Super structure destination (link)
	int superDest;
	bool doLink = false;
	if (strcmp(facetSuperDest->GetText(), "none") == 0 || strcmp(facetSuperDest->GetText(), "no") == 0 || strcmp(facetSuperDest->GetText(), "0") == 0) {
		doLink = true;
		superDest = 0;
	}
	else if (sscanf(facetSuperDest->GetText(), "%d", &superDest) > 0) {
		if (superDest == superStruct) {
			GLMessageBox::Display("Link and superstructure can't be the same", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return false;
		}
		else if (superDest < 0 || superDest > geom->GetNbStructure()) {
			GLMessageBox::Display("Link destination points to a structure that doesn't exist", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return false;
		}
		else
			doLink = true;
	}
	else if (strcmp(facetSuperDest->GetText(), "...") == 0) doLink = false;
	else {
		GLMessageBox::Display("Invalid superstructure destination", "Error", GLDLG_OK, GLDLG_ICONERROR);
		return false;
	}

	// teleport
	int teleport;
	bool doTeleport = false;

	if (facetTeleport->GetNumberInt(&teleport)) {
		if (teleport<-1 || teleport>geom->GetNbFacet()) {
			GLMessageBox::Display("Invalid teleport destination\n(If no teleport: set number to 0)", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return false;
		}
		else if (teleport > 0 && geom->GetFacet(teleport - 1)->selected) {
			char tmp[256];
			sprintf(tmp, "The teleport destination of facet #%d can't be itself!", teleport);
			GLMessageBox::Display(tmp, "Error", GLDLG_OK, GLDLG_ICONERROR);
			return false;
		}
		doTeleport = true;
	}
	else {
		if (strcmp(facetTeleport->GetText(), "...") == 0) doTeleport = false;
		else {
			GLMessageBox::Display("Invalid teleport destination\n(If no teleport: set number to 0)", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return false;
		}
	}

	// temp.accomodation factor
	double accfactor;
	bool doAccfactor = false;
	if (facetAccFactor->GetNumber(&accfactor)) {
		if (accfactor<0.0 || accfactor>1.0) {
			GLMessageBox::Display("Facet accomodation factor must be between 0 and 1", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return false;
		}
		doAccfactor = true;
	}
	else {
		if (strcmp(facetAccFactor->GetText(), "...") == 0) doAccfactor = false;
		else {
			GLMessageBox::Display("Invalid accomodation factor number", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return false;
		}
	}

	// Use desorption map
	int useMapA = 0;
	bool doUseMapA = false;
	if (strcmp(facetUseDesFile->GetSelectedValue(), "...") == 0) doUseMapA = false;
	else {
		useMapA = (facetUseDesFile->GetSelectedIndex() == 1);
		bool missingMap = false;
		int missingMapId;
		if (useMapA) {
			for (int i = 0; i < geom->GetNbFacet(); i++) {
				if (geom->GetFacet(i)->selected && !geom->GetFacet(i)->hasOutgassingFile) {
					missingMap = true;
					missingMapId = i;
				}
			}
		}
		if (missingMap) {
			char tmp[256];
			sprintf(tmp, "%d is selected but doesn't have any outgassing map loaded.", missingMapId + 1);
			GLMessageBox::Display(tmp, "Can't use map on all facets", GLDLG_OK, GLDLG_ICONERROR);
			return false;
		}
		doUseMapA = true;
	}

	// sojourn time coefficient 1
	double sojF;
	bool doSojF = false;

	if (sojournFreq->GetNumber(&sojF)) {
		if (sojF <= 0.0) {
			GLMessageBox::Display("Wall sojourn time frequency has to be positive", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return false;
		}
		doSojF = true;
	}
	else {
		if (enableSojournTime->GetState() == 0 || strcmp(sojournFreq->GetText(), "...") == 0) doSojF = false;
		else {
			GLMessageBox::Display("Invalid wall sojourn time frequency", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return false;
		}
	}

	// sojourn time coefficient 2
	double sojE;
	bool doSojE = false;

	if (sojournE->GetNumber(&sojE)) {
		if (sojE <= 0.0) {
			GLMessageBox::Display("Wall sojourn time second coefficient (Energy) has to be positive", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return false;
		}
		doSojE = true;
	}
	else {
		if (enableSojournTime->GetState() == 0 || strcmp(sojournE->GetText(), "...") == 0) doSojE = false;
		else {
			GLMessageBox::Display("Invalid wall sojourn time second coefficient (Energy)", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return false;
		}
	}

	// angle map width
	int angleMapWidth;
	bool doAngleMapWidth = false;

	if (angleMapRecordCheckbox->GetState() == 0 || strcmp(angleMapWidthBox->GetText(), "...") == 0) doAngleMapWidth = false;
	else if (angleMapWidthBox->GetNumberInt(&angleMapWidth)) {
		if (angleMapWidth <= 0) {
			GLMessageBox::Display("Angle map width has to be positive", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return false;
		}
		doAngleMapWidth = true;
	}
	else {
			GLMessageBox::Display("Invalid angle map width", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return false;
	}
	
	// angle map height
	int angleMapHeight;
	bool doAngleMapHeight = false;

	if (angleMapRecordCheckbox->GetState() == 0 || strcmp(angleMapHeightBox->GetText(), "...") == 0) doAngleMapHeight = false;
	else if (angleMapHeightBox->GetNumberInt(&angleMapHeight)) {
		if (angleMapHeight <= 0) {
			GLMessageBox::Display("Angle map height has to be positive", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return false;
		}
		doAngleMapHeight = true;
	}
	else {
		GLMessageBox::Display("Invalid angle map height", "Error", GLDLG_OK, GLDLG_ICONERROR);
		return false;
	}
	

	// Reflection type
	int reflType = facetReflType->GetSelectedIndex();

	//Check complete, let's apply
	//First applying angle map recording before a reset is done
	int angleMapState = angleMapRecordCheckbox->GetState();
	if (angleMapState < 2) {
		for (auto sel:selectedFacets) {
			geom->GetFacet(sel)->sh.recordAngleMap=angleMapState;
		}
	}
	
	if (!mApp->AskToReset(worker)) return false;
	progressDlg = new GLProgress("Applying facet parameters", "Please wait");
	progressDlg->SetVisible(true);
	progressDlg->SetProgress(0.0);
	int count = 0;
	for (auto sel:selectedFacets) {
		Facet *f = geom->GetFacet(sel);
		/*
		bool hadAnyTexture = f->sh.countDes || f->sh.countAbs || f->sh.countRefl || f->sh.countTrans || f->sh.countACD || f->sh.countDirection;
		bool hadDirCount = f->sh.countDirection;
		
		if (enableBtn->GetState() == 0 || ratio == 0.0) {
			//Let the user disable textures with the main switch or by typing 0 as resolution
			f->sh.countDes = f->sh.countAbs = f->sh.countRefl = f->sh.countTrans = f->sh.countACD = f->sh.countDirection = false;
		}
		else {
			if (recordDesBtn->GetState() < 2) f->sh.countDes = recordDesBtn->GetState();
			if (recordAbsBtn->GetState() < 2) f->sh.countAbs = recordAbsBtn->GetState();
			if (recordReflBtn->GetState() < 2) f->sh.countRefl = recordReflBtn->GetState();
			if (recordTransBtn->GetState() < 2) f->sh.countTrans = recordTransBtn->GetState();
			if (recordACBtn->GetState() < 2) f->sh.countACD = recordACBtn->GetState();
			if (recordDirBtn->GetState() < 2) f->sh.countDirection = recordDirBtn->GetState();
		}
		
		bool hasAnyTexture = f->sh.countDes || f->sh.countAbs || f->sh.countRefl || f->sh.countTrans || f->sh.countACD || f->sh.countDirection;
		*/

		if (doTeleport) f->sh.teleportDest = teleport;
		if (doAccfactor) f->sh.accomodationFactor = accfactor;
		if (reflType >= 0) f->sh.reflectType = reflType;
		if (doSuperStruct) {
			if (f->sh.superIdx != (superStruct - 1)) {
				f->sh.superIdx = superStruct - 1;
				structChanged = true;
			}
		}
		if (doLink) {
			f->sh.superDest = superDest;
			if (superDest) f->sh.opacity = 1.0; // Force opacity for link facet
		}
		//Moving or not
		if (facetMovingToggle->GetState() < 2) f->sh.isMoving = facetMovingToggle->GetState();
		if (enableSojournTime->GetState() < 2) f->sh.enableSojournTime = enableSojournTime->GetState();
		if (doSojF) f->sh.sojournFreq = sojF;
		if (doSojE) f->sh.sojournE = sojE;

		if (doUseMapA) {
			f->sh.useOutgassingFile = useMapA;
		}

		if (doAngleMapWidth) {
			if (angleMapWidth != f->sh.angleMapPhiWidth) {
				SAFE_FREE(f->angleMapCache);
				f->sh.hasRecordedAngleMap = false;
				f->sh.angleMapPhiWidth = angleMapWidth;
			}
		}
		if (doAngleMapHeight) {
			if (angleMapHeight != f->sh.angleMapThetaHeight) {
				SAFE_FREE(f->angleMapCache);
				f->sh.hasRecordedAngleMap = false;
				f->sh.angleMapThetaHeight = angleMapHeight;
			}
		}

		/*
		//set textures
		try {
			bool needsRemeshing = (hadAnyTexture != hasAnyTexture) || (hadDirCount != f->sh.countDirection) || (doRatio && (!IsEqual(geom->GetFacet(sel)->tRatio , ratio)));
			if (needsRemeshing) geom->SetFacetTexture(sel, hasAnyTexture ? ratio : 0.0, hasAnyTexture ? boundMap : false);
		}
		catch (Error &e) {
			GLMessageBox::Display((char *)e.GetMsg(), "Error", GLDLG_OK, GLDLG_ICONWARNING);
			progressDlg->SetVisible(false);
			SAFE_DELETE(progressDlg);
			return false;
		}
		catch (...) {
			GLMessageBox::Display("Unexpected error while setting textures", "Error", GLDLG_OK, GLDLG_ICONWARNING);
			progressDlg->SetVisible(false);
			SAFE_DELETE(progressDlg);
			return false;
		}
		*/
		//if (angleMapRecordCheckbox->GetState() < 2) f->sh.recordAngleMap = angleMapRecordCheckbox->GetState();
		if (showTexture->GetState() < 2) f->textureVisible = showTexture->GetState();
		if (showVolume->GetState() < 2) f->volumeVisible = showVolume->GetState();
		nbPerformed++;
		progressDlg->SetProgress((double)nbPerformed / (double)selectedFacets.size());
	} //main cycle end
	if (structChanged) geom->BuildGLList(); //Re-render facets

	if (progressDlg) progressDlg->SetVisible(false);
	SAFE_DELETE(progressDlg);

	return ApplyTexture(); //Finally, apply textures
}


void FacetAdvParams::ApplyDrawSettings() {
	//Apply view settings without stopping the simulation


	double nbSelected = (double)geom->GetNbSelectedFacets();
	double nbPerformed = 0.0;

	for (int i = 0; i < geom->GetNbFacet(); i++) {

		Facet *f = geom->GetFacet(i);
		if (f->selected) {

			if (showTexture->GetState() < 2) f->textureVisible = showTexture->GetState();
			if (showVolume->GetState() < 2) f->volumeVisible = showVolume->GetState();

			nbPerformed += 1.0;
			progressDlg->SetProgress(nbPerformed / nbSelected);
		}

	}
	geom->BuildGLList(); //Re-render facets
}


void FacetAdvParams::UpdateToggle(GLComponent *src) {

	/*if (src==boundaryBtn) {
		recordACBtn->SetState(false);
		} else if(src==enableBtn) {
		//boundaryBtn->SetState(enableBtn->GetState());
		} else */

	if (src == recordDesBtn) {
		enableBtn->SetState(true);
		//boundaryBtn->SetState(true);
		recordACBtn->SetState(false);
	}
	else if (src == recordAbsBtn) {
		enableBtn->SetState(true);
		//boundaryBtn->SetState(true);
		recordACBtn->SetState(false);
	}
	else if (src == recordReflBtn) {
		enableBtn->SetState(true);
		//boundaryBtn->SetState(true);
		recordACBtn->SetState(false);
	}
	else if (src == recordTransBtn) {
		enableBtn->SetState(true);
		//boundaryBtn->SetState(true);
		recordACBtn->SetState(false);
	}
	else if (src == recordDirBtn) {
		enableBtn->SetState(true);
		//boundaryBtn->SetState(true);
		recordACBtn->SetState(false);
	}
	else if (src == recordACBtn) {
		if (recordACBtn->GetState()) {
			enableBtn->SetState(true);
			//boundaryBtn->SetState(true);
			recordDesBtn->SetState(false);
			recordAbsBtn->SetState(false);
			recordReflBtn->SetState(false);
			recordTransBtn->SetState(false);
			recordDirBtn->SetState(false);
		}
	}
	else if (src == enableSojournTime) {
		sojournFreq->SetEditable(enableSojournTime->GetState());
		sojournE->SetEditable(enableSojournTime->GetState());
		if (enableSojournTime->GetState() == 0) {
			enableSojournTime->SetText("Wall sojourn time");
		}
		else {
			CalcSojournTime();
		}
	}
	else if (src == angleMapRecordCheckbox) {

	}

	//UpdateSizeForRatio();
}



void FacetAdvParams::ProcessMessage(GLComponent *src, int message) {

	switch (message) {

		// -------------------------------------------------------------
	case MSG_BUTTON:
		if (src == quickApply) {

			progressDlg = new GLProgress("Applying view settings", "Please wait");
			progressDlg->SetVisible(true);
			progressDlg->SetProgress(0.5);

			ApplyDrawSettings();

			progressDlg->SetVisible(false);
			SAFE_DELETE(progressDlg);

		}
		
		else if (src==SojournInfoButton) {

			char tmp[] = "f: Molecule's surface oscillation frequency [Hz]\n"
				"E: Adsorption energy [J/mole]\n"
				"A: Escape probability per oscillation:\n"
				"A = exp(-E/(R*T))\n\n"
				"Probability of sojourn time t:\n"
				"p(t)= A*f*exp(-A*f*t)\n\n"
				"Mean sojourn time:\n"
				"mean= 1/(A*f) = 1/f*exp(E/(kT))\n";
			GLMessageBox::Display(tmp, "Wall sojourn time", GLDLG_OK, GLDLG_ICONINFO);

		}
		else if (src == angleMapExportButton) {
			mApp->ExportAngleMaps();
		}
		else if (src == angleMapReleaseButton) {
			mApp->ClearAngleMapsOnSelection();
		}
		else if (src == remeshButton) {
			ApplyTexture(true);
		}
		
		break;

		// -------------------------------------------------------------
	case MSG_TEXT_UPD:
		mApp->facetApplyBtn->SetEnabled(true);
		if (src == resolutionText) {
			enableBtn->SetState(true);
			UpdateSizeForRatio();
			double res;
			if (resolutionText->GetNumber(&res) && res != 0.0)
				lengthText->SetText(1.0 / res);
			else
				lengthText->SetText("");
		}
		else if (src == lengthText) {
			enableBtn->SetState(true);
			double length;
			if (lengthText->GetNumber(&length) && length != 0.0) {
				resolutionText->SetText(1.0 / length);
				UpdateSizeForRatio();
			}
			else
				resolutionText->SetText("");
		}
		else if (src == sojournFreq || src == sojournE) {
			CalcSojournTime();
		}

		break;

		// -------------------------------------------------------------
	case MSG_TOGGLE:
		UpdateToggle(src);
		mApp->facetApplyBtn->SetEnabled(true);
		break;
	case MSG_TEXT:
		/*if (
			src == resolutionText
			|| src == lengthText
			|| src == facetAccFactor
			|| src == facetTeleport
			|| src == facetSuperDest
			|| src == facetStructure
			|| src == sojournFreq
			|| src == sojournE
			) {*/
			mApp->ApplyFacetParams();
		//}
		break;
	case MSG_COMBO:
		if (src == facetReflType) {
			mApp->facetApplyBtn->SetEnabled(true);
		}
		else if (src == facetUseDesFile) {
			mApp->facetApplyBtn->SetEnabled(true);
			if (facetUseDesFile->GetSelectedIndex() == 0) {
				//User values
				mApp->facetFlow->SetEditable(true);
				mApp->facetFlowArea->SetEditable(true);
			}
			else { //use desorption file
				mApp->facetFlow->SetEditable(false);
				mApp->facetFlowArea->SetEditable(false);
				//Values from last Refresh();
				char tmp[64];
				sprintf(tmp, "%.2E", sumOutgassing);
				mApp->facetFlow->SetText(tmp);
				sprintf(tmp, "%.2E", sumOutgassing / sumArea);
				mApp->facetFlowArea->SetText(tmp);
			}
		}
		break;
	}

	GLWindow::ProcessMessage(src, message);
}

void FacetAdvParams::CalcSojournTime() {
	double sojF,sojE,facetT;
	if (enableSojournTime->GetState() == 0
		|| !(sojournFreq->GetNumber(&sojF))
		|| !(sojournE->GetNumber(&sojE))
		|| !(mApp->facetTemperature->GetNumber(&facetT))) {
		enableSojournTime->SetText("Wall sojourn time");
		return;
	}
	std::ostringstream tmp;
	tmp<< "Wall sojourn time (mean=" << 1.0/(sojF*exp(-sojE /(8.31* facetT))) << " s)";
	enableSojournTime->SetText(tmp.str());
}