/*
File:        PressureEvolution.cpp
Description: Pressure Evolution plotter window
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

#include "PressureEvolution.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLMessageBox.h"
#include "Utils.h"
#include <math.h>
#include "Molflow.h"

extern GLApplication *theApp;

extern double gasMass;
extern double totalOutgassing;
extern double totalInFlux;

static const char*profType[] = { "None", "Pressure \201 [mbar]", "Pressure \202 [mbar]", "Angle", "Velocity" };

PressureEvolution::PressureEvolution() :GLWindow() {

	int wD = 750;
	int hD = 425;

	SetTitle("Pressure evolution plotter");
	SetIconfiable(TRUE);
	nbView = 0;
	worker = NULL;
	lastUpdate = 0.0f;

	nbColors = 8;
	colors[0] = new GLCColor(); colors[0]->r = 255; colors[0]->g = 000; colors[0]->b = 055; //red
	colors[1] = new GLCColor(); colors[1]->r = 000; colors[1]->g = 000; colors[1]->b = 255; //blue
	colors[2] = new GLCColor(); colors[2]->r = 000; colors[2]->g = 204; colors[2]->b = 051; //green
	colors[3] = new GLCColor(); colors[3]->r = 000; colors[3]->g = 000; colors[3]->b = 000; //black
	colors[4] = new GLCColor(); colors[4]->r = 255; colors[4]->g = 153; colors[4]->b = 051; //orange
	colors[5] = new GLCColor(); colors[5]->r = 153; colors[5]->g = 204; colors[5]->b = 255; //light blue
	colors[6] = new GLCColor(); colors[6]->r = 153; colors[6]->g = 000; colors[6]->b = 102; //violet
	colors[7] = new GLCColor(); colors[7]->r = 255; colors[7]->g = 230; colors[7]->b = 005; //yellow

	selectedSlice = 49;

	chart = new GLChart(0);
	chart->SetBorder(BORDER_BEVEL_IN);
	chart->GetY1Axis()->SetGridVisible(TRUE);
	chart->GetXAxis()->SetGridVisible(TRUE);
	chart->GetXAxis()->SetName("Time (s)");
	chart->GetY1Axis()->SetAutoScale(TRUE);
	chart->GetY2Axis()->SetAutoScale(TRUE);
	chart->GetY1Axis()->SetAnnotation(VALUE_ANNO);
	chart->GetXAxis()->SetAnnotation(VALUE_ANNO);
	Add(chart);

	dismissButton = new GLButton(0, "Dismiss");
	Add(dismissButton);

	selButton = new GLButton(0, "Show Facet");
	Add(selButton);

	addButton = new GLButton(0, "Add curve");
	Add(addButton);

	removeButton = new GLButton(0, "Remove curve");
	Add(removeButton);

	resetButton = new GLButton(0, "Remove all");
	Add(resetButton);

	profCombo = new GLCombo(0);
	profCombo->SetEditable(TRUE);
	Add(profCombo);

	logXToggle = new GLToggle(0, "Log X");
	Add(logXToggle);

	logYToggle = new GLToggle(0, "Log Y");
	Add(logYToggle);

	normLabel = new GLLabel("Normalize");
	Add(normLabel);
	//qLabel = new GLLabel("Q=");
	//Add(qLabel);
	//unitLabel = new GLLabel("units*l/s");
	//Add(unitLabel);

	normCombo = new GLCombo(0);
	normCombo->SetEditable(TRUE);
	normCombo->SetSize(8);
	normCombo->SetValueAt(0, "None");
	normCombo->SetValueAt(1, "Absorption (total)");
	normCombo->SetValueAt(2, "Desorption (total)");
	normCombo->SetValueAt(3, "Hit           (total)");
	normCombo->SetValueAt(4, "Absorption (local)");
	normCombo->SetValueAt(5, "Desorption (local)");
	normCombo->SetValueAt(6, "Hit           (local)");
	normCombo->SetValueAt(7, "Pressure [mbar]");
	normCombo->SetSelectedIndex(7);
	Add(normCombo);

	modeCombo = new GLCombo(0);
	modeCombo->SetEditable(TRUE);
	modeCombo->SetSize(2);
	modeCombo->SetValueAt(0, "Average");
	modeCombo->SetValueAt(1, "Slice #");
	modeCombo->SetSelectedIndex(0);
	Add(modeCombo);

	selectedSliceText = new GLTextField(0, "50");
	selectedSliceText->SetEditable(FALSE);
	Add(selectedSliceText);

	label1 = new GLLabel("Show evolution of:");
	Add(label1);

	formulaText = new GLTextField(0, "");
	formulaText->SetEditable(TRUE);
	Add(formulaText);

	formulaBtn = new GLButton(0, "-> Plot");
	Add(formulaBtn);

	// Center dialog
	int wS, hS;
	GLToolkit::GetScreenSize(&wS, &hS);
	int xD = (wS - wD) / 2;
	int yD = (hS - hD) / 2;
	SetBounds(xD, yD, wD, hD);
	SetResizable(TRUE);
	SetMinimumSize(wD, 220);

	RestoreDeviceObjects();

}

void PressureEvolution::SetBounds(int x, int y, int w, int h) {

	chart->SetBounds(7, 5, w - 15, h - 110);
	profCombo->SetBounds(7, h - 95, 117, 19);
	selButton->SetBounds(130, h - 95, 80, 19);
	addButton->SetBounds(215, h - 95, 80, 19);
	removeButton->SetBounds(300, h - 95, 80, 19);
	resetButton->SetBounds(385, h - 95, 80, 19);

	logXToggle->SetBounds(w - 110, h - 95, 40, 19);
	logYToggle->SetBounds(w - 55, h - 95, 40, 19);
	label1->SetBounds(w - 265, h - 70, 90, 19);
	normLabel->SetBounds(7, h - 68, 50, 19);
	normCombo->SetBounds(60, h - 70, 105, 19);

	modeCombo->SetBounds(w - 170, h - 70, 105, 19);
	selectedSliceText->SetBounds(w - 60, h - 70, 50, 19);

	formulaText->SetBounds(7, h - 45, 525, 19);
	formulaBtn->SetBounds(537, h - 45, 80, 19);
	dismissButton->SetBounds(w - 100, h - 45, 90, 19);

	GLWindow::SetBounds(x, y, w, h);

}

void PressureEvolution::Refresh() {

	if (!worker) return;

	Geometry *geom = worker->GetGeometry();
	int nb = geom->GetNbFacet();
	int nbProf = 0;
	for (int i = 0; i < nb; i++)
	if (geom->GetFacet(i)->sh.isProfile) nbProf++;
	profCombo->Clear();
	if (nbProf) profCombo->SetSize(nbProf);
	nbProf = 0;
	for (int i = 0; i < nb; i++) {
		Facet *f = geom->GetFacet(i);
		if (f->sh.isProfile) {
			char tmp[128];
			sprintf(tmp, "F#%d %s", i + 1, profType[f->sh.profileType]);
			profCombo->SetValueAt(nbProf, tmp, i);
			profCombo->SetSelectedIndex(0);
			nbProf++;
		}
	}
	//Remove profiles that aren't present anymore
	for (int v = 0; v < nbView; v++)
	if (views[v]->userData >= geom->GetNbFacet() || !geom->GetFacet(views[v]->userData)->sh.isProfile) {
		chart->GetY1Axis()->RemoveDataView(views[v]);
		SAFE_DELETE(views[v]);
		for (int j = v; j < nbView - 1; j++) views[j] = views[j + 1];
		nbView--;
	}
	refreshViews();

}

void PressureEvolution::Display(Worker *w) {

	worker = w;
	Refresh();
	SetVisible(TRUE);

}

void PressureEvolution::Update(float appTime, BOOL force) {

	if (!IsVisible() || IsIconic()) return;

	if (force) {
		refreshViews();
		lastUpdate = appTime;
		return;
	}

	if ((appTime - lastUpdate > 1.0f || force) && nbView) {
		if (worker->running) refreshViews();
		lastUpdate = appTime;
	}

}

void PressureEvolution::plot() {

	GLParser *parser = new GLParser();
	parser->SetExpression(formulaText->GetText());
	if (!parser->Parse()) {
		GLMessageBox::Display(parser->GetErrorMsg(), "Error", GLDLG_OK, GLDLG_ICONERROR);
		SAFE_DELETE(parser);
		return;
	}

	int nbVar = parser->GetNbVariable();
	if (nbVar == 0) {
		GLMessageBox::Display("Variable 'x' not found", "Error", GLDLG_OK, GLDLG_ICONERROR);
		SAFE_DELETE(parser);
		return;
	}
	if (nbVar > 1) {
		GLMessageBox::Display("Too much variables or unknown constant", "Error", GLDLG_OK, GLDLG_ICONERROR);
		SAFE_DELETE(parser);
		return;
	}
	VLIST *var = parser->GetVariableAt(0);
	if (_stricmp(var->name, "x") != 0) {
		GLMessageBox::Display("Variable 'x' not found", "Error", GLDLG_OK, GLDLG_ICONERROR);
		SAFE_DELETE(parser);
		return;
	}

	Geometry *geom = worker->GetGeometry();
	GLDataView *v;

	// Check that view is not already added
	BOOL found = FALSE;
	int i = 0;
	while (i < nbView && !found) {
		found = (views[i]->userData == -1);
		if (!found) i++;
	}

	if (found) {
		v = views[i];
		v->SetName(formulaText->GetText());
		v->Reset();
	}
	else {
		if (nbView < 32) {
			v = new GLDataView();
			v->SetName(formulaText->GetText());
			v->userData = -1;
			chart->GetY1Axis()->AddDataView(v);
			views[nbView] = v;
			nbView++;
		}
	}

	// Plot
	for (int p = 0; p < 1000; p++) {
		double x = (double)p;
		double y;
		var->value = x;
		parser->Evaluate(&y);
		v->Add(x, y, FALSE);
	}
	v->CommitChange();

	delete parser;

}

void PressureEvolution::refreshViews() {

	// Lock during update
	BYTE *buffer = worker->GetHits();
	int normalize = normCombo->GetSelectedIndex();
	if (!buffer) return;

	Geometry *geom = worker->GetGeometry();
	SHGHITS *gHits = (SHGHITS *)buffer;
	double nbAbs = (double)gHits->total.hit.nbAbsorbed;
	double nbDes = (double)gHits->total.hit.nbDesorbed;
	double nbHit = (double)gHits->total.hit.nbHit;

	double scale;

	for (int i = 0; i < nbView; i++) {

		GLDataView *v = views[i];
		if (v->userData >= 0 && v->userData < geom->GetNbFacet()) {
			Facet *f = geom->GetFacet(v->userData);
			SHHITS *fCount = (SHHITS *)(buffer + f->sh.hitOffset);
			double fnbAbs = (double)fCount->hit.nbAbsorbed;
			double fnbDes = (double)fCount->hit.nbDesorbed;
			double fnbHit = (double)fCount->hit.nbHit;
			//double q;
			v->Reset();
			for (int m = 1; m <= MIN((int)worker->moments.size(), 10000); m++) { //max 10000 points
				APROFILE *profilePtr = (APROFILE *)(buffer + f->sh.hitOffset + sizeof(SHHITS)+m*sizeof(APROFILE)*PROFILE_SIZE);
				double val;
				if (modeCombo->GetSelectedIndex() == 1) //plot one slice
					val = (double)profilePtr[selectedSlice].sum_v_ort;
				else { //plot sum/average
					double valLLong = 0.0;
					for (int j = 0; j < PROFILE_SIZE; j++)
						valLLong += profilePtr[j].sum_v_ort;
					val = (double)valLLong;
					if (normCombo->GetSelectedIndex() == 7) //pressure, calculate average instead of sum
						val /= (double)PROFILE_SIZE;
				}

				switch (normalize) {
				case 0:
					v->Add(worker->moments[m - 1], val, FALSE);
					break;
				case 1:
					v->Add(worker->moments[m - 1], val / nbAbs, FALSE);
					break;
				case 2:
					v->Add(worker->moments[m - 1], val / nbDes, FALSE);
					break;
				case 3:
					v->Add(worker->moments[m - 1], val / nbHit, FALSE);
					break;
				case 4:
					v->Add(worker->moments[m - 1], val / fnbAbs, FALSE);
					break;
				case 5:
					v->Add(worker->moments[m - 1], val / fnbDes, FALSE);
					break;
				case 6:
					v->Add(worker->moments[m - 1], val / fnbHit, FALSE);
					break;
				case 7: //Pressure
					scale = totalInFlux / nbDes / (f->sh.area / (double)PROFILE_SIZE*1E-4)* gasMass / 1000 / 6E23 *0.0100; //0.01: Pa->mbar
					scale *= (worker->desorptionStopTime - worker->desorptionStartTime)
						/ worker->timeWindowSize; //correction for time window length
					if (f->sh.is2sided) scale *= 0.5;
					v->Add(worker->moments[m - 1], val*scale, FALSE);
					break;
				}
			}
			v->CommitChange();
		}
		else {
			if (v->userData == -2 && nbDes != 0.0) {

				// Volatile profile
				v->Reset();
				int nb = geom->GetNbFacet();
				for (int j = 0; j < nb; j++) {
					Facet *f = geom->GetFacet(j);
					if (f->sh.isVolatile) {
						SHHITS *fCount = (SHHITS *)(buffer + f->sh.hitOffset);
						double z = geom->GetVertex(f->indices[0])->z;
						v->Add(z, (double)(fCount->hit.nbAbsorbed) / nbDes, FALSE);
					}
				}
				// Last
				Facet *f = geom->GetFacet(28);
				SHHITS *fCount = (SHHITS *)(buffer + f->sh.hitOffset);
				double fnbAbs = (double)fCount->hit.nbAbsorbed;
				v->Add(1000.0, fnbAbs / nbDes, FALSE);
				v->CommitChange();

				//v->Reset();
				//for(int j=0;j<BOUNCEMAX && nbAbs;j++)
				//  v->Add((double)j,(double)gHits->wallHits[j]/nbAbs);

			}
		}

	}

	worker->ReleaseHits();

}


void PressureEvolution::addView(int facet) {

	char tmp[128];
	Geometry *geom = worker->GetGeometry();

	// Check that view is not already added
	BOOL found = FALSE;
	int i = 0;
	while (i<nbView && !found) {
		found = (views[i]->userData == facet);
		if (!found) i++;
	}
	if (worker->moments.size()>10000) {
		GLMessageBox::Display("Only the first 10000 moments will be plotted", "Error", GLDLG_OK, GLDLG_ICONWARNING);
	}
	if (found) {
		GLMessageBox::Display("Profile already plotted", "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}
	if (nbView < 32) {
		Facet *f = geom->GetFacet(facet);
		GLDataView *v = new GLDataView();
		sprintf(tmp, "F#%d %s", facet + 1, profType[f->sh.profileType]);
		v->SetName(tmp);
		v->SetViewType(TYPE_BAR);
		v->SetMarker(MARKER_DOT);
		v->SetColor(*colors[nbView%nbColors]);
		v->SetMarkerColor(*colors[nbView%nbColors]);
		v->userData = facet;
		chart->GetY1Axis()->AddDataView(v);
		views[nbView] = v;
		nbView++;
	}

}

void PressureEvolution::remView(int facet) {

	Geometry *geom = worker->GetGeometry();

	BOOL found = FALSE;
	int i = 0;
	while (i < nbView && !found) {
		found = (views[i]->userData == facet);
		if (!found) i++;
	}
	if (!found) {
		GLMessageBox::Display("Profile not plotted", "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}
	chart->GetY1Axis()->RemoveDataView(views[i]);
	SAFE_DELETE(views[i]);
	for (int j = i; j < nbView - 1; j++) views[j] = views[j + 1];
	nbView--;

}

void PressureEvolution::Reset() {

	chart->GetY1Axis()->ClearDataView();
	for (int i = 0; i < nbView; i++) SAFE_DELETE(views[i]);
	nbView = 0;

}

void PressureEvolution::ProcessMessage(GLComponent *src, int message) {
	Geometry *geom = worker->GetGeometry();
	MolFlow *mApp = (MolFlow *)theApp;
	switch (message) {
	case MSG_BUTTON:
		if (src == dismissButton) {
			SetVisible(FALSE);
		}
		else if (src == selButton) {
			int idx = profCombo->GetSelectedIndex();
			geom->UnSelectAll();
			geom->GetFacet(profCombo->GetUserValueAt(idx))->selected = TRUE;
			geom->UpdateSelection();
			mApp->UpdateFacetParams(TRUE);
			mApp->facetList->SetSelectedRow(profCombo->GetUserValueAt(idx));
			mApp->facetList->ScrollToVisible(profCombo->GetUserValueAt(idx), 1, TRUE);
		}
		else if (src == addButton) {
			int idx = profCombo->GetSelectedIndex();
			if (idx >= 0) addView(profCombo->GetUserValueAt(idx));
			refreshViews();
		}
		else if (src == removeButton) {
			int idx = profCombo->GetSelectedIndex();
			if (idx >= 0) remView(profCombo->GetUserValueAt(idx));
			refreshViews();
		}
		else if (src == resetButton) {
			Reset();
		}
		else if (src == formulaBtn) {
			plot();
		}
		break;
	case MSG_COMBO:
		if (src == normCombo) {
			if (normCombo->GetSelectedIndex() == 7) {//pressure selected
				modeCombo->SetValueAt(0, "Average");
				if (modeCombo->GetSelectedIndex() == 0) modeCombo->SetSelectedIndex(0); //re-select to update text
			}
			else { //something else selected
				modeCombo->SetValueAt(0, "Sum");
				if (modeCombo->GetSelectedIndex() == 0) modeCombo->SetSelectedIndex(0); //re-select to update text
			}
		}
		else if (src == modeCombo) {
			selectedSliceText->SetEditable(modeCombo->GetSelectedIndex() == 1);
		}
		refreshViews();
		break;
	case MSG_TEXT: //enter pressed
		if (src == selectedSliceText) {
			int sliceNum;
			if (!selectedSliceText->GetNumberInt(&sliceNum)) {
				GLMessageBox::Display("Invalid slice number", "Error", GLDLG_OK, GLDLG_ICONERROR);
				return;
			}
			if (!(sliceNum > 0 && sliceNum <= 100)) {
				GLMessageBox::Display("Choose a slice between 1 and 100", "Error", GLDLG_OK, GLDLG_ICONERROR);
				return;
			}
			//valid input
			selectedSlice = sliceNum - 1;
			refreshViews();
		}
		break;
	case MSG_TOGGLE:
		if (src == logXToggle) {
			chart->GetXAxis()->SetScale(logXToggle->IsChecked());
		}
		else if (src == logYToggle) {
			chart->GetY1Axis()->SetScale(logYToggle->IsChecked());
		}
		break;
	}

	GLWindow::ProcessMessage(src, message);

}

