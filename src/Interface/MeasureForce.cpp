/*
Program:     MolFlow+ / Synrad+
Description: Monte Carlo simulator for ultra-high vacuum and synchrotron radiation
Authors:     Jean-Luc PONS / Roberto KERSEVAN / Marton ADY / Pascal BAEHR
Copyright:   E.S.R.F / CERN
Website:     https://cern.ch/molflow

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

Full license text: https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html
*/

#include "MeasureForce.h"
#include "GLApp/GLTitledPanel.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLMessageBox.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLTextField.h"
#include "GLApp/GLLabel.h"
#include "GLApp/GLToggle.h"
#include "GLApp/GLTitledPanel.h"
#include "Geometry_shared.h"
#include "MolFlow.h"
#include "Facet_shared.h"
#include <sstream>

extern MolFlow* mApp;

/**
* \brief Constructor with basic initialisation for the measure forces window (Tools/Measure forces)
* \param g pointer to InterfaceGeometry
* \param w pointer to Worker handler
*/
MeasureForce::MeasureForce(InterfaceGeometry *g,Worker *w):GLWindow() {

	int wD = 380;
	int hD = 160;

	enableMeasureCheckbox = new GLToggle(0, "Enable force measurement (has performance impact)");
	enableMeasureCheckbox->SetBounds(5, 5, 350, 20);
	Add(enableMeasureCheckbox);
	
	GLTitledPanel* groupBox1 = new GLTitledPanel("Torque relative to...");
	groupBox1->SetBounds(5, 25, wD-10, 75);
	Add(groupBox1);

	auto label3 = new GLLabel("mx0");
	groupBox1->SetCompBounds(label3, 5, 18, 20, 20);
	groupBox1->Add(label3);

	x0Text = new GLTextField(0, "0");
	groupBox1->SetCompBounds(x0Text, 30, 17, 80, 20);
	groupBox1->Add(x0Text);

	auto label4 = new GLLabel("my0");
	groupBox1->SetCompBounds(label4, 130, 18, 20, 20);
	groupBox1->Add(label4);

	y0Text = new GLTextField(0, "0");
	groupBox1->SetCompBounds(y0Text, 155, 17, 80, 20);
	groupBox1->Add(y0Text);

	auto label5 = new GLLabel("mz0");
	groupBox1->SetCompBounds(label5, 255, 18, 20, 20);
	groupBox1->Add(label5);

	z0Text = new GLTextField(0, "0");
	groupBox1->SetCompBounds(z0Text, 280, 17, 80, 20);
	groupBox1->Add(z0Text);

	useVertexButton = new GLButton(0, "Selected vertex");
	groupBox1->SetCompBounds(useVertexButton, 30, 45, 140, 20);
	groupBox1->Add(useVertexButton);

	centerOfFacetButton = new GLButton(0, "Center of selected facet");
	groupBox1->SetCompBounds(centerOfFacetButton, 175, 45, 140, 20);
	groupBox1->Add(centerOfFacetButton);

	applyButton = new GLButton(0, "Apply");
	applyButton->SetBounds(100, 110, 75, 21);
	Add(applyButton);

	dismissButton = new GLButton(0, "Dismiss");
	dismissButton->SetBounds(180, 110, 75, 21);
	Add(dismissButton);

	SetTitle("Measure forces");
	// Center dialog
	int wS, hS;
	GLToolkit::GetScreenSize(&wS, &hS);
	int xD = (wS - wD) / 2;
	int yD = (hS - hD) / 2;
	SetBounds(xD, yD, wD, hD);

	guiGeom = g;
	work = w;
}

/**
* \brief Function for processing various inputs (button, check boxes etc.)
* \param src Exact source of the call
* \param message Type of the source (button)
*/
void MeasureForce::ProcessMessage(GLComponent* src, int message) {

	switch (message) {


	case MSG_BUTTON:

		if (src == dismissButton) { //cancel

			GLWindow::ProcessMessage(NULL, MSG_CLOSE);

		}
		else if (src == applyButton) { //Apply
			double x0, y0, z0;
			bool enabled = enableMeasureCheckbox->GetState();

			if (enabled) {
				if (!(x0Text->GetNumber(&x0))) {
					GLMessageBox::Display("Invalid mx0 coordinate", "Error", GLDLG_OK, GLDLG_ICONERROR);
					return;
				}
				if (!(y0Text->GetNumber(&y0))) {
					GLMessageBox::Display("Invalid my0 coordinate", "Error", GLDLG_OK, GLDLG_ICONERROR);
					return;
				}
				if (!(z0Text->GetNumber(&z0))) {
					GLMessageBox::Display("Invalid mz0 coordinate", "Error", GLDLG_OK, GLDLG_ICONERROR);
					return;
				}
			}

			if (mApp->AskToReset()) {

				work->model->wp.enableForceMeasurement = enabled;
				if (enabled) {
					work->model->wp.torqueRefPoint = Vector3d(x0, y0, z0);
				}
				work->MarkToReload();
				mApp->changedSinceSave = true;
				return;
			}
		}
		else if (src == useVertexButton) { //Use selected vertex as base
			size_t nbs = guiGeom->GetNbSelectedVertex();
			if (nbs != 1) {
				std::ostringstream strstr;
				auto msg = fmt::format("Exactly one vertex needs to be selected.\n(You have selected {}.)", nbs);
				GLMessageBox::Display(msg.c_str(), "Can't use vertex as base", GLDLG_OK, GLDLG_ICONWARNING);
				return;
			}
			else {
				for (size_t i = 0; i < guiGeom->GetNbVertex(); i++) {
					auto v = guiGeom->GetVertex(i);
					if (v->selected) {
						x0Text->SetText(v->x);
						y0Text->SetText(v->y);
						z0Text->SetText(v->z);
						break;
					}
				}
			}
		}
		else if (src == centerOfFacetButton) { //Use center of selected facet as base
			size_t nbs = guiGeom->GetNbSelectedFacets();
			if (nbs != 1) {
				auto msg = fmt::format("Exactly one vertex needs to be selected.\n(You have selected {}.)", nbs);
				GLMessageBox::Display(msg.c_str(), "Can't use facet's center", GLDLG_OK, GLDLG_ICONWARNING);
				return;
			}
			else {
				for (size_t i = 0; i < guiGeom->GetNbFacet(); i++) {
					auto f = guiGeom->GetFacet(i);
					if (f->selected) {
						x0Text->SetText(f->sh.center.x);
						y0Text->SetText(f->sh.center.y);
						z0Text->SetText(f->sh.center.z);
						break;
					}
				}
			}
		}
		break;
	}

	GLWindow::ProcessMessage(src, message);
}

/**
* \brief Updating GUI when checkbox (for force measurement) is toggled
* \param src Exact source of the call (checkbox)
*/
void MeasureForce::UpdateToggle(GLComponent *src) {
	
	if (src == enableMeasureCheckbox) {
		bool enabled = enableMeasureCheckbox->GetState();
		x0Text->SetEditable(enabled);
		y0Text->SetEditable(enabled);
		z0Text->SetEditable(enabled);
	}
}

/**
* \brief Updates values in the GUI from worker for the Meassure Force window
*/
void MeasureForce::Update() {
	
	enableMeasureCheckbox->SetState(work->model->wp.enableForceMeasurement);
	x0Text->SetText(work->model->wp.torqueRefPoint.x);
	y0Text->SetText(work->model->wp.torqueRefPoint.y);
	z0Text->SetText(work->model->wp.torqueRefPoint.z);
}