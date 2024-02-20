/*
Program:     MolFlow+ / Synrad+
Description: Monte Carlo simulator for ultra-high vacuum and synchrotron radiation
Authors:     Jean-Luc PONS / Roberto KERSEVAN / Marton ADY
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


#include "BackgroundGas.h"
#include "GLApp/GLTitledPanel.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLMessageBox.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLTextField.h"
#include "GLApp/GLLabel.h"
#include "GLApp/GLToggle.h"
#include "Geometry_shared.h"
#include "MolFlow.h"
#include <sstream>

extern MolFlow *mApp;

/**
* \brief Constructor with basic initialisation for the movement window (Tools/moving parts)
* \param g pointer to InterfaceGeometry
* \param w pointer to Worker handler
*/
BackgroundGas::BackgroundGas(InterfaceGeometry *g,Worker *w):GLWindow() {
	
	interfGeom = g;
	work = w;
	
	int wD = 300;
	int hD = 215;
	
	GLTitledPanel* groupBox1 = new GLTitledPanel("Background gas");
	groupBox1->SetBounds(15, 30, wD-20, hD-85);
	Add(groupBox1);

	GLLabel* label1 = new GLLabel("Enables collision with a static background gas");
	label1->SetBounds(12, 8, 261, 22);
	Add(label1);

	enableCheckbox = new GLToggle(0, "Enable collisions");
	groupBox1->SetCompBounds(enableCheckbox, 6, 19, 103, 17);
	groupBox1->Add(enableCheckbox);

	GLLabel* label2 = new GLLabel("Mean free path (cm):");
	groupBox1->SetCompBounds(label2, 22, 50, 107, 13);
	groupBox1->Add(label2);

	mfpTextbox = new GLTextField(0,"");
	SetCompBoundsRelativeTo(label2, mfpTextbox, 110, 0, 100, 20);
	groupBox1->Add(mfpTextbox);


	GLLabel* label3 = new GLLabel("Mass ratio:");
	groupBox1->SetCompBounds(label3, 22, 75, 107, 13);
	groupBox1->Add(label3);

	massTextbox = new GLTextField(0, "");
	SetCompBoundsRelativeTo(label3, massTextbox, 110, 0, 100, 20);
	groupBox1->Add(massTextbox);

	GLLabel* label4 = new GLLabel("Rho:");
	groupBox1->SetCompBounds(label4, 22, 100, 107, 13);
	groupBox1->Add(label4);

	rhoTextbox = new GLTextField(0, "");
	SetCompBoundsRelativeTo(label4, rhoTextbox, 110, 0, 100, 20);
	groupBox1->Add(rhoTextbox);

	updateButton = new GLButton(0, "Update");
	updateButton->SetBounds(40, 170, 75, 21);
	Add(updateButton);

	applyButton = new GLButton(0, "Apply");
	applyButton->SetBounds(130, 170, 75, 21);
	Add(applyButton);

	
	SetTitle("Collisions with background gas");
	// Center dialog
	int wS, hS;
	GLToolkit::GetScreenSize(&wS, &hS);
	int xD = (wS - wD) / 2;
	int yD = (hS - hD) / 2;
	SetBounds(xD, yD, wD, hD);

	Update();
}

/**
* \brief Function for processing various inputs (button, check boxes etc.)
* \param src Exact source of the call
* \param message Type of the source (button)
*/
void BackgroundGas::ProcessMessage(GLComponent *src,int message) {

	switch(message) {

	case MSG_BUTTON:

		if(src==updateButton) {
			Update();

		}
		else if (src == applyButton) { //Apply
			double mfp_cm, massRatio, rho;

			if (!(mfpTextbox->GetNumber(&mfp_cm))) {
				GLMessageBox::Display("Invalid mean free path", "Error", GLDLG_OK, GLDLG_ICONERROR);
				return;
			}
			if (!(massTextbox->GetNumber(&massRatio))) {
				GLMessageBox::Display("Invalid mass ratio", "Error", GLDLG_OK, GLDLG_ICONERROR);
				return;
			}
			if (!(rhoTextbox->GetNumber(&rho))) {
				GLMessageBox::Display("Invalid rho", "Error", GLDLG_OK, GLDLG_ICONERROR);
				return;
			}
			
			work->model->sp.scattering.enabled = (bool)enableCheckbox->GetState();
			work->model->sp.scattering.meanFreePath_cm = mfp_cm;
			work->model->sp.scattering.massRatio = massRatio;
			work->model->sp.scattering.rho = rho;
		}
		break;
	case MSG_TEXT_UPD:
		
		break;
	}
	GLWindow::ProcessMessage(src,message);
}


/**
* \brief Updates values in the GUI for the Movements window
*/
void BackgroundGas::Update() {
	
	enableCheckbox->SetState(work->model->sp.scattering.enabled);
	mfpTextbox->SetText(work->model->sp.scattering.meanFreePath_cm);
	rhoTextbox->SetText(work->model->sp.scattering.rho);
	massTextbox->SetText(work->model->sp.scattering.massRatio);

}