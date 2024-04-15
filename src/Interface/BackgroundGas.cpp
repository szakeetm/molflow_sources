
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
	
	int wD = 260;
	int hD = 225;
	
	GLTitledPanel* groupBox1 = new GLTitledPanel("Background gas");
	groupBox1->SetBounds(15, 30, wD-25, hD-90);
	Add(groupBox1);

	GLLabel* label1 = new GLLabel("Enables collision with a static background gas");
	label1->SetBounds(12, 8, 261, 22);
	Add(label1);

	enableCheckbox = new GLToggle(0, "Enable collisions");
	groupBox1->SetCompBounds(enableCheckbox, 6, 19, 103, 17);
	groupBox1->Add(enableCheckbox);

	GLLabel* label2 = new GLLabel("Mean free path (cm):");
	groupBox1->SetCompBounds(label2, 10, 50, 107, 13);
	groupBox1->Add(label2);

	mfpTextbox = new GLTextField(0,"");
	SetCompBoundsRelativeTo(label2, mfpTextbox, 110, 0, 100, 20);
	groupBox1->Add(mfpTextbox);


	GLLabel* label3 = new GLLabel("Mass ratio (simu/bg):");
	groupBox1->SetCompBounds(label3, 10, 75, 107, 13);
	groupBox1->Add(label3);

	massTextbox = new GLTextField(0, "");
	SetCompBoundsRelativeTo(label3, massTextbox, 110, 0, 100, 20);
	groupBox1->Add(massTextbox);

	cutoffToggle = new GLToggle(0, "Cutoff below (m/s):");
	SetCompBoundsRelativeTo(label3, cutoffToggle, 0, 25, 100, 20);
	groupBox1->Add(cutoffToggle);

	cutoffSpeedTextbox = new GLTextField(0, "");
	SetCompBoundsRelativeTo(massTextbox, cutoffSpeedTextbox, 0, 25, 100, 20);
	groupBox1->Add(cutoffSpeedTextbox);

	applyButton = new GLButton(0, "Apply");
	applyButton->SetBounds(wD/2-40, hD-50, 80, 21);
	Add(applyButton);

	
	SetTitle("Collisions with background gas");
	SetBounds(17, 30, wD, hD);

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

		if (src == applyButton) { //Apply

			double mfp_cm, massRatio, cutoffSpeed;

			if (!(mfpTextbox->GetNumber(&mfp_cm))) {
				GLMessageBox::Display("Invalid mean free path", "Error", GLDLG_OK, GLDLG_ICONERROR);
				return;
			}
			if (!(massTextbox->GetNumber(&massRatio))) {
				GLMessageBox::Display("Invalid mass ratio", "Error", GLDLG_OK, GLDLG_ICONERROR);
				return;
			}

			bool cutoffEnabled = cutoffToggle->GetState() == 1;
			if (cutoffEnabled) {
				if (!(cutoffSpeedTextbox->GetNumber(&cutoffSpeed))) {
					GLMessageBox::Display("Invalid cutoff speed", "Error", GLDLG_OK, GLDLG_ICONERROR);
					return;
				}
			}
			
			if (mApp->AskToReset()) {
				work->model->sp.scattering.enabled = (bool)enableCheckbox->GetState();
				work->model->sp.scattering.meanFreePath_cm = mfp_cm;
				work->model->sp.scattering.massRatio = massRatio;
				work->model->sp.scattering.enableCutoff = cutoffEnabled;
				if (cutoffEnabled) work->model->sp.scattering.cutoffSpeed = cutoffSpeed;

				work->MarkToReload();
				mApp->changedSinceSave = true;
				return;
			}
		}
		break;
	case MSG_TEXT_UPD:
		/*
		if (src == densityTextbox) {
			double density,sigma;
			if (densityTextbox->GetNumber(&density) && sigmaTextbox->GetNumber(&sigma)) {
				mfpTextbox->SetText(100.0 / (density * sigma)); //m->cm
			}
		}
		else if (src == sigmaTextbox) {
			double density, sigma;
			if (sigmaTextbox->GetNumber(&sigma) && densityTextbox->GetNumber(&density)) {
				mfpTextbox->SetText(100.0 / (density * sigma)); //m->cm
			}
		}
		else if (src == mfpTextbox) {
			double mfp_cm, sigma; //let's keep sigma fixed, random choice
			if (mfpTextbox->GetNumber(&mfp_cm) && sigmaTextbox->GetNumber(&sigma)) {
				densityTextbox->SetText(1.0 / (mfp_cm * 0.01 * sigma));
			}
		}
		else if (src == massRatioTextbox) {
			double massRatio;
			if (massRatioTextbox->GetNumber(&massRatio)) {
				massTextbox->SetText(work->model->sp.gasMass / massRatio);
			}
		}
		else if (src == massTextbox) {
			double mass;
			if (massTextbox->GetNumber(&mass)) {
				massRatioTextbox->SetText(work->model->sp.gasMass / mass);
			}
		}
		*/
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
	massTextbox->SetText(work->model->sp.scattering.massRatio);
	cutoffToggle->SetState((int)work->model->sp.scattering.enableCutoff);
	cutoffSpeedTextbox->SetText(work->model->sp.scattering.cutoffSpeed);
}