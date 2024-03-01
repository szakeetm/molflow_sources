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
#include "TextureScaling.h"
#include "Facet_shared.h"
#include "MolflowGeometry.h"
#include "Worker.h"
#include "GLApp/GLMessageBox.h"
#include "Helper/MathTools.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLTextField.h"
#include "GLApp/GLLabel.h"
#include "GLApp/GLToggle.h"
#include "GLApp/GLTitledPanel.h"
#include "GLApp/GLGradient.h"
#include "GLApp/GLCombo.h"

/**
* \brief Constructor with initialisation for Texture scaling window (Tools/Texture Scaling)
*/
TextureScaling::TextureScaling(Worker* worker_, GeometryViewer** viewers_):GLWindow(),
	worker(worker_), interfGeom(worker_->GetMolflowGeometry()), viewers(viewers_) {

	int wD = 500;
	int hD = 225;

	SetTitle("Texture Scaling");
	SetIconfiable(true);

	GLTitledPanel *panel = new GLTitledPanel("Texture Range");
	panel->SetBounds(5,2,365,98);
	Add(panel);

	GLLabel *l1 = new GLLabel("Min");
	l1->SetBounds(10,20,30,18);
	Add(l1);

	manualScaleMinText = new GLTextField(0,"");
	manualScaleMinText->SetBounds(40,20,85,19);
	manualScaleMinText->SetEditable(true);
	Add(manualScaleMinText);

	GLLabel *l2 = new GLLabel("Max");
	l2->SetBounds(10,45,30,18);
	Add(l2);

	manualScaleMaxText = new GLTextField(0,"");
	manualScaleMaxText->SetBounds(40,45,85,19);
	manualScaleMaxText->SetEditable(true);
	Add(manualScaleMaxText);

	setToCurrentButton = new GLButton(0,"Set to geom.");
	setToCurrentButton->SetBounds(10,70,90,19);
	Add(setToCurrentButton);
	
	applyButton = new GLButton(0,"Apply min/max");
	applyButton->SetBounds(105,70,90,19);
	Add(applyButton);
	
	autoScaleToggle = new GLToggle(0,"Autoscale");
	autoScaleToggle->SetBounds(130,20,80,19);
	Add(autoScaleToggle);

    autoscaleTimedepModeCombo = new GLCombo(0);
    autoscaleTimedepModeCombo->SetSize(3);
    autoscaleTimedepModeCombo->SetValueAt(0,"Only moments");
    autoscaleTimedepModeCombo->SetValueAt(1,"Moments+const. flow");
    autoscaleTimedepModeCombo->SetValueAt(2,"Only constant flow");
    autoscaleTimedepModeCombo->SetBounds(130,45,125,19);
    autoscaleTimedepModeCombo->SetSelectedIndex(1);
    Add(autoscaleTimedepModeCombo);

	useColorToggle = new GLToggle(0,"Use colors");
	useColorToggle->SetBounds(260,20,85,19);
	Add(useColorToggle);

	logarithmicToggle = new GLToggle(0,"Logarithmic scale");
	logarithmicToggle->SetBounds(260,45,80,19);
	Add(logarithmicToggle);

	GLLabel *l3 = new GLLabel("Swap");
	l3->SetBounds(275,70,30,18);
	Add(l3);

	swapText = new GLTextField(0,"");
	swapText->SetEditable(false);
	swapText->SetBounds(305,70,55,18);
	Add(swapText);

	GLTitledPanel *panel2 = new GLTitledPanel("Geometry:");
	panel2->SetBounds(375,2,120,98);
	Add(panel2);

	GLLabel *l4 = new GLLabel("Min:");
	l4->SetBounds(391,31,20,19);
	Add(l4);

	geomMinLabel = new GLLabel("");
	geomMinLabel->SetBounds(420,30,70,19);
	Add(geomMinLabel);

	GLLabel *l5 = new GLLabel("Max:");
	l5->SetBounds(391,66,20,19);
	Add(l5);

	geomMaxLabel = new GLLabel("");
	geomMaxLabel->SetBounds(420,65,70,19);
	Add(geomMaxLabel);

	GLTitledPanel *panel3 = new GLTitledPanel("Gradient");
	panel3->SetBounds(5,102,490,65);
	Add(panel3);

	gradient = new GLGradient(0);
	gradient->SetMouseCursor(true);
	gradient->SetBounds(10,117,470,40);
	Add(gradient);

	GLLabel* displayLabel = new GLLabel("Physics mode:");
	displayLabel->SetBounds(160,179,35,20);
	Add(displayLabel);

	physicsModeCombo = new GLCombo(0);
	physicsModeCombo->SetSize(3);
	physicsModeCombo->SetValueAt(0,"Pressure [mbar]");
	physicsModeCombo->SetValueAt(1,"Impingement rate [1/sec/m\262]");
	physicsModeCombo->SetValueAt(2,"Particle density [1/m\263]");
	physicsModeCombo->SetBounds(200,178,150,20);
	physicsModeCombo->SetSelectedIndex(0);
	Add(physicsModeCombo);

	SetBounds(8,30,wD,hD);

	RestoreDeviceObjects();

}

/**
* \brief Updates text for memory requirement size (Swap)
*/
void TextureScaling::RecalcSwapSize() {

	size_t swap = 0;
	size_t nbFacet = interfGeom->GetNbFacet();
	for(size_t i=0;i<nbFacet;i++) {
		InterfaceFacet *f = interfGeom->GetFacet(i);
		if(f->sh.isTextured) {
			swap += f->GetTexSwapSize(useColorToggle->GetState());
		}
	}
	swapText->SetText(FormatMemory(swap));

}

/**
* \brief Updates all components in the window e.g. text labels and gradient
*/
void TextureScaling::Update() {

	char tmp[128];

	if (!IsVisible() || IsIconic()) return;

	UpdateAutoScaleLimits();

	//Set manual min/max text fields
	sprintf(tmp,"%g",interfGeom->texture_limits[interfGeom->textureMode].manual.min.steady_state);
	manualScaleMinText->SetText(tmp);
	manualScaleMinText->SetEditable(!interfGeom->texAutoScale);
	sprintf(tmp,"%g",interfGeom->texture_limits[interfGeom->textureMode].manual.max.steady_state);
	manualScaleMaxText->SetText(tmp);
	manualScaleMaxText->SetEditable(!interfGeom->texAutoScale);
	applyButton->SetEnabled(!interfGeom->texAutoScale);
	autoScaleToggle->SetState(interfGeom->texAutoScale);
	autoscaleTimedepModeCombo->SetVisible(interfGeom->texAutoScale);
	autoscaleTimedepModeCombo->SetSelectedIndex(static_cast<int>(interfGeom->texAutoScaleMode));
	logarithmicToggle->SetState(interfGeom->texLogScale);
	useColorToggle->SetState(interfGeom->texColormap);
	gradient->SetType(interfGeom->texColormap /*viewers[0]->showColormap*/ ? GRADIENT_COLOR : GRADIENT_BW);
	physicsModeCombo->SetSelectedIndex(interfGeom->textureMode);
	RecalcSwapSize();

}

void TextureScaling::UpdateAutoScaleLimits() {
		//Set autoscale minimum label
	auto [autoscaleMin, autoscaleMax] = interfGeom->GetTextureAutoscaleMinMax();

	//Set current geometry limits
	char tmp[128];
	sprintf(tmp, "%.3E", autoscaleMin);
	geomMinLabel->SetText(tmp);
	sprintf(tmp, "%.3E", autoscaleMax);
	geomMaxLabel->SetText(tmp);
		gradient->SetScale(interfGeom->texLogScale ? LOG_SCALE : LINEAR_SCALE);
	if (!interfGeom->texAutoScale) { // Set manual texture scaling
		//In case of manual scaling, "steady state" variable used always
		gradient->SetMinMax(
			interfGeom->texture_limits[interfGeom->textureMode].manual.min.steady_state,
			interfGeom->texture_limits[interfGeom->textureMode].manual.max.steady_state
		);
	}
	else { //Set auto texture scaling
		gradient->SetMinMax(autoscaleMin, autoscaleMax);
	}
}

/**
* \brief Displays the window
* \param w Worker handle
* \param v handle for the GeometryViewer (TODO: needed?)
*/
void TextureScaling::Display() {

	if(!interfGeom->IsLoaded()) {
		GLMessageBox::Display("No geometry loaded.","No geometry",GLDLG_OK,GLDLG_ICONERROR);
		return;
	}

	SetVisible(true);
	Update();
	char tmp[64];
	sprintf(tmp, "%g", interfGeom->texture_limits[interfGeom->textureMode].manual.min.steady_state);
	manualScaleMinText->SetText(tmp);
	sprintf(tmp, "%g", interfGeom->texture_limits[interfGeom->textureMode].manual.max.steady_state);
	manualScaleMaxText->SetText(tmp);
}

/**
* \brief Function for processing various inputs (button, check boxes etc.)
* \param src Exact source of the call
* \param message Type of the source (button)
*/
void TextureScaling::ProcessMessage(GLComponent *src,int message) {

	switch(message) {
	case MSG_BUTTON:

		if (src==applyButton) {

			double min,max;

			if( !manualScaleMinText->GetNumber(&min) ) {
				GLMessageBox::Display("Invalid minimum value","Error",GLDLG_OK,GLDLG_ICONERROR);
				Update();
				return;
			}
			if( !manualScaleMaxText->GetNumber(&max) ) {
				GLMessageBox::Display("Invalid maximum value","Error",GLDLG_OK,GLDLG_ICONERROR);
				Update();
				return;
			}
			if( min>=max ) {
				GLMessageBox::Display("min must be lower than max","Error",GLDLG_OK,GLDLG_ICONERROR);
				Update();
				return;
			}

			interfGeom->texture_limits[interfGeom->textureMode].manual.min.steady_state = min;
			interfGeom->texture_limits[interfGeom->textureMode].manual.max.steady_state = max;
			
			try {
				worker->Update(0.0f);
			} catch (const std::exception &e) {
				GLMessageBox::Display(e.what(),"Error (Worker::Update)",GLDLG_OK,GLDLG_ICONERROR);
			}
			Update();

		} else if (src==setToCurrentButton) {

			auto [autoscaleMin, autoscaleMax] = interfGeom->GetTextureAutoscaleMinMax();

			//Apply to geometry
			interfGeom->texture_limits[interfGeom->textureMode].manual.min.steady_state = autoscaleMin;
			interfGeom->texture_limits[interfGeom->textureMode].manual.max.steady_state = autoscaleMax;
			interfGeom->texAutoScale=false;
			try {
				worker->Update(0.0f);
			} catch (const std::exception &e) {
				GLMessageBox::Display(e.what(),"Error (Worker::Update)",GLDLG_OK,GLDLG_ICONERROR);
			}
			Update();
		}
		break;

	case MSG_TOGGLE:
		if (src==useColorToggle) {
			//for(int i=0;i<MAX_VIEWER;i++) viewers[i]->showColormap = useColorToggle->GetState();
			interfGeom->texColormap = useColorToggle->GetState();
			worker->Update(0.0f);
			Update();
		} else if (src==autoScaleToggle) {
			interfGeom->texAutoScale = autoScaleToggle->GetState();
			worker->Update(0.0f);
			Update();
		} else if (src==logarithmicToggle) {
			interfGeom->texLogScale = logarithmicToggle->GetState();
			worker->Update(0.0f);
			Update();
		}
		break;

	case MSG_TEXT:
		ProcessMessage(applyButton,MSG_BUTTON);
		break;

	case MSG_COMBO:
		if(src==physicsModeCombo) {
			interfGeom->textureMode=physicsModeCombo->GetSelectedIndex();
			try {
				worker->Update(0.0f);
			} catch (const std::exception &e) {
				GLMessageBox::Display(e.what(),"Error (Worker::Update)",GLDLG_OK,GLDLG_ICONERROR);
			}
			char tmp[256];
			Update();
		}
		else if (src == this->autoscaleTimedepModeCombo) {
			interfGeom->texAutoScaleMode = static_cast<AutoScaleMode>(autoscaleTimedepModeCombo->GetSelectedIndex());
			worker->Update(0.0f);
			Update();
		}
		break;
		
	}

	GLWindow::ProcessMessage(src,message);
}