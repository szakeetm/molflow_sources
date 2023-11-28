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
TextureScaling::TextureScaling():GLWindow() {

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

	texMinText = new GLTextField(0,"");
	texMinText->SetBounds(40,20,85,19);
	texMinText->SetEditable(true);
	Add(texMinText);

	GLLabel *l2 = new GLLabel("Max");
	l2->SetBounds(10,45,30,18);
	Add(l2);

	texMaxText = new GLTextField(0,"");
	texMaxText->SetBounds(40,45,85,19);
	texMaxText->SetEditable(true);
	Add(texMaxText);

	setCurrentButton = new GLButton(0,"Set to current");
	setCurrentButton->SetBounds(40,70,90,19);
	Add(setCurrentButton);
	
	updateButton = new GLButton(0,"Apply");
	updateButton->SetBounds(135,70,90,19);
	Add(updateButton);
	
	texAutoScale = new GLToggle(0,"Autoscale");
	texAutoScale->SetBounds(130,20,80,19);
	Add(texAutoScale);

	/*includeConstantFlow = new GLToggle(0,"Include constant flow");
	includeConstantFlow->SetBounds(130,45,80,19);
	Add(includeConstantFlow);*/

    includeConstantFlow = new GLCombo(0);
    includeConstantFlow->SetSize(3);
    includeConstantFlow->SetValueAt(0,"Only moments");
    includeConstantFlow->SetValueAt(1,"Include constant flow");
    includeConstantFlow->SetValueAt(2,"Only constant flow");
    includeConstantFlow->SetBounds(130,45,125,19);
    includeConstantFlow->SetSelectedIndex(1);
    Add(includeConstantFlow);

	colormapBtn = new GLToggle(0,"Use colors");
	colormapBtn->SetBounds(260,20,85,19);
	Add(colormapBtn);

	logBtn = new GLToggle(0,"Logarithmic scale");
	logBtn->SetBounds(260,45,80,19);
	Add(logBtn);

	GLLabel *l3 = new GLLabel("Swap");
	l3->SetBounds(275,70,30,18);
	Add(l3);

	swapText = new GLTextField(0,"");
	swapText->SetEditable(false);
	swapText->SetBounds(305,70,55,18);
	Add(swapText);

	GLTitledPanel *panel2 = new GLTitledPanel("Current");
	panel2->SetBounds(375,2,120,98);
	Add(panel2);

	GLLabel *l4 = new GLLabel("Min:");
	l4->SetBounds(390,30,20,19);
	Add(l4);

	texCMinText = new GLLabel("");
	texCMinText->SetBounds(420,30,70,19);
	Add(texCMinText);

	GLLabel *l5 = new GLLabel("Max:");
	l5->SetBounds(390,65,20,19);
	Add(l5);

	texCMaxText = new GLLabel("");
	texCMaxText->SetBounds(420,65,70,19);
	Add(texCMaxText);

	

	GLTitledPanel *panel3 = new GLTitledPanel("Gradient");
	panel3->SetBounds(5,102,490,65);
	Add(panel3);

	gradient = new GLGradient(0);
	gradient->SetMouseCursor(true);
	gradient->SetBounds(10,117,470,40);
	Add(gradient);

	GLLabel* displayLabel = new GLLabel("Show:");
	displayLabel->SetBounds(160,179,35,20);
	Add(displayLabel);

	modeCombo = new GLCombo(0);
	modeCombo->SetSize(3);
	modeCombo->SetValueAt(0,"Pressure [mbar]");
	modeCombo->SetValueAt(1,"Impingement rate [1/sec/m\262]");
	modeCombo->SetValueAt(2,"Particle density [1/m3]");
	modeCombo->SetBounds(200,178,150,20);
	modeCombo->SetSelectedIndex(0);
	Add(modeCombo);

	SetBounds(8,30,wD,hD);

	RestoreDeviceObjects();

	interfGeom = NULL;

}

/**
* \brief Updates text for memory requirement size (Swap)
*/
void TextureScaling::UpdateSize() {

	size_t swap = 0;
	size_t nbFacet = interfGeom->GetNbFacet();
	for(size_t i=0;i<nbFacet;i++) {
		InterfaceFacet *f = interfGeom->GetFacet(i);
		if(f->sh.isTextured) {
			swap += f->GetTexSwapSize(colormapBtn->GetState());
		}
	}
	swapText->SetText(FormatMemory(swap));

}

/**
* \brief Updates all components in the window e.g. text labels and gradient
*/
void TextureScaling::Update() {

	if(!IsVisible() || IsIconic()) return;  

	char tmp[128];
	/*
	//Manual minimum label
	sprintf(tmp, "%.3E", interfGeom->texture_limits[interfGeom->textureMode].manual.min.all);
	texMinText->SetText(tmp);
	//Manual maximum label
	sprintf(tmp, "%.3E", interfGeom->texture_limits[interfGeom->textureMode].manual.max.all);
	texMaxText->SetText(tmp);
	*/

	//Set autoscale minimum label
	double minScale = 0.0;
    double maxScale = 0.0;
	if (interfGeom->texAutoScaleIncludeConstantFlow == 0) { //moments only
		minScale =
			interfGeom->texture_limits[interfGeom->textureMode].autoscale.min.moments_only;
		maxScale =
			interfGeom->texture_limits[interfGeom->textureMode].autoscale.max.moments_only;
	}
	else if (interfGeom->texAutoScaleIncludeConstantFlow == 1) { // moments and const flow
		minScale =
			std::min(interfGeom->texture_limits[interfGeom->textureMode].autoscale.min.steady_state, interfGeom->texture_limits[interfGeom->textureMode].autoscale.min.moments_only);
		maxScale =
			std::max(interfGeom->texture_limits[interfGeom->textureMode].autoscale.max.steady_state, interfGeom->texture_limits[interfGeom->textureMode].autoscale.max.moments_only);
	}
	else { // const flow only
		minScale =
			interfGeom->texture_limits[interfGeom->textureMode].autoscale.min.steady_state;
		maxScale =
			interfGeom->texture_limits[interfGeom->textureMode].autoscale.max.steady_state;
	}

	sprintf(tmp,"%.3E",minScale);
	texCMinText->SetText(tmp);
	//Set autoscale maximum label
	sprintf(tmp,"%.3E",maxScale);
	texCMaxText->SetText(tmp);

	texAutoScale->SetState(interfGeom->texAutoScale);
	includeConstantFlow->SetVisible(interfGeom->texAutoScale);
	includeConstantFlow->SetSelectedIndex(interfGeom->texAutoScaleIncludeConstantFlow);
	logBtn->SetState(interfGeom->texLogScale);
	gradient->SetScale(interfGeom->texLogScale?LOG_SCALE:LINEAR_SCALE);
	if( !interfGeom->texAutoScale ) { // Set manual texture scaling
		gradient->SetMinMax(
			interfGeom->texture_limits[interfGeom->textureMode].manual.min.steady_state,
			interfGeom->texture_limits[interfGeom->textureMode].manual.max.steady_state
			);
	} else { //Set auto texture scaling
	    // Has been previously calculated
        /*double minScale = 0.0;
        double maxScale = 0.0;
        if(interfGeom->texAutoScaleIncludeConstantFlow == 0){
            minScale=
                    interfGeom->texture_limits[interfGeom->textureMode].autoscale.min.moments_only;
            maxScale=
                    interfGeom->texture_limits[interfGeom->textureMode].autoscale.max.moments_only;
        } else if(interfGeom->texAutoScaleIncludeConstantFlow == 1){
            minScale=
                    std::min(interfGeom->texture_limits[interfGeom->textureMode].autoscale.min.all, interfGeom->texture_limits[interfGeom->textureMode].autoscale.min.moments_only);
            maxScale=
                    std::max(interfGeom->texture_limits[interfGeom->textureMode].autoscale.max.all, interfGeom->texture_limits[interfGeom->textureMode].autoscale.max.moments_only);
        } else { // == 2
            minScale=
                    interfGeom->texture_limits[interfGeom->textureMode].autoscale.min.all;
            maxScale=
                    interfGeom->texture_limits[interfGeom->textureMode].autoscale.max.all;
        }*/
		gradient->SetMinMax(minScale, maxScale);
	}
	//colormapBtn->SetThreadStates(viewers[0]->showColormap);
	colormapBtn->SetState(interfGeom->texColormap);
	gradient->SetType(interfGeom->texColormap /*viewers[0]->showColormap*/?GRADIENT_COLOR:GRADIENT_BW );
	modeCombo->SetSelectedIndex(interfGeom->textureMode);
	UpdateSize();

}

/**
* \brief Displays the window
* \param w Worker handle
* \param v handle for the GeometryViewer (TODO: needed?)
*/
void TextureScaling::Display(Worker *w,GeometryViewer **v) {

	worker = w;
	interfGeom = w->GetMolflowGeometry();
	viewers = v;
	if(!interfGeom->IsLoaded()) {
		GLMessageBox::Display("No geometry loaded.","No geometry",GLDLG_OK,GLDLG_ICONERROR);
		return;
	}

	SetVisible(true);
	Update();
	char tmp[64];
	sprintf(tmp, "%g", interfGeom->texture_limits[interfGeom->textureMode].manual.min.steady_state);
	texMinText->SetText(tmp);
	sprintf(tmp, "%g", interfGeom->texture_limits[interfGeom->textureMode].manual.max.steady_state);
	texMaxText->SetText(tmp);
}

/**
* \brief Function for processing various inputs (button, check boxes etc.)
* \param src Exact source of the call
* \param message Type of the source (button)
*/
void TextureScaling::ProcessMessage(GLComponent *src,int message) {

	switch(message) {
	case MSG_BUTTON:

		if (src==updateButton) {

			double min,max;

			if( !texMinText->GetNumber(&min) ) {
				GLMessageBox::Display("Invalid minimum value","Error",GLDLG_OK,GLDLG_ICONERROR);
				Update();
				return;
			}
			if( !texMaxText->GetNumber(&max) ) {
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
			interfGeom->texAutoScale = texAutoScale->GetState();
			interfGeom->texAutoScaleIncludeConstantFlow = includeConstantFlow->GetSelectedIndex();
			try {
				worker->Update(0.0f);
			} catch (const std::exception &e) {
				GLMessageBox::Display(e.what(),"Error (Worker::Update)",GLDLG_OK,GLDLG_ICONERROR);
			}
			Update();

		} else if (src==setCurrentButton) {
			if (interfGeom->texAutoScaleIncludeConstantFlow == 0) { //moments only
				interfGeom->texture_limits[interfGeom->textureMode].manual.min.steady_state =
					interfGeom->texture_limits[interfGeom->textureMode].autoscale.min.moments_only;
				interfGeom->texture_limits[interfGeom->textureMode].manual.max.steady_state =
					interfGeom->texture_limits[interfGeom->textureMode].autoscale.max.moments_only;
			}
			else if (interfGeom->texAutoScaleIncludeConstantFlow == 1) { // moments and const. flow
				interfGeom->texture_limits[interfGeom->textureMode].manual.min.steady_state =
					std::min(interfGeom->texture_limits[interfGeom->textureMode].autoscale.min.steady_state, interfGeom->texture_limits[interfGeom->textureMode].autoscale.min.moments_only);
				interfGeom->texture_limits[interfGeom->textureMode].manual.max.steady_state =
					std::max(interfGeom->texture_limits[interfGeom->textureMode].autoscale.max.steady_state, interfGeom->texture_limits[interfGeom->textureMode].autoscale.max.moments_only);
			}
			else { // const. flow only
				interfGeom->texture_limits[interfGeom->textureMode].manual.min.steady_state =
					interfGeom->texture_limits[interfGeom->textureMode].autoscale.min.steady_state;
				interfGeom->texture_limits[interfGeom->textureMode].manual.max.steady_state =
					interfGeom->texture_limits[interfGeom->textureMode].autoscale.max.steady_state;
			}

			texMinText->SetText(texCMinText->GetText());
			texMaxText->SetText(texCMaxText->GetText());
			texAutoScale->SetState(false);
			includeConstantFlow->SetVisible(false);
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
		if (src==colormapBtn) {
			//for(int i=0;i<MAX_VIEWER;i++) viewers[i]->showColormap = colormapBtn->GetState();
			interfGeom->texColormap = colormapBtn->GetState();
			worker->Update(0.0f);
			Update();
		} else if (src==texAutoScale) {
			interfGeom->texAutoScale = texAutoScale->GetState();
			worker->Update(0.0f);
			Update();
		} else if (src==this->includeConstantFlow) {
			interfGeom->texAutoScaleIncludeConstantFlow = includeConstantFlow->GetSelectedIndex();
			worker->Update(0.0f);
			Update();
		} else if (src==logBtn) {
			interfGeom->texLogScale = logBtn->GetState();
			gradient->SetScale(interfGeom->texLogScale?LOG_SCALE:LINEAR_SCALE);
			worker->Update(0.0f);
			Update();
		}
		break;

	case MSG_TEXT:
		ProcessMessage(updateButton,MSG_BUTTON);
		break;

	case MSG_COMBO:
		if(src==modeCombo) {
			interfGeom->textureMode=modeCombo->GetSelectedIndex();
			try {
				worker->Update(0.0f);
			} catch (const std::exception &e) {
				GLMessageBox::Display(e.what(),"Error (Worker::Update)",GLDLG_OK,GLDLG_ICONERROR);
			}
			char tmp[256];
			sprintf(tmp,"%g",interfGeom->texture_limits[interfGeom->textureMode].manual.min.steady_state);
			texMinText->SetText(tmp);
			sprintf(tmp,"%g",interfGeom->texture_limits[interfGeom->textureMode].manual.max.steady_state);
			texMaxText->SetText(tmp);
			Update();
		}
		break;
		
	}

	GLWindow::ProcessMessage(src,message);
}
