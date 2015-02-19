/*
File:        TextureSettings.cpp
Description: Texture settings dialog (min,max,autoscale,gradient)
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

#include "TextureSettings.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLMessageBox.h"

TextureSettings::TextureSettings():GLWindow() {

	int wD = 500;
	int hD = 225;

	SetTitle("Texture Scaling");
	SetIconfiable(TRUE);

	GLTitledPanel *panel = new GLTitledPanel("Texture Range");
	panel->SetBounds(5,2,365,98);
	Add(panel);

	GLLabel *l1 = new GLLabel("Min");
	l1->SetBounds(10,20,30,18);
	Add(l1);

	texMinText = new GLTextField(0,"");
	texMinText->SetBounds(40,20,85,19);
	texMinText->SetEditable(TRUE);
	Add(texMinText);

	GLLabel *l2 = new GLLabel("Max");
	l2->SetBounds(10,45,30,18);
	Add(l2);

	texMaxText = new GLTextField(0,"");
	texMaxText->SetBounds(40,45,85,19);
	texMaxText->SetEditable(TRUE);
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

	includeConstantFlow = new GLToggle(0,"Include constant flow");
	includeConstantFlow->SetBounds(130,45,80,19);
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
	swapText->SetEditable(FALSE);
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

	// ---------------------------------------------------

	GLTitledPanel *panel3 = new GLTitledPanel("Gradient");
	panel3->SetBounds(5,102,490,65);
	Add(panel3);

	gradient = new GLGradient(0);
	gradient->SetMouseCursor(TRUE);
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

	geom = NULL;

}

// --------------------------------------------------------------------

void TextureSettings::UpdateSize() {

	DWORD swap = 0;
	int nbFacet = geom->GetNbFacet();
	for(int i=0;i<nbFacet;i++) {
		Facet *f = geom->GetFacet(i);
		if(f->sh.isTextured) {
			swap += f->GetTexSwapSize(colormapBtn->GetState());
		}
	}
	swapText->SetText(FormatMemory(swap));

}

// --------------------------------------------------------------------

void TextureSettings::Update() {

	if(!IsVisible() || IsIconic()) return;  

	char tmp[128];
	//Set autoscale minimum label
	sprintf(tmp,"%.3E",(geom->texAutoScaleIncludeConstantFlow?
		geom->texture_limits[geom->textureMode].autoscale.min.all
		:geom->texture_limits[geom->textureMode].autoscale.min.moments_only));
	texCMinText->SetText(tmp);
	//Set autoscale maximum label
	sprintf(tmp,"%.3E",(geom->texAutoScaleIncludeConstantFlow?
		geom->texture_limits[geom->textureMode].autoscale.max.all
		:geom->texture_limits[geom->textureMode].autoscale.max.moments_only));
	texCMaxText->SetText(tmp);
	texAutoScale->SetState(geom->texAutoScale);
	includeConstantFlow->SetVisible(geom->texAutoScale);
	includeConstantFlow->SetState(geom->texAutoScaleIncludeConstantFlow);
	logBtn->SetState(geom->texLogScale);
	gradient->SetScale(geom->texLogScale?LOG_SCALE:LINEAR_SCALE);
	if( !geom->texAutoScale ) { // Set manual texture scaling
		gradient->SetMinMax(
			geom->texture_limits[geom->textureMode].manual.min.all,
			geom->texture_limits[geom->textureMode].manual.max.all
			);
	} else { //Set auto texture scaling
		gradient->SetMinMax(
			geom->texAutoScaleIncludeConstantFlow?
			geom->texture_limits[geom->textureMode].autoscale.min.all
			:geom->texture_limits[geom->textureMode].autoscale.min.moments_only,
			geom->texAutoScaleIncludeConstantFlow?
			geom->texture_limits[geom->textureMode].autoscale.max.all
			:geom->texture_limits[geom->textureMode].autoscale.max.moments_only
			);
	}
	colormapBtn->SetState(viewers[0]->showColormap);
	gradient->SetType( viewers[0]->showColormap?GRADIENT_COLOR:GRADIENT_BW );
	UpdateSize();

}

void TextureSettings::Display(Worker *w,GeometryViewer **v) {

	worker = w;
	geom = w->GetGeometry();
	viewers = v;
	if(!geom->IsLoaded()) {
		GLMessageBox::Display("No geometry loaded.","No geometry",GLDLG_OK,GLDLG_ICONERROR);
		return;
	}
	char tmp[128];
	sprintf(tmp,"%.3E",geom->texture_limits[geom->textureMode].manual.min.all);
	texMinText->SetText(tmp);
	sprintf(tmp,"%.3E",geom->texture_limits[geom->textureMode].manual.max.all);
	texMaxText->SetText(tmp);

	SetVisible(TRUE);
	Update();

}

// --------------------------------------------------------------------

void TextureSettings::ProcessMessage(GLComponent *src,int message) {

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
			geom->texture_limits[geom->textureMode].manual.min.all = min;
			geom->texture_limits[geom->textureMode].manual.max.all = max;
			geom->texAutoScale = texAutoScale->GetState();
			geom->texAutoScaleIncludeConstantFlow = includeConstantFlow->GetState();
			try {
				worker->Update(0.0f);
			} catch(Error &e) {
				GLMessageBox::Display((char *)e.GetMsg(),"Error (Worker::Update)",GLDLG_OK,GLDLG_ICONERROR);
			}
			Update();

		} else if (src==setCurrentButton) {
			geom->texture_limits[geom->textureMode].manual.min.all=
				geom->texAutoScaleIncludeConstantFlow?
				geom->texture_limits[geom->textureMode].autoscale.min.all
				:geom->texture_limits[geom->textureMode].autoscale.min.moments_only;
			geom->texture_limits[geom->textureMode].manual.max.all=
				geom->texAutoScaleIncludeConstantFlow?
				geom->texture_limits[geom->textureMode].autoscale.max.all
				:geom->texture_limits[geom->textureMode].autoscale.max.moments_only;texMinText->SetText(texCMinText->GetText());
			texMinText->SetText(texCMinText->GetText());
			texMaxText->SetText(texCMaxText->GetText());
			texAutoScale->SetState(FALSE);
			includeConstantFlow->SetVisible(FALSE);
			geom->texAutoScale=false;
			try {
				worker->Update(0.0f);
			} catch(Error &e) {
				GLMessageBox::Display((char *)e.GetMsg(),"Error (Worker::Update)",GLDLG_OK,GLDLG_ICONERROR);
			}
			Update();
		}
		break;

	case MSG_TOGGLE:
		if (src==colormapBtn) {
			for(int i=0;i<MAX_VIEWER;i++) viewers[i]->showColormap = colormapBtn->GetState();
			geom->texColormap = colormapBtn->GetState();
			worker->Update(0.0f);
			Update();
		} else if (src==texAutoScale) {
			geom->texAutoScale = texAutoScale->GetState();
			worker->Update(0.0f);
			Update();
		} else if (src==this->includeConstantFlow) {
			geom->texAutoScaleIncludeConstantFlow = includeConstantFlow->GetState();
			worker->Update(0.0f);
			Update();
		} else if (src==logBtn) {
			geom->texLogScale = logBtn->GetState();
			gradient->SetScale(geom->texLogScale?LOG_SCALE:LINEAR_SCALE);
			worker->Update(0.0f);
			Update();
		}
		break;

	case MSG_TEXT:
		ProcessMessage(updateButton,MSG_BUTTON);
		break;

	case MSG_COMBO:
		if(src==modeCombo) {
			geom->textureMode=modeCombo->GetSelectedIndex();
			try {
				worker->Update(0.0f);
			} catch(Error &e) {
				GLMessageBox::Display((char *)e.GetMsg(),"Error (Worker::Update)",GLDLG_OK,GLDLG_ICONERROR);
			}
			char tmp[256];
			sprintf(tmp,"%g",geom->texture_limits[geom->textureMode].manual.min.all);
			texMinText->SetText(tmp);
			sprintf(tmp,"%g",geom->texture_limits[geom->textureMode].manual.max.all);
			texMaxText->SetText(tmp);
			Update();
		}
		break;
		
	}

	GLWindow::ProcessMessage(src,message);
}
