/*
File:        ProfilePlotter.cpp
Description: Profile plotter window
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

#include "ProfilePlotter.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLMessageBox.h"
#include "Utils.h"
#include <math.h>
#include "Molflow.h"

extern MolFlow *mApp;

static const char* profType[] = {
	"None",
	"Pressure \201 [mbar]",
	"Pressure \202 [mbar]",
	"Incident angle [deg]",
	"Speed [m/s]",
	"Ort. velocity [m/s]"};

ProfilePlotter::ProfilePlotter():GLWindow() {

	int wD = 650;
	int hD = 400;

	SetTitle("Profile plotter");
	SetIconfiable(TRUE);
	nbView = 0;
	worker = NULL;
	lastUpdate = 0.0f;

	nbColors=8;
	colors[0]=new GLCColor();colors[0]->r=255;colors[0]->g=000;colors[0]->b=055; //red
	colors[1]=new GLCColor();colors[1]->r=000;colors[1]->g=000;colors[1]->b=255; //blue
	colors[2]=new GLCColor();colors[2]->r=000;colors[2]->g=204;colors[2]->b=051; //green
	colors[3]=new GLCColor();colors[3]->r=000;colors[3]->g=000;colors[3]->b=000; //black
	colors[4]=new GLCColor();colors[4]->r=255;colors[4]->g=153;colors[4]->b=051; //orange
	colors[5]=new GLCColor();colors[5]->r=153;colors[5]->g=204;colors[5]->b=255; //light blue
	colors[6]=new GLCColor();colors[6]->r=153;colors[6]->g=000;colors[6]->b=102; //violet
	colors[7]=new GLCColor();colors[7]->r=255;colors[7]->g=230;colors[7]->b=005; //yellow

	chart = new GLChart(0);
	chart->SetBorder(BORDER_BEVEL_IN);
	chart->GetY1Axis()->SetGridVisible(TRUE);
	chart->GetXAxis()->SetGridVisible(TRUE);
	chart->GetY1Axis()->SetAutoScale(TRUE);
	chart->GetY2Axis()->SetAutoScale(TRUE);
	chart->GetY1Axis()->SetAnnotation(VALUE_ANNO);
	chart->GetXAxis()->SetAnnotation(VALUE_ANNO);
	Add(chart);

	dismissButton = new GLButton(0,"Dismiss");
	Add(dismissButton);

	selButton = new GLButton(0,"Show Facet");
	Add(selButton);

	addButton = new GLButton(0,"Add curve");
	Add(addButton);

	removeButton = new GLButton(0,"Remove curve");
	Add(removeButton);

	resetButton = new GLButton(0,"Remove all");
	Add(resetButton);

	profCombo = new GLCombo(0);
	profCombo->SetEditable(TRUE);
	Add(profCombo);

	normLabel = new GLLabel("Normalize:");
	Add(normLabel);
	//qLabel = new GLLabel("Q=");
	//Add(qLabel);
	//unitLabel = new GLLabel("units*l/s");
	//Add(unitLabel);

	normCombo = new GLCombo(0);
	normCombo->SetEditable(TRUE);
	normCombo->SetSize(5);
	normCombo->SetValueAt(0,"None (raw data)");
	normCombo->SetValueAt(1,"Pressure (mbar)");
	normCombo->SetValueAt(2,"Speed (m/s)");
	normCombo->SetValueAt(3,"Angle (deg)");
	normCombo->SetValueAt(4,"Normalize to 1");

	normCombo->SetSelectedIndex(1);
	Add(normCombo);

	logYToggle = new GLToggle(0,"Log Y");
	Add(logYToggle);

	correctForGas = new GLToggle(0, "Surface->Volume conversion");
	correctForGas->SetVisible(FALSE);
	Add(correctForGas);

	formulaText = new GLTextField(0,"");
	formulaText->SetEditable(TRUE);
	Add(formulaText);

	formulaBtn = new GLButton(0,"-> Plot");
	Add(formulaBtn);

	// Center dialog
	int wS,hS;
	GLToolkit::GetScreenSize(&wS,&hS);
	int xD = (wS-wD)/2;
	int yD = (hS-hD)/2;
	SetBounds(xD,yD,wD,hD);
	SetResizable(TRUE);
	SetMinimumSize(wD,220);

	RestoreDeviceObjects();

}

void ProfilePlotter::SetBounds(int x,int y,int w,int h) {

	chart->SetBounds(7,5,w-15,h-110);
	profCombo->SetBounds(7,h-95,180,19);
	selButton->SetBounds(190,h-95,80,19);
	addButton->SetBounds(275,h-95,80,19);
	removeButton->SetBounds(360,h-95,80,19);
	resetButton->SetBounds(445,h-95,80,19);
	logYToggle->SetBounds(190,h-70,40,19);
	correctForGas->SetBounds(240, h - 70, 80, 19);
	normLabel->SetBounds(7,h-68,50,19);
	normCombo->SetBounds(61,h-70,125,19);
	formulaText->SetBounds(7,h-45,350,19);
	formulaBtn->SetBounds(360,h-45,80,19);;
	dismissButton->SetBounds(w-100,h-45,90,19);

	GLWindow::SetBounds(x,y,w,h);
}

void ProfilePlotter::Refresh() {

	if(!worker) return;

	Geometry *geom = worker->GetGeometry();
	int nb = geom->GetNbFacet();
	int nbProf = 0;
	for(int i=0;i<nb;i++)
		if(geom->GetFacet(i)->sh.isProfile) nbProf++;
	profCombo->Clear();profCombo->SetSelectedIndex(0);
	if(nbProf) profCombo->SetSize(nbProf);
	nbProf=0;
	for(int i=0;i<nb;i++) {
		Facet *f = geom->GetFacet(i);
		if(f->sh.isProfile) {
			char tmp[128];
			sprintf(tmp,"F#%d %s",i+1,profType[f->sh.profileType]);
			profCombo->SetValueAt(nbProf,tmp,i);
			profCombo->SetSelectedIndex(0);
			nbProf++;
		}
	}
	//Remove profiles that aren't present anymore
	for (int v=0;v<nbView;v++) 
		if (views[v]->userData>=geom->GetNbFacet() || !geom->GetFacet(views[v]->userData)->sh.isProfile) {
			chart->GetY1Axis()->RemoveDataView(views[v]);
			SAFE_DELETE(views[v]);
			for(int j=v;j<nbView-1;j++) views[j] = views[j+1];
			nbView--;
		}

	refreshViews();

}

void ProfilePlotter::Display(Worker *w) {

	worker = w;
	Refresh();
	SetVisible(TRUE);

}

void ProfilePlotter::Update(float appTime,BOOL force) {

	if(!IsVisible() || IsIconic()) return;  

	if(force) {
		refreshViews();
		lastUpdate = appTime;
		return;
	}

	if( (appTime-lastUpdate>1.0f || force) && nbView ) {
		if(worker->running) refreshViews();
		lastUpdate = appTime;
	}

}

void ProfilePlotter::plot() {

	GLParser *parser = new GLParser();
	parser->SetExpression( formulaText->GetText() );
	if( !parser->Parse() ) {
		GLMessageBox::Display(parser->GetErrorMsg(),"Error",GLDLG_OK,GLDLG_ICONERROR);
		SAFE_DELETE(parser);
		return;
	}

	int nbVar = parser->GetNbVariable();
	if( nbVar==0 ) {
		GLMessageBox::Display("Variable 'x' not found","Error",GLDLG_OK,GLDLG_ICONERROR);
		SAFE_DELETE(parser);
		return;
	}
	if( nbVar>1 ) {
		GLMessageBox::Display("Too much variables or unknown constant","Error",GLDLG_OK,GLDLG_ICONERROR);
		SAFE_DELETE(parser);
		return;
	}
	VLIST *var = parser->GetVariableAt(0);
	if(_stricmp(var->name,"x")!=0) {
		GLMessageBox::Display("Variable 'x' not found","Error",GLDLG_OK,GLDLG_ICONERROR);
		SAFE_DELETE(parser);
		return;
	}

	Geometry *geom = worker->GetGeometry();
	GLDataView *v;

	// Check that view is not already added
	BOOL found = FALSE;
	int i = 0; 
	while(i<nbView && !found) {
		found = (views[i]->userData == -1);
		if(!found) i++;
	}

	if( found ) {
		v = views[i];
		v->SetName(formulaText->GetText());
		v->Reset();
	} else {
		if(nbView<50) {
			v = new GLDataView();
			v->SetName(formulaText->GetText());
			v->userData = -1;
			chart->GetY1Axis()->AddDataView(v);
			views[nbView] = v;
			nbView++;
		}
	}

	// Plot
	for(int i=0;i<1000;i++) {
		double x=(double)i;
		double y;
		var->value = x;
		parser->Evaluate(&y);
		v->Add(x,y,FALSE);
	}
	v->CommitChange();

	delete parser;

}

void ProfilePlotter::refreshViews() {

	// Lock during update
	BYTE *buffer = worker->GetHits();
	int displayMode = normCombo->GetSelectedIndex();
	if(!buffer) return;

	Geometry *geom = worker->GetGeometry();
	SHGHITS *gHits = (SHGHITS *)buffer;
	double nbDes = (double)gHits->total.hit.nbDesorbed;

	double scaleX,scaleY;

	for(int i=0;i<nbView;i++) {

		GLDataView *v = views[i];
		if( v->userData>=0 && v->userData<geom->GetNbFacet()) {
			Facet *f = geom->GetFacet(v->userData);
			v->Reset();
			APROFILE *profilePtr = (APROFILE *)(buffer + f->sh.hitOffset + sizeof(SHHITS)+worker->displayedMoment*sizeof(APROFILE)*PROFILE_SIZE);
			SHHITS *fCount = (SHHITS *)(buffer + f->sh.hitOffset);
			double fnbHit = (double)fCount->hit.nbHit;
			if (fnbHit == 0.0) fnbHit = 1.0;
			if (nbDes>0){
				switch (displayMode) {
				case 0: //Raw data
					for (int j = 0; j < PROFILE_SIZE; j++)
						v->Add((double)j, (double)profilePtr[j].count, FALSE);
					break;

				case 1: //Pressure
					scaleY = 1.0 / nbDes / (f->sh.area / (double)PROFILE_SIZE*1E-4)* worker->gasMass / 1000 / 6E23 * 0.0100; //0.01: Pa->mbar
					/*scaleY *= ((worker->displayedMoment == 0) ? 1.0 : ((worker->desorptionStopTime - worker->desorptionStartTime)
						/ worker->timeWindowSize)); //correction for time window length*/
					scaleY *= ((worker->displayedMoment == 0) ? worker->finalOutgassingRate : (worker->totalDesorbedMolecules
						/ worker->timeWindowSize));
					if (f->sh.is2sided) scaleY *= 0.5;
					//if(f->sh.opacity>0.0) scaleY *= f->sh.opacity;
					//if(IS_ZERO(f->sh.opacity)) scaleY*=2; //transparent profiles are profiled only once...

					for (int j = 0; j < PROFILE_SIZE; j++)
						v->Add((double)j, profilePtr[j].sum_v_ort*scaleY, FALSE);
					break;
				case 2: //Velocity
					scaleX = f->sh.maxSpeed / (double)PROFILE_SIZE;
					for (int j = 0; j < PROFILE_SIZE; j++) {
						if (!correctForGas->GetState()) {
							v->Add((double)j*scaleX, (double)profilePtr[j].count / fnbHit, FALSE);
						}
						else {
							v->Add((double)j*scaleX, (double)profilePtr[j].count / (((double)j+0.5)*scaleX) / fnbHit, FALSE);
						}
					}
					break;
				case 3: //Angle (deg)
					scaleX = 90.0 / (double)PROFILE_SIZE;
					for (int j = 0; j < PROFILE_SIZE ; j++) {
						if (!correctForGas->GetState()) {
							v->Add((double)j*scaleX, (double)profilePtr[j].count / fnbHit, FALSE);
						}
						else {
							v->Add((double)j*scaleX, (double)profilePtr[j].count / fnbHit / sin(((double)j + 0.5)*PI / 2.0 / (double)PROFILE_SIZE), FALSE);
						}
					}

					break;
				case 4: //To 1 (max value)
					llong max = 1;
					for (int j = 0; j<PROFILE_SIZE; j++)
					{
						if (profilePtr[j].count>max) max = profilePtr[j].count;
					}
					scaleY = 1.0 / (double)max;
					for (int j = 0; j < PROFILE_SIZE; j++)
						v->Add((double)j, (double)profilePtr[j].count*scaleY, FALSE);
					break;
				}
			}
			v->CommitChange();
		} else {
			if( v->userData==-2 && nbDes!=0.0 ) {

				// Volatile profile
				v->Reset();
				int nb = geom->GetNbFacet();
				for(int j=0;j<nb;j++) {
					Facet *f = geom->GetFacet(j);
					if( f->sh.isVolatile ) {
						SHHITS *fCount = (SHHITS *)(buffer + f->sh.hitOffset);
						double z = geom->GetVertex(f->indices[0])->z;
						v->Add(z,(double)(fCount->hit.nbAbsorbed)/nbDes,FALSE);
					}
				}
				// Last
				Facet *f = geom->GetFacet(28);
				SHHITS *fCount = (SHHITS *)(buffer + f->sh.hitOffset);
				double fnbAbs = (double)fCount->hit.nbAbsorbed;
				v->Add(1000.0,fnbAbs/nbDes,FALSE);
				v->CommitChange();

				//v->Reset();
				//for(int j=0;j<BOUNCEMAX && nbAbs;j++)
				//  v->Add((double)j,(double)gHits->wallHits[j]/nbAbs);

			}
		}

		}

		worker->ReleaseHits();

	}


	void ProfilePlotter::addView(int facet) {

		char tmp[128];
		Geometry *geom = worker->GetGeometry();

		// Check that view is not already added
		BOOL found = FALSE;
		int i = 0; 
		while(i<nbView && !found) {
			found = (views[i]->userData == facet);
			if(!found) i++;
		}
		if( found ) {
			GLMessageBox::Display("Profile already plotted","Error",GLDLG_OK,GLDLG_ICONERROR);
			return;
		}
		if(nbView<50) {
			Facet *f = geom->GetFacet(facet);
			GLDataView *v = new GLDataView();
			sprintf(tmp,"F#%d %s",facet+1,profType[f->sh.profileType]);
			v->SetName(tmp);
			v->SetColor(*colors[nbView%nbColors]);
			v->SetMarkerColor(*colors[nbView%nbColors]);
			v->SetLineWidth(2);
			v->userData = facet;
			chart->GetY1Axis()->AddDataView(v);
			views[nbView] = v;
			nbView++;
		}

	}

	void ProfilePlotter::remView(int facet) {

		Geometry *geom = worker->GetGeometry();

		BOOL found = FALSE;
		int i = 0; 
		while(i<nbView && !found) {
			found = (views[i]->userData == facet);
			if(!found) i++;
		}
		if( !found ) {
			GLMessageBox::Display("Profile not plotted","Error",GLDLG_OK,GLDLG_ICONERROR);
			return;
		}
		chart->GetY1Axis()->RemoveDataView(views[i]);
		SAFE_DELETE(views[i]);
		for(int j=i;j<nbView-1;j++) views[j] = views[j+1];
		nbView--;

	}

	void ProfilePlotter::Reset() {

		chart->GetY1Axis()->ClearDataView();
		for(int i=0;i<nbView;i++) SAFE_DELETE(views[i]);
		nbView=0;

	}

	void ProfilePlotter::ProcessMessage(GLComponent *src,int message) {
		Geometry *geom = worker->GetGeometry();
		switch(message) {
		case MSG_BUTTON:
			if(src==dismissButton) {
				SetVisible(FALSE);
			} else if(src==selButton) {
				int idx = profCombo->GetSelectedIndex();
				geom->UnSelectAll();
				geom->GetFacet(profCombo->GetUserValueAt(idx))->selected = TRUE;
				geom->UpdateSelection();
				mApp->UpdateFacetParams(TRUE);
				mApp->facetList->SetSelectedRow(profCombo->GetUserValueAt(idx));
				mApp->facetList->ScrollToVisible(profCombo->GetUserValueAt(idx),1,TRUE);
			} else if(src==addButton) {
				int idx = profCombo->GetSelectedIndex();
				if(idx>=0) addView(profCombo->GetUserValueAt(idx));
				refreshViews();
			} else if(src==removeButton) {
				int idx = profCombo->GetSelectedIndex();
				if(idx>=0) remView(profCombo->GetUserValueAt(idx));
				refreshViews();
			} else if(src==resetButton) {
				Reset();
			} else if(src==formulaBtn) {
				plot();
			}
			break;
		case MSG_COMBO:
			if( src==normCombo ) {
				int normMode=normCombo->GetSelectedIndex();
				correctForGas->SetVisible(normMode == 2 || normMode == 3);
				refreshViews();
			}
			break;
		case MSG_TOGGLE:
			if (src==logYToggle) {
				chart->GetY1Axis()->SetScale(logYToggle->GetState());
			}
			else if (src == correctForGas) {
				refreshViews();
			}
			break;
		}

		GLWindow::ProcessMessage(src,message);

	}

