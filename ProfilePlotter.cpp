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

extern GLApplication *theApp;

extern double gasMass;
extern double totalOutgassing;

static const char*profType[] = {"None","Pressure \201","Pressure \202","Angle"};

ProfilePlotter::ProfilePlotter():GLWindow() {

	int wD = 750;
	int hD = 400;

	SetTitle("Profile plotter");
	SetIconfiable(TRUE);
	nbView = 0;
	worker = NULL;
	lastUpdate = 0.0f;

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

	normLabel = new GLLabel("Normalize");
	Add(normLabel);
	//qLabel = new GLLabel("Q=");
	//Add(qLabel);
	//unitLabel = new GLLabel("units*l/s");
	//Add(unitLabel);

	normCombo = new GLCombo(0);
	normCombo->SetEditable(TRUE);
	normCombo->SetSize(8);
	normCombo->SetValueAt(0,"None");
	normCombo->SetValueAt(1,"Absorption (total)");
	normCombo->SetValueAt(2,"Desorption (total)");
	normCombo->SetValueAt(3,"Hit           (total)");
	normCombo->SetValueAt(4,"Absorption (local)");
	normCombo->SetValueAt(5,"Desorption (local)");
	normCombo->SetValueAt(6,"Hit           (local)");
	normCombo->SetValueAt(7,"Pressure");
	normCombo->SetSelectedIndex(7);
	Add(normCombo);

	logYToggle = new GLToggle(0,"Log Y");
	Add(logYToggle);

	//showAllMoments = new GLToggle(0,"Show moments");
	//Add(showAllMoments);

	formulaText = new GLTextField(0,"");
	formulaText->SetEditable(TRUE);
	Add(formulaText);
	//qText = new GLTextField(0,"1");
	//qText->SetEditable(TRUE);
	//Add(qText);


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

	chart->SetBounds(7,5,w-15,h-85);
	profCombo->SetBounds(7,h-70,117,19);
	selButton->SetBounds(130,h-70,80,19);
	addButton->SetBounds(215,h-70,80,19);
	removeButton->SetBounds(300,h-70,80,19);
	resetButton->SetBounds(385,h-70,80,19);
	logYToggle->SetBounds(w-55,h-70,40,19);
	normLabel->SetBounds(470,h-68,50,19);
	normCombo->SetBounds(520,h-70,105,19);
	//showAllMoments->SetBounds(640,h-70,100,19);
	formulaText->SetBounds(7,h-45,525,19);
	formulaBtn->SetBounds(537,h-45,80,19);
	//qLabel->SetBounds(437+105+5,h-70,10,19);
	//unitLabel->SetBounds(437+105+25+30+5,h-70,20,19);
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

	/*
	if( nbView==0 ) {
	GLDataView *v = new GLDataView();
	v->SetName("Transmission Prob.");
	v->userData = -2;
	GLCColor c;
	c.r=0;c.g=255;c.b=0;
	v->SetColor(c);
	chart->GetY1Axis()->AddDataView(v);
	views[nbView] = v;
	nbView++;
	}
	*/

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
		if(nbView<32) {
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
	int normalize = normCombo->GetSelectedIndex();
	//bool isPressureDisplayed=(normalize==7);
	//qText->SetVisible(isPressureDisplayed);
	//qLabel->SetVisible(isPressureDisplayed);
	//unitLabel->SetVisible(isPressureDisplayed);
	if(!buffer) return;

	Geometry *geom = worker->GetGeometry();
	SHGHITS *gHits = (SHGHITS *)buffer;
	double nbAbs = (double)gHits->total.hit.nbAbsorbed;
	double nbDes = (double)gHits->total.hit.nbDesorbed;
	double nbHit = (double)gHits->total.hit.nbHit;

	double scale;

	for(int i=0;i<nbView;i++) {

		GLDataView *v = views[i];
		if( v->userData>=0 && v->userData<geom->GetNbFacet()) {
			Facet *f = geom->GetFacet(v->userData);
			SHHITS *fCount = (SHHITS *)(buffer + f->sh.hitOffset);
			double fnbAbs = (double)fCount->hit.nbAbsorbed;
			double fnbDes = (double)fCount->hit.nbDesorbed;
			double fnbHit = (double)fCount->hit.nbHit;
			//double q;
			v->Reset();
			llong *profilePtr = (llong *)(buffer + f->sh.hitOffset + sizeof(SHHITS)+worker->displayedMoment*sizeof(llong)*PROFILE_SIZE);

			switch(normalize) {
			case 0:
				for(int j=0;j<PROFILE_SIZE;j++)
					v->Add((double)j,(double)profilePtr[j],FALSE);
				break;
			case 1:
				for(int j=0;j<PROFILE_SIZE && nbAbs!=0.0;j++)
					v->Add((double)j,(double)profilePtr[j]/nbAbs,FALSE);
				break;
			case 2:
				for(int j=0;j<PROFILE_SIZE && nbDes!=0.0;j++)
					v->Add((double)j,(double)profilePtr[j]/nbDes,FALSE);
				break;
			case 3:
				for(int j=0;j<PROFILE_SIZE && nbHit!=0.0;j++)
					v->Add((double)j,(double)profilePtr[j]/nbHit,FALSE);
				break;
			case 4:
				for(int j=0;j<PROFILE_SIZE && fnbAbs!=0.0;j++)
					v->Add((double)j,(double)profilePtr[j]/fnbAbs,FALSE);
				break;
			case 5:
				for(int j=0;j<PROFILE_SIZE && fnbDes!=0.0;j++)
					v->Add((double)j,(double)profilePtr[j]/fnbDes,FALSE);
				break;
			case 6:
				for(int j=0;j<PROFILE_SIZE && fnbHit!=0.0;j++)
					v->Add((double)j,(double)profilePtr[j]/fnbHit,FALSE);
				break;
			case 7: //Pressure
				scale = 40.0*totalOutgassing*((worker->displayedMoment==0)?1.0:((worker->desorptionStopTime-worker->desorptionStartTime)
					/worker->timeWindowSize))/(145.469*sqrt(f->sh.temperature
					/gasMass)*(f->sh.area/PROFILE_SIZE)*nbDes);
				if(f->sh.is2sided) scale = scale / 2.0;
				if(f->sh.opacity>0.0) scale = scale * f->sh.opacity;
				for(int j=0;j<PROFILE_SIZE && fnbHit!=0.0;j++)
					v->Add((double)j,(double)profilePtr[j]*scale,FALSE);
				break;
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
		if(nbView<32) {
			Facet *f = geom->GetFacet(facet);
			GLDataView *v = new GLDataView();
			sprintf(tmp,"F#%d %s",facet+1,profType[f->sh.profileType]);
			v->SetName(tmp);
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
		MolFlow *mApp = (MolFlow *)theApp;
		switch(message) {
		case MSG_BUTTON:
			if(src==dismissButton) {
				SetVisible(FALSE);
			} else if(src==selButton) {
				int idx = profCombo->GetSelectedIndex();
				geom->UnSelectAll();
				geom->GetFacet(profCombo->GetUserValueAt(idx))->selected = TRUE;
				mApp->UpdateFacetParams(TRUE);
				geom->UpdateSelection();
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
				refreshViews();
			}
			break;
		case MSG_TOGGLE:
			if (src==logYToggle) {
				chart->GetY1Axis()->SetScale(logYToggle->IsChecked());
			}
			break;
		}

		GLWindow::ProcessMessage(src,message);

	}

