/*
File:        TimewisePlotter.cpp
Description: Timewise Profile plotter window
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

#include "TimewisePlotter.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLMessageBox.h"
#include "Utils.h"
#include <math.h>
#include "Molflow.h"

extern GLApplication *theApp;

/*extern double gasMass;
extern double totalOutgassing;
extern double totalInFlux;*/

static const char*profType[] = {"None","Pressure \201 [mbar]","Pressure \202 [mbar]","Angle","Velocity","Ort.velocity"};

TimewisePlotter::TimewisePlotter():GLWindow() {

	int wD = 750;
	int hD = 400;

	SetTitle("Timewise plotter");
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

	profCombo = new GLCombo(0);
	profCombo->SetEditable(TRUE);
	Add(profCombo);

	normLabel = new GLLabel("Normalize");
	Add(normLabel);

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
	normCombo->SetValueAt(7,"Pressure [mbar]");
	normCombo->SetSelectedIndex(7);
	Add(normCombo);

	constantFlowToggle = new GLToggle(0,"Display constant flow");
	Add(constantFlowToggle);

	formulaText = new GLTextField(0,"");
	formulaText->SetEditable(TRUE);
	Add(formulaText);

	logYToggle = new GLToggle(0,"Log Y");
	Add(logYToggle);

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

void TimewisePlotter::SetBounds(int x,int y,int w,int h) {

	chart->SetBounds(7,5,w-15,h-85);
	profCombo->SetBounds(7,h-70,117,19);
	selButton->SetBounds(130,h-70,80,19);
	/*addButton->SetBounds(215,h-70,80,19);
	removeButton->SetBounds(300,h-70,80,19);
	resetButton->SetBounds(385,h-70,80,19);*/
	normLabel->SetBounds(300,h-68,50,19);
	normCombo->SetBounds(355,h-70,105,19);

	constantFlowToggle->SetBounds(w-180,h-70,120,19);
	logYToggle->SetBounds(w-55,h-70,40,19);
	
	//showAllMoments->SetBounds(640,h-70,100,19);
	formulaText->SetBounds(7,h-45,525,19);
	formulaBtn->SetBounds(537,h-45,80,19);
	
	//qLabel->SetBounds(437+105+5,h-70,10,19);
	//unitLabel->SetBounds(437+105+25+30+5,h-70,20,19);
	dismissButton->SetBounds(w-100,h-45,90,19);

	GLWindow::SetBounds(x,y,w,h);

}

void TimewisePlotter::Refresh() {
	Reset();
	if(!worker) return;

	Geometry *geom = worker->GetGeometry();
	int nb = geom->GetNbFacet();
	int nbProf = 0;
	for(int i=0;i<nb;i++)
		if(geom->GetFacet(i)->sh.isProfile) nbProf++;
	profCombo->Clear();
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
	if (nbProf>0 && nbView==0) addView(profCombo->GetUserValueAt(0));
	//Remove profiles that aren't present anymore
	if (nbView>0) 
		if (views[0]->userData>=geom->GetNbFacet() || !geom->GetFacet(views[0]->userData)->sh.isProfile) {
			Reset();
		}
	refreshViews();
}

void TimewisePlotter::Display(Worker *w) {

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

void TimewisePlotter::Update(float appTime,BOOL force) {

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

void TimewisePlotter::plot() {

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
	for(int p=0;p<1000;p++) {
		double x=(double)p;
		double y;
		var->value = x;
		parser->Evaluate(&y);
		v->Add(x,y,FALSE);
	}
	v->CommitChange();

	delete parser;

}

void TimewisePlotter::refreshViews() {

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

	for(int m=0;m<nbView;m++) {


		GLDataView *v = views[m];
		UpdateMoment();
		if( v->userData>=0 && v->userData<geom->GetNbFacet()) {
			Facet *f = geom->GetFacet(v->userData);
			SHHITS *fCount = (SHHITS *)(buffer + f->sh.hitOffset);
			double fnbAbs = (double)fCount->hit.nbAbsorbed;
			double fnbDes = (double)fCount->hit.nbDesorbed;
			double fnbHit = (double)fCount->hit.nbHit;
			//double q;
			v->Reset();
			int momentIndex;
			if (m==(nbView-1) && constantFlowToggle->IsChecked()) momentIndex=0; //Constant flow
			else momentIndex=m+1; //any other 'normal' moment
			APROFILE *profilePtr = (APROFILE *)(buffer + f->sh.hitOffset + sizeof(SHHITS)+momentIndex*sizeof(APROFILE)*PROFILE_SIZE);

			switch(normalize) {
			case 0:
				for(int j=0;j<PROFILE_SIZE;j++)
					v->Add((double)j,(double)profilePtr[j].count,FALSE);
				break;
			case 1:
				for(int j=0;j<PROFILE_SIZE && nbAbs!=0.0;j++)
					v->Add((double)j, (double)profilePtr[j].count / nbAbs, FALSE);
				break;
			case 2:
				for(int j=0;j<PROFILE_SIZE && nbDes!=0.0;j++)
					v->Add((double)j, (double)profilePtr[j].count / nbDes, FALSE);
				break;
			case 3:
				for(int j=0;j<PROFILE_SIZE && nbHit!=0.0;j++)
					v->Add((double)j, (double)profilePtr[j].count / nbHit, FALSE);
				break;
			case 4:
				for(int j=0;j<PROFILE_SIZE && fnbAbs!=0.0;j++)
					v->Add((double)j, (double)profilePtr[j].count / fnbAbs, FALSE);
				break;
			case 5:
				for(int j=0;j<PROFILE_SIZE && fnbDes!=0.0;j++)
					v->Add((double)j, (double)profilePtr[j].count / fnbDes, FALSE);
				break;
			case 6:
				for(int j=0;j<PROFILE_SIZE && fnbHit!=0.0;j++)
					v->Add((double)j, (double)profilePtr[j].count / fnbHit, FALSE);
				break;
			case 7: //Pressure
				scale = /*totalInFlux*/ 1.0 / nbDes / (f->sh.area / (double)PROFILE_SIZE*1E-4)* worker->gasMass / 1000 / 6E23 * 0.0100; //0.01: Pa->mbar
				/*scale *= ((momentIndex == 0) ? 1.0 : ((worker->desorptionStopTime - worker->desorptionStartTime)
					/ worker->timeWindowSize)); //correction for time window length*/
				//CHECK
				scale *= ((momentIndex == 0) ? worker->finalOutgassingRate : (worker->totalDesorbedMolecules / worker->timeWindowSize));
				if (f->sh.is2sided) scale *= 0.5;
				//if (f->sh.opacity>0.0) scale *= f->sh.opacity;
				for(int j=0;j<PROFILE_SIZE && fnbHit!=0.0;j++)
					v->Add((double)j,(double)profilePtr[j].sum_v_ort*scale,FALSE);
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


	void TimewisePlotter::addView(int facet) {

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
			GLMessageBox::Display("Profile already plotted","Error",GLDLG_OK,GLDLG_ICONWARNING);
			return;
		}
		
		Facet *f = geom->GetFacet(facet);
		
		if(nbView<49) {
			
			for (size_t m=0;m<MIN(worker->moments.size(),30);m++) {
				GLDataView *v = new GLDataView();
				sprintf(tmp,"F#%d %s - t=%gs",facet+1,profType[f->sh.profileType],worker->moments[m]);
				v->SetName(tmp);
				v->userData = facet;
				v->SetStyle(STYLE_DOT);
				chart->GetY1Axis()->AddDataView(v);
				views[nbView] = v;
				nbView++;
			}
		}

		if (constantFlowToggle->IsChecked()) { //add constant flow
			GLDataView *v = new GLDataView();
			sprintf(tmp,"F#%d %s - Constant Flow",facet+1,profType[f->sh.profileType]);
			v->SetName(tmp);
			v->userData = facet;
			v->SetStyle(STYLE_DOT);
			//v->SetLineWidth(2);
			chart->GetY1Axis()->AddDataView(v);
			views[nbView] = v;
			nbView++;
		}

	}

	void TimewisePlotter::remView(int facet) {

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

	void TimewisePlotter::Reset() {

		chart->GetY1Axis()->ClearDataView();
		for(int i=0;i<nbView;i++) SAFE_DELETE(views[i]);
		nbView=0;

	}

	void TimewisePlotter::ProcessMessage(GLComponent *src,int message) {
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
				geom->UpdateSelection();
				mApp->UpdateFacetParams(TRUE);
				mApp->facetList->SetSelectedRow(profCombo->GetUserValueAt(idx));
				mApp->facetList->ScrollToVisible(profCombo->GetUserValueAt(idx),1,TRUE);
			} /*else if(src==addButton) {
				int idx = profCombo->GetSelectedIndex();
				if(idx>=0) addView(profCombo->GetUserValueAt(idx));
				refreshViews();
			} else if(src==removeButton) {
				int idx = profCombo->GetSelectedIndex();
				if(idx>=0) remView(profCombo->GetUserValueAt(idx));
				refreshViews();
			} else if(src==resetButton) {
				Reset();
			}*/ else if(src==formulaBtn) {
				plot();
			}
		case MSG_COMBO:
			if( src==normCombo ) {
				refreshViews();
			} else if(src==profCombo) {
				Reset();
				int idx = profCombo->GetSelectedIndex();
				if(idx>=0) {
					addView(profCombo->GetUserValueAt(idx));
				}
				refreshViews();
			}
			/*
			case MSG_TEXT: //enter pressed
			if( src==qText ) {
			refreshViews();
			}
			break;
			*/
		case MSG_TOGGLE:
			if (src==logYToggle) {
				chart->GetY1Axis()->SetScale(logYToggle->IsChecked());
			} else if (src==constantFlowToggle) {
				Reset();
				int idx = profCombo->GetSelectedIndex();
				if(idx>=0) {
					addView(profCombo->GetUserValueAt(idx));
				}
				refreshViews();
			}
			break;
		
		}

		GLWindow::ProcessMessage(src,message);

	}

void TimewisePlotter::UpdateMoment() {
	for(int i=0;i<nbView;i++) {

		GLDataView *v = views[i];
		if ((i==(nbView-1) && constantFlowToggle->IsChecked() && worker->displayedMoment==0) || i==(worker->displayedMoment-1)) {
			v->SetStyle(STYLE_SOLID);
			v->SetLineWidth(2);
		} else {
			v->SetStyle(STYLE_DOT);
			v->SetLineWidth(1);
		}
	}
}