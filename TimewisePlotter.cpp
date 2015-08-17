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

extern MolFlow *mApp;

static const char*profType[] = { "None", "Pressure \201 [mbar]", "Pressure \202 [mbar]", "Angle", "Velocity", "Ort.velocity" };

TimewisePlotter::TimewisePlotter() :GLWindow() {

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

	dismissButton = new GLButton(0, "Dismiss");
	Add(dismissButton);

	selButton = new GLButton(0, "Show Facet");
	Add(selButton);

	profCombo = new GLCombo(0);
	profCombo->SetEditable(TRUE);
	Add(profCombo);

	normLabel = new GLLabel("Normalize");
	Add(normLabel);

	normCombo = new GLCombo(0);
	normCombo->SetEditable(TRUE);
	normCombo->SetSize(6);
	normCombo->SetValueAt(0, "None (raw data)");
	normCombo->SetValueAt(1, "Pressure (mbar)");
	normCombo->SetValueAt(2, "Density (1/m3)");
	normCombo->SetValueAt(3, "Speed (m/s)");
	normCombo->SetValueAt(4, "Angle (deg)");
	normCombo->SetValueAt(5, "Normalize to 1");

	normCombo->SetSelectedIndex(1);
	Add(normCombo);

	momLabel = new GLLabel("Displayed moments:");
	Add(momLabel);

	momentsText = new GLTextField(0, "1,1,32");
	momentsText->SetEditable(TRUE);
	Add(momentsText);

	momentsLabel = new GLLabel("32 moments");
	Add(momentsLabel);

	correctForGas = new GLToggle(0, "Surface->Volume conversion");
	correctForGas->SetVisible(FALSE);
	Add(correctForGas);

	constantFlowToggle = new GLToggle(0, "Display constant flow");
	Add(constantFlowToggle);

	logYToggle = new GLToggle(0, "Log Y");
	Add(logYToggle);

	

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

void TimewisePlotter::SetBounds(int x, int y, int w, int h) {

	chart->SetBounds(7, 5, w - 15, h - 85);
	profCombo->SetBounds(7, h - 70, 117, 19);
	selButton->SetBounds(295, h - 70, 80, 19);
	normLabel->SetBounds(130, h - 68, 50, 19);
	normCombo->SetBounds(185, h - 70, 105, 19);
	correctForGas->SetBounds(w - 340, h - 70, 150, 19);

	constantFlowToggle->SetBounds(w - 180, h - 70, 120, 19);
	logYToggle->SetBounds(w - 55, h - 70, 40, 19);
	momLabel->SetBounds(30, h - 45, 117, 19);
	momentsText->SetBounds(130, h - 45, 180, 19);
	momentsLabel->SetBounds(315, h - 45, 60, 19);
	dismissButton->SetBounds(w - 100, h - 45, 90, 19);

	GLWindow::SetBounds(x, y, w, h);

}

void TimewisePlotter::Refresh() {
	Reset();
	if (!worker) return;
	if (!ParseMoments()) return;
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
	if (nbProf>0 && nbView == 0) addView(profCombo->GetUserValueAt(0));
	//Remove profiles that aren't present anymore
	if (nbView>0)
		if (views[0]->userData >= geom->GetNbFacet() || !geom->GetFacet(views[0]->userData)->sh.isProfile) {
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

void TimewisePlotter::Update(float appTime, BOOL force) {

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


void TimewisePlotter::refreshViews() {

	// Lock during update
	BYTE *buffer = worker->GetHits();
	int displayMode = normCombo->GetSelectedIndex();
	if (!buffer) return;

	Geometry *geom = worker->GetGeometry();
	SHGHITS *gHits = (SHGHITS *)buffer;
	double nbDes = (double)gHits->total.hit.nbDesorbed;

	double scaleY;

	for (size_t i = 0; i < nbView; i++) {

		GLDataView *v = views[i];
		UpdateMoment();
		if (v->userData<0 || v->userData>worker->moments.size()) continue; //invalid moment
		int idx = profCombo->GetSelectedIndex();
		if (idx < 0) return;
		Facet *f = geom->GetFacet(profCombo->GetUserValueAt(idx));
		v->Reset();
		SHHITS *fCount = (SHHITS *)(buffer + f->sh.hitOffset);
		double fnbHit = (double)fCount->hit.nbHit;
		/*int momentIndex;
		if (m==(nbView-1) && constantFlowToggle->GetState()) momentIndex=0; //Constant flow
		else momentIndex=m+1; //any other 'normal' moment*/
		APROFILE *profilePtr = (APROFILE *)(buffer + f->sh.hitOffset + sizeof(SHHITS) + v->userData*sizeof(APROFILE)*PROFILE_SIZE);

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
		case 2: //Particle density
			scaleY = 1.0 / nbDes / (f->sh.area / (double)PROFILE_SIZE*1E-4);
			scaleY *= ((worker->displayedMoment == 0) ? worker->finalOutgassingRate : (worker->totalDesorbedMolecules
				/ worker->timeWindowSize));
			if (f->sh.is2sided) scaleY *= 0.5;
			//Correction for double-density effect (measuring density on desorbing/absorbing facets):
			if (f->sh.counter.hit.nbHit>0 || f->sh.counter.hit.nbDesorbed>0)
				if (f->sh.counter.hit.nbAbsorbed >0 || f->sh.counter.hit.nbDesorbed>0) //otherwise save calculation time
				scaleY *= 1.0 - ((double)f->sh.counter.hit.nbAbsorbed + (double)f->sh.counter.hit.nbDesorbed) / ((double)f->sh.counter.hit.nbHit + (double)f->sh.counter.hit.nbDesorbed) / 2.0;
			for (int j = 0; j < PROFILE_SIZE; j++)
				v->Add((double)j, profilePtr[j].sum_1_per_ort_velocity*scaleY, FALSE);
			break;
		case 3: {//Velocity
			double sum = 0.0;
			double val;
			double scaleX = f->sh.maxSpeed / (double)PROFILE_SIZE;
			std::vector<double> values;
			values.reserve(PROFILE_SIZE);
			for (int j = 0; j < PROFILE_SIZE; j++) {//count distribution sum
				if (!correctForGas->GetState())
					val = (double)profilePtr[j].count;
				else
					val = (double)profilePtr[j].count / (((double)j + 0.5)*scaleX); //fnbhit not needed, sum will take care of normalization
				sum += val;
				values.push_back(val);
			}
			for (int j = 0; j < PROFILE_SIZE; j++)
				v->Add((double)j*scaleX, values[j] / sum, FALSE);
			break; }
		case 4: {//Angle
			double sum = 0.0;
			double val;
			double scaleX = 90.0 / (double)PROFILE_SIZE;
			std::vector<double> values;
			values.reserve(PROFILE_SIZE);
			for (int j = 0; j < PROFILE_SIZE; j++) {//count distribution sum
				if (!correctForGas->GetState())
					val = (double)profilePtr[j].count;
				else
					val = (double)profilePtr[j].count / sin(((double)j + 0.5)*PI / 2.0 / (double)PROFILE_SIZE); //fnbhit not needed, sum will take care of normalization
				sum += val;
				values.push_back(val);
			}
			for (int j = 0; j < PROFILE_SIZE; j++)
				v->Add((double)j*scaleX, values[j] / sum, FALSE);
			break; }
		case 5: //To 1 (max value)
			llong max = 1;
			for (int j = 0; j < PROFILE_SIZE; j++)
			{
				if (profilePtr[j].count > max) max = profilePtr[j].count;
			}
			scaleY = 1.0 / (double)max;
			for (int j = 0; j < PROFILE_SIZE; j++)
				v->Add((double)j, (double)profilePtr[j].count*scaleY, FALSE);
			break;
		}
		v->CommitChange();


	}

	worker->ReleaseHits();

}


void TimewisePlotter::addView(int facet) {

	char tmp[128];
	Geometry *geom = worker->GetGeometry();

	Facet *f = geom->GetFacet(facet);

	if (constantFlowToggle->GetState()) { //add constant flow
		GLDataView *v = new GLDataView();
		sprintf(tmp, "Moment0 (Constant Flow)", facet + 1, profType[f->sh.profileType]);
		v->SetName(tmp);
		v->userData = 0;
		v->SetStyle(STYLE_DOT);
		//v->SetLineWidth(2);
		chart->GetY1Axis()->AddDataView(v);
		views[nbView] = v;
		nbView++;
	}



	for (size_t index : displayedMoments) {
		if (nbView < 49) {
			GLDataView *v = new GLDataView();
			if (index<1 || (index) > worker->moments.size()) continue; //invalid moment
			sprintf(tmp, "Moment%d (t=%gs)", (int)index, worker->moments[index - 1]);
			v->SetName(tmp);
			v->userData = (int)index;
			v->SetStyle(STYLE_DOT);
			chart->GetY1Axis()->AddDataView(v);
			views[nbView] = v;
			nbView++;
		}
	}



}

void TimewisePlotter::remView(int facet) {

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

void TimewisePlotter::Reset() {

	chart->GetY1Axis()->ClearDataView();
	for (int i = 0; i < nbView; i++) SAFE_DELETE(views[i]);
	nbView = 0;

}

void TimewisePlotter::ProcessMessage(GLComponent *src, int message) {
	Geometry *geom = worker->GetGeometry();
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
			}*/
	case MSG_COMBO:
		if (src == normCombo) {
			int normMode = normCombo->GetSelectedIndex();
			correctForGas->SetVisible(normMode == 3 || normMode == 4);
			refreshViews();
		}
		else if (src == profCombo) {
			Reset();
			int idx = profCombo->GetSelectedIndex();
			if (idx >= 0) {
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
		if (src == logYToggle) {
			chart->GetY1Axis()->SetScale(logYToggle->GetState());
		}
		else if (src == constantFlowToggle) {
			Reset();
			int idx = profCombo->GetSelectedIndex();
			if (idx >= 0) {
				addView(profCombo->GetUserValueAt(idx));
			}
			FormatParseText();
			refreshViews();
		}
		else if (src == correctForGas) {
			refreshViews();
		}
		break;
	case MSG_TEXT:
		if (src == momentsText) {
			if (!ParseMoments()) return;
			Reset();
			int idx = profCombo->GetSelectedIndex();
			if (idx >= 0) {
				addView(profCombo->GetUserValueAt(idx));
			}
			FormatParseText();
			refreshViews();
		}
		break;
	}

	GLWindow::ProcessMessage(src, message);

}

void TimewisePlotter::UpdateMoment() {
	for (int i = 0; i < nbView; i++) {

		GLDataView *v = views[i];
		if (worker->displayedMoment == v->userData) {
			v->SetStyle(STYLE_SOLID);
			v->SetLineWidth(2);
		}
		else {
			v->SetStyle(STYLE_DOT);
			v->SetLineWidth(1);
		}
	}
}

BOOL TimewisePlotter::ParseMoments(){
	//Quick string parsing from http://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c
	std::string s = momentsText->GetText();
	std::string delimiter = ";";
	size_t pos = 0;
	std::string token;
	displayedMoments = std::vector<size_t>();
	while ((pos = s.find(delimiter)) != std::string::npos) {
		token = s.substr(0, pos);
		//process tokens
		ParseToken(token);
		s.erase(0, pos + delimiter.length());
	}
	//last token
	ParseToken(s);
	return TRUE;
}

void TimewisePlotter::ParseToken(std::string token) {
	int begin, interval, end;
	int nb = sscanf(token.c_str(), "%d,%d,%d", &begin, &interval, &end);
	if (nb == 1 && (begin > 0)) {
		//One moment index
		displayedMoments.push_back(begin);
	}
	else if (nb == 3 && (begin > 0) && (end > begin) && (interval < (end - begin))) {
		//Range
		for (int time = begin; time <= end; time += interval)
			displayedMoments.push_back(time);
	}
}

void TimewisePlotter::FormatParseText() {
	std::ostringstream tmp;
	int dispNum = (int)displayedMoments.size() + constantFlowToggle->GetState();
	tmp << dispNum << " moments";
	if (dispNum > 50) {
		tmp << ", displaying first 50";
		momentsLabel->SetTextColor(255, 0, 0);
	}
	else {
		momentsLabel->SetTextColor(0, 0, 0);
	}
	momentsLabel->SetText(tmp.str().c_str());
}