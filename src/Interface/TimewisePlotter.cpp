
#include "TimewisePlotter.h"
#include "ProfileModes.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLMessageBox.h"
#include "GLApp/GLToggle.h"
#include "Helper/MathTools.h"
#include "Helper/StringHelper.h"
#include "GLApp/GLList.h"
#include "GLApp/GLChart/GLChart.h"
#include "GLApp/GLLabel.h"
#include "GLApp/GLCombo.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLFormula.h"
#include "GLApp/GLTextField.h"
#include "Geometry_shared.h"
#include "Facet_shared.h"
#include <math.h>
#if defined(MOLFLOW)
#include "MolFlow.h"
#endif

#if defined(SYNRAD)
#include "SynRad.h"
#endif

#if defined(MOLFLOW)
extern MolFlow *mApp;
#endif

#if defined(SYNRAD)
extern SynRad*mApp;
#endif

/**
* \brief Constructor with initialisation for Time settings window (Time/Timewise plotter)
*/
TimewisePlotter::TimewisePlotter() :GLWindow() {

	int wD = 750;
	int hD = 400;

	SetTitle("Timewise plotter");
	SetIconfiable(true);
	nbView = 0;
	worker = NULL;
	lastUpdate = 0.0f;

	chart = new GLChart(0);
	chart->SetBorder(BORDER_BEVEL_IN);
	chart->GetY1Axis()->SetGridVisible(true);
	chart->GetXAxis()->SetGridVisible(true);
	chart->GetY1Axis()->SetAutoScale(true);
	chart->GetY2Axis()->SetAutoScale(true);
	chart->GetY1Axis()->SetAnnotation(VALUE_ANNO);
	chart->GetXAxis()->SetAnnotation(VALUE_ANNO);
	Add(chart);

	selButton = new GLButton(0, "Show Facet");
	Add(selButton);

	profCombo = new GLCombo(0);
	profCombo->SetEditable(true);
	Add(profCombo);

	normLabel = new GLLabel("Normalize");
	Add(normLabel);

	displayModeCombo = new GLCombo(0);
	displayModeCombo->SetEditable(true);
	size_t nbDisplayModes = (size_t)ProfileDisplayModes::NUMITEMS;
	displayModeCombo->SetSize(nbDisplayModes);
	for (size_t i = 0;i<nbDisplayModes;i++) {
		displayModeCombo->SetValueAt(i, profileDisplayModeDescriptions[(ProfileDisplayModes)i]);
	}
	displayModeCombo->SetSelectedIndex(1);
	Add(displayModeCombo);

	momLabel = new GLLabel("Displayed moments:");
	Add(momLabel);

	momentsText = new GLTextField(0, "1,1,50");
	momentsText->SetEditable(true);
	Add(momentsText);

	momentsLabel = new GLLabel("32 moments");
	Add(momentsLabel);

	correctForGas = new GLToggle(0, "Surface->Volume conversion");
	correctForGas->SetVisible(false);
	Add(correctForGas);

	constantFlowToggle = new GLToggle(0, "Display constant flow");
	Add(constantFlowToggle);

	logYToggle = new GLToggle(0, "Log Y");
	Add(logYToggle);

	warningLabel = new GLLabel("Profiles can only be used on rectangular facets.");
	Add(warningLabel);

	dismissButton = new GLButton(0, "Dismiss");
	Add(dismissButton);

	// Center dialog
	int wS, hS;
	GLToolkit::GetScreenSize(&wS, &hS);
	int xD = (wS - wD) / 2;
	int yD = (hS - hD) / 2;
	SetBounds(xD, yD, wD, hD);
	SetResizable(true);
	SetMinimumSize(wD, 220);

	RestoreDeviceObjects();

}

/**
* \brief Sets positions and sizes of the window
* \param x x-coordinate of the element
* \param y y-coordinate of the element
* \param w width of the element
* \param h height of the element
*/
void TimewisePlotter::SetBounds(int x, int y, int w, int h) {

	chart->SetBounds(7, 5, w - 15, h - 85);
	profCombo->SetBounds(7, h - 70, 117, 19);
	selButton->SetBounds(295, h - 70, 80, 19);
	normLabel->SetBounds(130, h - 68, 50, 19);
	displayModeCombo->SetBounds(185, h - 70, 105, 19);
	correctForGas->SetBounds(w - 340, h - 70, 150, 19);

	constantFlowToggle->SetBounds(w - 180, h - 70, 120, 19);
	logYToggle->SetBounds(w - 55, h - 70, 40, 19);
	momLabel->SetBounds(30, h - 45, 117, 19);
	momentsText->SetBounds(130, h - 45, 180, 19);
	momentsLabel->SetBounds(315, h - 45, 60, 19);

	warningLabel->SetBounds(w - 340, h - 45, 235, 19);
	dismissButton->SetBounds(w - 100, h - 45, 90, 19);

	GLWindow::SetBounds(x, y, w, h);

}

/**
* \brief Refresh everything inside the window
*/
void TimewisePlotter::Refresh() {
	Reset();
	if (!worker) return;
	if (!ParseMoments()) return;
	InterfaceGeometry *interfGeom = worker->GetGeometry();
	size_t nb = interfGeom->GetNbFacet();
	size_t nbProf = 0;
	for (size_t i = 0; i < nb; i++)
		if (interfGeom->GetFacet(i)->sh.isProfile) nbProf++;
	profCombo->Clear();
	if (nbProf) profCombo->SetSize(nbProf);
	nbProf = 0;
    for (size_t i = 0; i < nb; i++) {
		InterfaceFacet *f = interfGeom->GetFacet(i);
		if (f->sh.isProfile) {
			std::ostringstream tmp;
			tmp << "F#" << (i + 1) << " " << profileRecordModeDescriptions[(ProfileRecordModes)f->sh.profileType].second; //short description
			profCombo->SetValueAt(nbProf, tmp.str(), (int)i);
			nbProf++;
		}
	}
	profCombo->SetSelectedIndex(nbProf ? 0 : -1);
	if (nbProf>0 && nbView == 0) addView(profCombo->GetUserValueAt(0));
	//Remove profiles that aren't present anymore
	if (nbView>0)
		if (views[0]->userData1 >= interfGeom->GetNbFacet() || !interfGeom->GetFacet(views[0]->userData1)->sh.isProfile) {
			Reset();
		}
	refreshViews();
}

/**
* \brief Displays window with refreshed values
* \param w worker handle
*/
void TimewisePlotter::Display(Worker *w) {

	/*
	if( nbView==0 ) {
	GLDataView *v = new GLDataView();
	v->SetName("Transmission Prob.");
	v->userData1 = -2;
	GLColor c;
	c.r=0;c.g=255;c.b=0;
	v->SetColor(c);
	chart->GetY1Axis()->AddDataView(v);
	views[nbView] = v;
	nbView++;
	}
	*/

	worker = w;
	Refresh();
	SetVisible(true);

}

/**
* \brief Refreshes the view if needed
* \param appTime current time of the applicaiton
* \param force if view should be refreshed no matter what
*/
void TimewisePlotter::Update(float appTime, bool force) {

	if (!IsVisible() || IsIconic()) return;

	if (force) {
		refreshViews();
		lastUpdate = appTime;
		return;
	}
}

/**
* \brief Refreshes view by updating the data for the plot depending on the selected normalisation mode
*/
void TimewisePlotter::refreshViews() {

    if(!nbView) return;

	// Lock during update
	if (!worker->ReloadIfNeeded()) return;
	auto lock = GetHitLock(worker->globalState.get(), 10000);
	if (!lock) return;
	ProfileDisplayModes displayMode = (ProfileDisplayModes)displayModeCombo->GetSelectedIndex(); //Choosing by index is error-prone


	InterfaceGeometry *interfGeom = worker->GetGeometry();

	double scaleY;

	size_t facetHitsSize = (1 + worker->interfaceMomentCache.size()) * sizeof(FacetHitBuffer);
	for (size_t i = 0; i < nbView; i++) {

		GLDataView *v = views[i];
		UpdateMoment();
		if (v->userData1<0 || v->userData1>worker->interfaceMomentCache.size()) continue; //invalid moment
		int idx = profCombo->GetSelectedIndex();
		if (idx < 0) return;
		InterfaceFacet *f = interfGeom->GetFacet(profCombo->GetUserValueAt(idx));
		v->Reset();
		//FacetHitBuffer *fCount = (FacetHitBuffer *)(buffer + f->sp.hitOffset);
		//double fnbHit = (double)fCount->hit.nbMCHit;
		/*int momentIndex;
		if (m==(nbView-1) && constantFlowToggle->GetState()) momentIndex=0; //Constant flow
		else momentIndex=m+1; //any other 'normal' moment*/
		const std::vector<ProfileSlice>& profile = worker->globalState->facetStates[profCombo->GetUserValueAt(idx)].momentResults[v->userData1].profile;
		//ProfileSlice *profilePtr = (ProfileSlice *)(buffer + f->sh.hitOffset + facetHitsSize + v->userData1*sizeof(ProfileSlice)*PROFILE_SIZE);
			if (worker->globalStatCache.globalHits.nbDesorbed > 0){

				if (displayMode == ProfileDisplayModes::Raw) {
					for (int j = 0; j < PROFILE_SIZE; j++)
						v->Add((double)j, profile[j].countEquiv, false);
				}
				else if (displayMode == ProfileDisplayModes::Pressure) {
					scaleY = 1.0 / (f->GetArea() * 1E-4 / (double)PROFILE_SIZE)* worker->model->sp.gasMass / 1000 / 6E23 * 0.0100; //0.01: Pa->mbar
					scaleY *= worker->GetMoleculesPerTP(v->userData1);

					for (int j = 0; j < PROFILE_SIZE; j++)
						v->Add((double)j, profile[j].sum_v_ort*scaleY, false);
				}
				else if (displayMode == ProfileDisplayModes::ImpRate) {

					scaleY = 1.0 / (f->GetArea() * 1E-4 / (double)PROFILE_SIZE);
					scaleY *= worker->GetMoleculesPerTP(v->userData1);

					for (int j = 0; j < PROFILE_SIZE; j++)
						v->Add((double)j, profile[j].countEquiv * scaleY, false);
				}
				else if (displayMode == ProfileDisplayModes::Density) {
					scaleY = 1.0 / ((f->GetArea() * 1E-4) / (double)PROFILE_SIZE);
					scaleY *= worker->GetMoleculesPerTP(v->userData1) * f->DensityCorrection();
					
					for (int j = 0; j < PROFILE_SIZE; j++)
						v->Add((double)j, profile[j].sum_1_per_ort_velocity*scaleY, false);
				}
				else if (displayMode == ProfileDisplayModes::Speed) {
					double sum = 0.0;
					double val;
					double scaleX = f->sh.maxSpeed / (double)PROFILE_SIZE;
					std::vector<double> values;
					values.reserve(PROFILE_SIZE);
					for (int j = 0; j < PROFILE_SIZE; j++) {//count distribution sum
						if (!correctForGas->GetState())
							val = profile[j].countEquiv;
						else
							val = profile[j].countEquiv / (((double)j + 0.5)*scaleX); //fnbhit not needed, sum will take care of normalization
						sum += val;
						values.push_back(val);
					}

					for (int j = 0; j < PROFILE_SIZE; j++)
						v->Add((double)j*scaleX, values[j] / sum, false);
				}
				else if (displayMode == ProfileDisplayModes::Angle) {
					double sum = 0.0;
					double val;
					double scaleX = 90.0 / (double)PROFILE_SIZE;
					std::vector<double> values;
					values.reserve(PROFILE_SIZE);
					for (int j = 0; j < PROFILE_SIZE; j++) {//count distribution sum
						if (!correctForGas->GetState())
							val = profile[j].countEquiv;
						else
							val = profile[j].countEquiv / sin(((double)j + 0.5)*PI / 2.0 / (double)PROFILE_SIZE); //fnbhit not needed, sum will take care of normalization
						sum += val;
						values.push_back(val);
					}

					for (int j = 0; j < PROFILE_SIZE; j++)
						v->Add((double)j*scaleX, values[j] / sum, false);
					break;
				}
				else if (displayMode == ProfileDisplayModes::NormalizeTo1) {
                    double max = 1.0;

                    for (int j = 0; j < PROFILE_SIZE; j++) {
                        max = std::max(max,profile[j].countEquiv);
                    }
                    scaleY = 1.0 / (double) max;

                    for (int j = 0; j < PROFILE_SIZE; j++)
                        v->Add((double) j, profile[j].countEquiv * scaleY, false);
                    break;
                }
				else{
                    // Unknown display mode, reset to RAW data
                    displayModeCombo->SetSelectedIndex(0);
                    break;
				}

			}
			v->CommitChange();
	}
}

/**
* \brief Uses a specified facet to the plot
* \param facet ID of the facet that should be added
*/
void TimewisePlotter::addView(int facet) {

	char tmp[128];
	InterfaceGeometry *interfGeom = worker->GetGeometry();

	InterfaceFacet *f = interfGeom->GetFacet(facet);

	if (constantFlowToggle->GetState()) { //add constant flow
		GLDataView *v = new GLDataView();
		sprintf(tmp, "Moment0 (Constant Flow)"/*, facet + 1, profType[f->sp.profileType]*/);
		v->SetName(tmp);
		v->userData1 = 0;
		v->SetStyle(STYLE_DOT);
        v->SetColor(ColorSchemes::defaultCol[0]);
        //v->SetLineWidth(2);
		chart->GetY1Axis()->AddDataView(v);
		views[nbView] = v;
		nbView++;
	}

	for (size_t index : displayedMoments) {
		if (nbView < 49) {
			GLDataView *v = new GLDataView();
			if (index<1 || (index) > worker->interfaceMomentCache.size()) continue; //invalid moment
			//sprintf(tmp, "Moment%d (t=%gs)", (int)index, worker->moments[index - 1]);
			sprintf(tmp, "t=%gs", worker->interfaceMomentCache[index - 1].time);
			v->SetName(tmp);
			v->userData1 = (int)index;
			v->SetStyle(STYLE_DOT);
            v->SetColor(ColorSchemes::defaultCol[0]);
            chart->GetY1Axis()->AddDataView(v);
			views[nbView] = v;
			nbView++;
		}
	}

}

/**
* \brief Removes a specified facet from the plot (TODO: check if funcion is needed)
* \param facet ID of the facet that should be removed
*/
void TimewisePlotter::remView(int facet) {

	InterfaceGeometry *interfGeom = worker->GetGeometry();

	bool found = false;
	int i = 0;
	while (i < nbView && !found) {
		found = (views[i]->userData1 == facet);
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

/**
* \brief Clears data inside the data view
*/
void TimewisePlotter::Reset() {

	chart->GetY1Axis()->ClearDataView();
	for (int i = 0; i < nbView; i++) SAFE_DELETE(views[i]);
	nbView = 0;

}

/**
* \brief Function for processing various inputs (button, check boxes etc.)
* \param src Exact source of the call
* \param message Type of the source (button)
*/
void TimewisePlotter::ProcessMessage(GLComponent *src, int message) {
	InterfaceGeometry *interfGeom = worker->GetGeometry();
	switch (message) {
	case MSG_BUTTON:
		if (src == dismissButton) {
			SetVisible(false);
		}
		else if (src == selButton) {
			int idx = profCombo->GetSelectedIndex();
			if (idx >= 0 && idx < interfGeom->GetNbFacet()) {
				interfGeom->UnselectAll();
				interfGeom->GetFacet(profCombo->GetUserValueAt(idx))->selected = true;
				interfGeom->UpdateSelection();
				mApp->UpdateFacetParams(true);
				mApp->facetList->SetSelectedRow(profCombo->GetUserValueAt(idx));
				mApp->facetList->ScrollToVisible(profCombo->GetUserValueAt(idx), 1, true);
			}
		} /*else if(src==addButton) {
			int idx = profCombo->GetSelectedIndex();
			if(idx>=0) addView(profCombo->GetUserValueAt(idx));
			refreshViews();
			} else if(src==removeButton) {
			int idx = profCombo->GetSelectedIndex();
			if(idx>=0) remView(profCombo->GetUserValueAt(idx));
			refreshViews();
			} else if(src==removeAllButton) {
			Reset();
			}*/
	case MSG_COMBO:
		if (src == displayModeCombo) {
			ProfileDisplayModes normMode = (ProfileDisplayModes)displayModeCombo->GetSelectedIndex();
			correctForGas->SetVisible(normMode == ProfileDisplayModes::Speed || normMode == ProfileDisplayModes::Angle);
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

/**
* \brief Changes line style to highlight the selected moment
*/
void TimewisePlotter::UpdateMoment() {
	for (int i = 0; i < nbView; i++) {

		GLDataView *v = views[i];
		if (worker->displayedMoment == v->userData1) {
			v->SetStyle(STYLE_SOLID);
			v->SetLineWidth(2);
		}
		else {
			v->SetStyle(STYLE_DOT);
			v->SetLineWidth(1);
		}
	}
}

/**
* \brief Parses (sets of) moments with ; delimiter
* \return true if parsing successful
*/
bool TimewisePlotter::ParseMoments(){
	//Quick string parsing from http://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c
	std::string s = momentsText->GetText();
	std::vector<std::string> tokens = SplitString(s,';');
	displayedMoments.clear();
	for (const auto& token : tokens) {
		ParseToken(token);
	}
	return true;
}

/**
* \brief Evaluates and parses moments or intervals of moments of the form "begin,interval,end" e.g. 1,2,50
* \param token string of a moment or interval of moments
*/
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

/**
* \brief Formats the text for how many moments will be displayed
*/
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