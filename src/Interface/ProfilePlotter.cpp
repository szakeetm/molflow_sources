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
#include "ProfilePlotter.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLMessageBox.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLLabel.h"
#include "GLApp/GLCombo.h"
#include "GLApp/GLTextField.h"
#include "GLApp/GLToggle.h"
#include "Helper/MathTools.h"
#include "GLApp/GLList.h"
#include "GLApp/GLChart/GLChart.h"
#include "Geometry_shared.h"
#include "Facet_shared.h"
#include <cmath>
#include <Helper/StringHelper.h>

#if defined(MOLFLOW)
#include "MolFlow.h"
#endif
#if defined(SYNRAD)
#include "SynRad.h"
#endif

extern GLApplication *theApp;

#if defined(MOLFLOW)
extern MolFlow *mApp;
#endif

#if defined(SYNRAD)
extern SynRad*mApp;
#endif

const char* profType[] = {
	"None",
	"Pressure \201 [mbar]",
	"Pressure \202 [mbar]",
	"Incident angle [deg]",
	"Speed [m/s]",
	"Ort. velocity [m/s]",
	"Tan. velocity [m/s]"
};

/**
* \brief Constructor with initialisation for Profile plotter window (Tools/Profile Plotter)
*/
ProfilePlotter::ProfilePlotter() :GLWindow() , views{}{

	int wD = 650;
	int hD = 400;

	SetTitle("Profile plotter");
	SetIconfiable(true);
	nbView = 0;
	worker = nullptr;

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

	dismissButton = new GLButton(0, "Dismiss");
	Add(dismissButton);

	selButton = new GLButton(0, "Show Facet");
	Add(selButton);

	addButton = new GLButton(0, "Add curve");
	Add(addButton);

	removeButton = new GLButton(0, "Remove curve");
	Add(removeButton);

	removeAllButton = new GLButton(0, "Remove all");
	Add(removeAllButton);

	profCombo = new GLCombo(0);
	profCombo->SetEditable(true);
    Add(profCombo);

    selFacInput = new GLTextField(0, "...");
    selFacInput->SetEditable(true);

    Add(selFacInput);

	normLabel = new GLLabel("Normalize:");
	Add(normLabel);

	normCombo = new GLCombo(0);
	normCombo->SetEditable(true);
	normCombo->SetSize(6);
	normCombo->SetValueAt(0, "None (raw data)");
	normCombo->SetValueAt(1, "Pressure (mbar)");
	normCombo->SetValueAt(2, "Density (1/m3)");
	normCombo->SetValueAt(3, "Speed (m/s)");
	normCombo->SetValueAt(4, "Angle (deg)");
	normCombo->SetValueAt(5, "Normalize to 1");
	normCombo->SetSelectedIndex(1);
	Add(normCombo);

	logYToggle = new GLToggle(0, "Log Y");
	Add(logYToggle);

    colorToggle = new GLToggle(0, "Colorblind mode");
    Add(colorToggle);

    fixedLineWidthText = new GLLabel("Change linewidth:");
    Add(fixedLineWidthText);
    fixedLineWidthButton = new GLButton(0, "-> Apply linewidth");
    Add(fixedLineWidthButton);
    fixedLineWidthField = new GLTextField(0, "2");
    fixedLineWidthField->SetEditable(true);
    Add(fixedLineWidthField);

    useProfColToggle = new GLToggle(0, "Identify profiles in geometry");
    Add(useProfColToggle);

	selectPlottedButton = new GLButton(0, "Select plotted facets");
	Add(selectPlottedButton);

	warningLabel = new GLLabel("Profiles can only be used on rectangular facets.");
	Add(warningLabel);

	correctForGas = new GLToggle(0, "Surface->Volume conversion");
	correctForGas->SetVisible(false);
	Add(correctForGas);

	formulaText = new GLTextField(0, "");
	formulaText->SetEditable(true);
	Add(formulaText);

	formulaBtn = new GLButton(0, "-> Plot expression");
	Add(formulaBtn);

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
* \brief Sets positions and sizes of all UI elements
* \param x x-coordinate of the element
* \param y y-coordinate of the element
* \param w width of the element
* \param h height of the element
*/
void ProfilePlotter::SetBounds(int x, int y, int w, int h) {
	
	chart->SetBounds(7, 5, w - 15, h - 135);

	profCombo->SetBounds(7, h - 120, 160, 19);
    selFacInput->SetBounds(170, h - 120, 120, 19);
    selButton->SetBounds(295, h - 120, 80, 19);
	addButton->SetBounds(380, h - 120, 80, 19);
	removeButton->SetBounds(465, h - 120, 80, 19);
	removeAllButton->SetBounds(550, h - 120, 80, 19);

    logYToggle->SetBounds(190, h - 95, 40, 19);
	warningLabel->SetBounds(w-240,h-95,235,19);
	correctForGas->SetBounds(240, h - 95, 80, 19);

	normLabel->SetBounds(7, h - 93, 50, 19);
	normCombo->SetBounds(61, h - 95, 125, 19);

    colorToggle->SetBounds(7, h - 70, 105, 19);
    fixedLineWidthText->SetBounds(112, h - 70, 93, 19);
    fixedLineWidthField->SetBounds(206, h - 70, 30, 19);
    fixedLineWidthButton->SetBounds(240, h - 70, 100, 19);
    useProfColToggle->SetBounds(350, h - 70, 105, 19);
	selectPlottedButton->SetBounds(w-130,h-70,120,19);

    formulaText->SetBounds(7, h - 45, 350, 19);
	formulaBtn->SetBounds(360, h - 45, 120, 19);;
	dismissButton->SetBounds(w - 100, h - 45, 90, 19);

	GLWindow::SetBounds(x, y, w, h);

}

/**
* \brief Refreshes all window components (combo, chart, values)
*/
void ProfilePlotter::Refresh() {

	if (!worker) return;

	//Rebuild selection combo box
	Geometry *geom = worker->GetGeometry();
	size_t nb = geom->GetNbFacet();
	size_t nbProf = 1; // minimum 1 for custom input
	for (size_t i = 0; i < nb; i++)
		if (geom->GetFacet(i)->sh.isProfile) nbProf++;
	profCombo->Clear();
	if (nbProf) profCombo->SetSize(nbProf);
    nbProf = 0;
    profCombo->SetValueAt(nbProf, "Select [v] or type ->", (int)-1);

    std::stringstream facetInputText;
    if(!nb) facetInputText << "...";
    else if(nb==1) facetInputText << "1";
    else facetInputText << "1-" << nb;
    selFacInput->SetText(facetInputText.str());
    nbProf = 1;
    for (size_t i = 0; i < nb; i++) {
		InterfaceFacet *f = geom->GetFacet(i);
		if (f->sh.isProfile) {
			char tmp[128];
			sprintf(tmp, "F#%zd %s", i + 1, profType[f->sh.profileType]);
			profCombo->SetValueAt(nbProf, tmp, (int)i);
			nbProf++;
		}
	}
	profCombo->SetSelectedIndex(nbProf ? 0 : -1);
	//Remove profiles that aren't present anymore
	for (int v = nbView - 1; v >= 0; v--) { //int because it can be -1, nbView is also int
		if (views[v]->userData1 >= geom->GetNbFacet() || !geom->GetFacet(views[v]->userData1)->sh.isProfile) {
			chart->GetY1Axis()->RemoveDataView(views[v]);
			SAFE_DELETE(views[v]);
			for (size_t j = v; j < nbView - 1; j++) views[j] = views[j + 1];
			nbView--;
		}
	}

	//Update values
	refreshViews();
    ResetHighlighting();

}


/**
* \brief Displays window with refreshed values
* \param w worker handle
*/
void ProfilePlotter::Display(Worker *w) {

	
    SetWorker(w);
	Refresh();
	SetVisible(true);

}

/**
* \brief Refreshes the view if needed
* \param appTime current time of the application
* \param force if view should be refreshed no matter what
*/
void ProfilePlotter::Update(float appTime, bool force) {

	if (!IsVisible() || IsIconic()) return;

	if (force) {
		refreshViews();
		lastUpdate = appTime;
		return;
	}

}

/**
* \brief Creates a plot from the expression given in the textbox of the form f(x)=EXPRESSION (e.g. 2*x+50)
*/
void ProfilePlotter::plot() {

	GLParser *parser = new GLParser();
	parser->SetExpression(formulaText->GetText().c_str());
	if (!parser->Parse()) {
		GLMessageBox::Display(parser->GetErrorMsg(), "Error", GLDLG_OK, GLDLG_ICONERROR);
		SAFE_DELETE(parser);
		return;
	}

	int nbVar = parser->GetNbVariable();
	if (nbVar == 0) {
		GLMessageBox::Display("Variable 'x' not found", "Error", GLDLG_OK, GLDLG_ICONERROR);
		SAFE_DELETE(parser);
		return;
	}
	if (nbVar > 1) {
		GLMessageBox::Display("Too much variables or unknown constant", "Error", GLDLG_OK, GLDLG_ICONERROR);
		SAFE_DELETE(parser);
		return;
	}
	VLIST *var = parser->GetVariableAt(0);
	if (!iequals(var->name, "x")) {
		GLMessageBox::Display("Variable 'x' not found", "Error", GLDLG_OK, GLDLG_ICONERROR);
		SAFE_DELETE(parser);
		return;
	}

	GLDataView *v;

	// Check that view is not already added
	bool found = false;
	int i = 0;
	while (i < nbView && !found) {
		found = (views[i]->userData1 == -1);
		if (!found) i++;
	}

	if (found) {
		v = views[i];
		v->SetName(formulaText->GetText().c_str());
		v->Reset();
	}
	else {

		if (nbView < 50) {
			v = new GLDataView();
			v->SetName(formulaText->GetText().c_str());
			v->userData1 = -1;
			chart->GetY1Axis()->AddDataView(v);
			views[nbView] = v;
			nbView++;
		}
		else {
			return;
		}
	}

	// Plot
	for (i = 0; i < 1000; i++) {
		double x = (double)i;
		double y;
		var->value = x;
		parser->Evaluate(&y);
		v->Add(x, y, false);
	}
	v->CommitChange();

	delete parser;

}

/**
* \brief Refreshes view by updating the data for the plot
*/
void ProfilePlotter::refreshViews() {

    if(!nbView) return;

	// Lock during update
	bool buffer_old = worker->GetHits();
	if (!buffer_old) return;
	int displayMode = normCombo->GetSelectedIndex();

	Geometry *geom = worker->GetGeometry();

	double scaleY;

	size_t facetHitsSize = (1 + worker->moments.size()) * sizeof(FacetHitBuffer);
	for (int i = 0; i < nbView; i++) {

		GLDataView *v = views[i];
		if (v->userData1 >= 0 && v->userData1 < geom->GetNbFacet()) {
			InterfaceFacet *f = geom->GetFacet(v->userData1);

			v->Reset();
			const std::vector<ProfileSlice>& profile = worker->globState.facetStates[v->userData1].momentResults[worker->displayedMoment].profile;

			//FacetHitBuffer *fCount = (FacetHitBuffer *)(buffer + f->wp.hitOffset+ worker->displayedMoment*sizeof(FacetHitBuffer));
			//double fnbHit = (double)fCount->hit.nbMCHit;
			//if (fnbHit == 0.0) fnbHit = 1.0;
			if (worker->globalHitCache.globalHits.nbDesorbed > 0){

				switch (displayMode) {
				case 0: //Raw data
					for (int j = 0; j < PROFILE_SIZE; j++)
						v->Add((double)j, profile[j].countEquiv, false);

					break;

				case 1: //Pressure
					scaleY = 1.0 / (f->GetArea() * 1E-4 / (double)PROFILE_SIZE)* worker->model.wp.gasMass / 1000 / 6E23 * 0.0100; //0.01: Pa->mbar
					scaleY *= worker->GetMoleculesPerTP(worker->displayedMoment);

					for (int j = 0; j < PROFILE_SIZE; j++)
						v->Add((double)j, profile[j].sum_v_ort*scaleY, false);
					break;
				case 2: //Particle density
					scaleY = 1.0 / ((f->GetArea() * 1E-4) / (double)PROFILE_SIZE);
					scaleY *= worker->GetMoleculesPerTP(worker->displayedMoment) * f->DensityCorrection();
					
					for (int j = 0; j < PROFILE_SIZE; j++)
						v->Add((double)j, profile[j].sum_1_per_ort_velocity*scaleY, false);
					break;
				case 3: {//Velocity
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
					break; }
				case 4: {//Angle
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
				case 5: {//To 1 (max value)
                    double max = 1.0;

                    for (int j = 0; j < PROFILE_SIZE; j++) {
                        if (profile[j].countEquiv > max) max = profile[j].countEquiv;
                    }
                    scaleY = 1.0 / (double) max;

                    for (int j = 0; j < PROFILE_SIZE; j++)
                        v->Add((double) j, profile[j].countEquiv * scaleY, false);
                    break;
                }
                default:
                    // Unknown display mode, reset to RAW data
                    normCombo->SetSelectedIndex(0);
                    break;
				}

			}
			v->CommitChange();
		}
		else {

			if (v->userData1 == -2 && worker->globalHitCache.globalHits.nbDesorbed != 0.0) {

				// Volatile profile
				v->Reset();
				size_t nb = geom->GetNbFacet();
				for (size_t j = 0; j < nb; j++) {
					InterfaceFacet *f = geom->GetFacet(j);
					if (f->sh.isVolatile) {
					    const FacetHitBuffer& fCount = worker->globState.facetStates[j].momentResults[worker->displayedMoment].hits;
						double z = geom->GetVertex(f->indices[0])->z;
						v->Add(z,fCount.nbAbsEquiv / worker->globalHitCache.globalHits.nbDesorbed, false);
					}
				}
				// Last
				InterfaceFacet *f = geom->GetFacet(28);
                const FacetHitBuffer& fCount = worker->globState.facetStates[28].momentResults[worker->displayedMoment].hits;
				double fnbAbs = fCount.nbAbsEquiv;
				v->Add(1000.0, fnbAbs / worker->globalHitCache.globalHits.nbDesorbed, false);
				v->CommitChange();

				//v->Reset();
				//for(int j=0;j<BOUNCEMAX && nbAbs;j++)
				//  v->Add((double)j,(double)gHits->wallHits[j]/nbAbs);

			}
		}

	}
	worker->ReleaseHits();
}

/**
* \brief Adds view/plot for a specific facet
* \param facet specific facet ID
* \return 0 if okay, 1 if already plotted
*/
int ProfilePlotter::addView(int facet) {

	char tmp[128];

	// Check that view is not already added
	bool found = false;
	int i = 0;
	while (i < nbView && !found) {
		found = (views[i]->userData1 == facet);
		if (!found) i++;
	}
	if (found) {
		return 1;
	}
	if (nbView < MAX_VIEWS) {
		GLDataView *v = new GLDataView();
		sprintf(tmp, "F#%d", facet + 1);
		v->SetName(tmp);
		//Look for first available color
		GLColor col = chart->GetFirstAvailableColor();
        int lineStyle = chart->GetFirstAvailableLinestyle(col);
		v->SetColor(col);
		v->SetMarkerColor(col);
		v->SetStyle(lineStyle);
		v->SetLineWidth(2);
		v->userData1 = facet;

		chart->GetY1Axis()->AddDataView(v);
		views[nbView] = v;
		nbView++;

        plottedFacets.insert(std::make_pair(facet, col));
    }

	return 0;
}

/**
* \brief Removes view/plot for a specific facet
* \param facet specific facet ID
* \return 0 if okay, 1 if not plotted
*/
int ProfilePlotter::remView(int facet) {

	bool found = false;
	int i = 0;
	while (i < nbView && !found) {
		found = (views[i]->userData1 == facet);
		if (!found) i++;
	}
	if (!found) {
		return 1;
	}
	chart->GetY1Axis()->RemoveDataView(views[i]);
	SAFE_DELETE(views[i]);
	for (int j = i; j < nbView - 1; j++) views[j] = views[j + 1];
	nbView--;

    plottedFacets.erase(facet);
    return 0;
}

/**
* \brief Resets the whole chart
*/
void ProfilePlotter::Reset() {

	chart->GetY1Axis()->ClearDataView();
	for (int i = 0; i < nbView; i++) SAFE_DELETE(views[i]);
	nbView = 0;

	plottedFacets.clear();
    applyFacetHighlighting();
}

/**
* \brief Resets the selection highlighting colors
*/
void ProfilePlotter::ResetHighlighting() {

    plottedFacets.clear();
    for(int viewId = 0; viewId < nbView; viewId++){
        GLDataView *v = views[viewId];
        plottedFacets.insert(std::make_pair(v->userData1, v->GetColor()));

    }

    applyFacetHighlighting();
}

/**
* \brief Function for processing various inputs (button, check boxes etc.)
* \param src Exact source of the call
* \param message Type of the source (button)
*/
void ProfilePlotter::ProcessMessage(GLComponent *src, int message) {
	Geometry *geom = worker->GetGeometry();

	switch (message) {
	case MSG_BUTTON:
		if (src == dismissButton) {
			SetVisible(false);
		}
		else if (src == selButton) {
			int idx = profCombo->GetSelectedIndex();
			if(idx >= 0) {
                geom->UnselectAll();
                size_t facetRow;
                if (idx > 0) { //Something selected, facets start with idx==1, custom input is idx==0 (not -1)
                    int facetId = profCombo->GetUserValueAt(idx);
                    geom->GetFacet(facetId)->selected = true;
                    mApp->facetList->SetSelectedRow(profCombo->GetUserValueAt(idx));
                    facetRow = profCombo->GetUserValueAt(idx);
                } else {
                    std::vector<size_t> facetIds;
                    try {
                        splitFacetList(facetIds, selFacInput->GetText(), geom->GetNbFacet());
                        if(facetIds.empty())
                            return;
                    }
                    catch (std::exception &e) {
                        GLMessageBox::Display(e.what(), "Error", GLDLG_OK, GLDLG_ICONERROR);
                        return;
                    }
                    catch (...) {
                        GLMessageBox::Display("Unknown exception", "Error", GLDLG_OK, GLDLG_ICONERROR);
                        return;
                    }

                    for (const auto &facetId : facetIds) {
                        geom->GetFacet(facetId)->selected = Contains({selButton, addButton}, src);
                    }

                    facetRow = facetIds.back();
                }

                geom->UpdateSelection();
                mApp->UpdateFacetParams(true);
                mApp->UpdateFacetlistSelected();
                mApp->facetList->ScrollToVisible(facetRow, 1, true);
            }
		}
		else if (src == addButton) {
			int idx = profCombo->GetSelectedIndex();
			if (idx >= 0) { //Something selected (not -1)
			    if(idx > 0){
                    if(addView(profCombo->GetUserValueAt(idx)))
                        GLMessageBox::Display("Profile already plotted", "Info", GLDLG_OK, GLDLG_ICONINFO);
                }
			    else {
                    std::vector<size_t> facetIds;
                    try {
                        splitFacetList(facetIds, selFacInput->GetText(), geom->GetNbFacet());
                    }
                    catch (std::exception &e) {
                        GLMessageBox::Display(e.what(), "Error", GLDLG_OK, GLDLG_ICONERROR);
                        return;
                    }

                    for (const auto &facetId : facetIds) {
                        if(geom->GetFacet(facetId)->sh.isProfile) {
                            if(addView(facetId))
                                GLMessageBox::Display("Profile already plotted", "Info", GLDLG_OK, GLDLG_ICONINFO);
                        }
                    }
			    }
				refreshViews();
                applyFacetHighlighting();
            }
		}
		else if (src == removeButton) {

			int idx = profCombo->GetSelectedIndex();
            if (idx >= 0) { //Something selected (not -1)
                if(idx > 0){
                    if(remView(profCombo->GetUserValueAt(idx))){
                        GLMessageBox::Display("Profile not plotted", "Error", GLDLG_OK, GLDLG_ICONERROR);
                    }
                }
                else {
                    std::vector<size_t> facetIds;
                    try {
                        splitFacetList(facetIds, selFacInput->GetText(), geom->GetNbFacet());
                    }
                    catch (std::exception &e) {
                        GLMessageBox::Display(e.what(), "Error", GLDLG_OK, GLDLG_ICONERROR);
                        return;
                    }

                    for (const auto &facetId : facetIds) {
                        /*
						if(remView(facetId)){
                            GLMessageBox::Display("Profile not plotted", "Error", GLDLG_OK, GLDLG_ICONERROR);
                        }
						*/ //Commented out: one wrong click on empty plotter can cause infinite error messages
						remView(facetId); //Try to remove, fail silently
                    }
                }
                refreshViews();
                applyFacetHighlighting();
            }
		}
		else if (src == removeAllButton) {
			Reset();
            applyFacetHighlighting();
		}
		else if (src == formulaBtn) {

			plot();
		}
		else if(src == fixedLineWidthButton) {
            int linW;
            fixedLineWidthField->GetNumberInt(&linW);
            for(int viewId = 0; viewId < nbView; viewId++){
                GLDataView *v = views[viewId];
                v->SetLineWidth(linW);
            }
        }  else if (src == selectPlottedButton) {
			std::vector<size_t> plottedFacetIds;
			for(int viewId = 0; viewId < nbView; viewId++){
                GLDataView *v = views[viewId];
				plottedFacetIds.push_back(v->userData1);
			}
			geom->UnselectAll();
			geom->SetSelection(plottedFacetIds,false,false);
		}
		break;
	case MSG_COMBO:
		if (src == normCombo) {
			int normMode = normCombo->GetSelectedIndex();
			correctForGas->SetVisible(normMode == 3 || normMode == 4);
			refreshViews();
		}
		else if(src == profCombo){
            int profMode = profCombo->GetSelectedIndex();
            selFacInput->SetEditable(!profMode);
        }
		break;
	case MSG_TOGGLE:
		if (src == logYToggle) {
			chart->GetY1Axis()->SetScale(logYToggle->GetState());
		}
		else if (src == correctForGas) {

			refreshViews();

		}
		else if(src == colorToggle) {
		    if(!colorToggle->GetState())
                chart->SetColorSchemeDefault();
		    else
		        chart->SetColorSchemeColorblind();

		    const auto& colors = chart->GetColorScheme();
		    for(int viewId = 0; viewId < nbView; viewId++){
                GLDataView *v = views[viewId];
                auto col = colors[viewId%colors.size()];
                int lineStyle = chart->GetFirstAvailableLinestyle(col);
                v->SetColor(col);
                v->SetMarkerColor(col);
                v->SetStyle(lineStyle);

                std::string facId(v->GetName());
                facId = facId.substr(facId.find('#') + 1);
                plottedFacets.at(std::stoi(facId)-1) = col;
		    }
            applyFacetHighlighting();
		}
		else if(src == useProfColToggle) {
            applyFacetHighlighting();
        }
		break;

    default:
        break;
	}
    //this->worker->GetGeometry()->SetPlottedFacets(plottedFacets);
	GLWindow::ProcessMessage(src, message);

}

/**
* \brief Adds views to the plotter if loaded form a file (XML)
* \param updatedViews vector containing the ids of the views
*/
void ProfilePlotter::SetViews(const std::vector<int> &updatedViews) {
	Reset();
	for (int view : updatedViews)
		if (view<worker->GetGeometry()->GetNbFacet() && worker->GetGeometry()->GetFacet(view)->sh.isProfile)
			addView(view);
	//Refresh(); //Commented out: at this point, simulation results are not yet loaded
}

/**
* \brief Create and return a vector of view IDs
* \return vector containing the IDs of the views
*/
std::vector<int> ProfilePlotter::GetViews() {
	std::vector<int>v;
	v.reserve(nbView);
	for (size_t i = 0; i < nbView; i++)
		v.push_back(views[i]->userData1);
	return v;
}

/**
* \brief Returns bool for active/inactive logarithmic scaling
* \return bool that expresses if Y axis is logarithmically scaled
*/
bool ProfilePlotter::IsLogScaled() {
	return chart->GetY1Axis()->GetScale();
}

/**
* \brief Sets type of scale (logarithmic or not)
* \param bool that expresses if Y axis should be logarithmically scaled or not
*/
void ProfilePlotter::SetLogScaled(bool logScale){
	chart->GetY1Axis()->SetScale(logScale);
	logYToggle->SetState(logScale);
}

/**
* \brief Sets worker handle for loading views before the full geometry
* \param w worker handle
*/
void ProfilePlotter::SetWorker(Worker *w) { //for loading views before the full geometry

	worker = w;

}

/**
* \brief Returns a list of plotted facets and the coupled colors
* \return map containing plotted facet IDs and colors
*/
std::map<int,GLColor> ProfilePlotter::GetIDColorPairs() const {
    return plottedFacets;
}

/**
* \brief Applies or removes facet highlighting
*/
void ProfilePlotter::applyFacetHighlighting() const {
    auto geom = this->worker->GetGeometry();
    if(useProfColToggle->GetState()){
        geom->SetPlottedFacets(plottedFacets);
    }
    else {
        geom->SetPlottedFacets(std::map<int, GLColor>());
    }
    geom->UpdateSelection();
}
