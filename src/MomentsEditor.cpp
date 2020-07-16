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
#include "MomentsEditor.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLMessageBox.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLLabel.h"
#include "GLApp/GLToggle.h"
#include "GLApp/GLTitledPanel.h"
#include "GLApp/GLTextField.h"
#include "GLApp/GLList.h"
#include "MolFlow.h"
#include "TimeSettings.h"
#include "TimewisePlotter.h"

extern MolFlow *mApp;

static const int   flWidth[] = { 35,100,35,80 };
static const char *flName[] = { "#","Time (s)","Nb","Window (s)" };
static const int   flAligns[] = { ALIGN_CENTER,ALIGN_LEFT,ALIGN_CENTER,ALIGN_LEFT };
static const int   fEdits[] = { 0,EDIT_STRING,0, EDIT_STRING };

/**
* \brief Constructor for initial creation of the window
* \param w vector GUI worker in charge for this window
*/
MomentsEditor::MomentsEditor(Worker *w) :GLWindow() {

	int wD = 290;
	int hD = 401;

	work = w;

	SetTitle("Edit time moments");

	panel1 = new GLTitledPanel("Moment list");
	panel1->SetBounds(5, 5, wD - 10, 250);
	Add(panel1);

	momentsList = new GLList(0);
	momentsList->SetBounds(10, 22, wD - 20, 200);
	momentsList->SetColumnLabelVisible(true);
	momentsList->SetGrid(true);
	Add(momentsList);

	clearButton = new GLButton(0, "Clear list");
	clearButton->SetBounds(10, 229, 95, 20);
	Add(clearButton);

	pasteButton = new GLButton(0, "Paste clipboard");
	pasteButton->SetBounds(110, 229, 95, 20);
	pasteButton->SetEnabled(false);
	Add(pasteButton);

	//char tmp[128];

	panel2 = new GLTitledPanel("Time parameters");
	panel2->SetBounds(5, 260, wD - 10, 95);
	Add(panel2);

	/*GLLabel *startLabel = new GLLabel("Desorption starts at:                   s");
	startLabel->SetBounds(15,275,170,25);
	Add(startLabel);*/

	//sprintf(tmp,"%g",work->desorptionStartTime);
  /*  desStartText = new GLTextField(0,"");
	desStartText->SetBounds(120,275,60,20);
	Add(desStartText);*/

	/*GLLabel *stopLabel = new GLLabel("Desorption stops at:                   s");
	stopLabel->SetBounds(15,300,170,25);
	Add(stopLabel);*/

	//sprintf(tmp,"%g",work->desorptionStopTime);
	/*desStopText = new GLTextField(0,"");
	desStopText->SetBounds(120,300,60,20);
	Add(desStopText);*/

	GLLabel *windowLabel = new GLLabel("Default time window:                  s");
	windowLabel->SetBounds(15, 275, 170, 25);
	Add(windowLabel);

	//sprintf(tmp,"%g",work->wp.wp.timeWindowSize);
	windowSizeText = new GLTextField(0, "");
	windowSizeText->SetBounds(120, 275, 60, 20);
	Add(windowSizeText);

	useMaxwellToggle = new GLToggle(0, "Use Maxwell-B. speed distr.");
	useMaxwellToggle->SetBounds(15, 300, wD - 25, 20);
	//useMaxwellToggle->SetState(work->wp.useMaxwellDistribution);
	Add(useMaxwellToggle);

	calcConstantFlow = new GLToggle(0, "Calculate constant flow");
	calcConstantFlow->SetBounds(15, 325, wD - 25, 20);
	//useMaxwellToggle->SetState(work->wp.useMaxwellDistribution);
	Add(calcConstantFlow);

	/*GLLabel *valveLabel = new GLLabel("Facets 1,2 open at:                   s");
	valveLabel->SetBounds(15,400,170,25);
	Add(valveLabel);*/

	//sprintf(tmp,"%g",work->desorptionStopTime);
   /* valveText = new GLTextField(0,"");
	valveText->SetBounds(120,400,60,20);
	Add(valveText);*/

	//RebuildList();

	setButton = new GLButton(0, "Apply");
	setButton->SetBounds(wD - 165, hD - 44, 75, 20);
	Add(setButton);

	cancelButton = new GLButton(0, "Dismiss");
	cancelButton->SetBounds(wD - 85, hD - 44, 75, 20);
	Add(cancelButton);

	// Center dialog
	int wS, hS;
	GLToolkit::GetScreenSize(&wS, &hS);
	int xD = (wS - wD) / 2;
	int yD = (hS - hD) / 2;
	SetBounds(xD, yD, wD, hD);

	RestoreDeviceObjects();

}

/**
* \brief Function for processing various inputs (button, check boxes etc.)
* \param src Exact source of the call
* \param message Type of the source (button)
*/
void MomentsEditor::ProcessMessage(GLComponent *src, int message) {
	switch (message) {
	case MSG_BUTTON:

		if (src == cancelButton) {

			GLWindow::ProcessMessage(NULL, MSG_CLOSE);

		}
		else if (src == setButton) {
			//validate user input
			double window;
			if (!(windowSizeText->GetNumber(&window))) {
				GLMessageBox::Display("Invalid window length", "Error", GLDLG_OK, GLDLG_ICONERROR);
				return;
			}

			//apply settings
			if (mApp->AskToReset()) {
				std::vector<std::vector<Moment>> parsedMoments;
				for (size_t u = 0; u != userMoments.size(); u++) {
					parsedMoments.emplace_back(ParseMoment(userMoments[u].first, userMoments[u].second));
				}

				auto overlapPair = Worker::CheckIntervalOverlap(parsedMoments);
				if(overlapPair.first != 0 || overlapPair.second != 0){
                    char tmp[128];
                    sprintf(tmp, "Overlapping time window detected! Check lines %d and %d.", overlapPair.first+1,overlapPair.second+1);
                    GLMessageBox::Display(tmp, "Error", GLDLG_OK, GLDLG_ICONERROR);
                    return;
                }

                moments.clear();
                for(auto& newMoment : parsedMoments)
                    AddMoment(newMoment);

                work->moments = moments;
				work->userMoments = userMoments;
				work->wp.timeWindowSize = window;
				work->wp.useMaxwellDistribution = useMaxwellToggle->GetState();
				work->wp.calcConstantFlow = calcConstantFlow->GetState();

				work->Reload();
				if (mApp->timeSettings) mApp->timeSettings->RefreshMoments();
				if (mApp->timewisePlotter) {
					mApp->timewisePlotter->Reset();
					mApp->timewisePlotter->refreshViews();
				}
			}
		}
		else if (src == clearButton) {
			if (GLDLG_OK == GLMessageBox::Display("Clear list?", "Moments", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONINFO)) {
				userMoments = std::vector<UserMoment>();
				RebuildList();
			}
		}
		else if (src == pasteButton) {
			PasteClipboard();
		}
		break;
	case MSG_TEXT:
	case MSG_LIST:
		for (size_t row = 0; row < (momentsList->GetNbRow() - 1); row++) {

		    // Change in time expression
            // Or change in windows size
            const char *momentEntry = momentsList->GetValueAt(1, row);
            const char *windowEntry = momentsList->GetValueAt(3, row);
            bool change = false;
            if (strcmp(momentEntry, userMoments[row].first.c_str()) != 0) {
				if (momentEntry != nullptr){//update
                    userMoments[row].first.assign(momentEntry);

                    // Try to use default window if unset
                    if(userMoments[row].second < DBL_EPSILON){
                        double window = 0.0;
                        if (!(windowSizeText->GetNumber(&window))) {
                            window = 0.0;
                        }
                        userMoments[row].second = window;
                        char tmp[10];
                        sprintf(tmp, "%g", window);
                        momentsList->SetValueAt(3, row, tmp);
                    }
                }
				else
					userMoments.erase(userMoments.begin() + row); //erase

				RebuildList();
				change = true;
			}
            else if (std::abs(std::strtod(windowEntry,nullptr) - userMoments[row].second) > DBL_EPSILON) {
                if (windowEntry != nullptr){//update
                    userMoments[row].second = std::strtod(windowEntry,nullptr);
                }
                else
                    userMoments.erase(userMoments.begin() + row); //erase

                RebuildList();
                change = true;
            }
            if(change){
                break;
            }
		}
		if (momentsList->GetValueAt(1, momentsList->GetNbRow() - 1) != nullptr) { //last line
			if (*(momentsList->GetValueAt(1, momentsList->GetNbRow() - 1)) != 0) {

			    double window = 0.0;
                if(windowSizeText != nullptr){
                    windowSizeText->GetNumber(&window);
                }

				//Add new line
				userMoments.emplace_back(momentsList->GetValueAt(1, momentsList->GetNbRow() - 1), window);
                RebuildList();
			}
		}
		break;
	}

	GLWindow::ProcessMessage(src, message);
}

/**
* \brief Rebuilds the moment list
*/
void MomentsEditor::RebuildList() {

    double defaultTimeWindow = 0.0;
    windowSizeText->GetNumber(&defaultTimeWindow);

    for (int u = 0; u < userMoments.size(); u++) {
        if(userMoments[u].first.empty() && (userMoments[u].second == 0.0 || userMoments[u].second == defaultTimeWindow)){
            userMoments.erase(userMoments.begin()+u);
            --u;
            continue;
        }
    }

	momentsList->SetSize(4, userMoments.size() + 1);
	momentsList->SetColumnWidths((int*)flWidth);
	momentsList->SetColumnLabels(flName);
	momentsList->SetColumnAligns((int *)flAligns);
	momentsList->SetColumnEditable((int *)fEdits);

	char tmp[128];
	size_t u; double latest = 0.0;

	for (u = 0; u < userMoments.size(); u++) {
		sprintf(tmp, "%zd", u + 1);
		momentsList->SetValueAt(0, u, tmp);
		sprintf(tmp, "%s", userMoments[u].first.c_str());
		momentsList->SetValueAt(1, u, tmp);
		sprintf(tmp, "%zd", (ParseMoment(userMoments[u].first, userMoments[u].second)).size());
		momentsList->SetValueAt(2, u, tmp);
        sprintf(tmp, "%g", userMoments[u].second);
        momentsList->SetValueAt(3, u, tmp);
	}
    //last line, possibility to enter new value
	sprintf(tmp, "%zd", u + 1);
	momentsList->SetValueAt(0, u, tmp);

}

void MomentsEditor::Refresh() {
	userMoments = work->userMoments;
	moments = work->moments;
	char tmp[128];
	sprintf(tmp, "%g", work->wp.timeWindowSize);
	windowSizeText->SetText(tmp);
	useMaxwellToggle->SetState(work->wp.useMaxwellDistribution);
	calcConstantFlow->SetState(work->wp.calcConstantFlow);
	RebuildList();
}

/**
* \brief Adds a time series to moments and returns the number of elements
* \param newMoments vector of new moments that should be inserted
* \return amount of new moments
*/
int MomentsEditor::AddMoment(std::vector<Moment> newMoments) {
    /*if(Worker::CheckIntervalOverlap(moments, newMoments)){
        return -1; // error
    }*/
    int nb = newMoments.size();
    moments.insert(moments.end(),newMoments.begin(),newMoments.end());
    std::sort(moments.begin(),moments.end());
    return nb;
}

/**
* \brief Parses a user input and returns a vector of time moments
* \param userInput string of the form "%lf,%lf,%lf" for beginning, interval step, ending of the moment series
* \return message Type of the source (button)
*/
std::vector<Moment> MomentsEditor::ParseMoment(const std::string& userInput, double timeWindow) {
	std::vector<Moment> parsedResult;
	double begin, interval, end;

	int nb = sscanf(userInput.c_str(), "%lf,%lf,%lf", &begin, &interval, &end);
	if (nb == 1 && (begin >= 0.0)) {
		//One moment
		parsedResult.emplace_back(begin,timeWindow);
		//} else if (nb==3 && (begin>0.0) && (end>begin) && (interval<(end-begin)) && ((end-begin)/interval<300.0)) {
	}
	else if (nb == 3 && (begin >= 0.0) && (end > begin) && (interval < (end - begin))) {
		//Range
		for (double time = begin; time <= end; time += interval)
			parsedResult.emplace_back(time,timeWindow);
	}
	return parsedResult;
}

/**
* \brief Pasting into the moments table from clipboard
*/
void MomentsEditor::PasteClipboard() {

	momentsList->PasteClipboardText(true, false);
}