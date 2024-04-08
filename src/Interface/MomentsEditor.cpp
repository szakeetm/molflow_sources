

#include <cfloat> // DBL_EPSILON

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
#include "TimeMoments.h"
#include "HistogramSettings.h"

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

	//sprintf(tmp,"%g",work->model->sp.sp.timeWindowSize);
	windowSizeText = new GLTextField(0, "");
	windowSizeText->SetBounds(120, 275, 60, 20);
	Add(windowSizeText);

	useMaxwellToggle = new GLToggle(0, "Use Maxwell-B. speed distr.");
	useMaxwellToggle->SetBounds(15, 300, wD - 25, 20);
	//useMaxwellToggle->SetThreadStates(work->model->sp.useMaxwellDistribution);
	Add(useMaxwellToggle);

	calcConstantFlow = new GLToggle(0, "Calculate constant flow");
	calcConstantFlow->SetBounds(15, 325, wD - 25, 20);
	//useMaxwellToggle->SetThreadStates(work->model->sp.useMaxwellDistribution);
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
					parsedMoments.emplace_back(ParseUserMoment(userMoments[u]));
				}

				auto overlapPair = TimeMoments::HasIntervalOverlap(parsedMoments);
				if(overlapPair.has_value()){
                    char tmp[128];
                    if(overlapPair->second < 0)
                        sprintf(tmp, "Interval length and time window would create overlap! Check line %d.", overlapPair->first+1);
                    else
                        sprintf(tmp, "Overlapping time window detected! Check lines %d and %d.", overlapPair->first+1,overlapPair->second+1);
                    GLMessageBox::Display(tmp, "Error", GLDLG_OK, GLDLG_ICONERROR);
                    return;
                }

                moments.clear();
                for(auto& newMoment : parsedMoments)
                    AddMoment(newMoment);

                work->interfaceMomentCache = moments;
				work->userMoments = userMoments;
				work->model->sp.timeWindowSize = window;
				work->model->sp.useMaxwellDistribution = useMaxwellToggle->GetState();
				work->model->sp.calcConstantFlow = calcConstantFlow->GetState();

				work->MarkToReload();
				if (mApp->timeSettings) mApp->timeSettings->RefreshMoments();
				if (mApp->timewisePlotter) {
					mApp->timewisePlotter->Reset();
					mApp->timewisePlotter->refreshViews();
				}
				if (mApp->histogramSettings && mApp->histogramSettings->IsVisible()) {
					mApp->histogramSettings->UpdateMemoryEstimate_Current();
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
            if (strcmp(momentEntry, userMoments[row].content.c_str()) != 0) {
				if (momentEntry != nullptr){//update
                    userMoments[row].content.assign(momentEntry);

                    // Try to use default window if unset
                    if(userMoments[row].timeWindow < DBL_EPSILON){
                        double window = 0.0;
                        if (!(windowSizeText->GetNumber(&window))) {
                            window = 0.0;
                        }
                        userMoments[row].timeWindow = window;
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
            else if (std::abs(std::strtod(windowEntry,nullptr) - userMoments[row].timeWindow) > DBL_EPSILON) {
                if (windowEntry != nullptr){//update
                    userMoments[row].timeWindow = std::strtod(windowEntry,nullptr);
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
				UserMoment um;
				um.content=momentsList->GetValueAt(1, momentsList->GetNbRow() - 1);
				um.timeWindow=window;
				userMoments.emplace_back(um);
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
        if(userMoments[u].content.empty() && (userMoments[u].timeWindow == 0.0 || userMoments[u].timeWindow == defaultTimeWindow)){
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
		sprintf(tmp, "%s", userMoments[u].content.c_str());
		momentsList->SetValueAt(1, u, tmp);
		sprintf(tmp, "%zd", (ParseUserMoment(userMoments[u]).size()));
		momentsList->SetValueAt(2, u, tmp);
        sprintf(tmp, "%g", userMoments[u].timeWindow);
        momentsList->SetValueAt(3, u, tmp);
	}
    //last line, possibility to enter new value
	sprintf(tmp, "%zd", u + 1);
	momentsList->SetValueAt(0, u, tmp);

}

void MomentsEditor::Refresh() {
	userMoments = work->userMoments;
	moments = work->interfaceMomentCache;
	char tmp[128];
	sprintf(tmp, "%g", work->model->sp.timeWindowSize);
	windowSizeText->SetText(tmp);
	useMaxwellToggle->SetState(work->model->sp.useMaxwellDistribution);
	calcConstantFlow->SetState(work->model->sp.calcConstantFlow);
	RebuildList();
}

/**
* \brief Adds a time series to moments and returns the number of elements
* \param newMoments vector of new moments that should be inserted
* \return amount of new moments
*/
int MomentsEditor::AddMoment(std::vector<Moment> newMoments) {

    int nb = (int)newMoments.size();
    moments.insert(moments.end(),newMoments.begin(),newMoments.end());
    	std::sort(moments.begin(), moments.end(), [](const Moment& a, const Moment& b) {
		return a.time - .5 * a.window < b.time - .5 * b.window; // This will sort in ascending order based on the 'startTime' member
		});
    return nb;
}

/**
* \brief Parses a user input and returns a vector of time moments
* \return message Type of the source (button)
*/
std::vector<Moment> MomentsEditor::ParseUserMoment(const UserMoment& input) {
	std::vector<Moment> parsedResult;
	double begin, interval, end;

	int nb = sscanf(input.content.c_str(), "%lf,%lf,%lf", &begin, &interval, &end);
	if (nb == 1 && (begin >= 0.0)) {
		//One moment
		Moment m;
		m.time=begin;
		m.window=input.timeWindow;
		parsedResult.emplace_back(m);
	}
	else if (nb == 3 && (begin >= 0.0) && (end > begin) && (interval < (end - begin))) {
		//Range
		// First check for potential overlap due to interval<timeWindow
		if(!(interval<input.timeWindow)){
            for (double time = begin; time <= end; time += interval){
				Moment m;
				m.time=time;
				m.window=input.timeWindow;
				parsedResult.emplace_back(m);	
			}
        }
	}
	return parsedResult;
}

/**
* \brief Pasting into the moments table from clipboard
*/
void MomentsEditor::PasteClipboard() {

	momentsList->PasteClipboardText(true, false);
}