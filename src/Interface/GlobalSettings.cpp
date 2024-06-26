

//#include "GLApp/GLLabel.h"
#include "GlobalSettings.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLToggle.h"
#include "GLApp/GLMessageBox.h"
#include "GLApp/GLList.h"
#include "GLApp/GLInputBox.h"
#include "Facet_shared.h"
#include "Geometry_shared.h"
#include "Helper/MathTools.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLTextField.h"
#include "GLApp/GLLabel.h"
#include "GLApp/GLTitledPanel.h"
#include "Buffer_shared.h"
#include "AppUpdater.h"
#include "SMP.h"
#include "ProcessControl.h"

#ifdef _WIN32

#else
//getpid in linux
#include <sys/types.h>
#include <unistd.h>
#endif

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

static const int   plWidth[] = { 60,40,70,70,295 };
static const char *plName[] = { "#","PID","Mem Usage","Mem Peak",/*"CPU",*/"Status" };
static const int   plAligns[] = { ALIGN_LEFT,ALIGN_LEFT,ALIGN_LEFT,ALIGN_LEFT,ALIGN_LEFT };

/**
* \brief Constructor for the global settings window with initial setup.
* \param worker Worker that handles the window.
*/
GlobalSettings::GlobalSettings(Worker *w) :GlobalSettingsBase(w) {

	worker = w;
	int windowWidth = 580;
	int windowHeight = 575;
	SetMinimumSize(windowWidth, windowHeight);
	SetResizable(true);

	const int checkboxHeight = 25;

	SetTitle("Global Settings");
	SetIconfiable(true);

	auto *settingsPanel = new GLTitledPanel("Program settings");
	settingsPanel->SetBounds(5, 2, 270, 267);
	Add(settingsPanel);

	auto *asLabel = new GLLabel("Autosave frequency (minutes):");
	asLabel->SetBounds(16, 22, 80, 19);
	settingsPanel->Add(asLabel);

	autoSaveText = new GLTextField(0, "");
	autoSaveText->SetBounds(170, 20, 30, 19);
	settingsPanel->Add(autoSaveText);

	int programSettingsPosY = 47;
	chkSimuOnly = new GLToggle(0, "Autosave only when simulation is running");
	chkSimuOnly->SetBounds(15, programSettingsPosY, 160, 19);
	settingsPanel->Add(chkSimuOnly);
    programSettingsPosY+=checkboxHeight;

	chkCompressSavedFiles = new GLToggle(0, "Use .zip as default extension (otherwise .xml)");
	chkCompressSavedFiles->SetBounds(15, programSettingsPosY, 100, 19);
	settingsPanel->Add(chkCompressSavedFiles);
    programSettingsPosY+=checkboxHeight;

	chkCheckForUpdates = new GLToggle(0, "Check for updates at startup");
	chkCheckForUpdates->SetBounds(15, programSettingsPosY, 160, 19);
	settingsPanel->Add(chkCheckForUpdates);
    programSettingsPosY+=checkboxHeight;

	chkAntiAliasing = new GLToggle(0, "Anti-Aliasing");
	chkAntiAliasing->SetBounds(15, programSettingsPosY, 160, 19);
	settingsPanel->Add(chkAntiAliasing);
    programSettingsPosY+=checkboxHeight;

	chkWhiteBg = new GLToggle(0, "White Background");
	chkWhiteBg->SetBounds(15, programSettingsPosY, 160, 19);
	settingsPanel->Add(chkWhiteBg);
    programSettingsPosY+=checkboxHeight;

	leftHandedToggle = new GLToggle(0, "Left-handed coord. system");
	leftHandedToggle->SetBounds(15, programSettingsPosY, 160, 19);
	settingsPanel->Add( leftHandedToggle );
    programSettingsPosY+=checkboxHeight;

	highlightNonplanarToggle = new GLToggle(0, "Highlight non-planar facets");
	highlightNonplanarToggle->SetBounds(15, programSettingsPosY, 160, 19);
	settingsPanel->Add(highlightNonplanarToggle);
    programSettingsPosY+=checkboxHeight;

    highlightSelectionToggle = new GLToggle(0, "Highlight selected facets");
    highlightSelectionToggle->SetBounds(15, programSettingsPosY, 160, 19);
    settingsPanel->Add(highlightSelectionToggle);
    programSettingsPosY+=checkboxHeight;

    useOldXMLFormat = new GLToggle(0, "Use old XML format");
    useOldXMLFormat->SetBounds(15, programSettingsPosY, 160, 19);
    settingsPanel->Add(useOldXMLFormat);


	auto *simuSettingsPanel = new GLTitledPanel("Simulation settings (current file)");
	simuSettingsPanel->SetBounds(280, 2, 290, 267);
	Add(simuSettingsPanel);

	auto *massLabel = new GLLabel("Gas molecular mass (g/mol):");
	massLabel->SetBounds(290, 22, 150, 19);
	simuSettingsPanel->Add(massLabel);

	gasMassText = new GLTextField(0, "");
	gasMassText->SetBounds(460, 20, 100, 19);
	simuSettingsPanel->Add(gasMassText);

	enableDecay = new GLToggle(0, "Gas half life (s):");
	enableDecay->SetBounds(290, 47, 150, 19);
	simuSettingsPanel->Add(enableDecay);

	halfLifeText = new GLTextField(0, "");
	halfLifeText->SetBounds(460, 45, 100, 19);
	simuSettingsPanel->Add(halfLifeText);
	
	auto *outgassingLabel = new GLLabel("Final outgassing rate (mbar*l/sec):");
	outgassingLabel->SetBounds(290, 72, 150, 19);
	simuSettingsPanel->Add(outgassingLabel);

	outgassingGasRateText = new GLTextField(0, "");
	outgassingGasRateText->SetBounds(460, 70, 100, 19);
	outgassingGasRateText->SetEditable(false);
	simuSettingsPanel->Add(outgassingGasRateText);

	auto *outgassingLabel2 = new GLLabel("Final outgassing rate (1/sec):");
	outgassingLabel2->SetBounds(290, 97, 150, 19);
	simuSettingsPanel->Add(outgassingLabel2);

	outgassingMoleculeRateText = new GLTextField(0, "");
	outgassingMoleculeRateText->SetBounds(460, 95, 100, 19);
	outgassingMoleculeRateText->SetEditable(false);
	simuSettingsPanel->Add(outgassingMoleculeRateText);

	desorbedMoleculesLabel = new GLLabel("Total desorbed molecules:");
	desorbedMoleculesLabel->SetBounds(290, 122, 150, 19);
	simuSettingsPanel->Add(desorbedMoleculesLabel);

	desorbedMoleculesText = new GLTextField(0, "");
	desorbedMoleculesText->SetBounds(460, 120, 100, 19);
	desorbedMoleculesText->SetEditable(false);
	simuSettingsPanel->Add(desorbedMoleculesText);

	recalcButton = new GLButton(0, "Recalc. outgassing");
	recalcButton->SetBounds(460, 148, 100, 19);
	simuSettingsPanel->Add(recalcButton);

	lowFluxToggle = new GLToggle(0, "Enable low flux mode");
	lowFluxToggle->SetBounds(290, 175, 120, 19);
	simuSettingsPanel->Add(lowFluxToggle);

	lowFluxInfo = new GLButton(0, "?");
	lowFluxInfo->SetBounds(420, 175, 20, 19);
	simuSettingsPanel->Add(lowFluxInfo);

	auto *cutoffLabel = new GLLabel("Cutoff ratio:");
	cutoffLabel->SetBounds(310, 201, 80, 19);
	simuSettingsPanel->Add(cutoffLabel);

	cutoffText = new GLTextField(0, "");
	cutoffText->SetBounds(370, 200, 70, 19);
	cutoffText->SetEditable(false);
	simuSettingsPanel->Add(cutoffText);

	applyButton = new GLButton(0, "Apply above settings");
	applyButton->SetBounds(windowWidth / 2 - 65, 273, 130, 19);
	Add(applyButton);

	processPanel = new GLTitledPanel("Process control");
	Add(processPanel);

	processList = new GLList(0);
	processList->SetSize(5, MAX_PROCESS+2);
	processList->SetColumnWidths((int*)plWidth);
	processList->SetColumnLabels((const char **)plName);
	processList->SetColumnAligns((int *)plAligns);
	processList->SetColumnLabelVisible(true);
	processList->SetHScrollVisible(false);
	processPanel->Add(processList);

	char tmp[128];
	sprintf(tmp, "Number of CPU cores:     %zd", mApp->numCPU);
	coreLabel = new GLLabel(tmp);
	processPanel->Add(coreLabel);

	subProcLabel = new GLLabel("Number of subprocesses:");
	processPanel->Add(subProcLabel);

	nbProcText = new GLTextField(0, "");
	nbProcText->SetEditable(true);
	processPanel->Add(nbProcText);

	restartButton = new GLButton(0, "Apply and restart processes");
	processPanel->Add(restartButton);

	desLimitButton = new GLButton(0, "Desorption limit");
	processPanel->Add(desLimitButton);
	UpdateDesLimitButtonText();

	// Center dialog
	int screenWidth, screenHeight;
	GLToolkit::GetScreenSize(&screenWidth, &screenHeight);
	int topLeftX = (screenWidth - windowWidth) / 2;
	int topLeftY = (screenHeight - windowHeight) / 2;
	SetBounds(topLeftX, topLeftY, windowWidth, windowHeight);

    GLContainer::RestoreDeviceObjects();

	lastUpdate = 0;
}

/**
* \brief Function to change global settings values (may happen due to change or load of file etc.)
*/
void GlobalSettings::Update() {
	Update_shared();
	UpdateOutgassing();
	gasMassText->SetText(worker->model->sp.gasMass);
	enableDecay->SetState(worker->model->sp.enableDecay);
	halfLifeText->SetText(worker->model->sp.halfLife);
	halfLifeText->SetEditable(worker->model->sp.enableDecay);
}

/**
* \brief Function for processing various inputs (button, check boxes etc.)
* \param src Exact source of the call
* \param message Type of the source (button)
*/
void GlobalSettings::ProcessMessage(GLComponent *src, int message) {
	
	ProcessMessage_shared(src, message); //Common Molflow-Synrad elements

	//Molflow-specific
	switch (message) {
	case MSG_BUTTON:

		if (src == recalcButton) {
			if (mApp->AskToReset()) {
				try {
					worker->ReloadIfNeeded();
					UpdateOutgassing();
				}
				catch (const std::exception &e) {
					GLMessageBox::Display(e.what(), "Recalculation failed: Couldn't convert model to simulation", GLDLG_OK, GLDLG_ICONWARNING);
				}
			}		
		}
		else if (src == applyButton) {
			mApp->useOldXMLFormat = useOldXMLFormat->GetState();
            double gm;
			if (!gasMassText->GetNumber(&gm) || gm <= 0.0) {
				GLMessageBox::Display("Invalid gas mass", "Error", GLDLG_OK, GLDLG_ICONERROR);
				return;
			}
			if (std::abs(gm - worker->model->sp.gasMass) > 1e-7) {
				if (mApp->AskToReset()) {
					worker->needsReload = true;
					worker->model->sp.gasMass = gm;
					if (worker->GetGeometry()->IsLoaded()) { //check if there are pumps
						bool hasPump = false;
						size_t nbFacet = worker->GetGeometry()->GetNbFacet();
						for (size_t i = 0; (i<nbFacet) && (!hasPump); i++) {
							if (worker->GetGeometry()->GetFacet(i)->sh.sticking>0.0) {
								hasPump = true;
							}
						}
						if (hasPump) GLMessageBox::Display("Don't forget the pumps: update pumping speeds and/or recalculate sticking factors.", "You have changed the gas mass.", GLDLG_OK, GLDLG_ICONINFO);
					}
				}
			}

			double hl;
			if (enableDecay->GetState() && (!halfLifeText->GetNumber(&hl) || hl <= 0.0)) {
				GLMessageBox::Display("Invalid half life", "Error", GLDLG_OK, GLDLG_ICONERROR);
				return;
			}
			if ((enableDecay->GetState()==1) != worker->model->sp.enableDecay || ((enableDecay->GetState()==1) && IsEqual(hl, worker->model->sp.halfLife))) {
				if (mApp->AskToReset()) {
					worker->needsReload = true;
					worker->model->sp.enableDecay = enableDecay->GetState();
					if (worker->model->sp.enableDecay) worker->model->sp.halfLife = hl;
				}
			}
			return;
		}
		else if (src == lowFluxInfo) {
			GLMessageBox::Display("Low flux mode helps to gain more statistics on low pressure parts of the system, at the expense\n"
				"of higher pressure parts. If a traced particle reflects from a high sticking factor surface, regardless of that probability,\n"
				"a reflected test particle representing a reduced flux will still be traced. Therefore test particles can reach low flux areas more easily, but\n"
				"at the same time tracing a test particle takes longer. The cutoff ratio defines what ratio of the originally generated flux\n"
				"can be neglected. If, for example, it is 0.001, then, when after subsequent reflections the test particle carries less than 0.1%\n"
				"of the original flux, it will be eliminated. A good advice is that if you'd like to see pressure across N orders of magnitude, set it to 1E-N"
				, "Low flux mode", GLDLG_OK, GLDLG_ICONINFO);
			return;
		}
		break;

	case MSG_TEXT:
		ProcessMessage(applyButton, MSG_BUTTON);
		break;

	case MSG_TOGGLE:
		if (src == enableDecay) {
			halfLifeText->SetEditable(enableDecay->GetState());
		}
		break;

    default:
        break;
	}

	GLWindow::ProcessMessage(src, message);
}

/**
* \brief Updates the values regarding outgassing (from memory and by calculation)
*/
void GlobalSettings::UpdateOutgassing() {
	if (worker->needsReload) {
		std::string changedTxt("Model changed");
		outgassingGasRateText->SetText(changedTxt);
		outgassingMoleculeRateText->SetText(changedTxt);
		desorbedMoleculesLabel->SetText("Tot.des.molecules["+changedTxt+"]");
		desorbedMoleculesText->SetText(changedTxt);
	}
	else {
		char tmp[128];
		sprintf(tmp, "%g", worker->model->sp.finalOutgassingRate_Pa_m3_sec * PAM3S_TO_MBARLS);
		outgassingGasRateText->SetText(tmp);
		sprintf(tmp, "%g", worker->model->sp.finalOutgassingRate); //In molecules/sec
		outgassingMoleculeRateText->SetText(tmp);
		sprintf(tmp, "Tot.des. molecules [0 to %g s]:", worker->model->sp.latestMoment);
		desorbedMoleculesLabel->SetText(tmp);
		sprintf(tmp, "%.3E", worker->model->sp.totalDesorbedMolecules);
		desorbedMoleculesText->SetText(tmp);
	}
}

void GlobalSettings::ResizeProcessPanel(int windowWidth, int windowHeight) {

	int processPanelHeight = windowHeight - 332;
	processPanel->SetBounds(7, 305, windowWidth - 15, processPanelHeight);
	processPanel->SetCompBounds(processList, 4, 18, windowWidth - 20, processPanelHeight - 73);

	processPanel->SetCompBounds(coreLabel,10, processPanelHeight - 50, 120, 19);
	processPanel->SetCompBounds(subProcLabel,10, processPanelHeight - 25, 120, 19);
	processPanel->SetCompBounds(nbProcText,135, processPanelHeight - 27, 30, 19);
	processPanel->SetCompBounds(restartButton,170, processPanelHeight - 27, 150, 19);
	processPanel->SetCompBounds(desLimitButton,windowWidth - 195, processPanelHeight - 27, 175, 19);

	processList->SetColumnWidth(4, windowWidth - 287);
}

void GlobalSettings::SetBounds(int x, int y, int width, int height) {
	ResizeProcessPanel(width, height);
	GLWindow::SetBounds(x, y, width, height);
}