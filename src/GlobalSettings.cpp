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

#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)

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

static const int   plWidth[] = { 60,40,70,70,335 };
static const char *plName[] = { "#","PID","Mem Usage","Mem Peak",/*"CPU",*/"Status" };
static const int   plAligns[] = { ALIGN_LEFT,ALIGN_LEFT,ALIGN_LEFT,ALIGN_LEFT,ALIGN_LEFT };

/**
* \brief Constructor for the global settings window with initial setup.
* \param worker Worker that handles the window.
*/
GlobalSettings::GlobalSettings(Worker *w) :GLWindow() {

	worker = w;
	int wD = 580;
	int hD = 565;

	const int checkboxHeight = 25;

	SetTitle("Global Settings");
	SetIconfiable(true);

	auto *settingsPanel = new GLTitledPanel("Program settings");
	settingsPanel->SetBounds(5, 2, 270, 292);
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

	chkAutoUpdateFormulas = new GLToggle(0, "Auto refresh formulas");
	chkAutoUpdateFormulas->SetBounds(15, programSettingsPosY, 160, 19);
	settingsPanel->Add(chkAutoUpdateFormulas);
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


	auto *simuSettingsPanel = new GLTitledPanel("Simulation settings");
	simuSettingsPanel->SetBounds(280, 2, 290, 292);
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
	applyButton->SetBounds(wD / 2 - 65, 298, 130, 19);
	Add(applyButton);

	/*chkNonIsothermal = new GLToggle(0,"Non-isothermal system (textures only, experimental)");
	chkNonIsothermal->SetBounds(315,125,100,19);
	Add(chkNonIsothermal);*/

	const int procControlPosY= 334;
	auto *panel3 = new GLTitledPanel("Process control");
	panel3->SetBounds(5, procControlPosY, wD - 10, hD - 310);
	Add(panel3);

	processList = new GLList(0);
	processList->SetHScrollVisible(true);
	processList->SetSize(5, MAX_PROCESS+1);
	processList->SetColumnWidths((int*)plWidth);
	processList->SetColumnLabels((const char **)plName);
	processList->SetColumnAligns((int *)plAligns);
	processList->SetColumnLabelVisible(true);
	processList->SetBounds(10, procControlPosY+20, wD - 20, hD - 433);
	panel3->Add(processList);

	char tmp[128];
	sprintf(tmp, "Number of CPU cores:     %zd", mApp->numCPU);
	auto *coreLabel = new GLLabel(tmp);
	coreLabel->SetBounds(10, hD - 74, 120, 19);
	panel3->Add(coreLabel);

	auto *l1 = new GLLabel("Number of subprocesses:");
	l1->SetBounds(10, hD - 49, 120, 19);
	panel3->Add(l1);

	nbProcText = new GLTextField(0, "");
	nbProcText->SetEditable(true);
	nbProcText->SetBounds(135, hD - 51, 30, 19);
	panel3->Add(nbProcText);

	restartButton = new GLButton(0, "Apply and restart processes");
	restartButton->SetBounds(170, hD - 51, 150, 19);
	panel3->Add(restartButton);

	maxButton = new GLButton(0, "Change MAX desorbed molecules");
	maxButton->SetBounds(wD - 195, hD - 51, 180, 19);
	panel3->Add(maxButton);

	

	// Center dialog
	int wS, hS;
	GLToolkit::GetScreenSize(&wS, &hS);
	int xD = (wS - wD) / 2;
	int yD = (hS - hD) / 2;
	GLWindow::SetBounds(xD, yD, wD, hD);

    GLContainer::RestoreDeviceObjects();

	lastUpdate = 0;
	//for (size_t i = 0; i < MAX_PROCESS; i++) lastCPUTime[i] = -1.0f;
	//memset(lastCPULoad, 0, MAX_PROCESS*sizeof(float));
}

/**
* \brief Function to change global settings values (may happen due to change or load of file etc.)
*/
void GlobalSettings::Update() {

	char tmp[256];
    useOldXMLFormat->SetState(mApp->useOldXMLFormat);
    chkAntiAliasing->SetState(mApp->antiAliasing);
	chkWhiteBg->SetState(mApp->whiteBg);
    highlightNonplanarToggle->SetState(mApp->highlightNonplanarFacets);
    highlightSelectionToggle->SetState(mApp->highlightSelection);
    leftHandedToggle->SetState(mApp->leftHandedView);
	//chkNonIsothermal->SetState(nonIsothermal);
	UpdateOutgassing();

	gasMassText->SetText(worker->wp.gasMass);

	enableDecay->SetState(worker->wp.enableDecay);
	halfLifeText->SetText(worker->wp.halfLife);
	halfLifeText->SetEditable(worker->wp.enableDecay);

	cutoffText->SetText(worker->ontheflyParams.lowFluxCutoff);
	cutoffText->SetEditable(worker->ontheflyParams.lowFluxMode);
	lowFluxToggle->SetState(worker->ontheflyParams.lowFluxMode);

	autoSaveText->SetText(mApp->autoSaveFrequency);
	chkSimuOnly->SetState(mApp->autoSaveSimuOnly);
	if (mApp->appUpdater) { //Updater initialized
		chkCheckForUpdates->SetState(mApp->appUpdater->IsUpdateCheckAllowed());
	}
	else {
		chkCheckForUpdates->SetState(0);
		chkCheckForUpdates->SetEnabled(false);
	}
	chkAutoUpdateFormulas->SetState(mApp->autoUpdateFormulas);
	chkCompressSavedFiles->SetState(mApp->compressSavedFiles);

	size_t nb = worker->GetProcNumber();
	sprintf(tmp, "%zd", nb);
	nbProcText->SetText(tmp);
}

/**
* \brief Function to update the thread information table in the global settings window.
*/
void GlobalSettings::SMPUpdate() {
	int time = SDL_GetTicks();

	if (!IsVisible() || IsIconic()) return;
	size_t nb = worker->GetProcNumber();
	if (processList->GetNbRow() != (nb + 1)) processList->SetSize(5, nb + 1,true);

	if( time-lastUpdate>333 ) {

	char tmp[512];
	size_t  states[MAX_PROCESS];
	std::vector<std::string> statusStrings(MAX_PROCESS);

	memset(states, 0, MAX_PROCESS * sizeof(int));
	worker->GetProcStatus(states, statusStrings);

    std::vector<SubProcInfo> procInfo;
    worker->GetProcStatus(procInfo);

    processList->ResetValues();

	//Interface
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
	size_t currPid = GetCurrentProcessId();
    PROCESS_INFO parentInfo;
    GetProcInfo(currPid, &parentInfo);

    processList->SetValueAt(0, 0, "Interface");
	sprintf(tmp, "%zd", currPid);
	processList->SetValueAt(1, 0, tmp, currPid);
	sprintf(tmp, "%.0f MB", (double)parentInfo.mem_use / (1024.0*1024.0));
	processList->SetValueAt(2, 0, tmp);
	sprintf(tmp, "%.0f MB", (double)parentInfo.mem_peak / (1024.0*1024.0));
	processList->SetValueAt(3, 0, tmp);

#else
    size_t currPid = getpid();
    PROCESS_INFO parentInfo;
    GetProcInfo(currPid, &parentInfo);

    //GetProcInfo(currpid, &parentInfo);
    processList->SetValueAt(0, 0, "Interface");
    sprintf(tmp, "%zd", currPid);
    processList->SetValueAt(1, 0, tmp, (int)currPid);

    sprintf(tmp, "%.0f MB", (double)parentInfo.mem_use / (1024.0));
    processList->SetValueAt(2, 0, tmp);
    sprintf(tmp, "%.0f MB", (double)parentInfo.mem_peak / (1024.0));
    processList->SetValueAt(3, 0, tmp);
    //sprintf(tmp, "%d %%", (int)parentInfo.cpu_time);
    //processList->SetValueAt(4, 0, tmp);
#endif

    size_t i = 1;
	for (auto& proc : procInfo) {
		DWORD pid = proc.procId;
		sprintf(tmp, "Subproc.%lu", i);
		processList->SetValueAt(0, i, tmp);
		sprintf(tmp, "%d", pid);
		processList->SetValueAt(1, i, tmp);

#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
        PROCESS_INFO pInfo;
        if (!GetProcInfo(pid, &pInfo)) {
			processList->SetValueAt(2, i, "0 KB");
			processList->SetValueAt(3, i, "0 KB");
			//processList->SetValueAt(4,i,"0 %");
			processList->SetValueAt(4, i, "Dead");
		}
		else {
			sprintf(tmp, "%.0f MB", (double)pInfo.mem_use / (1024.0*1024.0));
			processList->SetValueAt(2, i, tmp);
			sprintf(tmp, "%.0f MB", (double)pInfo.mem_peak / (1024.0*1024.0));
			processList->SetValueAt(3, i, tmp);

			// State/Status
			std::stringstream tmp; tmp << "[" << prStates[states[i-1]] << "] " << statusStrings[i-1];
			processList->SetValueAt(4, i, tmp.str().c_str());
		}

#else
        if (pid == currPid) { // TODO: Check if this is wanted
            processList->SetValueAt(2, i, "0 KB");
            processList->SetValueAt(3, i, "0 KB");
            //processList->SetValueAt(4,i,"0 %");
            processList->SetValueAt(4, i, "Dead");
        }
        else {
            PROCESS_INFO pInfo = proc.runtimeInfo;
            //GetProcInfo(pid, &pInfo);


            sprintf(tmp, "%.0f MB", (double)pInfo.mem_use / (1024.0));
            processList->SetValueAt(2, i, tmp);
            sprintf(tmp, "%.0f MB", (double)pInfo.mem_peak / (1024.0));
            processList->SetValueAt(3, i, tmp);
            //sprintf(tmp, "%d %%", (int)pInfo.cpu_time);
            //processList->SetValueAt(4, i, tmp);

            // State/Status
            std::stringstream tmp_ss; tmp_ss << "[" << prStates[states[i-1]] << "] " << statusStrings[i-1];
            processList->SetValueAt(4, i, tmp_ss.str().c_str());
        }
#endif
++i;
	}
	lastUpdate = SDL_GetTicks();
	}

}

/**
* \brief Function to apply changes to the number of processes.
*/
void GlobalSettings::RestartProc() {

	int nbProc;
	if (!nbProcText->GetNumberInt(&nbProc)) {
		GLMessageBox::Display("Invalid process number", "Error", GLDLG_OK, GLDLG_ICONERROR);
	}
	else {
			if (nbProc <= 0 || nbProc > MAX_PROCESS) {
				GLMessageBox::Display("Invalid process number [1..64]", "Error", GLDLG_OK, GLDLG_ICONERROR);
			}
			else {
				try {
					worker->Stop_Public();
                    worker->SetProcNumber(nbProc);
					worker->RealReload(true);
					mApp->SaveConfig();
				}
				catch (Error &e) {
					GLMessageBox::Display(e.what(), "Error", GLDLG_OK, GLDLG_ICONERROR);
				}
			}
	}

}

/**
* \brief Function for processing various inputs (button, check boxes etc.)
* \param src Exact source of the call
* \param message Type of the source (button)
*/
void GlobalSettings::ProcessMessage(GLComponent *src, int message) {

	switch (message) {
	case MSG_BUTTON:

		if (src == recalcButton) {
			if (mApp->AskToReset()) {
				try {
					worker->RealReload();
				}
				catch (std::exception &e) {
					GLMessageBox::Display(e.what(), "Recalculation failed: Couldn't reload Worker", GLDLG_OK, GLDLG_ICONWARNING);
				}
			}
		}
		else if (src == restartButton) {
			RestartProc();
		}
		else if (src == maxButton) {
			if (worker->GetGeometry()->IsLoaded()) {
				char tmp[128];
				sprintf(tmp, "%zd", worker->ontheflyParams.desorptionLimit);
				char *val = GLInputBox::GetInput(tmp, "Desorption max (0=>endless)", "Edit MAX");
				if (val) {
                    char* endptr;
                    size_t maxDes = strtold(val,&endptr); // use double function to allow exponential format
					if (val==endptr) {
						GLMessageBox::Display("Invalid 'maximum desorption' number", "Error", GLDLG_OK, GLDLG_ICONERROR);
					}
					else {
                        worker->ontheflyParams.desorptionLimit = maxDes;
                        worker->ChangeSimuParams(); //Sync with subprocesses
                    }
				}
			}
			else {
				GLMessageBox::Display("No geometry loaded.", "No geometry", GLDLG_OK, GLDLG_ICONERROR);
			}
		}
		else if (src == applyButton) {
            mApp->useOldXMLFormat = useOldXMLFormat->GetState();
            mApp->antiAliasing = chkAntiAliasing->GetState();
			mApp->whiteBg = chkWhiteBg->GetState();
            mApp->highlightSelection = highlightSelectionToggle->GetState();
            mApp->highlightNonplanarFacets = highlightNonplanarToggle->GetState();
            mApp->leftHandedView = (bool)leftHandedToggle->GetState();
			for (auto & i : mApp->viewer) {
				i->UpdateMatrix();
				i->UpdateLabelColors();
			}
			mApp->wereEvents = true;
			bool updateCheckPreference = chkCheckForUpdates->GetState();
			if (mApp->appUpdater) {
				if (mApp->appUpdater->IsUpdateCheckAllowed() != updateCheckPreference) {
					mApp->appUpdater->SetUserUpdatePreference(updateCheckPreference);
				}
			}
			mApp->autoUpdateFormulas = chkAutoUpdateFormulas->GetState();
			mApp->compressSavedFiles = chkCompressSavedFiles->GetState();
			mApp->autoSaveSimuOnly = chkSimuOnly->GetState();
			double gm;
			if (!gasMassText->GetNumber(&gm) || gm <= 0.0) {
				GLMessageBox::Display("Invalid gas mass", "Error", GLDLG_OK, GLDLG_ICONERROR);
				return;
			}
			if (std::abs(gm - worker->wp.gasMass) > 1e-7) {
				if (mApp->AskToReset()) {
					worker->needsReload = true;
					worker->wp.gasMass = gm;
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
			if ((enableDecay->GetState()==1) != worker->wp.enableDecay || ((enableDecay->GetState()==1) && IsEqual(hl, worker->wp.halfLife))) {
				if (mApp->AskToReset()) {
					worker->needsReload = true;
					worker->wp.enableDecay = enableDecay->GetState();
					if (worker->wp.enableDecay) worker->wp.halfLife = hl;
				}
			}

			double cutoffnumber;
			if (!cutoffText->GetNumber(&cutoffnumber) || !(cutoffnumber>0.0 && cutoffnumber<1.0)) {
				GLMessageBox::Display("Invalid cutoff ratio, must be between 0 and 1", "Error", GLDLG_OK, GLDLG_ICONWARNING);
				return;
			}

			if (!IsEqual(worker->ontheflyParams.lowFluxCutoff, cutoffnumber) || (int)worker->ontheflyParams.lowFluxMode != lowFluxToggle->GetState()) {
				worker->ontheflyParams.lowFluxCutoff = cutoffnumber;
				worker->ontheflyParams.lowFluxMode = lowFluxToggle->GetState();
				worker->ChangeSimuParams();
			}

			double autosavefreq;
			if (!autoSaveText->GetNumber(&autosavefreq) || autosavefreq <= 0.0) {
				GLMessageBox::Display("Invalid autosave frequency", "Error", GLDLG_OK, GLDLG_ICONERROR);
				return;
			}
			mApp->autoSaveFrequency = autosavefreq;

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
		} else if (src == lowFluxToggle) {
			cutoffText->SetEditable(lowFluxToggle->GetState());
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
	char tmp[128];
	sprintf(tmp, "%g", worker->wp.gasMass);
	gasMassText->SetText(tmp);
	sprintf(tmp, "%g", worker->wp.finalOutgassingRate_Pa_m3_sec * 10.00); //10: conversion Pa*m3/sec -> mbar*l/s
	outgassingGasRateText->SetText(tmp);
	sprintf(tmp, "%g", worker->wp.finalOutgassingRate); //In molecules/sec
	outgassingMoleculeRateText->SetText(tmp);
	sprintf(tmp,"Tot.des. molecules [0 to %g s]:",worker->wp.latestMoment);
	desorbedMoleculesLabel->SetText(tmp);
	sprintf(tmp, "%.3E", worker->wp.totalDesorbedMolecules);
	desorbedMoleculesText->SetText(tmp);
}