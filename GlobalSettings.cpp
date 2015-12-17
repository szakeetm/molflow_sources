/*
  File:        GlobalSettings.cpp
  Description: Global settings dialog (antiAliasing,units)
  Program:     MolFlow

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  */

#include "GLApp/GLLabel.h"
#include "GlobalSettings.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLMessageBox.h"
#include "GLApp/GLInputBox.h"
#include "Utils.h"
#include "Molflow.h"

extern MolFlow *mApp;

static const int   plWidth[] = { 15, 40, 70, 70, 50, 330 };
static const char *plName[] = { "#", "PID", "Mem Usage", "Mem Peak", "CPU", "Status" };
static const int   plAligns[] = { ALIGN_LEFT, ALIGN_CENTER, ALIGN_CENTER, ALIGN_CENTER, ALIGN_CENTER, ALIGN_LEFT };

// --------------------------------------------------------------------

GlobalSettings::GlobalSettings(Worker *w) :GLWindow() {

	worker = w;
	int wD = 580;
	int hD = 500;

	SetTitle("Global Settings");
	SetIconfiable(TRUE);

	GLTitledPanel *settingsPanel = new GLTitledPanel("Program settings");
	settingsPanel->SetBounds(5, 2, 270, 203);
	Add(settingsPanel);

	GLLabel *asLabel = new GLLabel("Autosave frequency (minutes):");
	asLabel->SetBounds(16, 25, 80, 19);
	settingsPanel->Add(asLabel);

	autoSaveText = new GLTextField(0, "");
	autoSaveText->SetBounds(170, 20, 30, 19);
	settingsPanel->Add(autoSaveText);

	chkSimuOnly = new GLToggle(0, "Autosave only when simulation is running");
	chkSimuOnly->SetBounds(15, 50, 160, 19);
	settingsPanel->Add(chkSimuOnly);

	chkCompressSavedFiles = new GLToggle(0, "Use .zip as default extension (otherwise .xml)");
	chkCompressSavedFiles->SetBounds(15, 75, 100, 19);
	settingsPanel->Add(chkCompressSavedFiles);

	chkCheckForUpdates = new GLToggle(0, "Check for updates at startup");
	chkCheckForUpdates->SetBounds(15, 100, 160, 19);
	settingsPanel->Add(chkCheckForUpdates);

	chkAutoUpdateFormulas = new GLToggle(0, "Auto refresh formulas");
	chkAutoUpdateFormulas->SetBounds(15, 125, 160, 19);
	settingsPanel->Add(chkAutoUpdateFormulas);

	chkAntiAliasing = new GLToggle(0, "Anti-Aliasing");
	chkAntiAliasing->SetBounds(15, 150, 160, 19);
	settingsPanel->Add(chkAntiAliasing);

	chkWhiteBg = new GLToggle(0, "White Background");
	chkWhiteBg->SetBounds(15, 175, 160, 19);
	settingsPanel->Add(chkWhiteBg);

	GLTitledPanel *gasPanel = new GLTitledPanel("Gas settings");
	gasPanel->SetBounds(280, 2, 295, 203);
	Add(gasPanel);

	GLLabel *massLabel = new GLLabel("Gas molecular mass (g/mol):");
	massLabel->SetBounds(290, 25, 150, 19);
	gasPanel->Add(massLabel);

	gasmassText = new GLTextField(0, "");
	gasmassText->SetBounds(460, 20, 100, 19);
	gasPanel->Add(gasmassText);

		enableDecay = new GLToggle(0, "Gas half life (s):");
	enableDecay->SetBounds(290, 75, 150, 19);
	gasPanel->Add(enableDecay);

	halfLifeText = new GLTextField(0, "");
	halfLifeText->SetBounds(460, 80, 100, 19);
	gasPanel->Add(halfLifeText);
	
	GLLabel *outgassingLabel = new GLLabel("Final outgassing rate (mbar*l/sec):");
	outgassingLabel->SetBounds(290, 75, 150, 19);
	gasPanel->Add(outgassingLabel);

	outgassingText = new GLTextField(0, "");
	outgassingText->SetBounds(460, 70, 100, 19);
	outgassingText->SetEditable(FALSE);
	gasPanel->Add(outgassingText);

	GLLabel *influxLabel = new GLLabel("Total desorbed molecules:");
	influxLabel->SetBounds(290, 100, 150, 19);
	gasPanel->Add(influxLabel);

	influxText = new GLTextField(0, "");
	influxText->SetBounds(460, 95, 100, 19);
	influxText->SetEditable(FALSE);
	gasPanel->Add(influxText);

	recalcButton = new GLButton(0, "Recalc. outgassing");
	recalcButton->SetBounds(460, 123, 100, 19);
	gasPanel->Add(recalcButton);


	applyButton = new GLButton(0, "Apply above settings");
	applyButton->SetBounds(wD / 2 - 65, 210, 130, 19);
	Add(applyButton);

	/*chkNonIsothermal = new GLToggle(0,"Non-isothermal system (textures only, experimental)");
	chkNonIsothermal->SetBounds(315,125,100,19);
	Add(chkNonIsothermal);*/

	GLTitledPanel *panel3 = new GLTitledPanel("Subprocess control");
	panel3->SetBounds(5, 259, wD - 10, hD - 285);
	Add(panel3);


	processList = new GLList(0);
	processList->SetHScrollVisible(TRUE);
	processList->SetSize(6, MAX_PROCESS);
	processList->SetColumnWidths((int*)plWidth);
	processList->SetColumnLabels((char **)plName);
	processList->SetColumnAligns((int *)plAligns);
	processList->SetColumnLabelVisible(TRUE);
	processList->SetBounds(10, 278, wD - 20, hD - 355);
	panel3->Add(processList);

	char tmp[128];
	sprintf(tmp, "Number of CPU cores:     %d", mApp->numCPU);
	GLLabel *coreLabel = new GLLabel(tmp);
	coreLabel->SetBounds(10, hD - 74, 120, 19);
	panel3->Add(coreLabel);
	GLLabel *l1 = new GLLabel("Number of subprocesses:");
	l1->SetBounds(10, hD - 49, 120, 19);
	panel3->Add(l1);

	nbProcText = new GLTextField(0, "");
	nbProcText->SetEditable(TRUE);
	nbProcText->SetBounds(135, hD - 51, 30, 19);
	panel3->Add(nbProcText);

	restartButton = new GLButton(0, "Apply and restart processes");
	restartButton->SetBounds(170, hD - 51, 150, 19);
	panel3->Add(restartButton);

	maxButton = new GLButton(0, "Change MAX desorbed molecules");
	maxButton->SetBounds(wD - 195, hD - 51, 180, 19);
	panel3->Add(maxButton);

	// ---------------------------------------------------



	// Center dialog
	int wS, hS;
	GLToolkit::GetScreenSize(&wS, &hS);
	int xD = (wS - wD) / 2;
	int yD = (hS - hD) / 2;
	SetBounds(xD, yD, wD, hD);

	RestoreDeviceObjects();

	lastUpdate = 0.0f;
	for (int i = 0; i < MAX_PROCESS; i++) lastCPUTime[i] = -1.0f;
	memset(lastCPULoad, 0, MAX_PROCESS*sizeof(float));
}


// --------------------------------------------------------------------

void GlobalSettings::Update() {

	char tmp[256];
	chkAntiAliasing->SetState(mApp->antiAliasing);
	chkWhiteBg->SetState(mApp->whiteBg);
	//chkNonIsothermal->SetState(nonIsothermal);
	UpdateOutgassing();

	gasmassText->SetText(worker->gasMass);

	enableDecay->SetState(worker->enableDecay);
	halfLifeText->SetText(worker->halfLife);
	halfLifeText->SetEditable(worker->enableDecay);

	autoSaveText->SetText(mApp->autoSaveFrequency);
	chkSimuOnly->SetState(mApp->autoSaveSimuOnly);
	chkCheckForUpdates->SetState(mApp->checkForUpdates);
	chkAutoUpdateFormulas->SetState(mApp->autoUpdateFormulas);
	chkCompressSavedFiles->SetState(mApp->compressSavedFiles);

	int nb = worker->GetProcNumber();
	sprintf(tmp, "%d", nb);
	nbProcText->SetText(tmp);
}

// ----------------------------------------------------------------

void GlobalSettings::SMPUpdate(float appTime) {

	if (!IsVisible() || IsIconic()) return;
	int nb = worker->GetProcNumber();

	if (appTime - lastUpdate > 1.0 && nb > 0) {

		char tmp[512];
		PROCESS_INFO pInfo;
		int  states[MAX_PROCESS];
		char statusStr[MAX_PROCESS][64];

		memset(states, 0, MAX_PROCESS*sizeof(int));
		worker->GetProcStatus(states, (char **)statusStr);

		processList->ResetValues();
		for (int i = 0; i < nb; i++) {
			DWORD pid = worker->GetPID(i);
			sprintf(tmp, "%d", i + 1);
			processList->SetValueAt(0, i, tmp);
			sprintf(tmp, "%d", pid);
			processList->SetValueAt(1, i, tmp);
			if (!GetProcInfo(pid, &pInfo)) {
				processList->SetValueAt(2, i, "0 KB");
				processList->SetValueAt(3, i, "0 KB");
				processList->SetValueAt(4, i, "0 %");
				processList->SetValueAt(5, i, "Dead");
			}
			else {
				processList->SetValueAt(2, i, FormatMemory(pInfo.mem_use));
				processList->SetValueAt(3, i, FormatMemory(pInfo.mem_peak));

				// CPU usage
				if (lastCPUTime[i] != -1.0f) {
					float dTime = appTime - lastUpdate;
					float dCPUTime = (float)pInfo.cpu_time - lastCPUTime[i];
					float cpuLoad = dCPUTime / dTime;
					lastCPULoad[i] = 0.85f*cpuLoad + 0.15f*lastCPULoad[i];
					int percent = (int)(100.0f*lastCPULoad[i] + 0.5f);
					if (percent < 0) percent = 0;
					sprintf(tmp, "%d %%", percent);
					processList->SetValueAt(4, i, tmp);
				}
				else {
					processList->SetValueAt(4, i, "---");
				}
				lastCPUTime[i] = (float)pInfo.cpu_time;

				// State/Status
				char status[128];
				_snprintf(tmp, 127, "%s: %s", prStates[states[i]], statusStr[i]);
				status[127] = 0;
				processList->SetValueAt(5, i, tmp);

			}
		}

		lastUpdate = appTime;
	}

}
// ----------------------------------------------------------------

void GlobalSettings::RestartProc() {


	int nbProc;
	if (sscanf(nbProcText->GetText(), "%d", &nbProc) == 0) {
		GLMessageBox::Display("Invalid process number", "Error", GLDLG_OK, GLDLG_ICONERROR);
	}
	else {
		//char tmp[128];
		//sprintf(tmp,"Kill all running sub-process(es) and start %d new ones ?",nbProc);
		//int rep = GLMessageBox::Display(tmp,"Question",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONWARNING);

		if (mApp->AskToReset()) {
			if (nbProc <= 0 || nbProc > MAX_PROCESS) {
				GLMessageBox::Display("Invalid process number [1..32]", "Error", GLDLG_OK, GLDLG_ICONERROR);
			}
			else {
				try {
					worker->SetProcNumber(nbProc);
					worker->Reload();
					mApp->SaveConfig();
				}
				catch (Error &e) {
					GLMessageBox::Display((char *)e.GetMsg(), "Error", GLDLG_OK, GLDLG_ICONERROR);
				}
			}
		}
	}

}
// --------------------------------------------------------------------

void GlobalSettings::ProcessMessage(GLComponent *src, int message) {

	switch (message) {
	case MSG_BUTTON:

		if (src == recalcButton) {
			if (mApp->AskToReset()) {
				try {
					worker->RealReload();
				}
				catch (Error &e) {
					GLMessageBox::Display(e.GetMsg(), "Recalculation failed: Couldn't reload Worker", GLDLG_OK, GLDLG_ICONWARNING);
				}
			}
		}
		else if (src == restartButton) {
			RestartProc();
		}
		else if (src == maxButton) {
			if (worker->GetGeometry()->IsLoaded()) {
				char tmp[128];
				sprintf(tmp, "%I64d", worker->maxDesorption);
				char *val = GLInputBox::GetInput(tmp, "Desorption max (0=>endless)", "Edit MAX");
				if (val) {
					llong maxDes;
					if (sscanf(val, "%I64d", &maxDes) == 0) {
						GLMessageBox::Display("Invalid 'maximum desorption' number", "Error", GLDLG_OK, GLDLG_ICONERROR);
					}
					else {
						worker->SetMaxDesorption(maxDes);
					}
				}
			}
			else {
				GLMessageBox::Display("No geometry loaded.", "No geometry", GLDLG_OK, GLDLG_ICONERROR);
			}
		}
		else if (src == applyButton) {
			mApp->antiAliasing = chkAntiAliasing->GetState();
			mApp->whiteBg = chkWhiteBg->GetState();
			mApp->checkForUpdates = chkCheckForUpdates->GetState();
			mApp->autoUpdateFormulas = chkAutoUpdateFormulas->GetState();
			mApp->compressSavedFiles = chkCompressSavedFiles->GetState();
			mApp->autoSaveSimuOnly = chkSimuOnly->GetState();
			double gm;
			if (!gasmassText->GetNumber(&gm) || !(gm > 0.0)) {
				GLMessageBox::Display("Invalid gas mass", "Error", GLDLG_OK, GLDLG_ICONERROR);
				return;
			}
			if (abs(gm - worker->gasMass) > 1e-7) {
				if (mApp->AskToReset()) {
					worker->needsReload = TRUE;
					worker->gasMass = gm;
					if (worker->GetGeometry()->IsLoaded()) { //check if there are pumps
						BOOL hasPump = FALSE;
						int nbFacet = worker->GetGeometry()->GetNbFacet();
						for (int i = 0; (i<nbFacet) && (!hasPump); i++) {
							if (worker->GetGeometry()->GetFacet(i)->sh.sticking>0.0) {
								hasPump = TRUE;
							}
						}
						if (hasPump) GLMessageBox::Display("Don't forget the pumps: update pumping speeds and/or recalculate sticking factors.", "You have changed the gas mass.", GLDLG_OK, GLDLG_ICONINFO);
					}
				}
			}

			double hl;
			if (enableDecay->GetState() && (!halfLifeText->GetNumber(&hl) || !(hl > 0.0))) {
				GLMessageBox::Display("Invalid half life", "Error", GLDLG_OK, GLDLG_ICONERROR);
				return;
			}
			if (enableDecay->GetState()!=worker->enableDecay || abs(hl - worker->halfLife) > 1e-7) {
				if (mApp->AskToReset()) {
					worker->needsReload = TRUE;
					worker->enableDecay = enableDecay->GetState();
					if (worker->enableDecay) worker->halfLife = hl;
				}
			}



			double autosavefreq;
			if (!autoSaveText->GetNumber(&autosavefreq) || !(autosavefreq > 0.0)) {
				GLMessageBox::Display("Invalid autosave frequency", "Error", GLDLG_OK, GLDLG_ICONERROR);
				return;
			}
			mApp->autoSaveFrequency = autosavefreq;

			return;

		}
		break;

	case MSG_TEXT:
		ProcessMessage(applyButton, MSG_BUTTON);
		break;

	case MSG_TOGGLE:
		if (src == enableDecay)
			halfLifeText->SetEditable(enableDecay->GetState());
		break;
	}

	GLWindow::ProcessMessage(src, message);
}

void GlobalSettings::UpdateOutgassing() {
	char tmp[128];
	sprintf(tmp, "%g", worker->gasMass);
	gasmassText->SetText(tmp);
	sprintf(tmp, "%g", worker->finalOutgassingRate_Pa_m3_sec * 10.00); //10: conversion Pa*m3/sec -> mbar*l/s
	outgassingText->SetText(tmp);
	sprintf(tmp, "%.3E", worker->totalDesorbedMolecules);
	influxText->SetText(tmp);
}