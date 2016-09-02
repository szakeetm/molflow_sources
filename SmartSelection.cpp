/*
File:        SmartSelection.cpp
Description: Smart selection dialog
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

#include "SmartSelection.h"
#include "GLApp/GLTitledPanel.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLWindowManager.h"
#include "GLApp/GLMessageBox.h"
#include "Molflow.h"

extern MolFlow *mApp;

SmartSelection::SmartSelection(Geometry *g,Worker *w):GLWindow() {

	int wD = 270;
	int hD = 125;

	SetTitle("Smart selection");

	enableToggle = new GLToggle(0,"Enable smart selection");
	enableToggle->SetBounds(5,5,170,18);
	enableToggle->SetState(TRUE);
	Add(enableToggle);

	angleThreshold = new GLTextField(0,"30");
	angleThreshold->SetBounds(50,30,80,18);
	Add(angleThreshold);

	resultLabel = new GLLabel("");
	resultLabel->SetBounds(10,75,wD-20,hD-125);
	Add(resultLabel);

	analyzeButton = new GLButton(0,"Analyze");
	analyzeButton->SetBounds(wD-265,hD-44,85,21);
	Add(analyzeButton);

	// Center dialog
	/*int wS,hS;
	GLToolkit::GetScreenSize(&wS,&hS);
	int xD = (wS-wD)/2;
	int yD = (hS-hD)/2;*/
	SetBounds(10,30,wD,hD);

	RestoreDeviceObjects();

	geom = g;
	work = w;
	isRunning = FALSE;
}

BOOL SmartSelection::IsSmartSelection()
{
	if (!IsVisible()) return FALSE;
	else return (enableToggle->GetState()==1);
}

double SmartSelection::GetMaxAngle()
{
	if (!IsVisible()) return -1.0;
	double maxAngleDiff;
	if (!angleThreshold->GetNumber(&maxAngleDiff) || !(maxAngleDiff>=0.0)) {
		GLMessageBox::Display("Invalid angle threshold in Smart Selection dialog\nMust be a non-negative number.", "Error", GLDLG_OK, GLDLG_ICONERROR);
		return -1.0;
	}
	return maxAngleDiff/180.0*3.14159;
}

void SmartSelection::ProcessMessage(GLComponent *src,int message) {
	

	switch(message) {
	case MSG_BUTTON:

		if (src==analyzeButton) {
			if (!isRunning) {
				if (!geom->IsLoaded()) {
					GLMessageBox::Display("Geometry not loaded yet", "Error", GLDLG_OK, GLDLG_ICONERROR);
					return;
				}
				analyzeButton->SetText("Stop analyzing");
				isRunning = TRUE;
				GLProgress *progressDlg = new GLProgress("Analyzing facets", "Please wait");
				progressDlg->SetProgress(0.0);
				progressDlg->SetVisible(TRUE);

				geom->AnalyzeNeighbors(work,progressDlg);

				progressDlg->SetVisible(FALSE);
				SAFE_DELETE(progressDlg);
				analyzeButton->SetText("Analyze");
				isRunning = FALSE;
			}
			else {
				analyzeButton->SetText("Analyze");
				isRunning = FALSE;
				work->abortRequested = TRUE;
			}
		}
		break;
	}

	GLWindow::ProcessMessage(src,message);
}

