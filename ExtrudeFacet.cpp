/*
File:        MoveFacet.cpp
Description: Move facet by offset dialog
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

#include "ExtrudeFacet.h"
#include "GLApp/GLTitledPanel.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLWindowManager.h"
#include "GLApp/GLMessageBox.h"
#include "MolFlow.h"

extern MolFlow *mApp;

ExtrudeFacet::ExtrudeFacet(Geometry *g, Worker *w) :GLWindow() {

	int wD = 312;
	int hD = 240;
	groupBox1 = new GLTitledPanel("Use facet normal");
	groupBox1->SetBounds(5, 3, 285, 63);
	Add(groupBox1);
	groupBox2 = new GLTitledPanel("Define offset");
	groupBox2->SetBounds(5, 69, 285, 114);
	Add(groupBox2);
	offsetCheckbox = new GLToggle(0, "Use offset");
	groupBox2->SetCompBounds(offsetCheckbox, 10, 18, 74, 17);
	groupBox2->Add(offsetCheckbox);

	label3 = new GLLabel("cm");
	groupBox2->SetCompBounds(label3, 75, 45, 21, 13);
	groupBox2->Add(label3);

	dxText = new GLTextField(0, "");
	groupBox2->SetCompBounds(dxText, 35, 39, 40, 20);
	groupBox2->Add(dxText);

	label4 = new GLLabel("dX:");
	groupBox2->SetCompBounds(label4, 10, 45, 23, 13);
	groupBox2->Add(label4);

	label5 = new GLLabel("dY:");
	groupBox2->SetCompBounds(label5, 100, 45, 23, 13);
	groupBox2->Add(label5);

	label6 = new GLLabel("cm");
	groupBox2->SetCompBounds(label6, 165, 45, 21, 13);
	groupBox2->Add(label6);

	dyText = new GLTextField(0, "");
	groupBox2->SetCompBounds(dyText, 125, 39, 40, 20);
	groupBox2->Add(dyText);

	label7 = new GLLabel("dZ:");
	groupBox2->SetCompBounds(label7, 195, 45, 23, 13);
	groupBox2->Add(label7);

	label8 = new GLLabel("cm");
	groupBox2->SetCompBounds(label8, 260, 45, 21, 13);
	groupBox2->Add(label8);

	dzText = new GLTextField(0, "");
	groupBox2->SetCompBounds(dzText, 220, 39, 40, 20);
	groupBox2->Add(dzText);

	extrudeButton = new GLButton(0, "Extrude");
	extrudeButton->SetBounds(95, 189, 105, 21);
	Add(extrudeButton);

	getBaseButton = new GLButton(0, "Get Base Vertex");
	groupBox2->SetCompBounds(getBaseButton, 35, 69, 105, 21);
	groupBox2->Add(getBaseButton);

	getDirButton = new GLButton(0, "Get Dir. Vertex");
	groupBox2->SetCompBounds(getDirButton, 145, 69, 105, 21);
	groupBox2->Add(getDirButton);

	towardsNormalCheckbox = new GLToggle(0, "Towards normal");
	groupBox1->SetCompBounds(towardsNormalCheckbox, 5, 15, 101, 17);
	groupBox1->Add(towardsNormalCheckbox);

	againstNormalCheckbox = new GLToggle(0, "Against normal");
	groupBox1->SetCompBounds(againstNormalCheckbox, 5, 36, 95, 17);
	groupBox1->Add(againstNormalCheckbox);

	label1 = new GLLabel("extrude by:");
	groupBox1->SetCompBounds(label1, 160, 36, 59, 13);
	groupBox1->Add(label1);

	distanceTextbox = new GLTextField(0, "");
	groupBox1->SetCompBounds(distanceTextbox, 215, 30, 45, 20);
	groupBox1->Add(distanceTextbox);

	label2 = new GLLabel("cm");
	groupBox1->SetCompBounds(label2, 260, 36, 21, 13);
	groupBox1->Add(label2);

	dirLabel = new GLLabel("dirLabel");
	groupBox2->SetCompBounds(dirLabel, 175, 93, 44, 13);
	groupBox2->Add(dirLabel);

	baseLabel = new GLLabel("baseLabel");
	groupBox2->SetCompBounds(baseLabel, 55, 93, 56, 13);
	groupBox2->Add(baseLabel);

	SetTitle("Extrude Facet");
	// Center dialog
	int wS, hS;
	GLToolkit::GetScreenSize(&wS, &hS);
	int xD = (wS - wD) / 2;
	int yD = (hS - hD) / 2;
	SetBounds(xD, yD, wD, hD);

	baseLabel->SetText("");
	dirLabel->SetText("");
	dxText->SetEditable(FALSE);
	dyText->SetEditable(FALSE);
	dzText->SetEditable(FALSE);
	towardsNormalCheckbox->SetState(1);
	distanceTextbox->SetText("1");
	baseId = dirId = -1;

	RestoreDeviceObjects();

	geom = g;
	work = w;

}

void ExtrudeFacet::ProcessMessage(GLComponent *src, int message) {
	double dX,dY,dZ,dist;

	switch(message) {
	case MSG_BUTTON:

		if (src==extrudeButton) {
			if (geom->GetNbSelected()==0) {
				GLMessageBox::Display("No facets selected","Nothing to move",GLDLG_OK,GLDLG_ICONINFO);
				return;
			}
			else if (geom->GetNbSelected() > 1) {
				char warningMsg[512];
				sprintf(warningMsg, "Extrude %d facets at once?", geom->GetNbSelected());
				int rep = GLMessageBox::Display(warningMsg, "Extrusion of more than one facet", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONINFO);
				if (rep != GLDLG_OK) {
					return;
				}
			}

			if( offsetCheckbox->GetState()==1 && !dxText->GetNumber(&dX) ) {
				GLMessageBox::Display("Invalid X offset value","Error",GLDLG_OK,GLDLG_ICONERROR);
				return;
			}
			if (offsetCheckbox->GetState() == 1 && !dyText->GetNumber(&dY)) {
				GLMessageBox::Display("Invalid Y offset value","Error",GLDLG_OK,GLDLG_ICONERROR);
				return;
			}
			if (offsetCheckbox->GetState() == 1 && !dzText->GetNumber(&dZ)) {
				GLMessageBox::Display("Invalid Z offset value","Error",GLDLG_OK,GLDLG_ICONERROR);
				return;
			}
			if ( (towardsNormalCheckbox->GetState() == 1 || againstNormalCheckbox->GetState()==1) 
				&& !distanceTextbox->GetNumber(&dist)) {
				GLMessageBox::Display("Invalid 'extrude by' value", "Error", GLDLG_OK, GLDLG_ICONERROR);
				return;
			}
			if (mApp->AskToReset()){

				VERTEX3D dir;
				double sign = (againstNormalCheckbox->GetState() == 1) ? -1.0 : 1.0;
				dir.x = (offsetCheckbox->GetState()==1)?dX:0.0;
				dir.y = (offsetCheckbox->GetState() == 1) ? dY : 0.0;
				dir.z = (offsetCheckbox->GetState() == 1) ? dZ : 0.0;
				geom->Extrude(dir, sign*dist);
				//mApp->UpdateModelParams();
				work->Reload(); 
				mApp->changedSinceSave = TRUE;
				mApp->UpdateFacetlistSelected();	
				mApp->UpdateViewers();

				//GLWindowManager::FullRepaint();
			}
		}
		else if (src == getBaseButton) {
			int nbFound = 0; int foundId;
			for (int i = 0; i < geom->GetNbVertex() && nbFound<2; i++) {
				if (geom->GetVertex(i)->selected) {
					nbFound++;
					foundId = i;
				}
			}
			if (nbFound == 0) {
				GLMessageBox::Display("No vertex selected", "Can't define base", GLDLG_OK, GLDLG_ICONINFO);
				return;
			}
			else if (nbFound > 1) {
				GLMessageBox::Display("More than one vertex is selected", "Can't define base", GLDLG_OK, GLDLG_ICONINFO);
				return;
			}
			else {

				towardsNormalCheckbox->SetState(0);
				againstNormalCheckbox->SetState(0);
				offsetCheckbox->SetState(1);
				dxText->SetEditable(TRUE);
				dyText->SetEditable(TRUE);
				dzText->SetEditable(TRUE);
				distanceTextbox->SetEditable(FALSE);

				baseId = foundId;
				char tmp[32];
				sprintf(tmp, "Vertex %d",baseId+1);
				baseLabel->SetText(tmp);
				if (dirId>0 && dirId < geom->GetNbVertex()) {
					dxText->SetText(geom->GetVertex(dirId)->x - geom->GetVertex(baseId)->x);
					dyText->SetText(geom->GetVertex(dirId)->y - geom->GetVertex(baseId)->y);
					dzText->SetText(geom->GetVertex(dirId)->z - geom->GetVertex(baseId)->z);
				}
			}
		}
		else if (src == getDirButton) {
			int nbFound = 0; int foundId;
			for (int i = 0; i < geom->GetNbVertex() && nbFound<2; i++) {
				if (geom->GetVertex(i)->selected) {
					nbFound++;
					foundId = i;
				}
			}
			if (nbFound == 0) {
				GLMessageBox::Display("No vertex selected", "Can't define direction", GLDLG_OK, GLDLG_ICONINFO);
				return;
			}
			else if (nbFound > 1) {
				GLMessageBox::Display("More than one vertex is selected", "Can't define direction", GLDLG_OK, GLDLG_ICONINFO);
				return;
			}
			else {

				towardsNormalCheckbox->SetState(0);
				againstNormalCheckbox->SetState(0);
				offsetCheckbox->SetState(1);
				dxText->SetEditable(TRUE);
				dyText->SetEditable(TRUE);
				dzText->SetEditable(TRUE);
				distanceTextbox->SetEditable(FALSE);

				dirId = foundId;
				char tmp[32];
				sprintf(tmp, "Vertex %d", dirId + 1);
				dirLabel->SetText(tmp);
				if (baseId>0 && baseId < geom->GetNbVertex()) {
					dxText->SetText(geom->GetVertex(dirId)->x - geom->GetVertex(baseId)->x);
					dyText->SetText(geom->GetVertex(dirId)->y - geom->GetVertex(baseId)->y);
					dzText->SetText(geom->GetVertex(dirId)->z - geom->GetVertex(baseId)->z);
				}
			}
		}
		break;
	case MSG_TOGGLE:
		towardsNormalCheckbox->SetState(src == towardsNormalCheckbox);
		againstNormalCheckbox->SetState(src == againstNormalCheckbox);
		offsetCheckbox->SetState(src == offsetCheckbox);
		distanceTextbox->SetEditable(src != offsetCheckbox);
		dxText->SetEditable(src == offsetCheckbox);
		dyText->SetEditable(src == offsetCheckbox);
		dzText->SetEditable(src == offsetCheckbox);
		break;
	}

	GLWindow::ProcessMessage(src,message);
}

