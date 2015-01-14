/*
  File:        Viewer3DSettings.cpp
  Description: 3D viewer settings dialog
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

#include "Viewer3DSettings.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLMessageBox.h"

// --------------------------------------------------------------------

Viewer3DSettings::Viewer3DSettings():GLWindow() {

  int wD = 215;
  int hD = 425;

  SetTitle("3D Viewer Settings");

  panel = new GLTitledPanel("3D Viewer settings");
  panel->SetBounds(5,5,wD-10,270);
  Add(panel);

  GLLabel *l4 = new GLLabel("Show facet");
  l4->SetBounds(10,25,90,18);
  Add(l4);

  showMode = new GLCombo(0);
  showMode->SetEditable(TRUE);
  showMode->SetSize(3);
  showMode->SetValueAt(0,"Front & Back");
  showMode->SetValueAt(1,"Back");
  showMode->SetValueAt(2,"Front");
  showMode->SetBounds(100,25,100,19);
  Add(showMode);

  GLLabel *l5 = new GLLabel("Translation step");
  l5->SetBounds(10,50,90,18);
  Add(l5);

  traStepText = new GLTextField(0,"");
  traStepText->SetBounds(100,50,100,18);
  Add(traStepText);

  GLLabel *l6 = new GLLabel("Angle step");
  l6->SetBounds(10,75,90,18);
  Add(l6);

  angStepText = new GLTextField(0,"");
  angStepText->SetBounds(100,75,100,18);
  Add(angStepText);

  //Number of hits displayed
  GLLabel *l8 = new GLLabel("Number of lines");
  l8->SetBounds(10,100,90,18);
  Add(l8);

  dispNumHits = new GLTextField(0,"");
  dispNumHits->SetBounds(100,100,100,18);
  Add(dispNumHits);

    //Number of hits displayed
  GLLabel *l9 = new GLLabel("Number of leaks");
  l9->SetBounds(10,125,90,18);
  Add(l9);

  dispNumLeaks = new GLTextField(0,"");
  dispNumLeaks->SetBounds(100,125,100,18);
  Add(dispNumLeaks);

  hiddenEdge = new GLToggle(0,"Show hidden edges (selected facets)");
  hiddenEdge->SetBounds(10,150,50,18);
  Add(hiddenEdge);

  hiddenVertex = new GLToggle(0,"Show hidden vertex (if selected)");
  hiddenVertex->SetBounds(10,175,50,18);
  Add(hiddenVertex);

  showMesh = new GLToggle(0,"Show texture mesh (Slow!)");
  showMesh->SetBounds(10,200,50,18);
  Add(showMesh);

  bigDots = new GLToggle(0,"Larger dots for hits");
  bigDots->SetBounds(10,225,50,18);
  Add(bigDots);

  showTimeToggle =  new GLToggle(0,"Show time overlay");
  showTimeToggle->SetBounds(10,250,50,18);
  Add(showTimeToggle);

  dirShowdirToggle = new GLToggle(0,"Show direction");
  dirShowdirToggle->SetBounds(10,280,190,18);
  Add(dirShowdirToggle);

  GLTitledPanel *panel2 = new GLTitledPanel("Direction field");
  panel2->SetBounds(5,305,wD-10,70);
  Add(panel2);

  GLLabel *l7 = new GLLabel("Norme ratio");
  l7->SetBounds(10,325,90,18);
  Add(l7);

  dirNormeText = new GLTextField(0,"");
  dirNormeText->SetBounds(100,325,100,18);
  Add(dirNormeText);

  dirNormalizeToggle = new GLToggle(0,"Normalize");
  dirNormalizeToggle->SetBounds(10,350,100,18);
  Add(dirNormalizeToggle);

  dirCenterToggle = new GLToggle(0,"Center");
  dirCenterToggle->SetBounds(110,350,90,18);
  Add(dirCenterToggle);

  applyButton = new GLButton(0,"Apply");
  applyButton->SetBounds(wD-170,hD-43,80,19);
  Add(applyButton);

  cancelButton = new GLButton(0,"Dismiss");
  cancelButton->SetBounds(wD-85,hD-43,80,19);
  Add(cancelButton);

  // Center dialog
  int wS,hS;
  GLToolkit::GetScreenSize(&wS,&hS);
  int xD = (wS-wD)/2;
  int yD = (hS-hD)/2;
  SetBounds(xD,yD,wD,hD);

  RestoreDeviceObjects();

  geom = NULL;

}

// --------------------------------------------------------------------

void Viewer3DSettings::Display(Geometry *s,GeometryViewer *v) {

  char tmp[128];

  geom = s;
  viewer = v;
  showMode->SetSelectedIndex(viewer->showBack);
  hiddenEdge->SetState(viewer->showHidden);
  hiddenVertex->SetState(viewer->showHiddenVertex);
  showMesh->SetState(viewer->showMesh);
  showTimeToggle->SetState(viewer->showTime);

  bigDots->SetState(viewer->bigDots);
  dirShowdirToggle->SetState(viewer->showDir);
  sprintf(tmp,"%g",viewer->transStep);
  traStepText->SetText(tmp);
  sprintf(tmp,"%g",viewer->angleStep);
  angStepText->SetText(tmp);
  sprintf(tmp,"%g",(double) viewer->dispNumHits);
  dispNumHits->SetText(tmp);
  sprintf(tmp,"%g",(double) viewer->dispNumLeaks);
  dispNumLeaks->SetText(tmp);
  sprintf(tmp,"Viewer #%d",viewer->GetId()+1);
  panel->SetTitle(tmp);
  sprintf(tmp,"%g",geom->GetNormeRatio());
  dirNormeText->SetText(tmp);
  dirNormalizeToggle->SetState( geom->GetAutoNorme() );
  dirCenterToggle->SetState( geom->GetCenterNorme() );
  DoModal();

}

// --------------------------------------------------------------------

void Viewer3DSettings::ProcessMessage(GLComponent *src,int message) {

  switch(message) {
    case MSG_BUTTON:

    if(src==cancelButton) {

      GLWindow::ProcessMessage(NULL,MSG_CLOSE);

    } else if (src==applyButton) {

      double tstep,astep,nratio;
	  int dnh,dnl;
	  //int dnh;

      if( !traStepText->GetNumber(&tstep) ) {
        GLMessageBox::Display("Invalid translation step value","Error",GLDLG_OK,GLDLG_ICONERROR);
        return;
      }
      if( !angStepText->GetNumber(&astep) ) {
        GLMessageBox::Display("Invalid angle step value","Error",GLDLG_OK,GLDLG_ICONERROR);
        return;
      }
	  if(( !dispNumHits->GetNumberInt(&dnh)||dnh<1||dnh>2048 )) {
        GLMessageBox::Display("Invalid number of displayed hits.\nMust be between 1 and 2048.","Error",GLDLG_OK,GLDLG_ICONERROR);
        return;
      }
	  if(( !dispNumLeaks->GetNumberInt(&dnl)||dnl<1||dnl>2048 )) {
        GLMessageBox::Display("Invalid number of displayed leaks.\nMust be between 1 and 2048.","Error",GLDLG_OK,GLDLG_ICONERROR);
        return;
      }
      viewer->showBack=showMode->GetSelectedIndex();
      viewer->transStep = tstep;
      viewer->angleStep = astep;
	  viewer->dispNumHits = dnh;
      viewer->dispNumLeaks = dnl;
      viewer->showHidden=hiddenEdge->GetState();
	  viewer->showHiddenVertex=hiddenVertex->GetState();
      viewer->showMesh=showMesh->GetState();

	  viewer->bigDots=bigDots->GetState();
      viewer->showDir=dirShowdirToggle->GetState();
	  viewer->showTime=showTimeToggle->GetState();

      if( !dirNormeText->GetNumber(&nratio) ) {
        GLMessageBox::Display("Invalid norme ratio value","Error",GLDLG_OK,GLDLG_ICONERROR);
        return;
      }
      geom->SetNormeRatio((float)nratio);
      geom->SetAutoNorme(dirNormalizeToggle->GetState());
      geom->SetCenterNorme(dirCenterToggle->GetState());

      GLWindow::ProcessMessage(NULL,MSG_CLOSE);

    }
    break;

    case MSG_TOGGLE:
    break;
  }

  GLWindow::ProcessMessage(src,message);
}
