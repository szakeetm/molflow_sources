/*
  File:        ViewEditor.cpp
  Description: View editor window
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

#include "ViewEditor.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLMessageBox.h"
#include "Utils.h"

ViewEditor::ViewEditor():GLWindow() {

  int wD = 350;
  int hD = 400;

  SetTitle("View editor");

  dismissButton = new GLButton(0,"Dismiss");
  dismissButton->SetBounds(wD-95,hD-43,90,19);
  Add(dismissButton);

  // Center dialog
  int wS,hS;
  GLToolkit::GetScreenSize(&wS,&hS);
  int xD = (wS-wD)/2;
  int yD = (hS-hD)/2;
  SetBounds(xD,yD,wD,hD);

  RestoreDeviceObjects();

}

void ViewEditor::Display(Worker *w) {

  worker = w;
  DoModal();

}

void ViewEditor::ProcessMessage(GLComponent *src,int message) {

  switch(message) {
    case MSG_BUTTON:
      if(src==dismissButton) {
        SetVisible(FALSE);
      }
    break;
  }

  GLWindow::ProcessMessage(src,message);

}

