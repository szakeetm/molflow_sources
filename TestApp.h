/*
  File:        TestApp.cpp
  Description: Test application
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

#include "GLApp/GLApp.h"
#include "GLApp/GLLabel.h"
#include "GLApp/GLChart.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLLabel.h"
#include "GLApp/GLTextField.h"
#include "GLApp/GLSpinner.h"
#include "GLApp/GLList.h"
#include "DrawingArea.h"

class TestApp : public GLApplication
{
private:
    
  DrawingArea *drawingArea;
  GLList      *list;
  GLLabel     *console;

protected:

    int  OneTimeSceneInit();
    int  RestoreDeviceObjects();
    int  InvalidateDeviceObjects();
    int  OnExit();
    int  FrameMove();
    void ProcessMessage(GLComponent *src,int message);

public:
    TestApp();

};
