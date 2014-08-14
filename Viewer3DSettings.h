/*
  File:        Viewer3DSettings.h
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

#include "GLApp/GLWindow.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLTextField.h"
#include "GLApp/GLLabel.h"
#include "GLApp/GLToggle.h"
#include "GLApp/GLCombo.h"
#include "GLApp/GLTitledPanel.h"
#include "Geometry.h"
#include "GeometryViewer.h"

#ifndef _VIEWER3DSETTINGSH_
#define _VIEWER3DSETTINGSH_

class Viewer3DSettings : public GLWindow {

public:

  // Construction
  Viewer3DSettings();

  // Component methods
  void Display(Geometry *s,GeometryViewer *v);

  // Implementation
  void ProcessMessage(GLComponent *src,int message);

private:

  Geometry       *geom;
  GeometryViewer *viewer;

  GLTitledPanel *panel;
  GLCombo       *showMode;
  GLTextField   *traStepText;
  GLTextField   *angStepText;
  GLTextField   *dispNumHits;
  GLTextField   *dispNumLeaks;
  GLToggle      *hiddenEdge;
  GLToggle      *hiddenVertex;
  GLToggle      *showMesh;
  GLToggle      *showDirection;
  GLToggle      *showTimeToggle;
  GLToggle      *autoNorme;
  GLToggle      *centerNorme;
  GLToggle		*antiAliasing;
  GLToggle		*bigDots;
  GLTextField   *normeText;

  GLButton    *applyButton;
  GLButton    *cancelButton;

};

#endif /* _VIEWER3DSETTINGSH_ */
