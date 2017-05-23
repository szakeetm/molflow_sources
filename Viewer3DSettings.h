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
#ifndef _VIEWER3DSETTINGSH_
#define _VIEWER3DSETTINGSH_

#include "GLApp/GLWindow.h"
class GLButton;
class GLTextField;
class GLLabel;
class GLToggle;
class GLCombo;
class GLTitledPanel;
class Geometry;
class GeometryViewer;

class Viewer3DSettings : public GLWindow {

public:

  // Construction
  Viewer3DSettings();

  // Component methods
  void Refresh(Geometry *s,GeometryViewer *v);
  void Reposition(int wD = 0, int hD = 0);

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
  GLToggle      *dirShowdirToggle;
  GLToggle      *showTimeToggle;
  GLToggle      *dirNormalizeToggle;
  GLToggle      *dirCenterToggle;
  GLToggle		*antiAliasing;
  GLToggle		*bigDots;
  GLToggle		*hideLotselected;
  GLTextField   *dirNormeText;
  GLTextField   *hideLotText;

  GLButton    *applyButton;
  GLButton    *cancelButton;

};

#endif /* _VIEWER3DSETTINGSH_ */
