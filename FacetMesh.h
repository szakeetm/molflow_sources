/*
  File:        FacetMesh.h
  Description: Facet mesh configuration dialog
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
#include "GLApp/GLToggle.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLTextField.h"
#include "GLApp/GLTitledPanel.h"
#include "GLApp/GLProgress.h"
#include "Worker.h"

#ifndef _FACETMESHH_
#define _FACETMESHH_

class FacetMesh : public GLWindow {

public:

  // Construction
  FacetMesh();

  // Component method
  void EditFacet(Worker *w);

  // Implementation
  void ProcessMessage(GLComponent *src,int message);

private:

  BOOL Apply();
  void UpdateSize();
  void UpdateSizeForRatio();
  void UpdateToggle(GLComponent *src);
  void QuickApply(); //Apply View Settings without stopping the simulation

  Worker   *worker;
  Geometry *geom;
  int       fIdx;

  GLTitledPanel *iPanel;
  GLTitledPanel *vPanel;
  GLTextField   *uLength;
  GLTextField   *vLength;

  GLToggle      *enableBtn;
  GLToggle      *boundaryBtn;
  GLToggle      *recordDesBtn;
  GLToggle      *recordAbsBtn;
  GLToggle      *recordReflBtn;
  GLToggle      *recordTransBtn;
  GLToggle      *recordACBtn;
  GLToggle      *recordDirBtn;
  GLToggle      *showTexture;
  GLToggle      *showVolume;
  GLTextField   *resolutionText;

  GLTextField   *ramText;
  GLTextField   *cellText;
  GLButton      *updateButton;

  GLButton    *applyButton;
  GLButton	  *quickApply; //Apply View Settings without stopping the simulation
  GLButton    *cancelButton;

  GLProgress  *progressDlg;

};

#endif /* _FACETMESHH_ */
