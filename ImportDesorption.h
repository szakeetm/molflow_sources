/*
  File:        ImportDesorption.h
  Description: Collapse settings dialog
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
#include "GLApp/GLTitledPanel.h"
#include "GLApp/GLCombo.h"

#include "Geometry.h"
#include "Worker.h"

#ifndef _IMPORTDESORPTIONH_
#define _IMPORTDESORPTIONH_

class ImportDesorption : public GLWindow {

public:

  // Construction
  ImportDesorption();

  // Component methods
  void SetGeometry(Geometry *s,Worker *w);

  // Implementation
  void ProcessMessage(GLComponent *src,int message);

private:

	//double InterpolateY(double x,const std::vector<std::pair<double,double>>& table,BOOL limitToBounds=FALSE,BOOL logarithmic=FALSE);
	void    LoadConvFile(char* fileName);
	void   EnableDisableComponents();

  Geometry      *geom;
  Worker	    *work;

  GLTitledPanel *filePanel;
  GLTextField   *synFileName;
  GLButton      *loadSynButton,*reloadButton;
  GLButton      *useCurrentButton;
  GLLabel       *analysisResultLabel;

  GLTitledPanel *importPanel;
  GLCombo       *sourceSelectCombo;
  GLTextField   *timeField;

  GLToggle      *r1,*r2,*r3;
  GLTextField   *eta0Text,*alphaText,*cutoffText;
  GLTextField   *convFileName;
  GLButton      *loadConvButton;
  GLButton      *convInfoButton;
  GLLabel       *convAnalysisLabel;

  GLButton      *setButton;
  GLButton      *cancelButton;

  size_t doseSource;
  size_t mode;
  std::vector<std::pair<double,double>> convDistr;
};

#endif /* _IMPORTDESORPTIONH_ */

