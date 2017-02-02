/*
File:        TextureSettings.h
Description: Texture settings dialog (min,max,autoscale,gradient)
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
#include "GLApp/GLGradient.h"
#include "GeometryViewer.h"

#ifndef _TEXTURESETTINGSH_
#define _TEXTURESETTINGSH_

class MolflowGeometry;

class TextureSettings : public GLWindow {

public:

	// Construction
	TextureSettings();

	// Component methods
	void Display(Worker *w,GeometryViewer **v);
	void Update();

	// Implementation
	void ProcessMessage(GLComponent *src,int message);

private:

	void UpdateSize();

	Worker         *worker;
	MolflowGeometry       *geom;
	GeometryViewer **viewers;

	GLToggle      *texAutoScale;
	GLToggle      *includeConstantFlow;
	GLTextField   *texMinText;
	GLTextField   *texMaxText;
	GLLabel       *texCMinText;
	GLLabel       *texCMaxText;
	GLToggle      *colormapBtn;
	GLTextField   *swapText;
	GLGradient    *gradient;
	GLToggle      *logBtn;
	GLButton    *setCurrentButton;

	GLCombo       *modeCombo;
	GLButton    *updateButton;

};

#endif /* _TEXTURESETTINGSH_ */
