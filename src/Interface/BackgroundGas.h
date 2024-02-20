/*
Program:     MolFlow+ / Synrad+
Description: Monte Carlo simulator for ultra-high vacuum and synchrotron radiation
Authors:     Jean-Luc PONS / Roberto KERSEVAN / Marton ADY
Website:     https://cern.ch/molflow

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

Full license text: https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html
*/

/*
  File:        Movement.h
  Description: Define moving parts in system
*/
#pragma once

#include "GLApp/GLWindow.h"
#include <vector>
class GLButton;
class GLTextField;
class GLLabel;
class GLToggle;
class GLTitledPanel;

class InterfaceGeometry;
class Worker;

class BackgroundGas : public GLWindow {

public:
  // Construction
	BackgroundGas(InterfaceGeometry *interfGeom, Worker *work);
  void ProcessMessage(GLComponent *src,int message) override;
  void Update();

  // Implementation
private:
  
  GLToggle* enableCheckbox;
  
  
  GLTextField	*mfpTextbox;
  GLTextField	*massTextbox;
  GLTextField	*rhoTextbox;
  
  GLButton	*applyButton;

  InterfaceGeometry     *interfGeom;
  Worker	   *work;

};
