/*
Program:     MolFlow+ / Synrad+
Description: Monte Carlo simulator for ultra-high vacuum and synchrotron radiation
Authors:     Jean-Luc PONS / Roberto KERSEVAN / Marton ADY / Pascal BAEHR
Copyright:   E.S.R.F / CERN
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
  File:        MeasureForce.h
  Description: Enable/disable force measurement
*/
#pragma once

#include "GLApp/GLWindow.h"

class GLButton;
class GLTextField;
class GLLabel;
class GLTitledPanel;
class GLToggle;

class InterfaceGeometry;
class Worker;

class MeasureForce : public GLWindow {

public:
  // Construction
	MeasureForce(InterfaceGeometry *interfGeom, Worker *work);
	void ProcessMessage(GLComponent *src,int message) override;
	void UpdateToggle(GLComponent* src);
	void Update();

  // Implementation
private:

  
  GLButton *centerOfFacetButton, *useVertexButton, *applyButton, *dismissButton;
  GLTextField *x0Text, *y0Text, *z0Text;
  GLToggle* enableMeasureCheckbox;
  
  InterfaceGeometry     *interfGeom;
  Worker	   *work;

};
