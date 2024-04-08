

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
