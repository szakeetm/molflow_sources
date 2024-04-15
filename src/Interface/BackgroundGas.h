

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
  
  
  GLTextField*	mfpTextbox;
  GLTextField*	massTextbox;
  GLTextField*	cutoffSpeedTextbox;
  GLToggle* cutoffToggle;
  
  GLButton*		applyButton;

  InterfaceGeometry*	interfGeom;
  Worker*				work;

};
