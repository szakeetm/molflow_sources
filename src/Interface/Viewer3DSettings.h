
#pragma once

#include "GLApp/GLWindow.h"
class GLButton;
class GLTextField;
class GLLabel;
class GLToggle;
class GLCombo;
class GLTitledPanel;
class InterfaceGeometry;
class GeometryViewer;

class Viewer3DSettings : public GLWindow {

public:
  // Construction
  Viewer3DSettings();

  // Component methods
  void Refresh(InterfaceGeometry *s,GeometryViewer *v);
  void Reposition(int wD = 0, int hD = 0);

  // Implementation
  void ProcessMessage(GLComponent *src,int message) override;

private:

  InterfaceGeometry*	interfGeom;
  GeometryViewer*		viewer;

  GLTitledPanel*	panel;
  GLCombo*			volumeRenderModeCombo;
  GLTextField*		traStepText;
  GLTextField*		angStepText;
  GLTextField*		dispNumHits;
  GLTextField*		dispNumLeaks;
  GLToggle*			hiddenEdge;
  GLToggle*			hiddenVertex;
  GLToggle*			showMesh;
  GLToggle*			dirShowdirToggle;
  GLToggle*			showTPtoggle;
  GLToggle*			showTimeToggle;
  GLToggle*			dirNormalizeToggle;
  GLToggle*			dirCenterToggle;
  GLToggle*			antiAliasing;
  GLToggle*			bigDots;
  GLToggle*			hideLotselected;
  GLTextField*		dirNormeText;
  GLTextField*		hideLotText;
  GLButton*			crossSectionButton;

  GLButton    *applyButton;
  GLButton    *cancelButton;

};
