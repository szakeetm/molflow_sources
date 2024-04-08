
#pragma once

#include "GLApp/GLWindow.h"
class GLButton;
class GLCombo;
class GLTextField;
class GLLabel;
class GLToggle;
class GLTitledPanel;
class GLGradient;
class GeometryViewer;
class MolflowGeometry;
class Worker;

class TextureScaling : public GLWindow {

public:

	// Construction
	TextureScaling(Worker *worker_, GeometryViewer **viewers_);

	// Component methods
	void Display();
	void Update();
	void UpdateAutoScaleLimits();

	// Implementation
	void ProcessMessage(GLComponent *src,int message);

private:

	void RecalcSwapSize();

	Worker*				worker;
	MolflowGeometry*	interfGeom;
	GeometryViewer**	viewers;

	GLToggle      *autoScaleToggle;
    GLCombo       *autoscaleTimedepModeCombo;
	GLTextField   *manualScaleMinText;
	GLTextField   *manualScaleMaxText;
	GLLabel       *geomMinLabel;
	GLLabel       *geomMaxLabel;
	GLToggle      *useColorToggle;
	GLTextField   *swapText;
	GLGradient    *gradient;
	GLToggle      *logarithmicToggle;

	GLButton    *setToCurrentButton;
	GLCombo     *physicsModeCombo;
	GLButton    *applyButton;

};
