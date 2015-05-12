/*
  File:        ExtrudeFacet.h
  Description: Extrude facet by offset dialog
*/

#include "GLApp/GLWindow.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLTextField.h"
#include "GLApp/GLLabel.h"
#include "GLApp/GLToggle.h"
#include "GLApp/GLTitledPanel.h"

#include "Geometry.h"
#include "Worker.h"

#ifndef _EXTRUDEFACETH_
#define _EXTRUDEFACETH_

class ExtrudeFacet : public GLWindow {

public:

  // Construction
  ExtrudeFacet(Geometry *geom,Worker *work);

  // Implementation
  void ProcessMessage(GLComponent *src,int message);

private:

  Geometry     *geom;
  Worker	   *work;

  GLToggle	*offsetCheckbox;
  GLLabel	*label3;
  GLTextField	*dxText;
  GLLabel	*label4;
  GLLabel	*label5;
  GLLabel	*label6;
  GLTextField	*dyText;
  GLLabel	*label7;
  GLLabel	*label8;
  GLTextField	*dzText;
  GLButton	*extrudeButton;
  GLButton	*getBaseButton;
  GLButton	*getDirButton;
  GLToggle	*towardsNormalCheckbox;
  GLToggle	*againstNormalCheckbox;
  GLLabel	*label1;
  GLTextField	*distanceTextbox;
  GLLabel	*label2;
  GLTitledPanel	*groupBox1;
  GLTitledPanel	*groupBox2;
  GLLabel	*dirLabel;
  GLLabel	*baseLabel;
	
  int baseId, dirId;

};

#endif /* _EXTRUDEFACETH_ */
