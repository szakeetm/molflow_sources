
#ifndef _TEXTUREPLOTTERH_
#define _TEXTUREPLOTTERH_

#include "GLApp/GLWindow.h"
class GLButton;
class GLList;
class GLCombo;
class GLLabel;
class GLToggle;
class InterfaceGeometry;
class GeometryViewer;
class Worker;
class InterfaceFacet;

class TexturePlotter : public GLWindow {

public:

  // Construction
  TexturePlotter();

  // Component methods
  void Display(Worker *w);
  void Update(float appTime,bool force = false);

  // Implementation
  void ProcessMessage(GLComponent *src,int message) override;
  void SetBounds(int x,int y,int w,int h) override;

private:

  void GetSelected();
  void UpdateTable();
  void PlaceComponents();
  void Close();
  void SaveFile();

  Worker       *worker;
  InterfaceFacet*      selFacet;
  int         selFacetId;
  float        lastUpdate;
  double    	maxValue;
  size_t			maxX,maxY;

  GLList      *mapList;
  GLButton    *saveButton;
  GLButton    *sizeButton;
  GLLabel     *viewLabel;
  GLCombo     *viewCombo;
  GLToggle    *autoSizeOnUpdate;

  GLButton    *cancelButton;
  GLButton	  *maxButton;

};

#endif /* _TEXTUREPLOTTERH_ */
