
#ifndef _OutgassingMapH_
#define _OutgassingMapH_

#include "GLApp/GLWindow.h"
class GLButton;
class GLList;
class GLCombo;
class GLLabel;
class GLTextField;
class InterfaceGeometry;
class GeometryViewer;
class Worker;
class InterfaceFacet;

class OutgassingMapWindow : public GLWindow {

public:

  // Construction
  OutgassingMapWindow();

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
  InterfaceFacet        *selFacet;
  float        lastUpdate;
  float			maxValue;
  int			maxX,maxY;

  GLLabel     *desLabel;
  GLList      *mapList;
  GLButton    *saveButton;
  GLButton    *sizeButton;
  GLLabel     *viewLabel;
  GLCombo     *desCombo;
  GLButton    *pasteButton;

  GLButton    *explodeButton;
  GLButton    *cancelButton;
  GLButton	  *maxButton;

  GLTextField *exponentText;

};

#endif /* _OutgassingMapH_ */
