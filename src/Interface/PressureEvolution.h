

#pragma once
#include <vector>
#include "GLApp/GLWindow.h"
#include "GLApp/GLChart/GLChartConst.h"
class GLChart;
class GLLabel;
class GLCombo;
class GLButton;
class GLDataView;
class GLToggle;
class GLTextField;
class Worker;
class InterfaceGeometry;


class PressureEvolution : public GLWindow {

public:

  // Construction
  PressureEvolution(Worker *w);

  // Component method
  void Refresh();
  void Update(float appTime,bool force=false);
  void Reset();

  // Implementation
  void ProcessMessage(GLComponent *src,int message) override;
  void SetBounds(int x,int y,int w,int h) override;

private:

  void addView(size_t facetId);
  void remView(size_t viewId);
  void refreshChart();

  Worker      *worker;

  GLChart     *chart;
  GLLabel *label1, *normLabel;


  GLCombo     *profCombo;
  GLCombo     *yScaleCombo;
  GLButton    *selButton;
  GLButton    *addButton;
  GLButton    *removeButton;
  GLButton    *removeAllButton;

  GLToggle *logXToggle,*logYToggle;

  std::vector<GLDataView*>  views;
  std::vector<GLColor>    colors;
  float        lastUpdate;

};
