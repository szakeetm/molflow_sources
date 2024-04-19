
#ifndef _TIMEWISEPLOTTERH_
#define _TIMEWISEPLOTTERH_

#include "GLApp/GLWindow.h"
#include "GLApp/GLChart/GLChartConst.h"
#include <vector>
class GLChart;
class GLLabel;
class GLCombo;
class GLButton;
class GLFormula;
class GLDataView;
class GLToggle;
class GLTextField;
class Worker;
class InterfaceGeometry;

class TimewisePlotter : public GLWindow {

public:

  // Construction
  TimewisePlotter();

  // Component method
  void Display(Worker *w);
  void Refresh();
  void Update(float appTime,bool force=false);
  void UpdateMoment();
  void Reset();

  // Implementation
  void ProcessMessage(GLComponent *src,int message) override;
  void SetBounds(int x,int y,int w,int h) override;
  void refreshViews();

private:
  bool ParseMoments();
  void ParseToken(std::string token);
  void FormatParseText();
  void addView(int facet);
  void remView(int facet);
  
  Worker      *worker;
  GLButton    *dismissButton;
  GLLabel     *warningLabel;
  GLChart     *chart;
  GLCombo     *profCombo;
  GLLabel     *normLabel;
  GLCombo     *displayModeCombo;
  GLButton    *selButton;
  GLTextField *momentsText;
  GLLabel     *momLabel,*momentsLabel;
  GLToggle    *logYToggle,*constantFlowToggle,*correctForGas;

  std::vector<size_t> displayedMoments;

  GLDataView  *views[50];

  int          nbView;
  float        lastUpdate;

};

#endif /* _TIMEWISEPLOTTERH_ */
