/*
Program:     MolFlow+ / Synrad+
Description: Monte Carlo simulator for ultra-high vacuum and synchrotron radiation
Authors:     Jean-Luc PONS / Roberto KERSEVAN / Marton ADY
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
#ifndef _PROFILEPLOTTERH_
#define _PROFILEPLOTTERH_

#include "GLApp/GLWindow.h"
#include "GLApp/GLChart/GLChartConst.h"
#include <vector>
#include <map>

class GLChart;
class GLLabel;
class GLCombo;
class GLButton;
class GLParser;
class GLDataView;
class GLToggle;
class GLTextField;
class Worker;
class Geometry;

class ProfilePlotter : public GLWindow {

public:

  // Construction
  ProfilePlotter();

  // Component method
  void Display(Worker *w);
  void Refresh();
  void Update(float appTime,bool force=false);
  void Reset();
    void ResetHighlighting();

  // Implementation
  void ProcessMessage(GLComponent *src,int message) override;
  void SetBounds(int x,int y,int w,int h);
  int addView(int facet);
  std::vector<int> GetViews();
  void SetViews(const std::vector<int> &updatedViews);
  bool IsLogScaled();
  void SetLogScaled(bool logScale);
  void SetWorker(Worker *w);
    std::map<int,GLColor> GetIDColorPairs() const;

private:  
  int remView(int facet);
  void refreshViews();
  void plot();
  void applyFacetHighlighting() const;

  Worker      *worker;
  GLButton    *dismissButton;
  GLChart     *chart;
  GLCombo     *profCombo;
  GLTextField *selFacInput;
  GLLabel     *normLabel;
  GLLabel     *warningLabel;
  GLCombo     *normCombo;
  //GLToggle    *showAllMoments;

  GLButton    *selButton;
  GLButton    *addButton;
  GLButton    *removeButton;
  GLButton    *removeAllButton;
  GLTextField *formulaText;
  GLButton    *formulaBtn;
  GLToggle    *logYToggle;
  GLToggle    *correctForGas;

    GLToggle    *colorToggle;
    GLLabel    *fixedLineWidthText;
    GLButton    *fixedLineWidthButton;
    GLTextField    *fixedLineWidthField;
    GLToggle    *useProfColToggle;
    GLButton    *selectPlottedButton;

    GLDataView  *views[MAX_VIEWS];
  int          nbView;
  float        lastUpdate;
    std::map<int,GLColor> plottedFacets;
};

#endif /* _PROFILEPLOTTERH_ */
