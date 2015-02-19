/*
  File:        TimewisePlotter.h
  Description: Timewise Profile plotter window
  Program:     MolFlow
  Author:      R. KERSEVAN / J-L PONS / M ADY
  Copyright:   E.S.R.F / CERN

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#include "GLApp/GLWindow.h"
#include "GLApp/GLChart.h"
#include "GLApp/GLLabel.h"
#include "GLApp/GLCombo.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLParser.h"
#include "GLApp/GLTextField.h"
#include "Worker.h"
#include "Geometry.h"

#ifndef _TIMEWISEPLOTTERH_
#define _TIMEWISEPLOTTERH_

class TimewisePlotter : public GLWindow {

public:

  // Construction
  TimewisePlotter();

  // Component method
  void Display(Worker *w);
  void Refresh();
  void Update(float appTime,BOOL force=FALSE);
  void UpdateMoment();
  void Reset();

  // Implementation
  void ProcessMessage(GLComponent *src,int message);
  void SetBounds(int x,int y,int w,int h);
  void refreshViews();

private:
  BOOL ParseMoments();
  void ParseToken(std::string token);
  void addView(int facet);
  void remView(int facet);
  
  Worker      *worker;
  GLButton    *dismissButton;
  GLChart     *chart;
  GLCombo     *profCombo;
  GLLabel     *normLabel;
  GLCombo     *normCombo;
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
