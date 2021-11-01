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
#pragma once

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
  void SetBounds(int x,int y,int w,int h) override;
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
  GLCombo     *displayModeCombo;
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

enum class profileRecordModes {
	None,
	RegularU,
	RegularV,
	IncAngle,
	Speed,
	OrtSpeed,
	TanSpeed,
	NUMITEMS
};

std::map<profileRecordModes,std::pair<std::string,std::string>> profileRecordModeDescriptions = { //mode, long description, short description
	{profileRecordModes::None, {"None","None"}},
	{profileRecordModes::RegularU, {"Pressure/imp/density (\201)","along \201"}},
	{profileRecordModes::RegularV, {"Pressure/imp/density (\202)","along \202"}},
	{profileRecordModes::IncAngle, {"Incident angle","Inc. angle"}},
	{profileRecordModes::Speed, {"Speed distribution","Speed"}},
	{profileRecordModes::OrtSpeed,{"Orthogonal velocity","Ort.velocity"}},
	{profileRecordModes::TanSpeed,{"Tangential velocity","Tan.velocity"}}
};

enum class profileDisplayModes {
  Raw,
  Pressure,
  ImpRate,
  Density,
  Speed,
  Angle,
  NormalizeTo1,
  NUMITEMS
};

std::map<profileDisplayModes,std::string> profileDisplayModeDescriptions = { //mode, description
	{profileDisplayModes::Raw, "Raw"},
	{profileDisplayModes::Pressure, "Pressure (mbar)"},
	{profileDisplayModes::ImpRate, "Impingement rate (1/m\262/sec)"},
	{profileDisplayModes::Density, "Density (1/m3)"},
	{profileDisplayModes::Speed, "Speed (m/s)"},
	{profileDisplayModes::Angle,"Angle (deg)"},
	{profileDisplayModes::NormalizeTo1,"Normalize to 1"}
};