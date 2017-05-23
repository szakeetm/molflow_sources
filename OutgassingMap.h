/*
  File:        OutgassingMap.h
  Description: Outgassing Map dialog
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
#ifndef _OutgassingMapH_
#define _OutgassingMapH_

#include "GLApp/GLWindow.h"
class GLButton;
class GLList;
class GLCombo;
class GLLabel;
class GLTextField;
class Geometry;
class GeometryViewer;
class Worker;
class Facet;

class OutgassingMap : public GLWindow {

public:

  // Construction
  OutgassingMap();

  // Component methods
  void Display(Worker *w);
  void Update(float appTime,bool force = false);

  // Implementation
  void ProcessMessage(GLComponent *src,int message);
  void SetBounds(int x,int y,int w,int h);

private:

  void GetSelected();
  void UpdateTable();
  void PlaceComponents();
  void Close();
  void SaveFile();

  Worker       *worker;
  Facet        *selFacet;
  float        lastUpdate;
  float			maxValue;
  int			maxX,maxY;
  char         currentDir[512];

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
