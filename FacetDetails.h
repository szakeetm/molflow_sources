/*
  File:        FacetDetails.h
  Description: Facet details window
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

#ifndef _FACETDETAILSH_
#define _FACETDETAILSH_

#include "GLApp/GLWindow.h"
#include "GLApp/GLList.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLToggle.h"
#include "GLApp/GLTitledPanel.h"
#include "Worker.h"

class Facet;

#define NB_FDCOLUMN 26

class FacetDetails : public GLWindow {

public:

  // Construction
  FacetDetails();

  // Component method
  void Display(Worker *w);
  void Update();

  // Implementation
  void ProcessMessage(GLComponent *src,int message);
  void SetBounds(int x,int y,int w,int h);

private:

  char *GetCountStr(Facet *f);
  void UpdateTable();
  char *FormatCell(int idx,Facet *f,int mode);
  void PlaceComponents();

  Worker      *worker;
  GLList      *facetListD;

  GLTitledPanel *sPanel;          // Toggle panel
  GLToggle      *show[NB_FDCOLUMN];
  int            shown[NB_FDCOLUMN];

  GLButton    *checkAllButton;
  GLButton    *uncheckAllButton;
  GLButton    *dismissButton;
  GLButton	  *updateButton;

};

#endif /* _FACETDETAILSH_ */
