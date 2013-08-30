/*
  File:        ViewEditor.h
  Description: View editor window
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
#include "GLApp/GLList.h"
#include "GLApp/GLButton.h"
#include "Worker.h"

#ifndef _VIEWEDITORH_
#define _VIEWEDITORH_

class ViewEditor : public GLWindow {

public:

  // Construction
  ViewEditor();

  // Component method
  void Display(Worker *w);

  // Implementation
  void ProcessMessage(GLComponent *src,int message);

private:

  Worker      *worker;
  GLButton    *dismissButton;

};

#endif /* _VIEWEDITORH_ */
