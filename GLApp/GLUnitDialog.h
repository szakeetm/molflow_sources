/*
  File:        GLUnitDialog.h
  Description: STL import unit dialog

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/


#include <SDL_opengl.h>
#include "GLWindow.h"

#ifndef _GLUNITDIALOGH_
#define _GLUNITDIALGOH_

// Buttons
#define GLDLG_MM		0x0001
#define GLDLG_CM		0x0002
#define GLDLG_M			0x0003
#define GLDLG_INCH		0x0004
#define GLDLG_FOOT		0x0005
#define GLDLG_CANCEL_U  0x0006

// Icons
#define GLDLG_ICONNONE    0
#define GLDLG_ICONERROR   1
#define GLDLG_ICONWARNING 2
#define GLDLG_ICONINFO    3

class GLUnitDialog : private GLWindow {

public:
  // Display a modal dialog and return the code of the pressed button
  static int Display(char *message,char *title=NULL,int mode=GLDLG_MM|GLDLG_CM|GLDLG_M|GLDLG_INCH|GLDLG_FOOT|GLDLG_CANCEL_U,int icon=GLDLG_ICONNONE);

  int  rCode;

private:
  GLUnitDialog(char *message,char *title,int mode,int icon);
  void ProcessMessage(GLComponent *src,int message);

};

#endif /* _GLMESSAGEBOXH_ */
