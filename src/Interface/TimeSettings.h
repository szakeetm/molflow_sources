/*
Program:     MolFlow+ / Synrad+
Description: Monte Carlo simulator for ultra-high vacuum and synchrotron radiation
Authors:     Jean-Luc PONS / Roberto KERSEVAN / Marton ADY / Pascal BAEHR
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

/*
  File:        TimeSettings.h
  Description: Move vertex by offset dialog
*/
#ifndef _TIMESETTINGSH_
#define _TIMESETTINGSH_

#include "GLApp/GLWindow.h"
class GLButton;
class GLTextField;
class GLLabel;

class InterfaceGeometry;
class Worker;

class TimeSettings : public GLWindow {

public:

  // Construction
  TimeSettings(Worker *work);

  // Implementation
  void ProcessMessage(GLComponent *src,int message) override;
  void RefreshMoments();

private:

  Worker	   *work;

  GLButton    *setButton;
  GLButton    *previousButton,*ffBackButton;
  GLButton    *nextButton,*ffForwardButton;
  GLButton    *editButton;
  GLLabel     *l1;
  GLLabel     *timeLabel;
  GLTextField *timeId,*ffStep;

};

#endif /* _TIMESETTINGSH_ */
