/*
  File:        TimeSettings.h
  Description: Move vertex by offset dialog
*/

#include "GLApp/GLWindow.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLTextField.h"
#include "GLApp/GLLabel.h"

#include "Geometry.h"
#include "Worker.h"

#ifndef _TIMESETTINGSH_
#define _TIMESETTINGSH_

class TimeSettings : public GLWindow {

public:

  // Construction
  TimeSettings(Worker *work);

  // Implementation
  void ProcessMessage(GLComponent *src,int message);
  void RefreshMoments();

private:

  Worker	   *work;

  GLButton    *setButton;
  GLButton    *previousButton,*ffBackButton;
  GLButton    *nextButton,*ffForwardButton;
  //GLButton    *cancelButton;
  GLButton    *editButton;
  GLLabel     *l1;
  GLLabel     *timeLabel;
  //GLLabel     *l3;
  GLTextField *timeId,*ffStep;
  //GLTextField *yOffset;
  //GLTextField *zOffset;

};

#endif /* _TIMESETTINGSH_ */
