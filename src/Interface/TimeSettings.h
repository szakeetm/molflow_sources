

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
