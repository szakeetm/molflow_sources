

/*
  File:        MomentsEditor.h
  Description: Moments Editor
*/
#ifndef _MOMENTSEDITORH_
#define _MOMENTSEDITORH_

#include "GLApp/GLWindow.h"
class GLButton;
class GLLabel;
class GLTextField;
class GLToggle;
class GLTitledPanel;
class GLList;

#include "Worker.h"

class MomentsEditor : public GLWindow {

public:

  // Construction
  MomentsEditor(Worker *work);
  void RebuildList();
  void Refresh();

  // Implementation
  void ProcessMessage(GLComponent *src,int message) override;
  

private:

  Worker	   *work;

  GLButton    *setButton;
  GLButton    *cancelButton;
  GLButton    *clearButton;
  GLButton    *pasteButton;
  GLLabel     *l1;
  GLList      *momentsList;
  GLTextField *windowSizeText;
  GLToggle    *useMaxwellToggle;
  GLToggle    *calcConstantFlow;
  GLTitledPanel *panel1;
  GLTitledPanel *panel2;

  
  std::vector<Moment> moments;
  std::vector<UserMoment> userMoments;

  int AddMoment(std::vector<Moment> newMoments); //Adds a time serie to moments and returns the number of elements
  void PasteClipboard();
  std::vector<Moment> ParseUserMoment(const UserMoment& input); //Parses a user input and returns a vector of time moments

};

#endif /* _MOMENTSEDITORH_ */
