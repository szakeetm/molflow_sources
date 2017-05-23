/*
  File:        ParameterEditor.h
  Description: Moments Editor
*/
#ifndef _PARAMETEREDITORH_
#define _PARAMETEREDITORH_

#include "GLApp/GLWindow.h"
#include "Parameter.h"
class GLButton;
class GLCombo;
class GLChart;
class GLTextField;
class GLToggle;
class GLTitledPanel;
class GLList;
class GLDataView;

class Worker;

class ParameterEditor : public GLWindow {

public:

  // Construction
  ParameterEditor(Worker *work);
  void UpdateCombo();
  void RebuildList(bool autoSize=true, bool refillValues=true);
  // Implementation
  void ProcessMessage(GLComponent *src,int message);
  void Reset();
  

private:

  Worker	   *work;

  GLCombo *selectorCombo;
  GLButton *newButton, *deleteButton,/* *copyButton ,*/ *pasteButton, *loadCSVbutton, *plotButton, *applyButton;
  GLList *list;
  GLTextField *nameField;
  GLTitledPanel *editorPanel;
  GLChart *plotArea;
  GLDataView *dataView;
  GLToggle *logXtoggle;
  GLToggle *logYtoggle;

  Parameter tempParam;
  std::vector<std::pair<std::string, std::string>> userValues;
  
  void CopyToClipboard();
  void PasteFromClipboard();
  void LoadCSV();
  void Plot();
  void UpdateUserValues();
  bool ValidateInput();
  
};

#endif /* _PARAMETEREDITORH_ */
