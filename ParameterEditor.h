/*
  File:        ParameterEditor.h
  Description: Moments Editor
*/

#include "GLApp/GLWindow.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLCombo.h"
#include "GLApp/GLChart.h"
#include "GLApp\GLTextField.h"
#include "GLApp\GLToggle.h"
#include "GLApp\GLTitledPanel.h"
#include "GLApp/GLList.h"
#include "GLApp/GLChart/GLDataView.h"

#include "Worker.h"

#ifndef _PARAMETEREDITORH_
#define _PARAMETEREDITORH_

class ParameterEditor : public GLWindow {

public:

  // Construction
  ParameterEditor(Worker *work);
  void UpdateCombo();
  void RebuildList(BOOL autoSize=TRUE, BOOL refillValues=TRUE);
  // Implementation
  void ProcessMessage(GLComponent *src,int message);
  

private:

  Worker	   *work;

  GLCombo *selectorCombo;
  GLButton *newButton, *deleteButton,/* *copyButton ,*/ *pasteButton, *loadCSVbutton, *plotButton, *applyButton;
  GLList *list;
  GLTextField *nameField;
  GLTitledPanel *editorPanel;
  GLChart *plotArea;
  GLDataView dataView;
  GLToggle *logXtoggle;
  GLToggle *logYtoggle;

  Parameter tempParam;
  std::vector<std::pair<std::string, std::string>> userValues;
  
  void CopyToClipboard();
  void PasteFromClipboard();
  void LoadCSV();
  void Plot();
  void UpdateUserValues();
  BOOL ValidateInput();
  
};

#endif /* _PARAMETEREDITORH_ */
