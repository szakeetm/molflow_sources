/*
  File:        ParameterEditor.cpp
  Description: Parameter Editor
  Program:     MolFlow


  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#include "ParameterEditor.h"
#include "GLApp/GLTitledPanel.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLWindowManager.h"
#include "GLApp/GLMessageBox.h"
#include "MolFlow.h"
#include <sstream>
#include "GLApp/GLMessageBox.h"

extern MolFlow *mApp;

  static const int   flWidth[] = {107,107};
  static const char *flName[] = {"Time (s)","Value"};
  static const int   flAligns[] = { ALIGN_LEFT,ALIGN_LEFT };
  static const int   fEdits[] = { EDIT_STRING,EDIT_STRING };

ParameterEditor::ParameterEditor(Worker *w):GLWindow() {

  int wD = 700;
  int hD = 400;

  work=w;
  tempParam = Parameter();
  userValues = std::vector<std::pair<std::string, std::string>>();
  dataView = GLDataView();

  int hSpace = 10;
  int vSpace = 5;
  int col1 = 10;
  int col2 = 250;
  int cursorX = col1;
  int cursorY = 5;
  int buttonWidth = 110;
  int labelWidth = 30;
  int compHeight = 20;
  int panelHeight = 345;
  int listHeight = 265;

  SetTitle("Edit parameters");

  selectorCombo = new GLCombo(0);
  selectorCombo->SetBounds(cursorX, cursorY, col2 - cursorX - 5, compHeight);
  Add(selectorCombo);

  cursorX = col2;
  deleteButton = new GLButton(2, "Delete");
  deleteButton->SetBounds(cursorX, cursorY, buttonWidth, compHeight);
  Add(deleteButton);

  /*
  cursorX += buttonWidth+hSpace;
  newButton = new GLButton(1, "Create new");
  newButton->SetBounds(cursorX, cursorY, buttonWidth, compHeight);
  Add(newButton);*/

  cursorX = col1;
  cursorY += compHeight+vSpace;
  editorPanel=new GLTitledPanel("Editor");
  editorPanel->SetBounds(5, cursorY, wD - 10, panelHeight);
  Add(editorPanel);

  cursorY += compHeight;
  GLLabel *nameLabel = new GLLabel("Name:");
  nameLabel->SetBounds(cursorX, cursorY, labelWidth, compHeight);
  Add(nameLabel);

  cursorX += labelWidth + hSpace;
  nameField = new GLTextField(3, "Param1");
  nameField->SetBounds(cursorX, cursorY, col2 - col1 - labelWidth - 2 * hSpace, compHeight);
  Add(nameField);

  cursorX = col2;
  plotArea = new GLChart(0);
  plotArea->SetBounds(cursorX, cursorY, wD - col2 - hSpace, compHeight + vSpace + listHeight);
  plotArea->SetBorder(BORDER_BEVEL_IN);
  plotArea->GetY1Axis()->SetGridVisible(TRUE);
  plotArea->SetLabelVisible(FALSE);
  dataView.SetMarker(MARKER_DOT);
  dataView.SetLineWidth(2);
  dataView.SetMarkerSize(8);
  plotArea->GetY1Axis()->AddDataView(&dataView);
  plotArea->GetXAxis()->SetGridVisible(TRUE);
  plotArea->GetY1Axis()->SetAutoScale(TRUE);
  plotArea->GetY2Axis()->SetAutoScale(TRUE);
  plotArea->GetY1Axis()->SetAnnotation(VALUE_ANNO);
  plotArea->GetXAxis()->SetAnnotation(VALUE_ANNO);
  Add(plotArea);

  cursorX = col1;
  cursorY += compHeight + vSpace;
  list = new GLList(0);
  list->SetBounds(cursorX, cursorY, col2-col1-hSpace, listHeight);
  list->SetColumnLabelVisible(TRUE);
  list->SetGrid(TRUE);
  //list->SetSelectionMode(BOX_CELL);
  Add(list);

  cursorX = col1;
  cursorY += listHeight + vSpace;
  copyButton = new GLButton(0,"Copy to clipboard");
  copyButton->SetBounds(cursorX, cursorY, buttonWidth, compHeight);
  Add(copyButton);

  cursorX += buttonWidth + hSpace;
  pasteButton = new GLButton(0,"Paste from clipboard");
  pasteButton->SetBounds(cursorX,cursorY,buttonWidth,compHeight);
  pasteButton->SetEnabled(FALSE);
  Add(pasteButton);

  cursorX = wD-2*buttonWidth-2*hSpace;
  plotButton = new GLButton(0, "Plot");
  plotButton->SetBounds(cursorX, cursorY, buttonWidth, compHeight);
  
  Add(plotButton);

  cursorX += buttonWidth + hSpace;
  applyButton = new GLButton(0, "Apply");
  applyButton->SetBounds(cursorX, cursorY, buttonWidth, compHeight);
  
  Add(applyButton);

  UpdateCombo();
  RebuildList();
  
  // Center dialog
  int wS,hS;
  GLToolkit::GetScreenSize(&wS,&hS);
  int xD = (wS-wD)/2;
  int yD = (hS-hD)/2;
  SetBounds(xD,yD,wD,hD);
  

  RestoreDeviceObjects();
  
}



void ParameterEditor::ProcessMessage(GLComponent *src,int message) {
  switch(message) {
    case MSG_BUTTON:
		if (src==applyButton) {
			if (ValidateInput() && mApp->AskToReset()) {
				if (selectorCombo->GetSelectedIndex() >= 0 && selectorCombo->GetSelectedIndex() < selectorCombo->GetNbRow() - 1) {//existing param
					work->parameters[selectorCombo->GetSelectedIndex()] = tempParam;
					UpdateCombo();
					UpdateUserValues();
					Plot();
				} else if (selectorCombo->GetSelectedIndex() == selectorCombo->GetNbRow() - 1) { //new Param
					work->parameters.push_back(tempParam);
					UpdateCombo();
					selectorCombo->SetSelectedIndex(selectorCombo->GetNbRow() - 2);
					deleteButton->SetEnabled(TRUE);
					UpdateUserValues();
					Plot();
				}
				RebuildList();
				work->Reload();
			}
		} else if (src == deleteButton) {
			if (mApp->AskToReset()) {
				work->parameters.erase(work->parameters.begin() + selectorCombo->GetSelectedIndex());
				UpdateCombo();
				selectorCombo->SetSelectedIndex(0);
				UpdateUserValues();
				RebuildList();
				deleteButton->SetEnabled(FALSE);
			}
		} else if (src == plotButton) {
			ValidateInput();
			Plot();
		} else if (src==pasteButton) {
			PasteFromClipboard();
		}
		break;
	case MSG_TEXT:

	case MSG_LIST:
		for (int row=0;row<(list->GetNbRow()-1);row++) {
			if ((strcmp(list->GetValueAt(0, row), userValues[row].first.c_str()) != 0) ||
				(strcmp(list->GetValueAt(1, row), userValues[row].second.c_str()) != 0)) { //something changed
				if (*(list->GetValueAt(0, row)) != 0 ||
					*(list->GetValueAt(1, row)) != 0) //there is some content
					userValues[row]=std::make_pair(list->GetValueAt(0, row),list->GetValueAt(1, row)); //update pair
				else{ //no content
					userValues.erase(userValues.begin()+row); //erase
					RebuildList();
				}
				break;
			}
		}
		if ((list->GetValueAt(0,list->GetNbRow()-1)!=0 && (*(list->GetValueAt(0, list->GetNbRow() - 1)) != 0)) ||
			(list->GetValueAt(1,list->GetNbRow()-1)!=0 && (*(list->GetValueAt(1, list->GetNbRow() - 1)) != 0))) {
			//Add new line
			std::string first = (list->GetValueAt(0, list->GetNbRow() - 1) != 0) ? list->GetValueAt(0, list->GetNbRow() - 1) : "";
			std::string second = (list->GetValueAt(1, list->GetNbRow() - 1) != 0) ? list->GetValueAt(1, list->GetNbRow() - 1) : "";
			userValues.push_back(std::make_pair(first,second)); //add moment
			RebuildList();
		}
		break;
	case MSG_COMBO:
		if (selectorCombo->GetSelectedIndex() >= 0 && selectorCombo->GetSelectedIndex() < selectorCombo->GetNbRow() - 1) { //existing param
			UpdateUserValues();
			ValidateInput();
			Plot();
			deleteButton->SetEnabled(TRUE);
		}
		else {
			if (selectorCombo->GetSelectedIndex() != -1) { //New param
				userValues = std::vector<std::pair<std::string, std::string>>();
				char tmp[32];
				sprintf(tmp, "Param%d", selectorCombo->GetSelectedIndex() + 1);
				nameField->SetText(tmp);
				dataView.Reset();
			}
			deleteButton->SetEnabled(FALSE);
		}
		RebuildList();
		break;
}

  GLWindow::ProcessMessage(src,message);
}

void ParameterEditor::UpdateCombo() {
	selectorCombo->SetSize((int)work->parameters.size()+1);
	for (size_t i = 0; i < work->parameters.size(); i++)
		selectorCombo->SetValueAt((int)i, work->parameters[i].name.c_str());
	selectorCombo->SetValueAt((int)work->parameters.size(), "New...");
	if (selectorCombo->GetSelectedIndex() == -1 || selectorCombo->GetSelectedIndex() == (selectorCombo->GetNbRow() - 1)) selectorCombo->SetSelectedIndex(selectorCombo->GetNbRow() - 1);
}

void ParameterEditor::PasteFromClipboard() {
	
}

void ParameterEditor::CopyToClipboard() {

}

void ParameterEditor::Plot() {
	dataView.Reset();
	for (auto i : tempParam.values)
		dataView.Add(i.first, i.second);
	dataView.CommitChange();
}

void ParameterEditor::RebuildList() {

	list->SetSize(2, userValues.size()+1);
	list->SetColumnWidths((int*)flWidth);
	list->SetColumnLabels((char **)flName);
	list->SetColumnAligns((int *)flAligns);
	list->SetColumnEditable((int *)fEdits);
	//list->cEdits[0] = list->cEdits[1]= EDIT_NUMBER;

	//char tmp[128];
	for (int row = 0; row<(int)userValues.size(); row++){
		list->SetValueAt(0,row,userValues[row].first.c_str());
		list->SetValueAt(1,row,userValues[row].second.c_str());
	}
	//last line, possibility to enter new value
}

BOOL ParameterEditor::ValidateInput() {
	//Validate name
	std::string tempName = nameField->GetText();
	if (tempName.length() == 0) {
		GLMessageBox::Display("Parameter name can't be empty", "Invalid parameter definition", GLDLG_OK, GLDLG_ICONWARNING);
		return FALSE;
	}
	if (selectorCombo->GetSelectedIndex() == selectorCombo->GetNbRow() - 1) {
		for (auto p : work->parameters) {
			if (tempName.compare(p.name) == 0) {
				GLMessageBox::Display("This parameter name is already used", "Invalid parameter definition", GLDLG_OK, GLDLG_ICONWARNING);
				return FALSE;
			}
		}
	}
	if (!((tempName[0] >= 65 && tempName[0] <= 90) || (tempName[0] >= 97 && tempName[0] <= 122))) {
		GLMessageBox::Display("Parameter name must begin with a letter", "Invalid parameter definition", GLDLG_OK, GLDLG_ICONWARNING);
		return FALSE;
	}
	tempParam = Parameter();
	tempParam.name = tempName;

	//TODO check:
	//no two values for the same moment

	BOOL atLeastOne = FALSE;
	for (size_t row = 0; row < userValues.size(); row++) {
		double valueX, valueY;
		try {
			valueX = std::stod(userValues[row].first);
		} catch (std::exception err){
			char tmp[256];
			sprintf(tmp, "Can't parse value \"%s\" in row %d, first column:\n%s", userValues[row].first.c_str(), row+1, err.what());
			GLMessageBox::Display(tmp, "Invalid parameter definition", GLDLG_OK, GLDLG_ICONWARNING);
			return FALSE;
		}
		try {
			valueY = std::stod(userValues[row].second);
		}
		catch (std::exception err){
			char tmp[256];
			sprintf(tmp, "Can't parse value \"%s\" in row %d, second column:\n%s", userValues[row].second.c_str(), row+1, err.what());
			GLMessageBox::Display(tmp, "Invalid parameter definition", GLDLG_OK, GLDLG_ICONWARNING);
			return FALSE;
		}
		tempParam.AddValue(std::make_pair(valueX, valueY));
		atLeastOne = TRUE;
	}
	if (!atLeastOne) {
		GLMessageBox::Display("At least one value is required", "Invalid parameter definition", GLDLG_OK, GLDLG_ICONWARNING);
		return FALSE;
	}

	for (int i = 0; i < tempParam.values.size();i++) {
		for (int j = i+1; j < tempParam.values.size(); j++) {
			if (abs(tempParam.values[i].first - tempParam.values[j].first) < 1E-8) {
				std::stringstream msg;
				msg << "There are two values for t=" << tempParam.values[i].first << "s.";
				GLMessageBox::Display(msg.str().c_str(), "Invalid parameter definition", GLDLG_OK, GLDLG_ICONWARNING);
				return FALSE;
			}
		}
	}

	return TRUE;
}

void ParameterEditor::UpdateUserValues() {
	userValues = std::vector<std::pair<std::string, std::string>>();
	nameField->SetText("");
	if (selectorCombo->GetSelectedIndex() <(int)work->parameters.size()) {
		Parameter *getParam = &work->parameters[selectorCombo->GetSelectedIndex()];

		for (int row = 0; row < (int)getParam->values.size(); row++) {
			std::ostringstream str1, str2;
			str1 << getParam->values[row].first;
			str2 << getParam->values[row].second;
			userValues.push_back(std::make_pair(str1.str(), str2.str()));
		}
		nameField->SetText(getParam->name.c_str());
	}
}