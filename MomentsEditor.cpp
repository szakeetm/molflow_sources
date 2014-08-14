/*
  File:        MomentsEditor.cpp
  Description: Moments Editor
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

#include "MomentsEditor.h"
#include "GLApp/GLTitledPanel.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLWindowManager.h"
#include "GLApp/GLMessageBox.h"
#include "MolFlow.h"

extern MolFlow *mApp;

  static const int   flWidth[] = {35,100,35};
  static const char *flName[] = {"#","Time (s)","Nb"};
  static const int   flAligns[] = { ALIGN_CENTER,ALIGN_LEFT,ALIGN_CENTER };
  static const int   fEdits[] = { 0,EDIT_STRING,0 };

MomentsEditor::MomentsEditor(Worker *w):GLWindow() {

  int wD = 220;
  int hD = 401;

  work=w;

  SetTitle("Edit time moments");

  panel1=new GLTitledPanel("Moment list");
  panel1->SetBounds(5,5,wD-10,250);
  Add(panel1);

  momentsList = new GLList(0);
  momentsList->SetBounds(10,22,wD-20,200);
  momentsList->SetColumnLabelVisible(TRUE);
  momentsList->SetGrid(TRUE);
  Add(momentsList);

  clearButton = new GLButton(0,"Clear list");
  clearButton->SetBounds(10,229,95,20);
  Add(clearButton);

  pasteButton = new GLButton(0,"Paste clipboard");
  pasteButton->SetBounds(110,229,95,20);
  pasteButton->SetEnabled(FALSE);
  Add(pasteButton);

  //char tmp[128];

  panel2=new GLTitledPanel("Time parameters");
  panel2->SetBounds(5,260,wD-10,95);
  Add(panel2);

  /*GLLabel *startLabel = new GLLabel("Desorption starts at:                   s");
  startLabel->SetBounds(15,275,170,25);
  Add(startLabel);*/


  //sprintf(tmp,"%g",work->desorptionStartTime);
/*  desStartText = new GLTextField(0,"");
  desStartText->SetBounds(120,275,60,20);
  Add(desStartText);*/

  /*GLLabel *stopLabel = new GLLabel("Desorption stops at:                   s");
  stopLabel->SetBounds(15,300,170,25);
  Add(stopLabel);*/
  
  //sprintf(tmp,"%g",work->desorptionStopTime);
  /*desStopText = new GLTextField(0,"");
  desStopText->SetBounds(120,300,60,20);
  Add(desStopText);*/

  GLLabel *windowLabel = new GLLabel("Time window length:                  s");
  windowLabel->SetBounds(15,275,170,25);
  Add(windowLabel);

  //sprintf(tmp,"%g",work->timeWindowSize);
  windowSizeText = new GLTextField(0,"");
  windowSizeText->SetBounds(120,275,60,20);
  Add(windowSizeText);

  useMaxwellToggle = new GLToggle(0,"Use Maxwell-B. speed distr.");
  useMaxwellToggle->SetBounds(15,300,wD-25,20);
  //useMaxwellToggle->SetCheck(work->useMaxwellDistribution);
  Add(useMaxwellToggle);

  calcConstantFlow = new GLToggle(0,"Calculate constant flow");
  calcConstantFlow->SetBounds(15,325,wD-25,20);
  //useMaxwellToggle->SetCheck(work->useMaxwellDistribution);
  Add(calcConstantFlow);

  /*GLLabel *valveLabel = new GLLabel("Facets 1,2 open at:                   s");
  valveLabel->SetBounds(15,400,170,25);
  Add(valveLabel);*/
  
  //sprintf(tmp,"%g",work->desorptionStopTime);
 /* valveText = new GLTextField(0,"");
  valveText->SetBounds(120,400,60,20);
  Add(valveText);*/

  //RebuildList();

  setButton = new GLButton(0,"Apply");
  setButton->SetBounds(wD-165,hD-44,75,20);
  Add(setButton);

  cancelButton = new GLButton(0,"Dismiss");
  cancelButton->SetBounds(wD-85,hD-44,75,20);
  Add(cancelButton);

  // Center dialog
  int wS,hS;
  GLToolkit::GetScreenSize(&wS,&hS);
  int xD = (wS-wD)/2;
  int yD = (hS-hD)/2;
  SetBounds(xD,yD,wD,hD);

  RestoreDeviceObjects();

}

void MomentsEditor::ProcessMessage(GLComponent *src,int message) {
  switch(message) {
    case MSG_BUTTON:

		if(src==cancelButton) {
			
			GLWindow::ProcessMessage(NULL,MSG_CLOSE);

		} else if (src==setButton) {
			//validate user input
			double /*start,stop,*/window/*,valve*/;
			/*if (!(desStartText->GetNumber(&start))) {
				GLMessageBox::Display("Invalid desorption start","Error",GLDLG_OK,GLDLG_ICONERROR);
				return;
			}*/
			/*if (!(desStopText->GetNumber(&stop))) {
				GLMessageBox::Display("Invalid desorption stop","Error",GLDLG_OK,GLDLG_ICONERROR);
				return;
			}*/
			if (!(windowSizeText->GetNumber(&window))) {
				GLMessageBox::Display("Invalid window length","Error",GLDLG_OK,GLDLG_ICONERROR);
				return;
			}
			
			/*if (!(valveText->GetNumber(&valve))) {
				GLMessageBox::Display("Invalid valve open time","Error",GLDLG_OK,GLDLG_ICONERROR);
				return;
			}*/

			/*if (!(start>=0)) {
				GLMessageBox::Display("Start moment mustn't be negative","Error",GLDLG_OK,GLDLG_ICONERROR);
				return;
			}

			if (!(stop>start)) {
				GLMessageBox::Display("Desorption stop must be later than start","Error",GLDLG_OK,GLDLG_ICONERROR);
				return;
			}*/
			
			//apply settings
			if (mApp->AskToReset()) {
				moments=std::vector<double>();
				latest=0.0;

				for (size_t u=0;u!=userMoments.size();u++) {
					std::vector<double> newMoments=ParseMoment(userMoments[u]);
					for (size_t i=0;i!=newMoments.size();i++)
						if (newMoments[i]>latest) latest=newMoments[i];
					AddMoment(newMoments);
				}

				work->moments=moments;
				work->userMoments=userMoments;
				/*work->desorptionStartTime=start;
				work->desorptionStopTime=stop;*/
				work->timeWindowSize=window;
				work->useMaxwellDistribution=useMaxwellToggle->IsChecked();
				work->calcConstantFlow=calcConstantFlow->IsChecked();
				/*work->valveOpenMoment=valve;*/

				work->Reload();
				if (mApp->timeSettings) mApp->timeSettings->RefreshMoments();
				if (mApp->timewisePlotter) {
					mApp->timewisePlotter->Reset();
					mApp->timewisePlotter->refreshViews();
				}
			}
			//this->Refresh();
		} else if (src==clearButton) {
			if (GLDLG_OK == GLMessageBox::Display("Clear list?","Moments",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONINFO)){
				userMoments=std::vector<std::string>();
				RebuildList();
			}
		} else if (src==pasteButton) {
			PasteClipboard();
		}
		break;
	case MSG_TEXT:
	case MSG_LIST:
		for (int row=0;row<(momentsList->GetNbRow()-1);row++) {
			//char tmp=*(momentsList->GetValueAt(1,row));
			//if (tmp=='\0'){
			if (strcmp(momentsList->GetValueAt(1,row),userMoments[row].c_str())!=0) {
				if (*(momentsList->GetValueAt(1,row))!=0) 
					userMoments[row].assign(momentsList->GetValueAt(1,row)); //update
				else
					userMoments.erase(userMoments.begin()+row); //erase
				
				RebuildList();
				break;
			}
		}
		if (momentsList->GetValueAt(1,momentsList->GetNbRow()-1)!=0) {
			if (*(momentsList->GetValueAt(1,momentsList->GetNbRow()-1))!=0) {
			//Add new line
				userMoments.push_back(momentsList->GetValueAt(1,momentsList->GetNbRow()-1));
				RebuildList();
			}
		}
		break;
}

  GLWindow::ProcessMessage(src,message);
}

void MomentsEditor::RebuildList() {

	momentsList->SetSize(3,userMoments.size()+1);
	momentsList->SetColumnWidths((int*)flWidth);
	momentsList->SetColumnLabels((char **)flName);
	momentsList->SetColumnAligns((int *)flAligns);
	momentsList->SetColumnEditable((int *)fEdits);
  
  char tmp[128];
  size_t u;double latest=0.0;

  for (u=0;u<userMoments.size();u++){ 
	  sprintf(tmp,"%d",u+1);
	  momentsList->SetValueAt(0,u,tmp);
	  sprintf(tmp,"%s",userMoments[u].c_str());
	  momentsList->SetValueAt(1,u,tmp);
	  std::vector<double> newMoments=ParseMoment(userMoments[u]);
	  sprintf(tmp,"%d",newMoments.size());
	  momentsList->SetValueAt(2,u,tmp);
  }
  //last line, possibility to enter new value
  sprintf(tmp,"%d",u+1);
  momentsList->SetValueAt(0,u,tmp);

}

void MomentsEditor::Refresh() {
	userMoments=work->userMoments;
	moments=work->moments;
	char tmp[128];
	/*sprintf(tmp,"%g",work->desorptionStartTime);
	desStartText->SetText(tmp);  
	sprintf(tmp,"%g",work->desorptionStopTime);
	desStopText->SetText(tmp);*/
	sprintf(tmp,"%g",work->timeWindowSize);
	windowSizeText->SetText(tmp);
	/*sprintf(tmp,"%g",work->valveOpenMoment);
	valveText->SetText(tmp);*/
	useMaxwellToggle->SetCheck(work->useMaxwellDistribution);
	calcConstantFlow->SetCheck(work->calcConstantFlow);
	RebuildList();
}

int MomentsEditor::AddMoment(std::vector<double> newMoments) {
	int nb=(int)newMoments.size();
	for (int i=0;i<nb;i++)
		moments.push_back(newMoments[i]);
	return nb;
}

std::vector<double> MomentsEditor::ParseMoment(std::string userInput) {
	std::vector<double> parsedResult;
	double begin,interval,end;

	int nb=sscanf(userInput.c_str(),"%lf,%lf,%lf",&begin,&interval,&end);
	if (nb==1 && (begin>0.0)) {
		//One moment
		parsedResult.push_back(begin);
	//} else if (nb==3 && (begin>0.0) && (end>begin) && (interval<(end-begin)) && ((end-begin)/interval<300.0)) {
	} else if (nb==3 && (begin>0.0) && (end>begin) && (interval<(end-begin))) {
		//Range
		for (double time=begin;time<=end;time+=interval)
			parsedResult.push_back(time);
	}
	return parsedResult;
}

void MomentsEditor::PasteClipboard() {

	momentsList->PasteClipboardText();
}