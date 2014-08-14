/*
  File:        TimeSettings.cpp
  Description: Time Settings dialog
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

#include "TimeSettings.h"
#include "GLApp/GLTitledPanel.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLWindowManager.h"
#include "GLApp/GLMessageBox.h"
#include "MolFlow.h"

extern MolFlow *mApp;

TimeSettings::TimeSettings(Worker *w):GLWindow() {

  int wD = 170;
  int hD = 70;

  SetTitle("Time Settings");

  
  
  GLLabel *l1 = new GLLabel("Moment #");
  l1->SetBounds(5,5,50,18);
  Add(l1);

  timeId = new GLTextField(0,"0");
  timeId->SetBounds(55,4,30,18);
  Add(timeId);  


  timeLabel = new GLLabel("Constant Flow");
  timeLabel->SetBounds(95,5,70,18);
  Add(timeLabel);

  /*yOffset = new GLTextField(0,"0");
  yOffset->SetBounds(100,30,80,18);
  Add(yOffset);
  
  GLLabel *l3 = new GLLabel("dZ");
  l3->SetBounds(10,55,170,18);
  Add(l3);

  zOffset = new GLTextField(0,"0");
  zOffset->SetBounds(100,55,80,18);
  Add(zOffset);*/

  /*setButton = new GLButton(0,"Go");
  setButton->SetBounds(100,5,25,21);
  Add(setButton);*/

  /*cancelButton = new GLButton(0,"Dismiss");
  cancelButton->SetBounds(150,hD-44,65,21);
  Add(cancelButton);*/

  previousButton = new GLButton(0,"<");
  previousButton->SetBounds(5,hD-43,25,18);
  Add(previousButton);

  char tmp[128];
  sprintf(tmp,"%d moments",w->moments.size());
  editButton = new GLButton(0,tmp);
  editButton->SetBounds(35,hD-43,100,18);
  Add(editButton);
  
  nextButton = new GLButton(0,">");
  nextButton->SetBounds(140,hD-43,20,18);
  Add(nextButton);


  // Center dialog
  /*int wS,hS;
  GLToolkit::GetScreenSize(&wS,&hS);
  int xD = (wS-wD)/2;
  int yD = (hS-hD)/2;*/
  SetBounds(8,30,wD,hD);

  RestoreDeviceObjects();
  work = w;

}



void TimeSettings::ProcessMessage(GLComponent *src,int message) {
  int id;
  int nbMoments=(int)work->moments.size();

  switch(message) {
  case MSG_TEXT:  
  case MSG_BUTTON:

    /*if(src==cancelButton) {

      GLWindow::ProcessMessage(NULL,MSG_CLOSE);

    } else */if (src==timeId || src==setButton || src==nextButton || src==previousButton) {

		if( !timeId->GetNumberInt(&id) ) {
        GLMessageBox::Display("Can't interpret time Id","Error",GLDLG_OK,GLDLG_ICONERROR);
        return;
      }
		if( id>nbMoments || id<0 ) {
        GLMessageBox::Display("This time id doesn't exist","Error",GLDLG_OK,GLDLG_ICONERROR);
        return;
      }
			      
		if (src==nextButton) {
			id++;
			if (id>nbMoments) id=0;
		}

		if (src==previousButton) {
			id--;
			if (id<0) id=nbMoments;
		}

		work->displayedMoment=id;
		char tmp[64];
		sprintf(tmp,"%d",id);
		timeId->SetText(tmp);
		if (id==0)
			timeLabel->SetText("Constant Flow");
		else {
			sprintf(tmp,"t=%gs",work->moments[id-1]);
			timeLabel->SetText(tmp);
		}
		try {
		  work->Update(0.0f);//update displayed profiles and textures
	  } catch(Error &e) {
		  GLMessageBox::Display((char *)e.GetMsg(),"Error (Worker::Update)",GLDLG_OK,GLDLG_ICONERROR);
	  } 
		if(mApp->profilePlotter) mApp->profilePlotter->Update(0.0f,TRUE);
		//if(mApp->pressureEvolution) mApp->pressureEvolution->Update(0.0f,TRUE);
		if(mApp->timewisePlotter) mApp->timewisePlotter->UpdateMoment();
		if(mApp->texturePlotter) mApp->texturePlotter->Update(0.0f,TRUE); 
    } else if (src==editButton) {
		if( mApp->momentsEditor==NULL ) mApp->momentsEditor = new MomentsEditor(work);
		mApp->momentsEditor->Refresh();
		//momentsEditor->Display(work,momentsEditor->GetSelectedIndex());
		mApp->momentsEditor->SetVisible(TRUE);
		//momentsEditor->DoModal();
		//SAFE_DELETE(momentsEditor);
	}
    break;
  }

  GLWindow::ProcessMessage(src,message);
}

void TimeSettings::RefreshMoments() {
	timeLabel->SetText("Constant Flow");
	timeId->SetText("0");
	work->displayedMoment=0;
	int nbMoments=(int)work->moments.size();
	char tmp[128];
	sprintf(tmp,"%d moments",nbMoments);
	editButton->SetText(tmp);
}