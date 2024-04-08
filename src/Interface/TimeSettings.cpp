
#include "TimeSettings.h"
#include "GLApp/GLTitledPanel.h"
#include "GLApp/GLMessageBox.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLTextField.h"
#include "GLApp/GLLabel.h"

#include "MolFlow.h"
#include "MomentsEditor.h"
#include "FacetDetails.h"
#include "Interface/FormulaEditor.h"

extern MolFlow *mApp;

/**
* \brief Constructor with initialisation for Time settings window (Time/Time settings)
*/
TimeSettings::TimeSettings(Worker *w):GLWindow() {

  int wD = 170;
  int hD = 90;

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

  previousButton = new GLButton(0,"<");
  previousButton->SetBounds(5,25,25,18);
  Add(previousButton);

  char tmp[128];
  sprintf(tmp,"%zd moments",w->interfaceMomentCache.size());
  editButton = new GLButton(0,tmp);
  editButton->SetBounds(35,25,100,18);
  Add(editButton);
  
  nextButton = new GLButton(0,">");
  nextButton->SetBounds(140,25,25,18);
  Add(nextButton);

  ffBackButton= new GLButton(0, "<<");
  ffBackButton->SetBounds(5, 48, 25, 18);
  Add(ffBackButton);

  GLLabel *stepLabel = new GLLabel("Fast step:");
  stepLabel->SetBounds(45,48,40,18);
  Add(stepLabel);

  ffStep = new GLTextField(0, "10");
  ffStep->SetBounds(95, 48, 30, 18);
  Add(ffStep);

  ffForwardButton = new GLButton(0, ">>");
  ffForwardButton->SetBounds(140, 48, 25, 18);
  Add(ffForwardButton);

  SetBounds(8,30,wD,hD);

  RestoreDeviceObjects();
  work = w;

}

/**
* \brief Function for processing various inputs (button, check boxes etc.)
* \param src Exact source of the call
* \param message Type of the source (button)
*/
void TimeSettings::ProcessMessage(GLComponent *src,int message) {
  int id;
  int nbMoments=(int)work->interfaceMomentCache.size();
  int stepSize;

  switch(message) {
  case MSG_TEXT:  
  case MSG_BUTTON:

    /*if(src==cancelButton) {

      GLWindow::ProcessMessage(NULL,MSG_CLOSE);

	  } else */if (src == timeId || src == setButton || src == nextButton || src == previousButton || src == ffBackButton || src == ffForwardButton ) {

		if( !timeId->GetNumberInt(&id) ) {
        GLMessageBox::Display("Can't interpret time Id","Error",GLDLG_OK,GLDLG_ICONWARNING);
        return;
      }
		if( id>nbMoments || id<0 ) {
        GLMessageBox::Display("This time id doesn't exist","Error",GLDLG_OK,GLDLG_ICONERROR);
        return;
      }
			      
		if (src==nextButton) {
			id++;
			if (id>nbMoments) id=0;
		} else if (src==previousButton) {
			id--;
			if (id<0) id=nbMoments;
		}
		else if (src == ffForwardButton) {
			if (!ffStep->GetNumberInt(&stepSize) || id<0) {
				GLMessageBox::Display("Can't interpret step size", "Error", GLDLG_OK, GLDLG_ICONWARNING);
				return;
			}
			id=(id+stepSize)%(nbMoments+1);
		}
		else if (src == ffBackButton) {
			if ((!ffStep->GetNumberInt(&stepSize)) || !(id>0)) {
				GLMessageBox::Display("Can't interpret step size", "Error", GLDLG_OK, GLDLG_ICONWARNING);
				return;
			}
			id=(id-stepSize)%(nbMoments+1);
			if (id < 0) id = nbMoments + 1 + id; //To correct for the modulo operator implementation for negative dividends in C++
		}

		work->displayedMoment=id;
		char tmp[64];
		sprintf(tmp,"%d",id);
		timeId->SetText(tmp);
		if (id==0)
			timeLabel->SetText("Constant Flow");
		else {
			sprintf(tmp,"t=%gs",work->interfaceMomentCache[id-1].time);
			timeLabel->SetText(tmp);
		}
		try {
		  work->Update(0.0f);//update displayed profiles and textures, facet hits, etc
	  } catch (const std::exception &e) {
		  GLMessageBox::Display(e.what(),"Error (Worker::Update)",GLDLG_OK,GLDLG_ICONERROR);
	  } 
	  mApp->UpdatePlotters();
	  if (mApp->facetDetails) mApp->facetDetails->Update();
		//if (mApp->autoUpdateFormulas) mApp->UpdateFormula();
        if (mApp->autoUpdateFormulas && mApp->formulaEditor && mApp->formulaEditor->IsVisible()) {
          mApp->appFormulas->EvaluateFormulas(work->globalStatCache.globalHits.nbDesorbed);
          mApp->formulaEditor->UpdateValues();
        }
    } else if (src==editButton) {
		if( mApp->momentsEditor==nullptr ) mApp->momentsEditor = new MomentsEditor(work);
		mApp->momentsEditor->Refresh();
		mApp->momentsEditor->SetVisible(true);
		
	}
    break;
  }

  GLWindow::ProcessMessage(src,message);
}

/**
* \brief Refreshes the labels and button texts in the window
*/
void TimeSettings::RefreshMoments() {
	timeLabel->SetText("Constant Flow");
	timeId->SetText("0");
	work->displayedMoment=0;
	int nbMoments=(int)work->interfaceMomentCache.size();
	char tmp[128];
	sprintf(tmp,"%d moments",nbMoments);
	editButton->SetText(tmp);
}