/*
  File:        FacetDetails.cpp
  Description: Facet details window
  Program:     MolFlow
  Author:      R. KERSEVAN / J-L PONS / M ADY
  Copyright:   E.S.R.F / CERN

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#include "FacetDetails.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLMessageBox.h"
#include "Utils.h"
#ifdef MOLFLOW
#include "MolFlow.h"
#endif

#ifdef SYNRAD
#include "SynRad.h"
#endif



#ifdef MOLFLOW
extern MolFlow *mApp;
#endif

#ifdef SYNRAD
extern SynRad*mApp;
#endif

typedef struct {

  char *name;
  int   width;
  int   align;
  int   timeDepColor; //Can display time-dependent value, change color accordingly

} COLUMN;

static COLUMN allColumn[] = {
  {"#"             , 40 , ALIGN_CENTER, 0} ,
  {"Sticking"      , 80 , ALIGN_CENTER, 0 } ,
  {"Opacity"       , 80 , ALIGN_CENTER, 0 } ,
  {"Structure"     , 80 , ALIGN_CENTER, 0 } ,
  {"Link"          , 40 , ALIGN_CENTER, 0 } ,
  {"Desorption"    , 80 , ALIGN_CENTER, 0 } ,
  {"Reflection"    , 80 , ALIGN_CENTER, 0 } ,
  {"2 Sided"       , 60 , ALIGN_CENTER, 0 } ,
  {"Vertex"        , 80 , ALIGN_CENTER, 0 } ,
  {"Area"          , 80 , ALIGN_CENTER, 0 } ,
  {"Temperature (K)",80 , ALIGN_CENTER, 0 } ,
  {"Facet 2D Box"  , 100, ALIGN_CENTER, 0 } ,
  {"Texture (u,v)" , 120, ALIGN_CENTER, 0 } ,
  {"Mesh sample/cm", 80 , ALIGN_CENTER, 0 } ,
  {"Count"         , 80 , ALIGN_CENTER, 0 } ,
  {"Memory"        , 80 , ALIGN_CENTER, 0 } ,
  {"Planarity"     , 80 , ALIGN_CENTER, 0 } ,
  {"Profile"       , 80 , ALIGN_CENTER, 0 } ,
  {"Imping.rate"   , 80, ALIGN_CENTER, COLOR_BLUE },
  {"Density [1/m3]", 80, ALIGN_CENTER, COLOR_BLUE },
  {"Density [kg/m3]",80, ALIGN_CENTER, COLOR_BLUE },
  {"Pressure [mbar]",80 , ALIGN_CENTER, COLOR_BLUE } ,
  {"Av.mol.speed[m/s]",80, ALIGN_CENTER, COLOR_BLUE },
  {"Hits"           , 80 , ALIGN_CENTER, COLOR_BLUE } ,
  {"Des."           , 80 , ALIGN_CENTER, COLOR_BLUE } ,
  {"Abs."           , 80 , ALIGN_CENTER, COLOR_BLUE } ,
};

static const char *desStr[] = {
  "None",
  "Uniform",
  "Cosine",
  "Cosine^"
};

static const char *refStr[] = {
  "Diffuse",
  "Mirror",
  "Uniform"
};

static const char *profStr[] = {
  "None",
  "Pressure (\201)",
  "Pressure (\202)",
  "Angular",
  "Speed distr.",
  "Ort. velocity"
};

static const char *ynStr[] = {
  "No",
  "Yes"
};


// -----------------------------------------------------------------

FacetDetails::FacetDetails():GLWindow() {

  int wD = 502;
  int hD = 400;
  worker = NULL;

  SetTitle("Facets details");
  SetIconfiable(TRUE);
  SetResizable(TRUE);
  SetMinimumSize(502,200);

  checkAllButton = new GLButton(0,"Check All");
  Add(checkAllButton);
  uncheckAllButton = new GLButton(0,"Uncheck All");
  Add(uncheckAllButton);
  updateButton = new GLButton(0,"Update");
  Add(updateButton);
  dismissButton = new GLButton(0,"Dismiss");
  Add(dismissButton);

  facetListD = new GLList(0);
  facetListD->SetColumnLabelVisible(TRUE);
  facetListD->SetGrid(TRUE);
  //facetListD->Sortable=TRUE;
  Add(facetListD);

  sPanel = new GLTitledPanel("Show column");
  sPanel->SetClosable(TRUE);
  Add(sPanel);

  show[0] = new GLToggle(0,"#");
  show[0]->SetState(TRUE);  // Always visible (not displayed)

  show[1] = new GLToggle(1,"Sticking");
  show[1]->SetState(TRUE);
  sPanel->Add(show[1]);
  show[2] = new GLToggle(2,"Opacity");
  show[2]->SetState(TRUE);
  sPanel->Add(show[2]);
  show[3] = new GLToggle(3,"Structure");
  show[3]->SetState(TRUE);
  sPanel->Add(show[3]);
  show[4] = new GLToggle(4,"Link");
  show[4]->SetState(TRUE);
  sPanel->Add(show[4]);
  show[5] = new GLToggle(5,"Desorption");
  show[5]->SetState(TRUE);
  sPanel->Add(show[5]);
  show[6] = new GLToggle(6,"Reflection");
  show[6]->SetState(TRUE);
  sPanel->Add(show[6]);
  show[7] = new GLToggle(7,"2 Sided");
  show[7]->SetState(TRUE);
  sPanel->Add(show[7]);
  show[8] = new GLToggle(8,"Vertex nb");
  show[8]->SetState(TRUE);
  sPanel->Add(show[8]);
  show[9] = new GLToggle(9,"Area");
  show[9]->SetState(TRUE);
  sPanel->Add(show[9]);
  show[10] = new GLToggle(10, "Temperature");
  show[10]->SetState(TRUE);
  sPanel->Add(show[10]);
  show[11] = new GLToggle(11,"2D Box");
  show[11]->SetState(TRUE);
  sPanel->Add(show[11]);
  show[12] = new GLToggle(12,"Texture UV");
  show[12]->SetState(TRUE);
  sPanel->Add(show[12]);
  show[13] = new GLToggle(13,"Mesh sample/cm");
  show[13]->SetState(TRUE);
  sPanel->Add(show[13]);
  show[14] = new GLToggle(14,"Count mode");
  show[14]->SetState(TRUE);
  sPanel->Add(show[14]);
  show[15] = new GLToggle(15,"Memory");
  show[15]->SetState(TRUE);
  sPanel->Add(show[15]);
  show[16] = new GLToggle(16,"Planarity");
  show[16]->SetState(TRUE);
  sPanel->Add(show[16]);
  show[17] = new GLToggle(17,"Profile");
  show[17]->SetState(TRUE);
  sPanel->Add(show[17]);
  show[18] = new GLToggle(18, "Imping.rate");
  show[18]->SetState(TRUE);
  sPanel->Add(show[18]);
  show[19] = new GLToggle(19, "Density[1/m3]");
  show[19]->SetState(TRUE);
  sPanel->Add(show[19]);
  show[20] = new GLToggle(20, "Density[kg/m3]");
  show[20]->SetState(TRUE);
  sPanel->Add(show[20]);
  show[21] = new GLToggle(21, "Pressure");
  show[21]->SetState(TRUE);
  sPanel->Add(show[21]);
  show[22] = new GLToggle(22, "Mol.speed");
  show[22]->SetState(TRUE);
  sPanel->Add(show[22]);
  show[23] = new GLToggle(23,"Hits");
  show[23]->SetState(TRUE);
  sPanel->Add(show[23]);
  show[24] = new GLToggle(24,"Des.");
  show[24]->SetState(TRUE);
  sPanel->Add(show[24]);
  show[25] = new GLToggle(25,"Abs.");
  show[25]->SetState(TRUE);
  sPanel->Add(show[25]);

  // Center dialog
  int wS,hS;
  GLToolkit::GetScreenSize(&wS,&hS);
  int xD = (wS-wD)/2;
  int yD = (hS-hD)/2;
  SetBounds(xD,yD,wD,hD);

  RestoreDeviceObjects();

}

// -----------------------------------------------------------------

void FacetDetails::PlaceComponents() {

  // Show toggle panel
  int nbW = (width-20)/80;
  int nbL = (NB_FDCOLUMN-1)/nbW + (((NB_FDCOLUMN-1)%nbW)?1:0);
  int hp;
  if(!sPanel->IsClosed())
    hp = 20*(nbL+1);
  else
    hp = 20;
  sPanel->SetBounds(5,height-(hp+52),width-10,hp);
  for(int i=0;i<NB_FDCOLUMN;i++)
    sPanel->SetCompBounds(show[i],5+80*((i-1)%nbW),18+20*((i-1)/nbW),85,19);

  facetListD->SetBounds(5,5,width-10,height-(62+hp));

  checkAllButton->SetBounds(5,height-45,90,19);
  uncheckAllButton->SetBounds(100,height-45,90,19);
  updateButton->SetBounds(195,height-45,90,19);
  dismissButton->SetBounds(width-100,height-45,90,19);

}

// -----------------------------------------------------------------

void FacetDetails::SetBounds(int x,int y,int w,int h) {

  GLWindow::SetBounds(x,y,w,h);
  PlaceComponents();

}

// -----------------------------------------------------------------

char *FacetDetails::GetCountStr(Facet *f) {
  static char ret[128];
  strcpy(ret,"");
  if(f->sh.countDes) strcat(ret,"DES");
  if(f->sh.countAbs) if(strlen(ret)==0) strcat(ret,"ABS"); else strcat(ret,"+ABS");
  if(f->sh.countRefl) if(strlen(ret)==0) strcat(ret,"REFL"); else strcat(ret,"+REFL");
  if(f->sh.countTrans) if(strlen(ret)==0) strcat(ret,"TRANS"); else strcat(ret,"+TRANS");
  return ret;
}

// -----------------------------------------------------------------

char *FacetDetails::FormatCell(int idx,Facet *f,int mode) {
  static char ret[256];
  strcpy(ret,"");

  switch(mode) {
    case 0:
      sprintf(ret,"%d",idx+1);
      break;
    case 1:
      sprintf(ret,"%g",f->sh.sticking);
      break;
    case 2:
      sprintf(ret,"%g",f->sh.opacity);
      break;
    case 3:
      sprintf(ret,"%d",f->sh.superIdx+1);
      break;
    case 4:
      sprintf(ret,"%d",f->sh.superDest);
      break;
    case 5:
      sprintf(ret,"%s",desStr[f->sh.desorbType]);
	  if (f->sh.desorbType==DES_COSINE_N) sprintf(ret,"%s%g",ret,f->sh.desorbTypeN); //append exponent
      break;
    case 6:
      sprintf(ret,"%s",refStr[f->sh.reflectType]);
      break;
    case 7:
      sprintf(ret,"%s",ynStr[f->sh.is2sided]);      
      break;
    case 8:
      sprintf(ret,"%d",f->sh.nbIndex);
      break;
    case 9:
		if (f->sh.is2sided) sprintf(ret,"2*%g",f->sh.area);
		else sprintf(ret,"%g",f->sh.area);
      break;
	case 10:
		sprintf(ret, "%g", f->sh.temperature);
		break;
    case 11:
      sprintf(ret,"%g x %g",Norme(f->sh.U),Norme(f->sh.V));
      break;
    case 12:
      if( f->sh.isTextured ) {
        sprintf(ret,"%dx%d (%g x %g)",f->sh.texWidth,f->sh.texHeight,f->sh.texWidthD,f->sh.texHeightD);
      } else {
        sprintf(ret,"None");
      }
      break;
    case 13:
      sprintf(ret,"%g",f->tRatio);
      break;
    case 14:
      sprintf(ret,"%s",GetCountStr(f));
      break;
    case 15:
		sprintf(ret,"%s",FormatMemory(f->GetTexRamSize(1+worker->moments.size())));
      break;
    case 16:
      sprintf(ret,"%f",f->err);
      break;
    case 17:
		sprintf(ret,"%s",profStr[f->sh.profileType]);
		break;
	case 18: //imp.rate
	{
	double dCoef = 1E4 * worker->GetMoleculesPerTP();  //1E4 is conversion from m2 to cm2; 0.01 is Pa->mbar
	sprintf(ret, "%g", f->counterCache.hit.nbHit / f->GetArea()*dCoef);
	//11.77=sqrt(8*8.31*293.15/3.14/0.028)/4/10
	break; }
	case 19: //particle density
	{
	double dCoef = 1E4 * worker->GetMoleculesPerTP();  //1E4 is conversion from m2 to cm2; 0.01 is Pa->mbar	
	//Correction for double-density effect (measuring density on desorbing/absorbing facets):
	if (f->counterCache.hit.nbHit>0 || f->counterCache.hit.nbDesorbed>0)
		if (f->counterCache.hit.nbAbsorbed >0 || f->counterCache.hit.nbDesorbed>0) //otherwise save calculation time
		dCoef *= 1.0 - ((double)f->counterCache.hit.nbAbsorbed + (double)f->counterCache.hit.nbDesorbed) / ((double)f->counterCache.hit.nbHit+(double)f->counterCache.hit.nbDesorbed) / 2.0;
	
	sprintf(ret, "%g", f->counterCache.hit.sum_1_per_ort_velocity / f->GetArea()*dCoef);

	break; }
	case 20: //gas density
	{
	double dCoef =  1E4 * worker->GetMoleculesPerTP();  //1E4 is conversion from m2 to cm2; 0.01 is Pa->mbar
	
	//Correction for double-density effect (measuring density on desorbing/absorbing facets):
	if (f->counterCache.hit.nbHit>0 || f->counterCache.hit.nbDesorbed>0)
		if (f->counterCache.hit.nbAbsorbed >0 || f->counterCache.hit.nbDesorbed>0) //otherwise save calculation time
		dCoef *= 1.0 - ((double)f->counterCache.hit.nbAbsorbed + (double)f->counterCache.hit.nbDesorbed) / ((double)f->counterCache.hit.nbHit + (double)f->counterCache.hit.nbDesorbed) / 2.0;
	
	sprintf(ret, "%g", f->counterCache.hit.sum_1_per_ort_velocity / f->GetArea()*dCoef*mApp->worker.gasMass / 1000.0 / 6E23);
	break; }
	case 21: //avg.pressure
	{
	double dCoef = 1E4 * worker->GetMoleculesPerTP() * (worker->gasMass / 1000 / 6E23) * 0.0100;  //1E4 is conversion from m2 to cm2; 0.01 is Pa->mbar
	
	sprintf(ret, "%g", f->counterCache.hit.sum_v_ort*dCoef / f->GetArea());
	break; }
	case 22: //avg. gas speed (estimate)
		sprintf(ret, "%g", 4.0*(double)(f->counterCache.hit.nbHit+f->counterCache.hit.nbDesorbed) / f->counterCache.hit.sum_1_per_ort_velocity);
		//<v_surf>=2*<v_surf_ort>
		//<v_gas>=1/<1/v_surf>
		break;
	case 23:
		sprintf(ret,"%I64d",f->counterCache.hit.nbHit);
		break;
	case 24:
		sprintf(ret,"%I64d",f->counterCache.hit.nbDesorbed);
		break;
	case 25:
		sprintf(ret,"%I64d",f->counterCache.hit.nbAbsorbed);
		break;
  }

  return ret;

}

// -----------------------------------------------------------------

void FacetDetails::UpdateTable() {

  Geometry *s = worker->GetGeometry();
  int nbFacet = s->GetNbFacet();
  int nbS = s->GetNbSelected();
  static char ret[256];
  strcpy(ret,"");

  /*
  //SUM Counters
  double sumArea=0;
  DWORD sumMemory=0;
  */

  char *tmpName[NB_FDCOLUMN];
  int  tmpWidth[NB_FDCOLUMN];
  int  tmpAlign[NB_FDCOLUMN];
  int  tmpColor[NB_FDCOLUMN];

  int nbCol = 0;

  for(int i=0;i<NB_FDCOLUMN;i++) {
    if(i==0 || show[i]->GetState()) {
      tmpName[nbCol]  = allColumn[i].name;
      tmpWidth[nbCol] = allColumn[i].width;
      tmpAlign[nbCol] = allColumn[i].align;
	  tmpColor[nbCol] = allColumn[i].timeDepColor;
      shown[nbCol] = i;
      nbCol++;
    }
  }

  facetListD->SetSize(nbCol,nbS);
  facetListD->SetColumnWidths(tmpWidth);
  facetListD->SetColumnLabels(tmpName);
  facetListD->SetColumnAligns(tmpAlign);
  
  if (worker->displayedMoment == 0)
	  facetListD->SetAllColumnColors(COLOR_BLACK);
  else 
	facetListD->SetColumnColors(tmpColor);
  

  nbS = 0;
  for(int i=0;i<nbFacet;i++) {
    Facet *f = s->GetFacet(i);
    if(f->selected) {
      for(int j=0;j<nbCol;j++)
        facetListD->SetValueAt(j,nbS,FormatCell(i,f,shown[j]));
		
	  /*
		sumArea+=f->sh.area;
		sumMemory+=f->GetTexRamSize();
      */
		nbS++;
		}
	
	}
	/*
  //SUM values
  for(int j=0;j<nbCol;j++) {
	  switch (shown[j]) {
		case 0:
			facetListD->SetValueAt(j,nbS,"SUM");
			break;
		case 9:
			sprintf(ret,"%g",sumArea);
			facetListD->SetValueAt(j,nbS,ret);
			break;
		case 14:
			sprintf(ret,"%s",FormatMemory(sumMemory));
			facetListD->SetValueAt(j,nbS,ret);
			break;
		default:
			facetListD->SetValueAt(j,nbS,"");
	  }
	  
  }
  */
}

// -----------------------------------------------------------------

void FacetDetails::Update() {

  if(!worker) return;
  if(!IsVisible()) return;

  Geometry *s = worker->GetGeometry();
  int nbS = s->GetNbSelected();
  
  if(nbS==0) {
    facetListD->Clear();
    return;
  }

  UpdateTable();

}

// -----------------------------------------------------------------

void FacetDetails::Display(Worker *w) {

  worker = w;
  SetVisible(TRUE);
  Update();

}

// -----------------------------------------------------------------

void FacetDetails::ProcessMessage(GLComponent *src,int message) {

  switch(message) {

    case MSG_BUTTON:
      if(src==dismissButton) {
        SetVisible(FALSE);
      } else if (src==checkAllButton) {
        for(int i=0;i<NB_FDCOLUMN;i++) show[i]->SetState(TRUE);
        UpdateTable();
      } else if (src==uncheckAllButton) {
        for(int i=0;i<NB_FDCOLUMN;i++) show[i]->SetState(FALSE);
        UpdateTable();
	  }	else if (src==updateButton) {
        UpdateTable();
      }
      break;

    case MSG_TOGGLE:
      UpdateTable();
      break;

    case MSG_LIST_COL:
      if( src==facetListD ) {
        // Save column width
        int c = facetListD->GetDraggedCol();
        allColumn[shown[c]].width = facetListD->GetColWidth(c);
      }
      break;

    case MSG_PANELR:
      PlaceComponents();
      break;

  }

  GLWindow::ProcessMessage(src,message);

}

