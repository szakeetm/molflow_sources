/*
Program:     MolFlow+ / Synrad+
Description: Monte Carlo simulator for ultra-high vacuum and synchrotron radiation
Authors:     Jean-Luc PONS / Roberto KERSEVAN / Marton ADY / Pascal BAEHR
Copyright:   E.S.R.F / CERN
Website:     https://cern.ch/molflow

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

Full license text: https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html
*/
#include "FacetDetails.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLMessageBox.h"
#include "Helper/MathTools.h"
#include "Geometry_shared.h"
#include "Facet_shared.h"

#include "MolFlow.h"
extern MolFlow *mApp;

typedef struct {

  const char *name;
  int   width;
  int   align;
  int   timeDepColor; //Can display time-dependent value, change color accordingly

} COLUMN;

/**
* \brief Values for the table headers for when a facet is opened in the window
*/
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
  {"MC Hits"           , 80 , ALIGN_CENTER, COLOR_BLUE } ,
  {"Equiv.hits"           , 80 , ALIGN_CENTER, COLOR_BLUE } ,
  {"Des."           , 80 , ALIGN_CENTER, COLOR_BLUE } ,
  {"Equiv.abs."           , 80 , ALIGN_CENTER, COLOR_BLUE } ,
  {"Force"           , 240 , ALIGN_CENTER, COLOR_BLUE } ,
  {"Force^2"           , 240 , ALIGN_CENTER, COLOR_BLUE } ,
  {"Torque"           , 240 , ALIGN_CENTER, COLOR_BLUE } ,
};

static const char *desStr[] = {
  "None",
  "Uniform",
  "Cosine",
  "Cosine^"
};

static const char *profStr[] = {
  "None",
  "Pressure (\201)",
  "Pressure (\202)",
  "Angular",
  "Speed distr.",
  "Ort. velocity",
  "Tan. velocity"
};

static const char *ynStr[] = {
  "No",
  "Yes"
};

/**
* \brief Constructor for the FacetDetails window with default initialisation
*/
FacetDetails::FacetDetails():GLWindow() {

  int wD = 502;
  int hD = 400;
  worker = NULL;

  SetTitle("Facets details");
  SetIconfiable(true);
  SetResizable(true);
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
  facetListD->SetColumnLabelVisible(true);
  facetListD->SetGrid(true);
  //facetListD->Sortable=true;
  Add(facetListD);

  sPanel = new GLTitledPanel("Show column");
  sPanel->SetClosable(true);
  Add(sPanel);

  size_t index = 0;

  show[index] = new GLToggle(index,"#");
  show[index]->SetState(true);  // Always visible (not displayed)

  show[++index] = new GLToggle(index,"Sticking");
  show[index]->SetState(true);
  sPanel->Add(show[index]);
  show[++index] = new GLToggle(index,"Opacity");
  show[index]->SetState(true);
  sPanel->Add(show[index]);
  show[++index] = new GLToggle(index,"Structure");
  show[index]->SetState(true);
  sPanel->Add(show[index]);
  show[++index] = new GLToggle(index,"Link");
  show[index]->SetState(true);
  sPanel->Add(show[index]);
  show[++index] = new GLToggle(index,"Desorption");
  show[index]->SetState(true);
  sPanel->Add(show[index]);
  show[++index] = new GLToggle(index,"Reflection");
  show[index]->SetState(true);
  sPanel->Add(show[index]);
  show[++index] = new GLToggle(index,"2 Sided");
  show[index]->SetState(true);
  sPanel->Add(show[index]);
  show[++index] = new GLToggle(index,"Vertex nb");
  show[index]->SetState(true);
  sPanel->Add(show[index]);
  show[++index] = new GLToggle(index,"Area");
  show[index]->SetState(true);
  sPanel->Add(show[index]);
  show[++index] = new GLToggle(index, "Temperature");
  show[index]->SetState(true);
  sPanel->Add(show[index]);
  show[++index] = new GLToggle(index,"2D Box");
  show[index]->SetState(true);
  sPanel->Add(show[index]);
  show[++index] = new GLToggle(index,"Texture UV");
  show[index]->SetState(true);
  sPanel->Add(show[index]);
  show[++index] = new GLToggle(index,"Mesh sample/cm");
  show[index]->SetState(true);
  sPanel->Add(show[index]);
  show[++index] = new GLToggle(index,"Count mode");
  show[index]->SetState(true);
  sPanel->Add(show[index]);
  show[++index] = new GLToggle(index,"Memory");
  show[index]->SetState(true);
  sPanel->Add(show[index]);
  show[++index] = new GLToggle(index,"Planarity");
  show[index]->SetState(true);
  sPanel->Add(show[index]);
  show[++index] = new GLToggle(index,"Profile");
  show[index]->SetState(true);
  sPanel->Add(show[index]);
  show[++index] = new GLToggle(index, "Imping.rate");
  show[index]->SetState(true);
  sPanel->Add(show[index]);
  show[++index] = new GLToggle(index, "Density[1/m3]");
  show[index]->SetState(true);
  sPanel->Add(show[index]);
  show[++index] = new GLToggle(index, "Density[kg/m3]");
  show[index]->SetState(true);
  sPanel->Add(show[index]);
  show[++index] = new GLToggle(index, "Pressure");
  show[index]->SetState(true);
  sPanel->Add(show[index]);
  show[++index] = new GLToggle(index, "Mol.speed");
  show[index]->SetState(true);
  sPanel->Add(show[index]);
  show[++index] = new GLToggle(index,"MC Hits");
  show[index]->SetState(true);
  sPanel->Add(show[index]);
  show[++index] = new GLToggle(index, "Equiv.hits");
  show[index]->SetState(true);
  sPanel->Add(show[index]);
  show[++index] = new GLToggle(index, "Des.");
  show[index]->SetState(true);
  sPanel->Add(show[index]);
  show[++index] = new GLToggle(index,"Equiv.Abs.");
  show[index]->SetState(true);
  sPanel->Add(show[index]);
  show[++index] = new GLToggle(index, "Force");
  show[index]->SetState(true);
  sPanel->Add(show[index]);
  show[++index] = new GLToggle(index, "Force^2");
  show[index]->SetState(true);
  sPanel->Add(show[index]);
  show[++index] = new GLToggle(index, "Torque");
  show[index]->SetState(true);
  sPanel->Add(show[index]);

  // Center dialog
  int wS,hS;
  GLToolkit::GetScreenSize(&wS,&hS);
  int xD = (wS-wD)/2;
  int yD = (hS-hD)/2;
  SetBounds(xD,yD,wD,hD);

  RestoreDeviceObjects();

}

/**
* \brief Places the various components inside the window
*/
void FacetDetails::PlaceComponents() {

  // Show toggle panel
  int nbW = (_width - 20) / 80;
  int nbL = (NB_FDCOLUMN-1)/nbW + (((NB_FDCOLUMN-1)%nbW)?1:0);
  int hp;
  if(!sPanel->IsClosed())
    hp = 20*(nbL+1);
  else
    hp = 20;
  sPanel->SetBounds(5, _height - (hp + 52), _width - 10, hp);
  for(size_t i=0;i<NB_FDCOLUMN;i++)
    sPanel->SetCompBounds(show[i],5+80*(((int)i-1)%nbW),18+20*(((int)i-1)/nbW),85,19);

  facetListD->SetBounds(5, 5, _width - 10, _height - (62 + hp));

  checkAllButton->SetBounds(5, _height - 45, 90, 19);
  uncheckAllButton->SetBounds(100, _height - 45, 90, 19);
  updateButton->SetBounds(195, _height - 45, 90, 19);
  dismissButton->SetBounds(_width - 100, _height - 45, 90, 19);

}

/**
* \brief To set placement of window and its size and what's inside
* \param x x position of window
* \param y y position of window
* \param w width of window
* \param h height of window
*/
void FacetDetails::SetBounds(int x,int y,int w,int h) {

  GLWindow::SetBounds(x,y,w,h);
  PlaceComponents();

}

/**
* \brief Gives a string which counts values corresponding to the facet settings
* \param f Pointer to a facet
* \return char pointer taking a string with the count value(s)
*/
char *FacetDetails::GetCountStr(InterfaceFacet *f) {
  static char ret[128];
  strcpy(ret,"");
  if(f->sh.countDes) strcat(ret,"DES");
  if(f->sh.countAbs) { if (strlen(ret) == 0) { strcat(ret, "ABS"); } else { strcat(ret, "+ABS"); }}
  if(f->sh.countRefl) { if (strlen(ret) == 0) { strcat(ret, "REFL"); } else { strcat(ret, "+REFL"); }}
  if(f->sh.countTrans) { if (strlen(ret) == 0) { strcat(ret, "TRANS"); } else { strcat(ret, "+TRANS"); }}
  return ret;
}

/**
* \brief Prints table values inside the corresponding cell
* \param idx Facet ID (local for table)
* \param f Pointer to a facet
* \param mode which kind of value has to be evaluated and printed
* \return char pointer taking a string with the count value(s)
*/
char *FacetDetails::FormatCell(size_t idx, InterfaceFacet *f, size_t mode) {
  static char ret[512];
  strcpy(ret,"");

  switch(mode) {
    case 0: //index
      sprintf(ret,"%zd",idx+1);
      break;
    case 1: //sticking factor
      sprintf(ret,"%g",f->sh.sticking);
      break;
    case 2: //opacity
      sprintf(ret,"%g",f->sh.opacity);
      break;
    case 3: //Structure
	{
		std::ostringstream out;
		if (f->sh.superIdx == -1) out << "All";
		else out << (f->sh.superIdx + 1);
		sprintf(ret, "%s", out.str().c_str());
		break;
	}
    case 4: //Link destination
      sprintf(ret,"%zd",f->sh.superDest);
      break;
    case 5: //Desorption type
	  if (f->sh.desorbType == DES_COSINE_N)
	  {
		  sprintf(ret, "%s%g", desStr[f->sh.desorbType], f->sh.desorbTypeN); //append exponent
	  }
	  else
	  {
		  sprintf(ret, "%s", desStr[f->sh.desorbType]);
	  }
      break;
    case 6: //Reflection type
      sprintf(ret,"%g diff. %g spec. %g cos^%g",f->sh.reflection.diffusePart,f->sh.reflection.specularPart,1.0-f->sh.reflection.diffusePart-f->sh.reflection.specularPart,f->sh.reflection.cosineExponent);
      break;
    case 7: //2-sided
      sprintf(ret,"%s",ynStr[f->sh.is2sided]);      
      break;
    case 8: //Nb of vertex
      sprintf(ret,"%zd",f->sh.nbIndex);
      break;
    case 9: //Area
		if (f->sh.is2sided) sprintf(ret,"2*%g",f->sh.area);
		else sprintf(ret,"%g",f->sh.area);
      break;
	case 10: //Temperature
		sprintf(ret, "%g", f->sh.temperature);
		break;
    case 11: //2D box
      sprintf(ret,"%g x %g",f->sh.U.Norme(),f->sh.V.Norme());
      break;
    case 12: //Texture type
      if( f->sh.isTextured ) {
        sprintf(ret,"%zdx%zd (%g x %g)",f->sh.texWidth,f->sh.texHeight,f->sh.texWidth_precise,f->sh.texHeight_precise);
      } else {
        sprintf(ret,"None");
      }
      break;
    case 13: //Texture sample/cm
        if(IsEqual(f->tRatioU,f->tRatioV))
            sprintf(ret,"%g",f->tRatioU);
        else
            sprintf(ret,"%g x %g",f->tRatioU,f->tRatioV);
      break;
    case 14: //Texture record type
      sprintf(ret,"%s",GetCountStr(f));
      break;
    case 15: //Texture memory
		sprintf(ret,"%s",FormatMemory(f->GetTexRamSize(1+worker->moments.size())));
      break;
    case 16: //Planarity
      sprintf(ret,"%f",f->planarityError);
      break;
    case 17: //Profile type
		sprintf(ret,"%s",profStr[f->sh.profileType]);
		break;
	case 18: //imp.rate
	{
	double dCoef = 1E4 * worker->GetMoleculesPerTP(worker->displayedMoment);  //1E4 is conversion from m2 to cm2; 0.01 is Pa->mbar
	sprintf(ret, "%g", f->facetHitCache.nbHitEquiv / f->GetArea()*dCoef);
	//11.77=sqrt(8*8.31*293.15/3.14/0.028)/4/10
	break; }
	case 19: //particle density
	{
	double dCoef = 1E4 * worker->GetMoleculesPerTP(worker->displayedMoment)*f->DensityCorrection();  //1E4 is conversion from m2 to cm2; 0.01 is Pa->mbar	
	
	sprintf(ret, "%g", f->facetHitCache.sum_1_per_ort_velocity / f->GetArea()*dCoef);

	break; }
	case 20: //gas density
	{
	double dCoef =  1E4 * worker->GetMoleculesPerTP(worker->displayedMoment)*f->DensityCorrection();  //1E4 is conversion from m2 to cm2; 0.01 is Pa->mbar
	
	sprintf(ret, "%g", f->facetHitCache.sum_1_per_ort_velocity / f->GetArea()*dCoef*mApp->worker.model->wp.gasMass / 1000.0 / 6E23);
	break; }
	case 21: //avg.pressure
	{
	double dCoef = 1E4 * worker->GetMoleculesPerTP(worker->displayedMoment) * (worker->model->wp.gasMass / 1000 / 6E23) * 0.0100;  //1E4 is conversion from m2 to cm2; 0.01 is Pa->mbar
	
	sprintf(ret, "%g", f->facetHitCache.sum_v_ort*dCoef / f->GetArea());
	break; }
	case 22: //avg. gas speed (estimate)
		/*sprintf(ret, "%g", 4.0*(double)(f->facetHitCache.hit.nbMCHit+f->facetHitCache.hit.nbDesorbed) / f->facetHitCache.hit.sum_1_per_ort_velocity);*/
		sprintf(ret, "%g", (f->facetHitCache.nbHitEquiv + static_cast<double>(f->facetHitCache.nbDesorbed)) / f->facetHitCache.sum_1_per_velocity);
		//<v_surf>=2*<v_surf_ort>
		//<v_gas>=1/<1/v_surf>
		break;
	case 23: //MC Hits
		sprintf(ret,"%zd",f->facetHitCache.nbMCHit);
		break;
	case 24: //Equiv. hits (low-flux)
		sprintf(ret, "%g", f->facetHitCache.nbHitEquiv);
		break;
	case 25: //Des Abs.
		sprintf(ret,"%zd",f->facetHitCache.nbDesorbed);
		break;
	case 26: //MC Abs.
		sprintf(ret,"%g",f->facetHitCache.nbAbsEquiv);
		break;
    case 27: //Force
    {
        auto force = f->facetHitCache.impulse * worker->GetMoleculesPerTP(worker->displayedMoment) * (worker->model->wp.gasMass / 1000 / 6E23);
        strcpy(ret, fmt::format("{:.4g} N ({:.4g},{:.4g},{:.4g})", force.Norme(), force.x, force.y, force.z).c_str());
        break;
    }
    case 28: //Force^2
    {
        auto force_sqr = f->facetHitCache.impulse_square * worker->GetMoleculesPerTP(worker->displayedMoment) * Sqr(worker->model->wp.gasMass / 1000 / 6E23);
        strcpy(ret, fmt::format("{:.4g} N^2 ({:.4g},{:.4g},{:.4g})", force_sqr.Norme(), force_sqr.x, force_sqr.y, force_sqr.z).c_str());
        break;
    }
    case 29: //Torque
    {
        auto torque = f->facetHitCache.impulse_momentum * worker->GetMoleculesPerTP(worker->displayedMoment) * (worker->model->wp.gasMass / 1000 / 6E23) * 0.01; //0.01: N*cm to Nm
        strcpy(ret, fmt::format("{:.4g} Nm ({:.4g},{:.4g},{:.4g})", torque.Norme(), torque.x, torque.y, torque.z).c_str());
        break;
    }
  }

  return ret;

}

/**
* \brief Prints table with header values and facet values
*/
void FacetDetails::UpdateTable() {

  Geometry *geom = worker->GetGeometry();
  auto selectedFacets = geom->GetSelectedFacets();
  static char ret[256];
  strcpy(ret,"");

  const char *tmpName[NB_FDCOLUMN];
  int  tmpWidth[NB_FDCOLUMN];
  int  tmpAlign[NB_FDCOLUMN];
  int  tmpColor[NB_FDCOLUMN];

  size_t nbCol = 0;

  for(size_t i=0;i<NB_FDCOLUMN;i++) {
    if(i==0 || show[i]->GetState()) {
      tmpName[nbCol]  = allColumn[i].name;
      tmpWidth[nbCol] = allColumn[i].width;
      tmpAlign[nbCol] = allColumn[i].align;
	  tmpColor[nbCol] = allColumn[i].timeDepColor;
      shown[nbCol] = i;
      nbCol++;
    }
  }

  facetListD->SetSize(nbCol,selectedFacets.size());
  facetListD->SetColumnWidths(tmpWidth);
  facetListD->SetColumnLabels(tmpName);
  facetListD->SetColumnAligns(tmpAlign);
  
  if (worker->displayedMoment == 0)
	  facetListD->SetAllColumnColors(COLOR_BLACK);
  else 
	facetListD->SetColumnColors(tmpColor);
  

  size_t nbS = 0;
  for(auto& sel:selectedFacets) {
    InterfaceFacet *f = geom->GetFacet(sel);
    for(size_t j=0;j<nbCol;j++)
        facetListD->SetValueAt(j,nbS,FormatCell(sel,f,shown[j]));
	nbS++;
	}
}

/**
* \brief Initial update of the table if it should be displayed
*/
void FacetDetails::Update() {

  if(!worker) return;
  if(!IsVisible()) return;

  Geometry *s = worker->GetGeometry();
  size_t nbS = s->GetNbSelectedFacets();
  
  if(nbS==0) {
    facetListD->Clear();
    return;
  }

  UpdateTable();

}

/**
* \brief Initial display function that renders the table
* \param w Worker for this task
*/
void FacetDetails::Display(Worker *w) {

  worker = w;
  SetVisible(true);
  Update();

}

/**
* \brief Processes events like button clicks for the advanced facet details window.
* \param src the component that got used to call this event
* \param message the type that triggered change (button, text change etc.)
*/
void FacetDetails::ProcessMessage(GLComponent *src,int message) {

  switch(message) {

    case MSG_BUTTON:
      if(src==dismissButton) {
        SetVisible(false);
      } else if (src==checkAllButton) {
        for(size_t i=0;i<NB_FDCOLUMN;i++) show[i]->SetState(true);
        UpdateTable();
      } else if (src==uncheckAllButton) {
        for(size_t i=0;i<NB_FDCOLUMN;i++) show[i]->SetState(false);
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

