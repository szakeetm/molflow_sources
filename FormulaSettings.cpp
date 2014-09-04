/*
  File:        FormulaSettings.cpp
  Description: Formula edition dialog
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

#include "FormulaSettings.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLMessageBox.h"

FormulaSettings::FormulaSettings():GLWindow() {

  int wD = 500;
  int hD = 425;

  SetTitle("Formula Editor");

  nameL = new GLLabel("Name");
  nameL->SetBounds(5,5,50,18);
  Add(nameL);
  nameT = new GLTextField(0,"Formula name");
  nameT->SetBounds(60,5,wD-65,18);
  Add(nameT);

  exprL = new GLLabel("Expression");
  exprL->SetBounds(5,30,50,18);
  Add(exprL);
  exprT = new GLTextField(0,"");
  exprT->SetBounds(60,30,wD-65,18);
  Add(exprT);

  descL = new GLLabel(
    "MC Variables: An (Absorption on facet n), Dn (Desorption on facet n), Hn (Hit on facet n)\n"
    "                 SUMABS (total absorbed), SUMDES (total desorbed), SUMHIT (total hit)\n\n"
    "AC Variables: _An (Absorption on facet n), _Dn (Desorption on facet n), _Hn (Density on facet n)\n\n"
    "Area variables: ARn (Area of facet n), DESAR (total desorption area), ABSAR (total absorption area)\n\n"
	"Final (constant) outgassing rate [mbar*l/s]: QCONST\n"
	"Final (constant) outgassing rate [molecules/s]: QCONST_N\n"
	"Total desorbed molecules until last moment: [molecules]: NTOT\n\n"
	"Mean Pumping Path: MPP (average path of molecules in the system before absorption)\n"
	"Mean Free Path:      MFP (average path of molecules between two wall hits)\n\n"
    "Math functions: sin(), cos(), tan(), sinh(), cosh(), tanh()\n"
    "                   asin(), acos(), atan(), exp(), ln(), pow(x,y)\n"
    "                   log2(), log10(), inv(), sqrt(), abs()\n\n"
    //"Utils functions: ci95(p,N) 95% confidence interval (p=prob,N=count)\n"
    "Utils functions: sum(prefix,i,j) sum variables ex: sum(AR,1,10)=AR1+AR2+...+AR10 \n"
	"Selections:      sum(prefix,Sn) sum for facets in selection group n.\n"
	"                   Ex: SUM(AR,S2) is the total area of all facets in selection group 2\n\n"
    "Constants:  Kb (Boltzmann's constant [J.K\270\271]), R (Gas constant [J.K\270\271.mol\270\271])\n"
    "              Na (Avogadro's number [mol\270\271]), PI\n\n"
    );
  descL->SetBounds(5,55,wD-10,hD-100);
  Add(descL);

  createButton = new GLButton(0,"Create");
  createButton->SetBounds(wD-300,hD-43,95,19);
  Add(createButton);

  deleteButton = new GLButton(0,"Delete");
  deleteButton->SetEnabled(FALSE);
  deleteButton->SetBounds(wD-200,hD-43,95,19);
  Add(deleteButton);

  cancelButton = new GLButton(0,"Cancel");
  cancelButton->SetBounds(wD-100,hD-43,95,19);
  Add(cancelButton);

  // Center dialog
  int wS,hS;
  GLToolkit::GetScreenSize(&wS,&hS);
  int xD = (wS-wD)/2;
  int yD = (hS-hD)/2;
  SetBounds(xD,yD,wD,hD);

  RestoreDeviceObjects();

}

BOOL FormulaSettings::EditFormula(GLParser *f) {

  createButton->SetText("Apply");
  deleteButton->SetEnabled(TRUE);
  nameT->SetText(f->GetName());
  exprT->SetText(f->GetExpression());
  rCode = 0;
  DoModal();

  if( rCode==2 ) {
    // Delete
    return FALSE;
  } else if (rCode==1) {
    // Apply
    f->SetExpression(exprT->GetText());
    f->SetName(nameT->GetText());
    if( !f->Parse() ) 
      DisplayError(f);
  }

  return TRUE;
}

GLParser *FormulaSettings::NewFormula() {

  createButton->SetText("Create");
  deleteButton->SetEnabled(FALSE);
  nameT->SetText("Formula name");
  exprT->SetText("");
  rCode = 0;
  DoModal();

  if (rCode==1) {
    // Create
    GLParser *f = new GLParser();
    f->SetExpression(exprT->GetText());
    f->SetName(nameT->GetText());
    if( !f->Parse() )
      DisplayError(f);
    return f;
  }

  return NULL;

}

void FormulaSettings::DisplayError(GLParser *f) {

  char tmp[512];
  char tmp2[512];
  sprintf(tmp2,f->GetExpression());
  if(strlen(tmp2)) {
    int pos = f->GetCurrentPos();
    tmp2[pos] = 0;
    int ew = GLToolkit::GetDialogFont()->GetTextWidth(tmp2);
    tmp2[0] = ' ';
    tmp2[1] = 0;
    int sw = GLToolkit::GetDialogFont()->GetTextWidth(tmp2);
    int nbSpace = ew / sw;
    memset(tmp2,' ',512);
    tmp2[nbSpace]=0;
    sprintf(tmp,"%s\n%s^\n%s",f->GetExpression(),tmp2,f->GetErrorMsg());
  } else {
    GLMessageBox::Display(f->GetErrorMsg(),"Error",GLDLG_OK,GLDLG_ICONERROR);
  }

}

void FormulaSettings::ProcessMessage(GLComponent *src,int message) {

  switch(message) {
    case MSG_BUTTON:
    if(src==createButton) {
      rCode = 1;
      GLWindow::ProcessMessage(NULL,MSG_CLOSE);
    } else if(src==deleteButton) {
      rCode = 2;
      GLWindow::ProcessMessage(NULL,MSG_CLOSE);
    } else if(src==cancelButton) {
      rCode = 0;
      GLWindow::ProcessMessage(NULL,MSG_CLOSE);
    }
    break;
  }

  GLWindow::ProcessMessage(src,message);
}

