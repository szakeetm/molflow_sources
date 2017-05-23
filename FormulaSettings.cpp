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
#include "GLApp/GLLabel.h"
#include "GLApp/GLTextField.h"
#include "GLApp/GLButton.h"

FormulaSettings::FormulaSettings():GLWindow() {

	formulaId = -1;

  int wD = 500;
  int hD = 500;

  SetTitle("Formula Editor");

  exprL = new GLLabel("Expression");
  exprL->SetBounds(5,5,75,18);
  Add(exprL);
  exprT = new GLTextField(0,"");
  exprT->SetBounds(85,5,wD-95,18);
  Add(exprT);

  nameL = new GLLabel("Name (optional)");
  nameL->SetBounds(5,30,75,18);
  Add(nameL);
  nameT = new GLTextField(0,"");
  nameT->SetBounds(85,30,wD-95,18);
  Add(nameT);

  descL = new GLLabel(
    R"(MC Variables: An (Absorption on facet n), Dn (Desorption on facet n), Hn (Hit on facet n)
Pn (Pressure [mbar] on facet n), DENn (Density [1/m3] on facet n)
Zn (Imp. rate on facet n), Vn (avg. speed [m/s] on facet n), Tn (temp[K] of facet n)
SUMABS (total absorbed), SUMDES (total desorbed), SUMHIT (total hit)

SUM(H,3,8)    calculates the sum of hits on facets 3,4,... ...7,8. (Works with H,A,D,AR).
AVG(P,3,8)    calculates the average pressure (area-wise) on facets 3 to 8 (Works with P,D,Z)
SUM(H,S3)    calculates the sum of hits on selection group 3 (works with H,A,D,AR)
AVG(DEN,S2) calculates the average (area-wise) on facets belonging to sel. group 2
SUM(H,SEL)    calculates the sum of hits on the current selection. (Works with H,A,D,AR)
AVG(Z,SEL)    calculates the average impingement rate on the current selection
For the last two, might need to manually refresh formulas after you change the selection.

Area variables: ARn (Area of facet n), DESAR (total desorption area), ABSAR (total absorption area)

Final (constant) outgassing rate [mbar*l/s]: QCONST
Final (constant) outgassing rate [molecules/s]: QCONST_N
Total desorbed molecules until last moment: [molecules]: NTOT
Gas mass [g/mol]: GASMASS

Mean Pumping Path: MPP (average path of molecules in the system before absorption)
Mean Free Path:      MFP (average path of molecules between two wall hits)

Math functions: sin(), cos(), tan(), sinh(), cosh(), tanh(), asin(), acos(),
                     atan(), exp(), ln(), pow(x,y), log2(), log10(), inv(), sqrt(), abs()

Constants:  Kb (Boltzmann's constant), R (Gas constant), Na (Avogadro's number), PI
)");
  descL->SetBounds(5,55,wD-10,hD-100);
  Add(descL);

  applyButton = new GLButton(0,"Create");
  applyButton->SetBounds(wD-300,hD-43,95,19);
  Add(applyButton);

  deleteButton = new GLButton(0,"Delete");
  deleteButton->SetEnabled(false);
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