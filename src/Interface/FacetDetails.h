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
#pragma once

#include "GLApp/GLWindow.h"
#include "GLApp/GLList.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLToggle.h"
#include "GLApp/GLTitledPanel.h"
#include "Worker.h"

class InterfaceFacet;

struct ColumnData {

	std::string colName;
	int   width;
	int   align;
	int   timeDepColor; //Can display time-dependent value, change color accordingly

};

class FacetDetails : public GLWindow {

public:

	// Construction
	FacetDetails();

	// Component method
	void Display(Worker* w);
	void Update();

	// Implementation
	void ProcessMessage(GLComponent* src, int message) override;
	void SetBounds(int x, int y, int w, int h) override;

private:

	std::string GetCountStr(InterfaceFacet* f);
	void UpdateTable();
	std::string FormatCell(size_t idx, InterfaceFacet* f, size_t mode);
	void PlaceComponents();

	Worker* worker;
	GLList* facetListD;

	GLTitledPanel* sPanel;          // Toggle panel
	std::vector<std::unique_ptr<GLToggle>> colShowToggle;
	std::vector<size_t> shownColIds;
	std::vector<ColumnData> columnData;

	GLButton* checkAllButton;
	GLButton* uncheckAllButton;
	GLButton* dismissButton;
	GLButton* updateButton;

};
