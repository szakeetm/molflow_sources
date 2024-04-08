
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
