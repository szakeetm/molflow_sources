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
#include "TexturePlotter.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLMessageBox.h"
//#include "GLApp/GLFileBox.h"
#include "NativeFileDialog/molflow_wrapper/nfd_wrapper.h"
#include "GLApp/GLLabel.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLList.h"
#include "GLApp/GLCombo.h"
#include "GLApp/GLToggle.h"
#include "Geometry_shared.h"
#include "Facet_shared.h"
#if defined(MOLFLOW)
#include "MolFlow.h"
#endif

#if defined(SYNRAD)
#include "SynRad.h"
#endif

#if defined(MOLFLOW)
extern MolFlow* mApp;
#endif

#if defined(SYNRAD)
extern SynRad* mApp;
#endif

std::string fileFilters = "txt";

/**
* \brief Constructor with initialisation for Texture plotter window (Tools/Texture Plotter)
*/
TexturePlotter::TexturePlotter() :GLWindow() {

	int wD = 500;
	int hD = 300;
	lastUpdate = 0.0f;

	SetTitle("Texture plotter");
	SetResizable(true);
	SetIconfiable(true);
	SetMinimumSize(wD, hD);

	mapList = new GLList(0);
	mapList->SetColumnLabelVisible(true);
	mapList->SetRowLabelVisible(true);
	mapList->SetAutoColumnLabel(true);
	mapList->SetAutoRowLabel(true);
	mapList->SetRowLabelMargin(20);
	mapList->SetGrid(true);
	mapList->SetSelectionMode(BOX_CELL);
	mapList->SetCornerLabel("\202\\\201");
	Add(mapList);

	viewLabel = new GLLabel("View:");
	Add(viewLabel);
	viewCombo = new GLCombo(0);
	viewCombo->SetSize(9);
	viewCombo->SetValueAt(0, "Cell area (cm\262)");
	viewCombo->SetValueAt(1, "# of MC hits");
	viewCombo->SetValueAt(2, "Impingement rate [1/m\262/sec]");
	viewCombo->SetValueAt(3, "Particle density [1/m3]");
	viewCombo->SetValueAt(4, "Gas density [kg/m3]");
	viewCombo->SetValueAt(5, "Pressure [mbar]");
	viewCombo->SetValueAt(6, "Avg.speed estimate[m/s]");
	viewCombo->SetValueAt(7, "Incident velocity vector[m/s]");
	viewCombo->SetValueAt(8, "# of velocity vectors");

	viewCombo->SetSelectedIndex(5); //Pressure by default
	Add(viewCombo);

	saveButton = new GLButton(0, "Save");
	Add(saveButton);

	sizeButton = new GLButton(0, "Autosize");
	Add(sizeButton);

	maxButton = new GLButton(0, "Find Max.");
	Add(maxButton);

	cancelButton = new GLButton(0, "Dismiss");
	Add(cancelButton);

	autoSizeOnUpdate = new GLToggle(0, "Autosize on every update (disable for smooth scrolling)");
	autoSizeOnUpdate->SetState(true);
	Add(autoSizeOnUpdate);

	// Center dialog
	int wS, hS;
	GLToolkit::GetScreenSize(&wS, &hS);
	int xD = (wS - wD) / 2;
	int yD = (hS - hD) / 2;
	SetBounds(xD, yD, wD, hD);

	RestoreDeviceObjects();

	worker = NULL;

}

/**
* \brief Places all components (buttons, text etc.) at the right position inside the window
*/
void TexturePlotter::PlaceComponents() {

	mapList->SetBounds(5, 5, _width - 15, _height - 80);
	saveButton->SetBounds(10, _height - 70, 70, 19);
	sizeButton->SetBounds(10, _height - 45, 70, 19);
	autoSizeOnUpdate->SetBounds(90, _height - 45, 120, 19);
	maxButton->SetBounds(90, _height - 70, 70, 19);
	viewLabel->SetBounds(320, _height - 70, 30, 19);
	viewCombo->SetBounds(350, _height - 70, 130, 19);
	cancelButton->SetBounds(_width - 90, _height - 45, 80, 19);

}

/**
* \brief Sets positions and sizes of the window
* \param x x-coordinate of the element
* \param y y-coordinate of the element
* \param w width of the element
* \param h height of the element
*/
void TexturePlotter::SetBounds(int x, int y, int w, int h) {

	GLWindow::SetBounds(x, y, w, h);
	PlaceComponents();

}


/**
* \brief Sets selected facet ID in the title of the window
*/
void TexturePlotter::GetSelected() {

	if (!worker) return;

	Geometry *geom = worker->GetGeometry();
	selFacetId = -1;
	selFacet = NULL;
	int i = 0;
	size_t nb = geom->GetNbFacet();
	while (selFacetId==-1 && i < nb) {
		if (geom->GetFacet(i)->selected) {
			selFacetId = i;
			selFacet = geom->GetFacet(i);
		}
		else {
			i++;
		}
	}

	char tmp[64];
	sprintf(tmp, "Texture plotter [Facet #%d]", i + 1);
	SetTitle(tmp);

}

/**
* \brief Updates table values if necessary
* \param appTime curent time of the application
* \param force if update should be forced
*/
void TexturePlotter::Update(float appTime, bool force) {

	if (!IsVisible()) return;

	if (force) {
		UpdateTable();
		lastUpdate = appTime;
		return;
	}

}

/**
* \brief Update table values for selected facets (only shows facet with lowest ID) corresponding to the texture values
* \param appTime curent time of the application
* \param force if update should be forced
*/
void TexturePlotter::UpdateTable() {
	size_t nbMoments = mApp->worker.moments.size();
	size_t facetHitsSize = (1 + nbMoments) * sizeof(FacetHitBuffer);
	maxValue = 0.0f;
	//double scale;
	GetSelected();
	if (!selFacet || selFacet->cellPropertiesIds.empty()) {
		mapList->Clear();
		return;
	}

	//SHELEM *mesh = selFacet->mesh;
	if (!selFacet->cellPropertiesIds.empty()) {

		char tmp[256];
		size_t w = selFacet->sh.texWidth;
		size_t h = selFacet->sh.texHeight;
		mapList->SetSize(w, h);
		mapList->SetAllColumnAlign(ALIGN_CENTER);

		int mode = viewCombo->GetSelectedIndex();

		switch (mode) {

		case 0: {// Cell area
			for (size_t i = 0; i < w; i++) {
				for (size_t j = 0; j < h; j++) {
				    FacetMomentSnapshot* dummPointer = nullptr;
					double val = worker->GetGeometry()->GetPhysicalValue(selFacet, PhysicalMode::CellArea, 1.0, 1.0, 1.0, (int)(i + j*w), *dummPointer).value;
					sprintf(tmp, "%g", val);
					if (val > maxValue) {
						maxValue = val;
						maxX = i; maxY = j;
					}
					mapList->SetValueAt(i, j, tmp);
				}
			}
			break; }

		case 1: {// MC Hits

			bool buffer_old = worker->ReloadIfNeeded(); //Locks and returns handle
			if (!buffer_old) return;

			size_t profSize = (selFacet->sh.isProfile) ? (PROFILE_SIZE * sizeof(ProfileSlice)*(1 + nbMoments)) : 0;
			const auto& facetSnapshot = worker->globalState.facetStates[selFacetId].momentResults[mApp->worker.displayedMoment];
			//TextureCell *texture = (TextureCell *)((BYTE *)buffer + (selFacet->sh.hitOffset + facetHitsSize + profSize + mApp->worker.displayedMoment*w*h * sizeof(TextureCell)));
			for (size_t i = 0; i < w; i++) {
				for (size_t j = 0; j < h; j++) {
					//int tSize = selFacet->wp.texWidth*selFacet->wp.texHeight;

					PhysicalValue val = worker->GetGeometry()->GetPhysicalValue(selFacet, PhysicalMode::MCHits, 1.0, 1.0, 1.0, (int)(i + j*w), facetSnapshot);
					double realVal = val.value;
					if (realVal > maxValue) {
						maxValue = realVal;
						maxX = i; maxY = j;
					}
					sprintf(tmp, "%g", realVal);
					mapList->SetValueAt(i, j, tmp);
				}
			}
			break; }

		case 2: {// Impingement rate

			// Lock during update
			bool buffer_old = worker->ReloadIfNeeded();
			if (!buffer_old) return;
				
				size_t profSize = (selFacet->sh.isProfile) ? (PROFILE_SIZE * sizeof(ProfileSlice)*(1 + nbMoments)) : 0;
				const auto& facetSnapshot = worker->globalState.facetStates[selFacetId].momentResults[mApp->worker.displayedMoment];
				double moleculesPerTP = mApp->worker.GetMoleculesPerTP(worker->displayedMoment);
				for (size_t i = 0; i < w; i++) {
					for (size_t j = 0; j < h; j++) {
						PhysicalValue val = worker->GetGeometry()->GetPhysicalValue(selFacet, PhysicalMode::ImpingementRate, moleculesPerTP,1.0, worker->model->wp.gasMass, (int)(i + j * w), facetSnapshot);
						double realVal = val.value;
						if (realVal > maxValue) {
							maxValue = realVal;
							maxX = i; maxY = j;
						}
						sprintf(tmp, "%g", realVal);
						mapList->SetValueAt(i, j, tmp);
					}
				}				
			
			break; }

		case 3: {// Particle density [1/m3]

					 // Lock during update
			bool buffer_old = worker->ReloadIfNeeded();
			if (!buffer_old) return;
			size_t profSize = (selFacet->sh.isProfile) ? (PROFILE_SIZE * sizeof(ProfileSlice)*(1 + nbMoments)) : 0;
			const auto& facetSnapshot = worker->globalState.facetStates[selFacetId].momentResults[mApp->worker.displayedMoment];

			double moleculesPerTP = mApp->worker.GetMoleculesPerTP(worker->displayedMoment);
			double densityCorrection = selFacet->DensityCorrection();

			for (size_t i = 0; i < w; i++) {
				for (size_t j = 0; j < h; j++) {
					PhysicalValue val = worker->GetGeometry()->GetPhysicalValue(selFacet, PhysicalMode::ParticleDensity, moleculesPerTP, densityCorrection, worker->model->wp.gasMass, (int)(i + j*w), facetSnapshot);
					double rho = val.value;

					if (rho > maxValue) {
						maxValue = rho;
						maxX = i; maxY = j;
					}

					sprintf(tmp, "%g", rho);
					mapList->SetValueAt(i, j, tmp);
				}
			}

			break; }

		case 4: {// Gas density [kg/m3]

			// Lock during update
			bool buffer_old = worker->ReloadIfNeeded();
			if (!buffer_old) return;
			size_t profSize = (selFacet->sh.isProfile) ? (PROFILE_SIZE * sizeof(ProfileSlice)*(1 + nbMoments)) : 0;
			const auto& facetSnapshot = worker->globalState.facetStates[selFacetId].momentResults[mApp->worker.displayedMoment];

			double moleculesPerTP = mApp->worker.GetMoleculesPerTP(worker->displayedMoment);
			double densityCorrection = selFacet->DensityCorrection();

			for (size_t i = 0; i < w; i++) {
				for (size_t j = 0; j < h; j++) {
					PhysicalValue val = worker->GetGeometry()->GetPhysicalValue(selFacet, PhysicalMode::GasDensity, moleculesPerTP, densityCorrection, worker->model->wp.gasMass, (int)(i + j*w), facetSnapshot);
					double rho_mass = val.value;
					if (rho_mass > maxValue) {
						maxValue = rho_mass;
						maxX = i; maxY = j;
					}

					sprintf(tmp, "%g", rho_mass);
					mapList->SetValueAt(i, j, tmp);
				}
			}

			break; }

		case 5: {// Pressure

					 // Lock during update
			bool buffer_old = worker->ReloadIfNeeded();
			if (!buffer_old) return;
				
			size_t profSize = (selFacet->sh.isProfile) ? (PROFILE_SIZE * sizeof(ProfileSlice)*(1 + nbMoments)) : 0;
			const auto& facetSnapshot = worker->globalState.facetStates[selFacetId].momentResults[mApp->worker.displayedMoment];

			double moleculesPerTP = mApp->worker.GetMoleculesPerTP(worker->displayedMoment);

			for (size_t i = 0; i < w; i++) {
				for (size_t j = 0; j < h; j++) {

					PhysicalValue val = worker->GetGeometry()->GetPhysicalValue(selFacet, PhysicalMode::Pressure, moleculesPerTP, 1.0, worker->model->wp.gasMass, (int)(i + j*w), facetSnapshot);
					double p = val.value;

					if (p > maxValue) {
						maxValue = p;
						maxX = i; maxY = j;
					}

					sprintf(tmp, "%g", p);

					mapList->SetValueAt(i, j, tmp);
				}
			}
			break; }

		case 6: {// Average gas velocity [m/s]

					// Lock during update
			bool buffer_old = worker->ReloadIfNeeded();
			if (!buffer_old) return;
				
				size_t profSize = (selFacet->sh.isProfile) ? (PROFILE_SIZE * sizeof(ProfileSlice)*(1 + nbMoments)) : 0;
				const auto& facetSnapshot = worker->globalState.facetStates[selFacetId].momentResults[mApp->worker.displayedMoment];

                for (size_t i = 0; i < w; i++) {
					for (size_t j = 0; j < h; j++) {
						PhysicalValue val = worker->GetGeometry()->GetPhysicalValue(selFacet, PhysicalMode::AvgGasVelocity, 1.0, 1.0, 1.0, (int)(i + j*w), facetSnapshot);
						double realVal = val.value;

						if (realVal > maxValue) {
							maxValue = realVal;
							maxX = i; maxY = j;
						}
						sprintf(tmp, "%g", realVal);
						mapList->SetValueAt(i, j, tmp);
					}
				}
			break; }

		case 7: {// Gas velocity vector
			bool buffer_old = worker->ReloadIfNeeded();
			if (!buffer_old) return;
			size_t profSize = (selFacet->sh.isProfile) ? (PROFILE_SIZE * sizeof(ProfileSlice)*(1 + nbMoments)) : 0;
			size_t nbElem = selFacet->sh.texWidth*selFacet->sh.texHeight;
			size_t tSize = nbElem * sizeof(TextureCell);
			size_t dSize = nbElem * sizeof(DirectionCell);
			const auto& facetSnapshot = worker->globalState.facetStates[selFacetId].momentResults[mApp->worker.displayedMoment];
			for (size_t i = 0; i < w; i++) {
				for (size_t j = 0; j < h; j++) {
					if (selFacet->sh.countDirection) {
						PhysicalValue val = worker->GetGeometry()->GetPhysicalValue(selFacet, PhysicalMode::GasVelocityVector, 1.0, 1.0, 1.0, (int)(i + j*w), facetSnapshot);
						Vector3d v_vect = val.vect;

						sprintf(tmp, "%g,%g,%g",
							v_vect.x, v_vect.y, v_vect.z);
						double length = v_vect.Norme();
						if (length > maxValue) {
							maxValue = length;
							maxX = i; maxY = j;
						}
					}
					else {
						sprintf(tmp, "Direction not recorded");
					}
					mapList->SetValueAt(i, j, tmp);
				}
			}
			break; }

		case 8: {// Nb of velocity vectors
			bool buffer_old = worker->ReloadIfNeeded();
			if (!buffer_old) return;
			size_t profSize = (selFacet->sh.isProfile) ? (PROFILE_SIZE * sizeof(ProfileSlice)*(1 + nbMoments)) : 0;
			size_t nbElem = selFacet->sh.texWidth*selFacet->sh.texHeight;
			size_t tSize = nbElem * sizeof(TextureCell);
			size_t dSize = nbElem * sizeof(DirectionCell);
			const auto& facetSnapshot = worker->globalState.facetStates[selFacetId].momentResults[mApp->worker.displayedMoment];

				for (size_t i = 0; i < w; i++) {
					for (size_t j = 0; j < h; j++) {
						if (selFacet->sh.countDirection) {
							PhysicalValue val = worker->GetGeometry()->GetPhysicalValue(selFacet, PhysicalMode::NbVelocityVectors, 1.0, 1.0, 1.0, (int)(i + j*w), facetSnapshot);
							size_t count = val.count;

							sprintf(tmp, "%zd", count);
							double countEq = (double)count;
							if (countEq > maxValue) {
								maxValue = countEq;
								maxX = i; maxY = j;
							}
						}
						else {
							sprintf(tmp, "Direction not recorded");
						}
						mapList->SetValueAt(i, j, tmp);
					}
				}	
			break; }
		}

	}
	if (autoSizeOnUpdate->GetState()) mapList->AutoSizeColumn();
}

/**
* \brief Displays the window
* \param w Worker handle
*/
void TexturePlotter::Display(Worker *w) {

	worker = w;
	UpdateTable();
	SetVisible(true);

}

/**
* \brief Closes the window
*/
void TexturePlotter::Close() {
	worker = NULL;
	if (selFacet) selFacet->UnselectElem();
	mapList->Clear();
}

/**
* \brief Saves table values to a file
*/
void TexturePlotter::SaveFile() {

	if (!selFacet) return;

	std::string fn = NFD_SaveFile_Cpp(fileFilters, "");
	if (!fn.empty()) {

		size_t u, v, wu, wv;
		if (!mapList->GetSelectionBox(&u, &v, &wu, &wv)) {
			u = 0;
			v = 0;
			wu = mapList->GetNbRow();
			wv = mapList->GetNbColumn();
		}

		// Save tab separated text
		FILE *f = fopen(fn.c_str(), "w");

		if (f == NULL) {
			char errMsg[512];
			sprintf(errMsg, "Cannot open file\nFile:%s", fn.c_str());
			GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}

		for (size_t i = u; i < u + wu; i++) {
			for (size_t j = v; j < v + wv; j++) {
				char *str = mapList->GetValueAt(j, i);
				if (str) fprintf(f, "%s", str);
				if (j < v + wv - 1)
					fprintf(f, "\t");
			}
			fprintf(f, "\r\n");
		}
		fclose(f);

	}

}

/**
* \brief Function for processing various inputs (button, check boxes etc.)
* \param src Exact source of the call
* \param message Type of the source (button)
*/
void TexturePlotter::ProcessMessage(GLComponent *src, int message) {

	switch (message) {

	case MSG_CLOSE:
		Close();
		break;

	case MSG_BUTTON:
		if (src == cancelButton) {
			Close();
			GLWindow::ProcessMessage(NULL, MSG_CLOSE);
		}
		else if (src == sizeButton) {
			mapList->AutoSizeColumn();
		}
		else if (src == saveButton) {
			SaveFile();
		}
		else if (src == maxButton) {
			size_t u, v, wu, wv;
			mapList->SetSelectedCell(maxX, maxY);
			if (mapList->GetSelectionBox(&v, &u, &wv, &wu))
				selFacet->SelectElem(u, v, wu, wv);
		}
		break;

	case MSG_LIST:
		if (src == mapList) {
			size_t u, v, wu, wv;
			if (mapList->GetSelectionBox(&v, &u, &wv, &wu))
				selFacet->SelectElem(u, v, wu, wv);
		}
		break;

	case MSG_COMBO:
		if (src == viewCombo) {
			UpdateTable();
			maxButton->SetEnabled(true);
			//maxButton->SetEnabled(viewCombo->GetSelectedIndex()!=2);
		}
		break;

	}

	GLWindow::ProcessMessage(src, message);
}
