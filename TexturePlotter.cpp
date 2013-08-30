/*
File:        TexturePlotter.cpp
Description: Texture plotter dialog
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

#include "TexturePlotter.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLMessageBox.h"
#include "GLApp/GLFileBox.h"
#include "MolFlow.h"

extern MolFlow *theApp;

extern double totalOutgassing;
extern double gasMass;

static const char *fileFilters = "Text files\0*.txt";
static const int   nbFilter = sizeof(fileFilters) / (2*sizeof(char *));

// --------------------------------------------------------------------

TexturePlotter::TexturePlotter():GLWindow() {

	int wD = 500;
	int hD = 300;
	lastUpdate = 0.0f;
	strcpy(currentDir,".");

	SetTitle("Texture plotter");
	SetResizable(TRUE);
	SetIconfiable(TRUE);
	SetMinimumSize(wD,hD);

	mapList = new GLList(0);
	mapList->SetColumnLabelVisible(TRUE);
	mapList->SetRowLabelVisible(TRUE);
	mapList->SetAutoColumnLabel(TRUE);
	mapList->SetAutoRowLabel(TRUE);
	mapList->SetRowLabelMargin(20);
	mapList->SetGrid(TRUE);
	mapList->SetSelectionMode(BOX_CELL);
	mapList->SetCornerLabel("\202\\\201");
	Add(mapList);

	viewLabel = new GLLabel("View");
	Add(viewLabel);
	viewCombo = new GLCombo(0);
	viewCombo->SetSize(6);
	viewCombo->SetValueAt(0,"Cell area (cm\262)");
	viewCombo->SetValueAt(1,"Pressure");
	viewCombo->SetValueAt(2,"# of Hits");
	viewCombo->SetValueAt(3,"Avg. velocity (m/s)");
	viewCombo->SetValueAt(4,"Velocity vector (m/s)");
	viewCombo->SetValueAt(5,"# of velocity vectors");

	viewCombo->SetSelectedIndex(1);
	Add(viewCombo);

	saveButton = new GLButton(0,"Save");
	Add(saveButton);

	sizeButton = new GLButton(0,"Autosize");
	Add(sizeButton);

	maxButton = new GLButton(0,"Find Max.");
	Add(maxButton);

	cancelButton = new GLButton(0,"Dismiss");
	Add(cancelButton);

	autoSizeOnUpdate = new GLToggle(0,"Autosize on every update");
	autoSizeOnUpdate->SetCheck(TRUE);
	Add(autoSizeOnUpdate);

	// Center dialog
	int wS,hS;
	GLToolkit::GetScreenSize(&wS,&hS);
	int xD = (wS-wD)/2;
	int yD = (hS-hD)/2;
	SetBounds(xD,yD,wD,hD);

	RestoreDeviceObjects();

	worker = NULL;

}

// --------------------------------------------------------------------

void TexturePlotter::PlaceComponents() {

	mapList->SetBounds(5,5,width-15,height-80);
	saveButton->SetBounds(10,height-70,70,19);
	sizeButton->SetBounds(10,height-45,70,19);
	autoSizeOnUpdate->SetBounds(90,height-45,120,19);
	maxButton->SetBounds(90,height-70,70,19);
	viewLabel->SetBounds(330,height-70,30,19);
	viewCombo->SetBounds(360,height-70,120,19);
	cancelButton->SetBounds(width-90,height-45,80,19);

}

// -----------------------------------------------------------------

void TexturePlotter::SetBounds(int x,int y,int w,int h) {

	GLWindow::SetBounds(x,y,w,h);
	PlaceComponents();

}

// --------------------------------------------------------------------

void TexturePlotter::GetSelected() {

	if(!worker) return;

	Geometry *geom = worker->GetGeometry();
	selFacet = NULL;
	int i = 0;
	int nb = geom->GetNbFacet();
	while(!selFacet && i<nb) {
		if( geom->GetFacet(i)->selected ) selFacet = geom->GetFacet(i);
		if(!selFacet) i++;
	}

	char tmp[32];
	sprintf(tmp,"Texture plotter #%d",i+1);
	SetTitle(tmp);

}

// --------------------------------------------------------------------

void TexturePlotter::Update(float appTime,BOOL force) {

	if(!IsVisible()) return;

	if(force) {
		UpdateTable();
		lastUpdate = appTime;
		return;
	}

	if( (appTime-lastUpdate>1.0f) ) {
		if(worker->running) UpdateTable();
		lastUpdate = appTime;
	}

}

// --------------------------------------------------------------------

void TexturePlotter::UpdateTable() {

	MolFlow *mApp = (MolFlow *)theApp;
	int nbMoments=(int)mApp->worker.moments.size();
	maxValue=0.0f;
	double scale;
	GetSelected();
	if( !selFacet ) {
		mapList->Clear();
		return;
	}

	SHELEM *mesh = selFacet->mesh;
	if( mesh ) {

		char tmp[256];
		int w = selFacet->sh.texWidth;
		int h = selFacet->sh.texHeight;
		mapList->SetSize(w,h);
		mapList->SetColumnAlign(ALIGN_CENTER);


		int mode = viewCombo->GetSelectedIndex();

		switch(mode) {

		case 0: // Element area
			for(int i=0;i<w;i++) {
				for(int j=0;j<h;j++) {
					float val=selFacet->mesh[i+j*w].area;
					sprintf(tmp,"%g",val);
					if (val>maxValue) {
						maxValue=val;
						maxX=i;maxY=j;
					}
					mapList->SetValueAt(i,j,tmp);
				}
			}
			break;

		case 1:  {// Pressure

			// Lock during update
			BYTE *buffer = worker->GetHits();
			try {
				if(buffer) {
					SHGHITS *shGHit = (SHGHITS *)buffer;
					int profSize = (selFacet->sh.isProfile)?(PROFILE_SIZE*sizeof(llong)*(1+nbMoments)):0;
					AHIT *hits = (AHIT *)((BYTE *)buffer + (selFacet->sh.hitOffset + sizeof(SHHITS) + profSize + mApp->worker.displayedMoment*w*h*sizeof(AHIT)));
					float dCoef = (float)totalOutgassing ;
					if( shGHit->mode == MC_MODE ) dCoef *= ((mApp->worker.displayedMoment==0)?1.0f:(((float)worker->desorptionStopTime-(float)worker->desorptionStartTime)
						/(float)worker->timeWindowSize	))/ (float)shGHit->total.hit.nbDesorbed;
					for(int i=0;i<w;i++) {
						for(int j=0;j<h;j++) {
							float val=hits[i+j*w]*dCoef; //already divided by area
							if (val>maxValue) {
								maxValue=val;
								maxX=i;maxY=j;
							}
							sprintf(tmp,"%g",val);
							mapList->SetValueAt(i,j,tmp);
						}
					}
					worker->ReleaseHits();
				}
			} catch (...) { //incorrect hits reference
				worker->ReleaseHits();
			}
				 
			break;}

		case 2:  {// Hits

			// Lock during update
			BYTE *buffer = worker->GetHits();
			try {
			if(buffer) {
				SHGHITS *shGHit = (SHGHITS *)buffer;
				int profSize = (selFacet->sh.isProfile)?(PROFILE_SIZE*sizeof(llong)*(1+nbMoments)):0;
				AHIT *hits = (AHIT *)((BYTE *)buffer + (selFacet->sh.hitOffset + sizeof(SHHITS) + profSize + mApp->worker.displayedMoment*w*h*sizeof(AHIT) ));
				for(int i=0;i<w;i++) {
					for(int j=0;j<h;j++) {
						int tSize = selFacet->sh.texWidth*selFacet->sh.texHeight;

						scale = 0.27497/sqrt(selFacet->sh.temperature/gasMass); //0.27497 = 40 / 145.469

						if(selFacet->sh.opacity>0.0 && selFacet->sh.opacity!=1.0) scale = scale/selFacet->sh.opacity;
						if(selFacet->sh.is2sided) scale = scale / 2.0;
						float val=(float)((hits[i+j*w]/(float)scale*selFacet->mesh[i+j*w].area));
						if (val>maxValue) {
							maxValue=val;
							maxX=i;maxY=j;
						}
						sprintf(tmp,"%g",val);
						mapList->SetValueAt(i,j,tmp);
					}
				}
				worker->ReleaseHits();
			}
			} catch (...) {
				worker->ReleaseHits();
				
			}
				 
			break;}

		case 3: // Average velocity
			for(int i=0;i<w;i++) {
				for(int j=0;j<h;j++) {
					if( selFacet->dirCache ) {
						double v = selFacet->dirCache[i+j*w].sumSpeed/(double)selFacet->dirCache[i+j*w].count;

						if ((float)v>maxValue) {
							maxValue=(float)v;
							maxX=i;maxY=j;
						}

						sprintf(tmp,"%g",v);
					} else {
						sprintf(tmp,"None");
					}
					mapList->SetValueAt(i,j,tmp);
				}
			}
			break;

		case 4: // Velocity
			for(int i=0;i<w;i++) {
				for(int j=0;j<h;j++) {
					if( selFacet->dirCache ) {
						sprintf(tmp,"%g,%g,%g",
							selFacet->dirCache[i+j*w].dir.x,
							selFacet->dirCache[i+j*w].dir.y,
							selFacet->dirCache[i+j*w].dir.z);
					} else {
						sprintf(tmp,"None");
					}
					mapList->SetValueAt(i,j,tmp);
				}
			}
			break;

		case 5: // Count (velocity vectors)
			for(int i=0;i<w;i++) {
				for(int j=0;j<h;j++) {
					if( selFacet->dirCache ) {
						llong val=selFacet->dirCache[i+j*w].count;
						if ((float)val>maxValue) {
							maxValue=(float)val;
							maxX=i;maxY=j;
						}
						sprintf(tmp,"%I64d",val);
					} else {
						sprintf(tmp,"None");
					}
					mapList->SetValueAt(i,j,tmp);
				}
			}
			break;
		}

	}
	if (autoSizeOnUpdate->IsChecked()) mapList->AutoSizeColumn();
}

// --------------------------------------------------------------------

void TexturePlotter::Display(Worker *w) {

	worker = w;
	UpdateTable();
	SetVisible(TRUE);

}

// --------------------------------------------------------------------

void TexturePlotter::Close() {
	worker = NULL;
	if(selFacet) selFacet->UnselectElem();
	mapList->Clear();
}

// --------------------------------------------------------------------

void TexturePlotter::SaveFile() {

	if(!selFacet) return;

	FILENAME *fn = GLFileBox::SaveFile(currentDir,NULL,"Save File",fileFilters,nbFilter);

	if( fn ) {

		int u,v,wu,wv;
		if( !mapList->GetSelectionBox(&u,&v,&wu,&wv) ) {
			u=0;
			v=0;
			wu = mapList->GetNbRow();
			wv = mapList->GetNbColumn();
		}

		// Save tab separated text
		FILE *f = fopen(fn->fullName,"w");

		if( f==NULL ) {
			char errMsg[512];
			sprintf(errMsg,"Cannot open file\nFile:%s",fn->fullName);
			GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
			return;
		}

		for(int i=u;i<u+wu;i++) {
			for(int j=v;j<v+wv;j++) {
				char *str = mapList->GetValueAt(j,i);
				if( str ) fprintf(f,"%s",str);
				if( j<v+wv-1 ) 
					fprintf(f,"\t");
			}
			fprintf(f,"\r\n");
		}
		fclose(f);

	}

}

// --------------------------------------------------------------------

void TexturePlotter::ProcessMessage(GLComponent *src,int message) {

	switch(message) {

	case MSG_CLOSE:
		Close();
		break;

	case MSG_BUTTON:
		if(src==cancelButton) {
			Close();
			GLWindow::ProcessMessage(NULL,MSG_CLOSE);
		} else if (src==sizeButton) {
			mapList->AutoSizeColumn();
		} else if (src==saveButton) {
			SaveFile();
		} else if (src==maxButton) {
			int u,v,wu,wv;
			mapList->SetSelectedCell(maxX,maxY);
			if( mapList->GetSelectionBox(&v,&u,&wv,&wu) )
				selFacet->SelectElem(u,v,wu,wv);
		}
		break;

	case MSG_LIST:
		if(src==mapList) {
			int u,v,wu,wv;
			if( mapList->GetSelectionBox(&v,&u,&wv,&wu) )
				selFacet->SelectElem(u,v,wu,wv);
		}
		break;

	case MSG_COMBO:
		if(src==viewCombo) {
			UpdateTable();
			maxButton->SetEnabled(TRUE);
			//maxButton->SetEnabled(viewCombo->GetSelectedIndex()!=2);
		}
		break;

	}

	GLWindow::ProcessMessage(src,message);
}
