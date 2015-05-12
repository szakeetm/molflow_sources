/*
  File:        MolFlow.cpp
  Description: Main application class (GUI management)
  Program:     MolFlow+
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


#include <crtdbg.h> //To debug heap corruptions in memory

#include "GLApp/GLApp.h"
#include "GLApp/GLTextField.h"
#include "GLApp/GLToggle.h"
#include "GLApp/GLLabel.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLTitledPanel.h"
#include "GLApp/GLList.h"
#include "GLApp/GLCombo.h"
#include "GLApp/GLMenuBar.h"

#include "Worker.h"
#include "GeometryViewer.h"
#include "FormulaSettings.h"
#include "CollapseSettings.h"
#include "ImportDesorption.h"
#include "MoveVertex.h"
#include "TimeSettings.h"
#include "ScaleVertex.h"
#include "ScaleFacet.h"
#include "MoveFacet.h"
#include "ExtrudeFacet.h"
#include "MirrorFacet.h"
#include "RotateFacet.h"
#include "Movement.h"
#include "AlignFacet.h"
#include "AddVertex.h"
#include "FacetMesh.h"
#include "FacetCoordinates.h"
#include "VertexCoordinates.h"
#include "FacetDetails.h"
#include "Viewer3DSettings.h"
#include "TextureSettings.h"
#include "GlobalSettings.h"
#include "ProfilePlotter.h"
#include "PressureEvolution.h"
#include "TimewisePlotter.h"
#include "ViewEditor.h"
#include "TexturePlotter.h"
#include "OutgassingMap.h"
#include "SelectDialog.h"
#include "MomentsEditor.h"
#include "ParameterEditor.h"

#define MAX_FORMULA 10
#define MAX_VIEW    19
#define MAX_SELECTION 19
#define MAX_RECENT  10

typedef struct {
  GLLabel     *name;
  GLTextField *value;
  GLButton    *setBtn;
  GLParser    *parser;
} FORMULA;

class MolFlow : public GLApplication
{
public:
    MolFlow();

    // Worker handle
    Worker worker;

    // Current directory
    void UpdateCurrentDir(char *fileName);
    char currentDir[1024];
    void UpdateCurrentSelDir(char *fileName);
    char currentSelDir[1024];

    // Simulation state
    float    lastUpdate;   // Last 'hit update' time
    double   hps;          // Hit per second
    double   dps;          // Desorption (or new particle) per second
    double   lastHps;      // hps measurement
    double   lastDps;      // dps measurement
    llong    lastNbHit;    // measurement
    llong    lastNbDes;    // measurement
    llong    nbDesStart;   // measurement
    llong    nbHitStart;   // measurement
    int      nbProc;       // Temporary var (use Worker::GetProcNumber)

	float    lastAppTime;
	BOOL     antiAliasing;
	BOOL     whiteBg;
	float    lastMeasTime; // Last measurement time (for hps and dps)
	double   tolerance; //Select coplanar tolerance
	double   largeArea; //Selection filter
	double   planarityThreshold; //Planarity threshold
	int      checkForUpdates;
	int      autoUpdateFormulas;
	int      compressSavedFiles;
	int      autoSaveSimuOnly;
	int      numCPU;
	BOOL     changedSinceSave; //For saving and autosaving
	double   autoSaveFrequency; //autosave period, in minutes
	float    lastSaveTime;
	float    lastSaveTimeSimu;
	std::string autosaveFilename; //only delete files that this instance saved

	//Public textfields so we can disable them from "Advanced facet parameters":
	GLTextField   *facetFlow;
	GLTextField   *facetFlowArea;

	HANDLE compressProcessHandle;
	//float m_fTime;

    // Util functions
	//void SendHeartBeat(BOOL forced=FALSE);
    char *FormatInt(llong v,char *unit);
    char *FormatPS(double v,char *unit);
    char *FormatSize(DWORD size);
    char *FormatTime(float t);
    void PlaceComponents();
    void LoadSelection(char *fName=NULL);
    void SaveSelection();
    void LoadFile(char *fName=NULL);
	void LoadFileMemory();
	void InsertGeometry(BOOL newStr,char *fName=NULL);
	void SaveFile();
    void SaveFileAs();
    void ExportSelection();
	void ImportDesorption_DES();
	void ExportTextures(int mode);
	void ExportProfiles();
    void ClearFacetParams();
    void UpdateFacetParams(BOOL updateSelection=FALSE);
    void ApplyFacetParams();
    void UpdateModelParams();
	void UpdateViewerFlags();
    void StartStopSimulation();
    void ResetSimulation(BOOL askConfirm=TRUE);
    void SaveConfig();
    void LoadConfig();
    void UpdateStructMenu();
    void UpdateTitle();
    void UpdateFacetHits();
    void AnimateViewerChange(int next,BOOL init=FALSE);
    void UpdateViewerParams();
    void SelectViewer(int s);
    void Place3DViewer();
	void UpdateMeasurements();
	void QuickPipe();
	void UpdateFacetlistSelected();
	BOOL AskToSave();
	BOOL AskToReset(Worker *work=NULL);
	float GetAppTime();
	void ResetAutoSaveTimer();
	BOOL AutoSave(BOOL crashSave=FALSE);
	void CheckForRecovery();
	void UpdateViewers();
	void SetFacetSearchPrg(BOOL visible,char *text);
	void DisplayCollapseDialog();
	void RenumberSelections(int startFacetId);

    // Formula management
    int nbFormula;
	FORMULA formulas[MAX_FORMULA];
	void ProcessFormulaButtons(GLComponent *src);
    void UpdateFormula();
	void OffsetFormula(char* expression,int offset,int filter=0);
	void RenumberFormulas(int startId);
    void AddFormula(GLParser *f,BOOL doUpdate=TRUE);
	void AddFormula(const char *fName, const char *formula);
	void ClearFormula();
	void ClearParameters();

	//Flow/sticking coeff. conversion
	void calcFlow();
	void calcSticking();

    // Recent files
    char *recents[MAX_RECENT];
    int  nbRecent;
    void AddRecent(char *fileName);
    void RemoveRecent(char *fileName);

    // Components
    GLMenuBar     *menu;
    GeometryViewer *viewer[MAX_VIEWER];
    GLTextField   *geomNumber;
    GLToggle      *showNormal;
    GLToggle      *showRule;
    GLToggle      *showUV;
    GLToggle      *showLeak;
    GLToggle      *showHit;
    GLToggle      *showLine;
    GLToggle      *showVolume;
    GLToggle      *showTexture;
    GLToggle      *showFilter;
    GLToggle      *showIndex;
    GLToggle      *showVertex;
    GLButton      *showMoreBtn;
    GLLabel       *facetInfo;
    GLButton      *startSimu;
    GLButton      *resetSimu;
	GLButton      *texturePlotterShortcut;
	GLButton      *profilePlotterShortcut;
    //GLButton      *statusSimu;
    GLCombo       *modeCombo;
    GLTextField   *hitNumber;
    GLTextField   *desNumber;
    GLTextField   *leakNumber;
    GLTextField   *sTime;
    GLMenu        *facetMenu;
    GLTextField   *facetSticking;
    GLTextField   *facetOpacity;
    GLTextField   *facetTemperature;
	GLTextField   *facetArea;
	GLCombo       *facetSideType;
    GLCombo       *facetDesType;
	GLTextField   *facetDesTypeN;
    GLCombo       *facetRecType;
    GLButton      *facetApplyBtn;
    GLButton      *facetDetailsBtn;
    GLButton      *facetCoordBtn;
    GLButton      *facetMoreBtn;
    GLTitledPanel *facetPanel;
    GLList        *facetList;
    GLTitledPanel *togglePanel;
	GLLabel       *facetUseDesFileLabel;
	GLLabel       *modeLabel;
	GLLabel       *facetAreaLabel;
	GLLabel       *facetPumpingLabel;
	GLTextField   *facetPumping;
    GLButton      *compACBtn;
    GLButton      *singleACBtn;
    GLLabel       *hitLabel;
    GLLabel       *desLabel;
    GLLabel       *leakLabel;
    GLLabel       *sTimeLabel;
	GLTitledPanel *shortcutPanel;
	GLButton      *profilePlotterBtn;
	GLButton      *texturePlotterBtn;
	GLButton      *textureScalingBtn;
	GLButton      *globalSettingsBtn;
	GLTitledPanel *simuPanel;
    GLLabel       *facetSLabel;
	GLLabel       *facetSideLabel;
    GLLabel       *facetTLabel;
    GLLabel       *facetTempLabel;
    GLLabel       *facetDLabel;
    GLLabel       *facetReLabel;
    GLToggle       *facetFILabel;
	GLToggle      *facetFIAreaLabel;
    //GLLabel       *facetMLabel;
    GLMenu        *structMenu;
    GLMenu        *viewsMenu;
	GLMenu        *selectionsMenu;
	GLMenu        *memorizeSelectionsMenu;
	GLMenu        *memorizeViewsMenu;
	GLMenu        *clearSelectionsMenu;
	GLMenu        *clearViewsMenu;
	GLTitledPanel *inputPanel;
	GLTitledPanel *outputPanel;

    // Views
    void SelectView(int v);
	void AddView(char *selectionName,AVIEW v);
    void AddView();
    void ClearViewMenus();
	void ClearAllViews();
	void OverWriteView(int idOvr);
	void ClearView(int idClr);
	void RebuildViewMenus();
	
    // Selections
    void SelectSelection(int v);
    void AddSelection(char *selectionName,ASELECTION s);
    void AddSelection();
    void ClearSelectionMenus();
	void ClearAllSelections();
	void OverWriteSelection(int idOvr);
	void ClearSelection(int idClr);
	void RebuildSelectionMenus();

	void AddStruct();
	void DeleteStruct();

    AVIEW   views[MAX_VIEW];
    int     nbView;
    int     idView;
    int     curViewer;
    int     modeSolo;

	ASELECTION selections[MAX_SELECTION];
	int nbSelection;
	int idSelection;

    //Dialog
    FormulaSettings  *formulaSettings;
    CollapseSettings *collapseSettings;
	ImportDesorption *importDesorption;
	MoveVertex		 *moveVertex;
	TimeSettings     *timeSettings;
	ScaleFacet       *scaleFacet;
	ScaleVertex      *scaleVertex;
	SelectDialog     *selectDialog;
	ExtrudeFacet	 *extrudeFacet;
	MoveFacet		 *moveFacet;
	MirrorFacet	     *mirrorFacet;
	RotateFacet      *rotateFacet;
	Movement         *movement;
	AlignFacet       *alignFacet;
	AddVertex		 *addVertex;
    FacetMesh        *facetMesh;
    FacetDetails     *facetDetails;
    Viewer3DSettings *viewer3DSettings;
    TextureSettings  *textureSettings;
	GlobalSettings	 *globalSettings;
    FacetCoordinates *facetCoordinates;
	VertexCoordinates *vertexCoordinates;
    ProfilePlotter   *profilePlotter;
	PressureEvolution *pressureEvolution;
	TimewisePlotter  *timewisePlotter;
    ViewEditor       *viewEditor;
    TexturePlotter   *texturePlotter;
	OutgassingMap    *outgassingMap;
	MomentsEditor    *momentsEditor;
	ParameterEditor  *parameterEditor;
	char *nbF;

    // Testing
    //int     nbSt;
    //void LogProfile();
    void BuildPipe(double ratio,int steps=0);
    //void BuildPipeStick(double s);
	
	void CrashHandler(Error *e);

protected:

    int  OneTimeSceneInit();
    int  RestoreDeviceObjects();
    int  InvalidateDeviceObjects();
    int  OnExit();
    int  FrameMove();
    void ProcessMessage(GLComponent *src,int message);
    int  Resize(DWORD width, DWORD height, BOOL forceWindowed);

};
