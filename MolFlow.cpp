/*
File:        MolFlow.cpp
Description: Main application class (GUI management)
Program:     MolFlow+
Author:      Roberto KERSEVAN / J-L Pons / Marton ADY
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

#include <math.h>
#include <malloc.h>
#include "MolFlow.h"
#include "File.h"
#include "GLApp/GLMessageBox.h"
#include "GLApp/GLInputBox.h"
#include "GLApp/GLFileBox.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLWindowManager.h"
#include "RecoveryDialog.h"
#include "Utils.h" //for Remainder
#include "direct.h"
#include <vector>
#include <string>
#include <io.h>

/*
//Leak detection
#ifdef _DEBUG
#define _CRTDBG_MAP_ALLOC
#define DEBUG_NEW new(_NORMAL_BLOCK, __FILE__, __LINE__)
#define new DEBUG_NEW
#endif
*/

static const char *fileLFilters = "All MolFlow supported files\0*.txt;*.xml;*.zip;*.geo;*.geo7z;*.syn;*.syn7z;*.str;*.stl;*.ase\0"
"All files\0*.*\0";
static const char *fileInsFilters = "All insertable geometries\0*.txt;*.xml;*.zip;*.geo;*.geo7z;*.syn;*.syn7z;*.stl\0"
"All files\0*.*\0";
const char *fileSFilters = "MolFlow saveable files\0*.xml;*.zip;*.geo;*.geo7z;*.txt\0All files\0*.*\0";
static const char *fileDesFilters = "Desorption files\0*.des\0All files\0*.*\0";

int   cWidth[] = { 30, 56, 50, 50 };
char *cName[] = { "#", "Hits", "Des", "Abs" };

#ifdef _DEBUG
std::string appName = "MolFlow+ development version 64-bit (Compiled " __DATE__ " " __TIME__ ") DEBUG MODE";
#else
std::string appName = "Molflow+ 2.6.33 64-bit (" __DATE__ ")";
#endif

std::vector<string> formulaPrefixes = { "A","D","H","P","DEN","Z","V","T","AR","a","d","h","ar","," };

float m_fTime;
MolFlow *mApp;

//Menu elements, Molflow specific:
#define MENU_FILE_IMPORTDES_DES 120
#define MENU_FILE_IMPORTDES_SYN 121


#define MENU_FILE_EXPORTTEXTURE_AREA 151
#define MENU_FILE_EXPORTTEXTURE_MCHITS 152
#define MENU_FILE_EXPORTTEXTURE_IMPINGEMENT 153
#define MENU_FILE_EXPORTTEXTURE_PART_DENSITY 154
#define MENU_FILE_EXPORTTEXTURE_GAS_DENSITY 155
#define MENU_FILE_EXPORTTEXTURE_PRESSURE 156
#define MENU_FILE_EXPORTTEXTURE_AVG_V 157
#define MENU_FILE_EXPORTTEXTURE_V_VECTOR 158
#define MENU_FILE_EXPORTTEXTURE_N_VECTORS 159

#define MENU_FILE_EXPORTTEXTURE_AREA_COORD 1510
#define MENU_FILE_EXPORTTEXTURE_MCHITS_COORD  1520
#define MENU_FILE_EXPORTTEXTURE_IMPINGEMENT_COORD  1530
#define MENU_FILE_EXPORTTEXTURE_PART_DENSITY_COORD  1540
#define MENU_FILE_EXPORTTEXTURE_GAS_DENSITY_COORD  1550
#define MENU_FILE_EXPORTTEXTURE_PRESSURE_COORD  1560
#define MENU_FILE_EXPORTTEXTURE_AVG_V_COORD  1570
#define MENU_FILE_EXPORTTEXTURE_V_VECTOR_COORD  1580
#define MENU_FILE_EXPORTTEXTURE_N_VECTORS_COORD  1590


#define MENU_EDIT_MOVINGPARTS 26

#define MENU_FACET_MESH        306
#define MENU_SELECT_HASDESFILE 313
#define MENU_FACET_OUTGASSINGMAP 335

#define MENU_TIME_SETTINGS          50
#define MENU_TIMEWISE_PLOTTER       51
#define MENU_TIME_PRESSUREEVOLUTION 52
#define MENU_TIME_MOMENTS_EDITOR    53
#define MENU_TIME_PARAMETER_EDITOR  54


//-----------------------------------------------------------------------------
// Name: WinMain()
// Desc: Entry point to the program. Initializes everything, and goes into a
//       message-processing loop. Idle time is used to render the scene.
//-----------------------------------------------------------------------------

INT WINAPI WinMain(HINSTANCE hInst, HINSTANCE, LPSTR, INT)
{
	MolFlow *mApp = new MolFlow();

	if (!mApp->Create(1024, 800, FALSE)) {
		char *logs = GLToolkit::GetLogs();
#ifdef WIN
		if (logs) MessageBox(NULL, logs, "Molflow [Fatal error]", MB_OK);
#else
		if( logs ) {
			printf("Molflow [Fatal error]\n");
			printf(logs);
		}
#endif
		SAFE_FREE(logs);
		delete mApp;
		return -1;
	}
	try {
		mApp->Run();
	}
	catch (Error &e) {
		((MolFlow*)mApp)->CrashHandler(&e);
	}
	delete mApp;
	return 0;
}



//-----------------------------------------------------------------------------
// Name: MolFlow()
// Desc: Application constructor. Sets default attributes for the app.
//-----------------------------------------------------------------------------

MolFlow::MolFlow()
{
	mApp = this; //to refer to the app as extern variable

	//Different Molflow implementation:
	facetMesh = NULL;
	facetDetails = NULL;
	smartSelection = NULL;
	viewer3DSettings = NULL;
	textureSettings = NULL;
	globalSettings = NULL;
	profilePlotter = NULL;
	texturePlotter = NULL;
				 
	//Molflow only:
	movement = NULL;
	timewisePlotter = NULL;
	pressureEvolution = NULL;
	outgassingMap = NULL;
	momentsEditor = NULL;
	parameterEditor = NULL;
	importDesorption = NULL;
	timeSettings = NULL;
}

//-----------------------------------------------------------------------------
// Name: OneTimeSceneInit()
// Desc: Called during initial app startup, this function performs all the
//       permanent initialization.
//-----------------------------------------------------------------------------
int MolFlow::OneTimeSceneInit()
{
	/*
	//Enable memory check at EVERY malloc/free operation:
	int tmpFlag = _CrtSetDbgFlag( _CRTDBG_REPORT_FLAG );
	tmpFlag |= _CRTDBG_LEAK_CHECK_DF;
	_CrtSetDbgFlag( tmpFlag );
	*/

	OneTimeSceneInit_shared();

	menu->GetSubMenu("File")->Add("Export selected textures");
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->Add("Facet by facet");
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("Facet by facet")->Add("Cell Area (cm\262)", MENU_FILE_EXPORTTEXTURE_AREA);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("Facet by facet")->Add("# of MC Hits", MENU_FILE_EXPORTTEXTURE_MCHITS);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("Facet by facet")->Add("Impingement rate (1/s/m\262)", MENU_FILE_EXPORTTEXTURE_IMPINGEMENT);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("Facet by facet")->Add("Particle density (1/m\263)", MENU_FILE_EXPORTTEXTURE_PART_DENSITY);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("Facet by facet")->Add("Gas density (kg/m\263)", MENU_FILE_EXPORTTEXTURE_GAS_DENSITY);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("Facet by facet")->Add("Pressure (mbar)", MENU_FILE_EXPORTTEXTURE_PRESSURE);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("Facet by facet")->Add("Avg. Velocity (m/s)", MENU_FILE_EXPORTTEXTURE_AVG_V);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("Facet by facet")->Add("Velocity vector (m/s)", MENU_FILE_EXPORTTEXTURE_V_VECTOR);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("Facet by facet")->Add("# of velocity vectors", MENU_FILE_EXPORTTEXTURE_N_VECTORS);

	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->Add("By X,Y,Z coordinates");
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("By X,Y,Z coordinates")->Add("Cell Area (cm\262)", MENU_FILE_EXPORTTEXTURE_AREA_COORD);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("By X,Y,Z coordinates")->Add("# of MC Hits", MENU_FILE_EXPORTTEXTURE_MCHITS_COORD);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("By X,Y,Z coordinates")->Add("Impingement rate (1/s/m\262)", MENU_FILE_EXPORTTEXTURE_IMPINGEMENT_COORD);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("By X,Y,Z coordinates")->Add("Particle density (1/m\263)", MENU_FILE_EXPORTTEXTURE_PART_DENSITY_COORD);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("By X,Y,Z coordinates")->Add("Gas density (kg/m\263)", MENU_FILE_EXPORTTEXTURE_GAS_DENSITY_COORD);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("By X,Y,Z coordinates")->Add("Pressure (mbar)", MENU_FILE_EXPORTTEXTURE_PRESSURE_COORD);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("By X,Y,Z coordinates")->Add("Avg. Velocity (m/s)", MENU_FILE_EXPORTTEXTURE_AVG_V_COORD);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("By X,Y,Z coordinates")->Add("Velocity vector (m/s)", MENU_FILE_EXPORTTEXTURE_V_VECTOR_COORD);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("By X,Y,Z coordinates")->Add("# of velocity vectors", MENU_FILE_EXPORTTEXTURE_N_VECTORS_COORD);
	
	menu->GetSubMenu("File")->Add("Import desorption file");
	menu->GetSubMenu("File")->GetSubMenu("Import desorption file")->Add("SYN file", MENU_FILE_IMPORTDES_SYN);
	menu->GetSubMenu("File")->GetSubMenu("Import desorption file")->Add("DES file (deprecated)", MENU_FILE_IMPORTDES_DES);

	menu->GetSubMenu("File")->Add(NULL); // Separator
	menu->GetSubMenu("File")->Add("E&xit", MENU_FILE_EXIT); //Moved here from OnetimeSceneinit_shared to assert it's the last menu item

	menu->GetSubMenu("Selection")->Add(NULL); // Separator
	menu->GetSubMenu("Selection")->Add("Select Desorption", MENU_FACET_SELECTDES);
	menu->GetSubMenu("Selection")->Add("Select Outgassing Map", MENU_SELECT_HASDESFILE);
	menu->GetSubMenu("Selection")->Add("Select Reflective", MENU_FACET_SELECTREFL);
	menu->GetSubMenu("Selection")->Add("Select volatile facets", MENU_FACET_SELECTVOL);

	menu->GetSubMenu("Tools")->Add(NULL);
	menu->GetSubMenu("Tools")->Add("Moving parts...", MENU_EDIT_MOVINGPARTS);
	menu->GetSubMenu("Facet")->Add("Convert to outgassing map...", MENU_FACET_OUTGASSINGMAP);

	menu->Add("Time");
	menu->GetSubMenu("Time")->Add("Time settings...", MENU_TIME_SETTINGS, SDLK_i, ALT_MODIFIER);
	menu->GetSubMenu("Time")->Add("Edit moments...", MENU_TIME_MOMENTS_EDITOR);
	menu->GetSubMenu("Time")->Add("Edit parameters...", MENU_TIME_PARAMETER_EDITOR);
	menu->GetSubMenu("Time")->Add(NULL);
	menu->GetSubMenu("Time")->Add("Timewise plotter", MENU_TIMEWISE_PLOTTER);
	menu->GetSubMenu("Time")->Add("Pressure evolution", MENU_TIME_PRESSUREEVOLUTION);

	

	showFilter = new GLToggle(0, "Filtering");
	//togglePanel->Add(showFilter);

	showMoreBtn = new GLButton(0, "<< View");
	togglePanel->Add(showMoreBtn);

	shortcutPanel = new GLTitledPanel("Shortcuts");
	shortcutPanel->SetClosable(TRUE);
	shortcutPanel->Close();
	Add(shortcutPanel);

	profilePlotterBtn = new GLButton(0, "Profile pl.");
	shortcutPanel->Add(profilePlotterBtn);

	texturePlotterBtn = new GLButton(0, "Texture pl.");
	shortcutPanel->Add(texturePlotterBtn);

	textureScalingBtn = new GLButton(0, "Tex.scaling");
	shortcutPanel->Add(textureScalingBtn);

	

	globalSettingsBtn = new GLButton(0, "<< Sim");
	simuPanel->Add(globalSettingsBtn);

	/*
	statusSimu = new GLButton(0,"...");
	simuPanel->Add(statusSimu);
	*/

	modeLabel = new GLLabel("Mode");
	//simuPanel->Add(modeLabel);

	modeCombo = new GLCombo(0);
	modeCombo->SetEditable(TRUE);
	modeCombo->SetSize(2);
	modeCombo->SetValueAt(0, "Monte Carlo");
	modeCombo->SetValueAt(1, "Angular Coef");
	modeCombo->SetSelectedIndex(0);
	//simuPanel->Add(modeCombo);

	compACBtn = new GLButton(0, "Calc AC");
	compACBtn->SetEnabled(FALSE);
	//simuPanel->Add(compACBtn);

	singleACBtn = new GLButton(0, "1");
	singleACBtn->SetEnabled(FALSE);
	//simuPanel->Add(singleACBtn);

	

	inputPanel = new GLTitledPanel("Particles in");
	facetPanel->Add(inputPanel);

	facetDLabel = new GLLabel("Desorption");
	facetPanel->Add(facetDLabel);
	facetDesType = new GLCombo(0);
	facetDesType->SetSize(4);
	facetDesType->SetValueAt(0, "None");
	facetDesType->SetValueAt(1, "Uniform");
	facetDesType->SetValueAt(2, "Cosine");
	facetDesType->SetValueAt(3, "Cosine^N");
	inputPanel->Add(facetDesType);

	facetDesTypeN = new GLTextField(0, NULL);
	facetDesTypeN->SetEditable(FALSE);
	facetPanel->Add(facetDesTypeN);

	facetFILabel = new GLToggle(0, "Outgassing (mbar*l/s):");
	facetFILabel->SetEnabled(FALSE);
	facetFILabel->SetState(TRUE);
	inputPanel->Add(facetFILabel);
	facetFlow = new GLTextField(0, NULL);
	inputPanel->Add(facetFlow);

	facetFIAreaLabel = new GLToggle(1, "Outg/area(mbar*l/s/cm\262):");
	facetFIAreaLabel->SetEnabled(FALSE);
	inputPanel->Add(facetFIAreaLabel);
	facetFlowArea = new GLTextField(0, NULL);
	inputPanel->Add(facetFlowArea);

	/*facetUseDesFileLabel = new GLLabel("Desorp. file");
	facetPanel->Add(facetUseDesFileLabel);
	facetUseDesFile = new GLCombo(0);
	facetUseDesFile->SetSize(1);
	facetUseDesFile->SetValueAt(0,"No desorption map");
	inputPanel->Add(facetUseDesFile);*/

	outputPanel = new GLTitledPanel("Particles out");
	facetPanel->Add(outputPanel);

	facetSLabel = new GLLabel("Sticking factor:");
	outputPanel->Add(facetSLabel);
	facetSticking = new GLTextField(0, NULL);
	outputPanel->Add(facetSticking);

	facetPumpingLabel = new GLLabel("Pumping Speed (l/s):");
	outputPanel->Add(facetPumpingLabel);
	facetPumping = new GLTextField(0, NULL);
	outputPanel->Add(facetPumping);

	facetReLabel = new GLLabel("Profile:");
	facetPanel->Add(facetReLabel);
	facetRecType = new GLCombo(0);
	facetRecType->SetSize(6);
	facetRecType->SetValueAt(0, "None");
	facetRecType->SetValueAt(1, "Pressure, density (\201)");
	facetRecType->SetValueAt(2, "Pressure, density (\202)");
	facetRecType->SetValueAt(3, "Incident angle");
	facetRecType->SetValueAt(4, "Speed distribution");
	facetRecType->SetValueAt(5, "Orthogonal velocity");
	facetPanel->Add(facetRecType);

	facetMoreBtn = new GLButton(0, "<< Adv");
	facetPanel->Add(facetMoreBtn);

	facetList = new GLList(0);
	facetList->SetWorker(&worker);
	facetList->SetGrid(TRUE);
	facetList->SetSelectionMode(MULTIPLE_ROW);
	facetList->SetSize(4, 1);
	facetList->SetColumnWidths((int*)cWidth);
	facetList->SetColumnLabels((char **)cName);
	facetList->SetColumnLabelVisible(TRUE);
	facetList->Sortable = TRUE;
	Add(facetList);

	
	facetMesh = new FacetMesh(&worker); //To use its UpdatefacetParams() routines
	
	ClearFacetParams();
	LoadConfig();
	UpdateViewerParams();
	PlaceComponents();
	CheckNeedsTexture();

	//LoadFile();
	try {

		worker.SetProcNumber(nbProc);

	}
	catch (Error &e) {
		char errMsg[512];
		sprintf(errMsg, "Failed to start working sub-process(es), simulation not available\n%s", e.GetMsg());
		GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
	}

	//PlaceComponents(); //Why was it here?

	//SelectViewer(0);

	//viewer[0]->Paint();

	if (checkForUpdates) {
		//Launch updater tool
		char command[1024];
		char CWD[MAX_PATH];
		_getcwd(CWD, MAX_PATH);
		//sprintf(tmp5,"%s\\molflow_updater_tmp.exe",CWD);
		if (FileUtils::Exist("molflow_updater_tmp.exe")) { //rename after new installation
			sprintf(command, "move \"%s\\molflow_updater_tmp.exe\" \"%s\\molflow_updater.exe\"", CWD, CWD);
			system(command);
		}

		if (FileUtils::Exist("molflow_updater.exe"))
			StartProc("synrad_updater.exe", STARTPROC_BACKGROUND);
		else GLMessageBox::Display("molflow_updater.exe not found. You will not receive updates to Molflow."
			"\n(You can disable checking for updates in Tools/Global Settings)", "Updater module missing.", GLDLG_OK, GLDLG_ICONINFO);
	}
	return GL_OK;
}


void MolFlow::PlaceComponents() {

	int sx = m_screenWidth - 205;
	int sy = 30;

	Place3DViewer();

	geomNumber->SetBounds(sx, 3, 202, 18);

	// Viewer settings ----------------------------------------
	togglePanel->SetBounds(sx, sy, 202, 112);

	togglePanel->SetCompBounds(showRule, 5, 20, 60, 18);
	togglePanel->SetCompBounds(showNormal, 70, 20, 60, 18);
	togglePanel->SetCompBounds(showUV, 135, 20, 60, 18);

	togglePanel->SetCompBounds(showLine, 5, 42, 60, 18);
	togglePanel->SetCompBounds(showLeak, 70, 42, 60, 18);
	togglePanel->SetCompBounds(showHit, 135, 42, 60, 18);

	togglePanel->SetCompBounds(showVolume, 5, 64, 60, 18);
	togglePanel->SetCompBounds(showTexture, 70, 64, 60, 18);
	togglePanel->SetCompBounds(showFilter, 135, 64, 60, 18);

	togglePanel->SetCompBounds(showMoreBtn, 5, 86, 55, 18);
	togglePanel->SetCompBounds(showVertex, 70, 86, 60, 18);
	togglePanel->SetCompBounds(showIndex, 137, 86, 60, 18);

	sy += (togglePanel->GetHeight() + 5);

	// Selected facet -----------------------------------------
	facetPanel->SetBounds(sx, sy, 202, 330);

	facetPanel->SetCompBounds(inputPanel, 5, 16, 192, 90);

	int cursorY = 15;
	inputPanel->SetCompBounds(facetDLabel, 5, cursorY, 60, 18);
	inputPanel->SetCompBounds(facetDesType, 65, cursorY, 80, 18);
	inputPanel->SetCompBounds(facetDesTypeN, 150, cursorY, 30, 18);

	inputPanel->SetCompBounds(facetFILabel, 5, cursorY += 25, 110, 18);
	inputPanel->SetCompBounds(facetFlow, 140, cursorY, 45, 18);

	inputPanel->SetCompBounds(facetFIAreaLabel, 5, cursorY += 25, 110, 18);
	inputPanel->SetCompBounds(facetFlowArea, 140, cursorY, 45, 18);

	//inputPanel->SetCompBounds(facetUseDesFileLabel,5,90,60,18);
	//inputPanel->SetCompBounds(facetUseDesFile,65,90,120,18);

	facetPanel->SetCompBounds(outputPanel, 5, cursorY += 45, 192, 65);

	outputPanel->SetCompBounds(facetSLabel, 7, cursorY = 15, 100, 18);
	outputPanel->SetCompBounds(facetSticking, 140, cursorY, 45, 18);

	outputPanel->SetCompBounds(facetPumpingLabel, 7, cursorY += 25, 100, 18);
	outputPanel->SetCompBounds(facetPumping, 140, cursorY, 45, 18);

	facetPanel->SetCompBounds(facetSideLabel, 7, cursorY = 180, 50, 18);
	facetPanel->SetCompBounds(facetSideType, 65, cursorY, 130, 18);

	facetPanel->SetCompBounds(facetTLabel, 7, cursorY += 25, 100, 18);
	facetPanel->SetCompBounds(facetOpacity, 110, cursorY, 82, 18);

	facetPanel->SetCompBounds(facetTempLabel, 7, cursorY += 25, 100, 18);
	facetPanel->SetCompBounds(facetTemperature, 110, cursorY, 82, 18);

	facetPanel->SetCompBounds(facetAreaLabel, 7, cursorY += 25, 100, 18);
	facetPanel->SetCompBounds(facetArea, 110, cursorY, 82, 18);

	facetPanel->SetCompBounds(facetReLabel, 7, cursorY += 25, 60, 18);
	facetPanel->SetCompBounds(facetRecType, 65, cursorY, 130, 18);

	facetPanel->SetCompBounds(facetMoreBtn, 5, cursorY += 25, 48, 18);
	facetPanel->SetCompBounds(facetDetailsBtn, 56, cursorY, 45, 18);
	facetPanel->SetCompBounds(facetCoordBtn, 104, cursorY, 45, 18);
	facetPanel->SetCompBounds(facetApplyBtn, 153, cursorY, 44, 18);

	sy += facetPanel->GetHeight() + 5;

	shortcutPanel->SetBounds(sx, sy, 202, 40);
	shortcutPanel->SetCompBounds(profilePlotterBtn, 5, 15, 60, 18);
	shortcutPanel->SetCompBounds(texturePlotterBtn, 70, 15, 60, 18);
	shortcutPanel->SetCompBounds(textureScalingBtn, 135, 15, 60, 18);

	sy += shortcutPanel->GetHeight() + 5;

	// Simulation ---------------------------------------------
	simuPanel->SetBounds(sx, sy, 202, 169);

	simuPanel->SetCompBounds(globalSettingsBtn, 5, 20, 48, 19);
	simuPanel->SetCompBounds(startSimu, 58, 20, 66, 19);
	simuPanel->SetCompBounds(resetSimu, 128, 20, 66, 19);
	//simuPanel->SetCompBounds(statusSimu,175,20,20,19);
	simuPanel->SetCompBounds(modeLabel, 5, 45, 30, 18);
	simuPanel->SetCompBounds(modeCombo, 40, 45, 85, 18);
	simuPanel->SetCompBounds(compACBtn, 130, 45, 65, 19);
	simuPanel->SetCompBounds(autoFrameMoveToggle, 5, 45, 65, 19);
	simuPanel->SetCompBounds(forceFrameMoveButton, 128, 45, 66, 19);

	//simuPanel->SetCompBounds(compACBtn,123,45,52,19);
	//simuPanel->SetCompBounds(singleACBtn,178,45,20,19);
	simuPanel->SetCompBounds(hitLabel, 5, 70, 30, 18);
	simuPanel->SetCompBounds(hitNumber, 40, 70, 155, 18);
	simuPanel->SetCompBounds(desLabel, 5, 95, 30, 18);
	simuPanel->SetCompBounds(desNumber, 40, 95, 155, 18);
	simuPanel->SetCompBounds(leakLabel, 5, 120, 30, 18);
	simuPanel->SetCompBounds(leakNumber, 40, 120, 155, 18);
	simuPanel->SetCompBounds(sTimeLabel, 5, 145, 30, 18);


	simuPanel->SetCompBounds(sTime, 40, 145, 155, 18);

	sy += (simuPanel->GetHeight() + 5);

	// ---------------------------------------------------------
	int lg = m_screenHeight - (nbFormula * 25 + 23);

	facetList->SetBounds(sx, sy, 202, lg - sy);

	// ---------------------------------------------------------

	for (int i = 0; i < nbFormula; i++) {
		formulas[i].name->SetBounds(sx, lg + 5, 95, 18);
		formulas[i].value->SetBounds(sx + 90, lg + 5, 87, 18);
		formulas[i].setBtn->SetBounds(sx + 182, lg + 5, 20, 18);
		lg += 25;
	}

}

//-----------------------------------------------------------------------------
// Name: ClearFacetParams()
// Desc: Reset selected facet parameters.
//-----------------------------------------------------------------------------

void MolFlow::ClearFacetParams() {
	facetPanel->SetTitle("Selected Facet (none)");
	facetSticking->Clear();
	facetSticking->SetEditable(FALSE);
	facetFILabel->SetEnabled(FALSE);
	facetFIAreaLabel->SetEnabled(FALSE);
	facetFlow->Clear();
	facetFlow->SetEditable(FALSE);
	facetFlowArea->Clear();
	facetFlowArea->SetEditable(FALSE);
	facetArea->SetEditable(FALSE);
	facetArea->Clear();
	facetPumping->SetEditable(FALSE);
	facetPumping->Clear();
	facetOpacity->Clear();
	facetOpacity->SetEditable(FALSE);
	facetTemperature->Clear();
	facetTemperature->SetEditable(FALSE);
	facetSideType->SetSelectedValue("");
	facetSideType->SetEditable(FALSE);
	facetDesType->SetSelectedValue("");
	facetDesType->SetEditable(FALSE);
	facetDesTypeN->SetText("");
	facetDesTypeN->SetEditable(FALSE);
	facetRecType->SetSelectedValue("");
	facetRecType->SetEditable(FALSE);
}

//-----------------------------------------------------------------------------
// Name: ApplyFacetParams()
// Desc: Apply facet parameters.
//-----------------------------------------------------------------------------

void MolFlow::ApplyFacetParams() {


	Geometry *geom = worker.GetGeometry();
	int nbFacet = geom->GetNbFacet();
	if (!AskToReset()) return;
	changedSinceSave = TRUE;


	// Sticking
	double sticking;
	BOOL stickingNotNumber;
	BOOL doSticking = FALSE;
	if (facetSticking->GetNumber(&sticking)) {
		if (sticking<0.0 || sticking>1.0) {
			GLMessageBox::Display("Sticking must be in the range [0,1]", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}
		doSticking = TRUE;
		stickingNotNumber = FALSE;
	}
	else {
		if (strcmp(facetSticking->GetText(), "...") == 0) doSticking = FALSE;
		else {/*
			GLMessageBox::Display("Invalid sticking number","Error",GLDLG_OK,GLDLG_ICONERROR);
			UpdateFacetParams();
			return;*/
			doSticking = TRUE;
			stickingNotNumber = TRUE;
		}
	}

	// opacity
	double opacity;
	BOOL doOpacity = FALSE;
	BOOL opacityNotNumber;
	if (facetOpacity->GetNumber(&opacity)) {
		if (opacity<0.0 || opacity>1.0) {
			GLMessageBox::Display("Opacity must be in the range [0,1]", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}
		doOpacity = TRUE;
		opacityNotNumber = FALSE;
	}
	else {
		if (strcmp(facetOpacity->GetText(), "...") == 0) doOpacity = FALSE;
		else {/*
			GLMessageBox::Display("Invalid opacity number","Error",GLDLG_OK,GLDLG_ICONERROR);
			UpdateFacetParams();
			return;*/
			doOpacity = TRUE;
			opacityNotNumber = TRUE;
		}
	}

	// temperature
	double temperature;
	BOOL doTemperature = FALSE;
	if (facetTemperature->GetNumber(&temperature)) {
		if (temperature < 0.0) {
			GLMessageBox::Display("Temperature must be positive", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}
		doTemperature = TRUE;
	}
	else {
		if (strcmp(facetTemperature->GetText(), "...") == 0) doTemperature = FALSE;
		else {
			GLMessageBox::Display("Invalid temperature number", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}
	}

	// Outgassing
	double flow = 0.0;
	BOOL doFlow = FALSE;
	BOOL outgassingNotNumber;
	//Calculate flow
	if (facetFILabel->GetState() && strcmp(facetFlow->GetText(), "...") != 0 && facetDesType->GetSelectedIndex() != 0
		&& strcmp(facetDesType->GetSelectedValue(), "...") != 0 && facetFlow->IsEditable()) {  //We want outgassing
		if (facetFlow->GetNumber(&flow)) { //If we can parse the number
			if (!(flow > 0.0)) {
				GLMessageBox::Display("Outgassing must be positive", "Error", GLDLG_OK, GLDLG_ICONERROR);
				return;
			}
			doFlow = TRUE;
			outgassingNotNumber = FALSE;
		}
		else { //could not parse as number
			doFlow = TRUE;
			outgassingNotNumber = TRUE;
		}
	}

	// Outgassing per area
	double flowA = 0;
	BOOL doFlowA = FALSE;
	//Calculate flow

	if (facetFIAreaLabel->GetState() && strcmp(facetFlowArea->GetText(), "...") != 0
		&& facetDesType->GetSelectedIndex() != 0 && strcmp(facetDesType->GetSelectedValue(), "...") != 0 && facetFlowArea->IsEditable()) { //We want outgassing per area
		if (facetFlowArea->GetNumber(&flowA)) { //Can be parsed as number
			if (!(flowA > 0.0)) {
				GLMessageBox::Display("Outgassing per area must be positive", "Error", GLDLG_OK, GLDLG_ICONERROR);
				return;
			}
			doFlowA = TRUE;
		}
		else {
			GLMessageBox::Display("Invalid outgassing per area number", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}
	}

	// Desorption type
	int desorbType = facetDesType->GetSelectedIndex();

	double desorbTypeN;
	BOOL doDesorbTypeN = FALSE;
	if (desorbType == 3) {
		if (facetDesTypeN->GetNumber(&desorbTypeN)) {
			if (!(desorbTypeN > 0.0)) {
				GLMessageBox::Display("Desorption type exponent must be greater than 0.0", "Error", GLDLG_OK, GLDLG_ICONERROR);
				UpdateFacetParams();
				return;
			}
			doDesorbTypeN = TRUE;
		}
		else {
			if (strcmp(facetDesTypeN->GetText(), "...") == 0) doDesorbTypeN = FALSE;
			else {
				GLMessageBox::Display("Invalid desorption type exponent", "Error", GLDLG_OK, GLDLG_ICONERROR);
				return;
			}
		}
	}

	// Record (profile) type
	int rType = facetRecType->GetSelectedIndex(); // -1 if "..."

	// 2sided
	int is2Sided = facetSideType->GetSelectedIndex();


	// Update facets (local)
	for (int i = 0; i < nbFacet; i++) {
		Facet *f = geom->GetFacet(i);
		if (f->selected) {
			if (doSticking) {
				if (!stickingNotNumber) {
					f->sh.sticking = sticking;
					f->userSticking = "";
				}
				else {
					f->userSticking = facetSticking->GetText();
				}
			}

			if (doOpacity) {
				if (!opacityNotNumber) {
					f->sh.opacity = opacity;
					f->userOpacity = "";
				}
				else {
					f->userOpacity = facetOpacity->GetText();
				}
			}
			if (doTemperature) f->sh.temperature = temperature;
			if (doFlow) {
				if (!outgassingNotNumber) {
					f->sh.flow = flow*0.100; //0.1: mbar*l/s -> Pa*m3/s
					f->userOutgassing = "";
				}
				else {
					f->userOutgassing = facetFlow->GetText();
				}
			}
			if (doFlowA/* && !useMapA*/) f->sh.flow = flowA*f->sh.area*(f->sh.is2sided ? 2.0 : 1.0)*0.100;
			if (desorbType >= 0) {
				if (desorbType == 0) f->sh.flow = 0.0;
				if (desorbType != 3) f->sh.desorbTypeN = 0.0;
				f->sh.desorbType = desorbType;
				if (doDesorbTypeN) f->sh.desorbTypeN = desorbTypeN;
			}


			if (rType >= 0) {
				f->sh.profileType = rType;
				//f->sh.isProfile = (rType!=REC_NONE); //included below by f->UpdateFlags();
			}
			if (is2Sided >= 0) f->sh.is2sided = is2Sided;

			f->sh.maxSpeed = 4.0 * sqrt(2.0*8.31*f->sh.temperature / 0.001 / worker.gasMass);
			f->UpdateFlags();
		}
	}
	if (facetMesh && facetMesh->IsVisible()) {
		if (!facetMesh->Apply()) {
			return;
		}
	}

	// Mark "needsReload"
	try { worker.Reload(); }
	catch (Error &e) {
		GLMessageBox::Display((char *)e.GetMsg(), "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}
	geom->CalcTotalOutGassing();
	UpdateFacetParams();
	if (profilePlotter) profilePlotter->Refresh();
	if (pressureEvolution) pressureEvolution->Refresh();
	if (timewisePlotter) timewisePlotter->Refresh();
	//if (facetMesh) facetMesh->Refresh();
}



//-----------------------------------------------------------------------------
// Name: UpdateFacetParams()
// Desc: Update selected facet parameters.
//-----------------------------------------------------------------------------

void MolFlow::UpdateFacetParams(BOOL updateSelection) { //Calls facetMesh->Refresh()

	char tmp[256];

	// Update params
	Geometry *geom = worker.GetGeometry();
	int nbSel = geom->GetNbSelected();
	if (nbSel > 0) {

		Facet *f0;
		Facet *f;

		// Get list of selected facet
		int *selection = (int *)malloc(nbSel*sizeof(int));
		int count = 0;
		for (int i = 0; i < geom->GetNbFacet(); i++)
			if (geom->GetFacet(i)->selected)
				selection[count++] = i;

		f0 = geom->GetFacet(selection[0]);

		double f0Area = f0->sh.area*(f0->sh.is2sided ? 2.0 : 1.0);
		double area = f0Area; //sum facet area

		BOOL stickingE = TRUE;
		BOOL opacityE = TRUE;
		BOOL temperatureE = TRUE;
		BOOL flowE = TRUE;
		BOOL flowAreaE = TRUE;
		BOOL desorbTypeE = TRUE;
		BOOL desorbTypeNE = TRUE;
		BOOL recordE = TRUE;
		BOOL is2sidedE = TRUE;

		for (int i = 1; i < count; i++) {
			f = geom->GetFacet(selection[i]);
			double fArea = f->sh.area*(f->sh.is2sided ? 2.0 : 1.0);
			stickingE = stickingE && (f0->userSticking.compare(f->userSticking) == 0) && IsEqual(f0->sh.sticking, f->sh.sticking);
			opacityE = opacityE && (f0->userOpacity.compare(f->userOpacity) == 0) && IsEqual(f0->sh.opacity, f->sh.opacity);
			temperatureE = temperatureE && IsEqual(f0->sh.temperature, f->sh.temperature);
			flowE = flowE && f0->userOutgassing.compare(f->userOutgassing) == 0 && IsEqual(f0->sh.flow, f->sh.flow);
			flowAreaE = flowAreaE && IsEqual(f0->sh.flow / f0Area, f->sh.flow / fArea, 1e-20);
			is2sidedE = is2sidedE && (f0->sh.is2sided == f->sh.is2sided);
			desorbTypeE = desorbTypeE && (f0->sh.desorbType == f->sh.desorbType);
			desorbTypeNE = desorbTypeNE && IsEqual(f0->sh.desorbTypeN, f->sh.desorbTypeN);
			recordE = recordE && (f0->sh.profileType == f->sh.profileType);  //profiles
			area += fArea;
		}

		if (nbSel == 1)
			sprintf(tmp, "Selected Facet (#%d)", selection[0] + 1);
		else
			sprintf(tmp, "Selected Facet (%d selected)", count);

		// Old STR compatibility
		//if (stickingE && f0->sh.superDest) stickingE = FALSE;

		facetPanel->SetTitle(tmp);
		if (count > 1) facetAreaLabel->SetText("Sum Area (cm\262):");
		else facetAreaLabel->SetText("Area (cm\262):");
		facetArea->SetText(area);
		if (stickingE) {
			if (f0->userSticking.length() == 0)
				facetSticking->SetText(f0->sh.sticking);
			else facetSticking->SetText(f0->userSticking.c_str());
		}
		else facetSticking->SetText("...");

		if (opacityE) {
			if (f0->userOpacity.length() == 0)
				facetOpacity->SetText(f0->sh.opacity);
			else facetOpacity->SetText(f0->userOpacity.c_str());
		}
		else facetOpacity->SetText("...");

		if (temperatureE) facetTemperature->SetText(f0->sh.temperature); else facetTemperature->SetText("...");
		if (is2sidedE) facetSideType->SetSelectedIndex(f0->sh.is2sided); else facetSideType->SetSelectedValue("...");
		if (desorbTypeNE) facetDesTypeN->SetText(f0->sh.desorbTypeN); else facetDesTypeN->SetText("...");
		if (recordE) facetRecType->SetSelectedIndex(f0->sh.profileType); else facetRecType->SetSelectedValue("...");

		if (count == 1) {
			facetPumping->SetEditable(TRUE);
			calcFlow();
		}
		else{
			facetPumping->SetEditable(FALSE);
			facetPumping->SetText("...");
		}

		if (desorbTypeE) {
			facetDesType->SetSelectedIndex(f0->sh.desorbType);
			if (f0->sh.desorbType > DES_NONE) { //There is some desorption

				//facetFlow->SetEnabled(TRUE);
				facetFIAreaLabel->SetEnabled(TRUE);
				facetFlowArea->SetEditable(TRUE);
				facetFILabel->SetEnabled(TRUE);
				facetFlow->SetEditable(TRUE);
				if (flowE) {
					if (f0->userOutgassing.length() == 0)
						facetFlow->SetText(f0->sh.flow*10.00); //10: Pa*m3/sec -> mbar*l/s
					else facetFlow->SetText(f0->userOutgassing);
				}
				else facetFlow->SetText("...");
				if (flowAreaE) facetFlowArea->SetText(f0->sh.flow / f0Area*10.00); else facetFlowArea->SetText("...");
				if (f0->sh.desorbType == 3) {
					facetDesTypeN->SetEditable(TRUE);
					if (desorbTypeNE) facetDesTypeN->SetText(f0->sh.desorbTypeN); else facetDesTypeN->SetText("...");
				}
				else {
					facetDesTypeN->SetText("");
					facetDesTypeN->SetEditable(FALSE);
				};

			}
			else { //No desorption
				facetFILabel->SetEnabled(FALSE);
				facetFlow->SetEditable(FALSE);
				facetFIAreaLabel->SetEnabled(FALSE);
				facetFlowArea->SetEditable(FALSE);
				facetDesTypeN->SetText("");
				facetDesTypeN->SetEditable(FALSE);
				facetFlow->SetText("");
				facetFlowArea->SetText("");
			}
		}
		else { //Mixed state
			facetDesType->SetSelectedValue("...");
			facetFILabel->SetEnabled(FALSE);
			facetFlow->SetEditable(FALSE);
			facetFIAreaLabel->SetEnabled(FALSE);
			facetFlowArea->SetEditable(FALSE);
			facetDesTypeN->SetText("");
			facetDesTypeN->SetEditable(FALSE);
			facetFlow->SetText("");
			facetFlowArea->SetText("");
		}

		if (facetMesh) facetMesh->Refresh(count, selection); //Refresh advanced facet parameters panel
		if (updateSelection) {

			if (nbSel > 1000 || geom->GetNbFacet() > 50000) { //If it would take too much time to look up every selected facet in the list
				facetList->ReOrder();
				facetList->SetSelectedRows(selection, nbSel, FALSE);
			}
			else {
				facetList->SetSelectedRows(selection, nbSel, TRUE);
			}
			facetList->lastRowSel = -1;
		}

		free(selection);

		//Enabled->Editable
		facetSticking->SetEditable(TRUE);
		facetOpacity->SetEditable(TRUE);
		facetTemperature->SetEditable(TRUE);
		facetSideType->SetEditable(TRUE);
		facetDesType->SetEditable(TRUE);
		facetRecType->SetEditable(TRUE);
		facetApplyBtn->SetEnabled(FALSE);
	}
	else {
		ClearFacetParams();
		if (facetMesh) facetMesh->Refresh(0, NULL); //Clear
		if (updateSelection) facetList->ClearSelection();
	}

	if (facetDetails) facetDetails->Update();
	if (facetCoordinates) facetCoordinates->UpdateFromSelection();
	if (texturePlotter) texturePlotter->Update(m_fTime, TRUE);
	if (outgassingMap) outgassingMap->Update(m_fTime, TRUE);
}

/*
void MolFlow::LogProfile() {

Geometry *geom = worker.GetGeometry();
BYTE *buffer = worker.GetHits();
if(!buffer) return;

char filename[256];
sprintf(filename,"C:\\Temp\\dataS%d.txt",nbSt);
FILE *file = fopen(filename,"a");

SHGHITS *gHits = (SHGHITS *)buffer;
double nbAbs = (double)gHits->total.hit.nbAbsorbed;
double nbDes = (double)gHits->total.hit.nbDesorbed;
double nbHit = (double)gHits->total.hit.nbHit;

fprintf(file,"Time:%s Sticking=%g Des=%g\n",FormatTime(worker.simuTime),(double)nbSt/10.0,nbDes);

// Volatile profile
int nb = geom->GetNbFacet();
for(int j=0;j<nb;j++) {
Facet *f = geom->GetFacet(j);
if( f->sh.isVolatile ) {
SHHITS *fCount = (SHHITS *)(buffer + f->sh.hitOffset);
double z = geom->GetVertex(f->indices[0])->z;
fprintf(file,"%g %.10g\n",z,(double)(fCount->hit.nbAbsorbed)/nbDes);
}
}

// Last
Facet *f = geom->GetFacet(28);
SHHITS *fCount = (SHHITS *)(buffer + f->sh.hitOffset);
double fnbAbs = (double)fCount->hit.nbAbsorbed;
fprintf(file,"1000 %.10g\n",fnbAbs/nbDes);

fclose(file);

worker.ReleaseHits();

}
*/

// ----------------------------------------------------------------------------
void MolFlow::ResetAutoSaveTimer() {
	if (autoSaveSimuOnly) lastSaveTimeSimu = worker.simuTime + (m_fTime - worker.startTime);
	else lastSaveTime = m_fTime;
}

//-----------------------------------------------------------------------------
BOOL MolFlow::AutoSave(BOOL crashSave) {
	if (!changedSinceSave) return TRUE;
	GLProgress *progressDlg2 = new GLProgress("Peforming autosave...", "Please wait");
	progressDlg2->SetProgress(0.0);
	progressDlg2->SetVisible(TRUE);
	//GLWindowManager::Repaint();
	char CWD[MAX_PATH];
	_getcwd(CWD, MAX_PATH);

	std::string shortFn(worker.GetShortFileName());
	std::string newAutosaveFilename = "Molflow_Autosave";
	if (shortFn != "") newAutosaveFilename += "(" + shortFn + ")";
	newAutosaveFilename += ".zip";
	char fn[1024];
	strcpy(fn, newAutosaveFilename.c_str());
	try {
		worker.SaveGeometry(fn, progressDlg2, FALSE, FALSE, TRUE, crashSave);
		//Success:
		if (autosaveFilename != "" && autosaveFilename != newAutosaveFilename) remove(autosaveFilename.c_str());
		autosaveFilename = newAutosaveFilename;
		ResetAutoSaveTimer(); //deduct saving time from interval
	}
	catch (Error &e) {
		//delete fn;
		char errMsg[512];
		sprintf(errMsg, "%s\nFile:%s", e.GetMsg(), fn);
		GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
		progressDlg2->SetVisible(FALSE);
		SAFE_DELETE(progressDlg2);
		ResetAutoSaveTimer();
		return FALSE;
	}
	//lastSaveTime=(worker.simuTime+(m_fTime-worker.startTime));
	progressDlg2->SetVisible(FALSE);
	SAFE_DELETE(progressDlg2);
	return TRUE;
}

//-------------------------------------------------------------------------------
void MolFlow::CheckForRecovery() {
	// Check for autosave files in current dir.
	intptr_t file;
	_finddata_t filedata;
	file = _findfirst("Molflow_Autosave*.zip", &filedata);
	if (file != -1)
	{
		do
		{
			std::ostringstream msg;
			msg << "Autosave file found:\n" << filedata.name << "\n";
			int rep = RecoveryDialog::Display(msg.str().c_str(), "Autosave recovery", GLDLG_LOAD | GLDLG_SKIP, GLDLG_DELETE);
			if (rep == GLDLG_LOAD) {
				LoadFile(filedata.name);
				RemoveRecent(filedata.name);
			}
			else if (rep == GLDLG_CANCEL) return;
			else if (rep == GLDLG_SKIP) continue;
			else if (rep == GLDLG_DELETE) remove(filedata.name);
		} while (_findnext(file, &filedata) == 0);
	}
	_findclose(file);
}


//-----------------------------------------------------------------------------
// Name: FrameMove()
// Desc: Called once per frame, the call is the entry point for animating
//       the scene.
//-----------------------------------------------------------------------------
int MolFlow::FrameMove()
{
	char tmp[256];
	Geometry *geom = worker.GetGeometry();

	//Autosave routines
	BOOL timeForAutoSave = FALSE;
	if (geom->IsLoaded()) {
		if (autoSaveSimuOnly) {
			if (worker.running) {
				if (((worker.simuTime + (m_fTime - worker.startTime)) - lastSaveTimeSimu) >= (float)autoSaveFrequency*60.0f) {
					timeForAutoSave = TRUE;
				}
			}
		}
		else {
			if ((m_fTime - lastSaveTime) >= (float)autoSaveFrequency*60.0f) {
				timeForAutoSave = TRUE;
			}
		}
	}
	
	if (globalSettings) globalSettings->SMPUpdate(m_fTime);

	if (worker.running) {
		if (frameMoveRequested || autoFrameMove && (m_fTime - lastUpdate >= 1.0f)) {
			forceFrameMoveButton->SetEnabled(FALSE);
			forceFrameMoveButton->SetText("Updating...");
			//forceFrameMoveButton->Paint();
			GLWindowManager::Repaint();
			frameMoveRequested = FALSE;

			// Update hits
			try {
				worker.Update(m_fTime);
			}
			catch (Error &e) {
				GLMessageBox::Display((char *)e.GetMsg(), "Error (Stop)", GLDLG_OK, GLDLG_ICONERROR);
			}
			// Simulation monitoring
			if (profilePlotter) profilePlotter->Update(m_fTime);
			if (pressureEvolution) pressureEvolution->Update(m_fTime);
			if (timewisePlotter) timewisePlotter->Update(m_fTime);
			if (texturePlotter) texturePlotter->Update(m_fTime);
			//if(facetDetails) facetDetails->Update();
			if (textureSettings) textureSettings->Update();
			// Facet parameters and hits
			

			// Formulas
			if (autoUpdateFormulas) UpdateFormula();

			lastUpdate = GetTick(); //changed from m_fTime: include update duration

			// Update timing measurements
			if (worker.nbHit != lastNbHit || worker.nbDesorption != lastNbDes) {
				double dTime = (double)(m_fTime - lastMeasTime);
				hps = (double)(worker.nbHit - lastNbHit) / dTime;
				dps = (double)(worker.nbDesorption - lastNbDes) / dTime;
				if (lastHps != 0.0) {
					hps = 0.2*(hps)+0.8*lastHps;
					dps = 0.2*(dps)+0.8*lastDps;
				}
				lastHps = hps;
				lastDps = dps;
				lastNbHit = worker.nbHit;
				lastNbDes = worker.nbDesorption;
				lastMeasTime = m_fTime;
			}

		}
		if (worker.calcAC) {
			sprintf(tmp, "Calc AC: %s (%zd %%)", FormatTime(worker.simuTime + (m_fTime - worker.startTime)),
				worker.calcACprg);
		}
		else {
			sprintf(tmp, "Running: %s", FormatTime(worker.simuTime + (m_fTime - worker.startTime)));
		}
		sTime->SetText(tmp);

		forceFrameMoveButton->SetEnabled(!autoFrameMove);
		forceFrameMoveButton->SetText("Update");
	}
	else {
		if (worker.simuTime > 0.0) {
			hps = (double)(worker.nbHit - nbHitStart) / worker.simuTime;
			dps = (double)(worker.nbDesorption - nbDesStart) / worker.simuTime;
		}
		else {
			hps = 0.0;
			dps = 0.0;
		}
		sprintf(tmp, "Stopped: %s", FormatTime(worker.simuTime));
		sTime->SetText(tmp);
	}
	
	if (viewer[0]->SelectionChanged() ||
		viewer[1]->SelectionChanged() ||
		viewer[2]->SelectionChanged() ||
		viewer[3]->SelectionChanged()) {
		UpdateFacetParams(TRUE);
	}
	UpdateFacetHits();

	//Autosave
	if (timeForAutoSave) AutoSave();

	if ((m_fTime - worker.startTime <= 2.0f) && worker.running) {
		hitNumber->SetText("Starting...");
		desNumber->SetText("Starting...");

	}
	else {

		if (worker.mode != AC_MODE) {
			sprintf(tmp, "%s (%s)", FormatInt(worker.nbHit, "hit"), FormatPS(hps, "hit"));
			hitNumber->SetText(tmp);
		}
		else {
			hitNumber->SetText("");
		}
		sprintf(tmp, "%s (%s)", FormatInt(worker.nbDesorption, "des"), FormatPS(dps, "des"));
		desNumber->SetText(tmp);
	}


	if (worker.nbLeakTotal) {
		sprintf(tmp, "%g (%.4f%%)", (double)worker.nbLeakTotal, (double)(worker.nbLeakTotal)*100.0 / (double)worker.nbDesorption);
		leakNumber->SetText(tmp);
	}
	else {
		leakNumber->SetText("None");
	}

	resetSimu->SetEnabled(!worker.running&&worker.nbDesorption > 0);

	if (worker.running) {
		startSimu->SetText("Pause");
		//startSimu->SetFontColor(255, 204, 0);
	}
	else if (worker.nbHit > 0) {
		startSimu->SetText("Resume");
		//startSimu->SetFontColor(0, 140, 0);
	}
	else {
		startSimu->SetText("Begin");
		//startSimu->SetFontColor(0, 140, 0);
	}

	// Sleep a bit to avoid unwanted CPU load
	if (viewer[0]->IsDragging() ||
		viewer[1]->IsDragging() ||
		viewer[2]->IsDragging() ||
		viewer[3]->IsDragging() || !worker.running)
		SDL_Delay(32);
	else
		SDL_Delay(60);

	return GL_OK;
}

// ----------------------------------------------------------------

void MolFlow::UpdateFacetHits(BOOL allRows) {
	char tmp[256];
	Geometry *geom = worker.GetGeometry();

	try{
		// Facet list
		if (geom->IsLoaded()) {

			int sR, eR;
			if (allRows)
			{
				sR = 0;
				eR = facetList->GetNbRow() - 1;
			}
			else
			{
				facetList->GetVisibleRows(&sR, &eR);
			}

			
			if (worker.displayedMoment == 0) {
				int colors[] = { COLOR_BLACK, COLOR_BLACK, COLOR_BLACK, COLOR_BLACK };
				facetList->SetColumnColors(colors);
			}
			else
			{
				int colors[] = { COLOR_BLACK, COLOR_BLUE, COLOR_BLUE, COLOR_BLUE };
				facetList->SetColumnColors(colors);
			}
			

			for (int i = sR; i <= eR; i++) {
				int facetId = facetList->GetValueInt(i, 0) - 1;
				if (facetId == -2) facetId = i;
				if (i >= geom->GetNbFacet()) {
					char errMsg[512];
					sprintf(errMsg, "Molflow::UpdateFacetHits()\nError while updating facet hits. Was looking for facet #%d in list.\nMolflow will now autosave and crash.", i + 1);
					GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
					AutoSave();
				}
				Facet *f = geom->GetFacet(facetId);
				sprintf(tmp, "%d", facetId + 1);
				facetList->SetValueAt(0, i, tmp);
				switch (modeCombo->GetSelectedIndex()) {
				case MC_MODE:
					facetList->SetColumnLabel(1, "Hits");
					sprintf(tmp, "%I64d", f->counterCache.hit.nbHit);
					facetList->SetValueAt(1, i, tmp);
					sprintf(tmp, "%I64d", f->counterCache.hit.nbDesorbed);
					facetList->SetValueAt(2, i, tmp);
					sprintf(tmp, "%I64d", f->counterCache.hit.nbAbsorbed);
					facetList->SetValueAt(3, i, tmp);
					break;
				case AC_MODE:
					facetList->SetColumnLabel(1, "Density");
					sprintf(tmp, "%g", f->counterCache.density.value);
					facetList->SetValueAt(1, i, tmp);

					sprintf(tmp, "%g", f->counterCache.density.desorbed);
					facetList->SetValueAt(2, i, tmp);
					sprintf(tmp, "%g", f->counterCache.density.absorbed);
					facetList->SetValueAt(3, i, tmp);
					break;
				}
			}

		}
	}
	catch (Error &e) {
		char errMsg[512];
		sprintf(errMsg, "%s\nError while updating facet hits", e.GetMsg());
		GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
	}

}


//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Name: RestoreDeviceObjects()
// Desc: Initialize scene objects.
//-----------------------------------------------------------------------------
int MolFlow::RestoreDeviceObjects()
{
	RestoreDeviceObjects_shared();

	//Different Molflow implementations:
	RVALIDATE_DLG(facetMesh);
	RVALIDATE_DLG(facetDetails);
	RVALIDATE_DLG(smartSelection);
	RVALIDATE_DLG(viewer3DSettings);
	RVALIDATE_DLG(textureSettings);
	RVALIDATE_DLG(globalSettings);
	RVALIDATE_DLG(profilePlotter);
	RVALIDATE_DLG(texturePlotter);

	//Molflow only:
	RVALIDATE_DLG(importDesorption);
	RVALIDATE_DLG(timeSettings);
	RVALIDATE_DLG(movement);
	RVALIDATE_DLG(outgassingMap);
	RVALIDATE_DLG(parameterEditor);
	RVALIDATE_DLG(pressureEvolution);
	RVALIDATE_DLG(timewisePlotter);

	return GL_OK;
}

//-----------------------------------------------------------------------------
// Name: InvalidateDeviceObjects()
// Desc: Free all alocated resource
//-----------------------------------------------------------------------------

int MolFlow::InvalidateDeviceObjects()
{
	InvalidateDeviceObjects_shared();

	//Different Molflow implementations:
	IVALIDATE_DLG(facetMesh);
	IVALIDATE_DLG(facetDetails);
	IVALIDATE_DLG(smartSelection);
	IVALIDATE_DLG(viewer3DSettings);
	IVALIDATE_DLG(textureSettings);
	IVALIDATE_DLG(globalSettings);
	IVALIDATE_DLG(profilePlotter);
	IVALIDATE_DLG(texturePlotter);

	//Molflow only:
	IVALIDATE_DLG(importDesorption);
	IVALIDATE_DLG(timeSettings);
	IVALIDATE_DLG(movement);
	IVALIDATE_DLG(outgassingMap);
	IVALIDATE_DLG(parameterEditor);
	IVALIDATE_DLG(pressureEvolution);
	IVALIDATE_DLG(timewisePlotter);

	return GL_OK;
}



void MolFlow::SaveFileAs() {

	FILENAME *fn = GLFileBox::SaveFile(currentDir, worker.GetShortFileName(), "Save File", fileSFilters, 0);

	GLProgress *progressDlg2 = new GLProgress("Saving file...", "Please wait");
	progressDlg2->SetProgress(0.0);
	progressDlg2->SetVisible(TRUE);
	//GLWindowManager::Repaint();  
	if (fn) {

		try {

			worker.SaveGeometry(fn->fullName, progressDlg2);
			ResetAutoSaveTimer();
			changedSinceSave = FALSE;
			UpdateCurrentDir(worker.fullFileName);
			UpdateTitle();
			AddRecent(worker.fullFileName);
		}
		catch (Error &e) {
			char errMsg[512];
			sprintf(errMsg, "%s\nFile:%s", e.GetMsg(), fn->fullName);
			GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
			RemoveRecent(fn->fullName);
		}

	}

	progressDlg2->SetVisible(FALSE);
	SAFE_DELETE(progressDlg2);
}

void MolFlow::ExportTextures(int grouping, int mode) {

	Geometry *geom = worker.GetGeometry();
	if (geom->GetNbSelected() == 0) {
		GLMessageBox::Display("Empty selection", "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}
	if (!worker.IsDpInitialized()) {
		GLMessageBox::Display("Worker Dataport not initialized yet", "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}

	FILENAME *fn = GLFileBox::SaveFile(currentDir, NULL, "Save File", fileTexFilters, 0);

	if (fn) {

		try {
			worker.ExportTextures(fn->fullName, grouping, mode, TRUE, TRUE);
			//UpdateCurrentDir(fn->fullName);
			//UpdateTitle();
		}
		catch (Error &e) {
			char errMsg[512];
			sprintf(errMsg, "%s\nFile:%s", e.GetMsg(), fn->fullName);
			GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
		}

	}

}

void MolFlow::ExportProfiles() {

	Geometry *geom = worker.GetGeometry();
	if (geom->GetNbSelected() == 0) {
		GLMessageBox::Display("Empty selection", "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}
	if (!worker.IsDpInitialized()) {
		GLMessageBox::Display("Worker Dataport not initialized yet", "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}

	FILENAME *fn = GLFileBox::SaveFile(currentDir, NULL, "Save File", fileProfFilters, 0);

	if (fn) {

		try {
			worker.ExportProfiles(fn->fullName);
		}
		catch (Error &e) {
			char errMsg[512];
			sprintf(errMsg, "%s\nFile:%s", e.GetMsg(), fn->fullName);
			GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
		}

	}

}

void MolFlow::ImportDesorption_DES() {

	FILENAME *fn = GLFileBox::OpenFile(currentDir, NULL, "Import desorption File", fileDesFilters, 0);

	if (fn) {

		try {
			worker.ImportDesorption_DES(fn->fullName);
		}
		catch (Error &e) {
			char errMsg[512];
			sprintf(errMsg, "%s\nFile:%s", e.GetMsg(), fn->fullName);
			GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
		}
		worker.GetGeometry()->CalcTotalOutGassing();
		UpdateFacetParams();
	}

}

void MolFlow::SaveFile() {
	if (strlen(worker.fullFileName) > 0){

		GLProgress *progressDlg2 = new GLProgress("Saving...", "Please wait");
		progressDlg2->SetProgress(0.5);
		progressDlg2->SetVisible(TRUE);
		//GLWindowManager::Repaint();

		try {
			worker.SaveGeometry(worker.fullFileName, progressDlg2, FALSE);
			ResetAutoSaveTimer();
		}
		catch (Error &e) {
			char errMsg[512];
			sprintf(errMsg, "%s\nFile:%s", e.GetMsg(), worker.GetFileName());
			GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
		}
		progressDlg2->SetVisible(FALSE);
		SAFE_DELETE(progressDlg2);
		changedSinceSave = FALSE;

	}
	else SaveFileAs();
}

//-----------------------------------------------------------------------------



void MolFlow::LoadFile(char *fName) {


	char fullName[1024];
	char shortName[512];
	strcpy(fullName, "");

	if (fName == NULL) {
		FILENAME *fn = GLFileBox::OpenFile(currentDir, NULL, "Open file", fileLFilters, 2);
		if (fn)
			strcpy(fullName, fn->fullName);
	}
	else {
		strcpy(fullName, fName);
	}

	GLProgress *progressDlg2 = new GLProgress("Preparing to load file...", "Please wait");
	progressDlg2->SetVisible(TRUE);
	progressDlg2->SetProgress(0.0);
	//GLWindowManager::Repaint();

	if (strlen(fullName) == 0) {
		progressDlg2->SetVisible(FALSE);
		SAFE_DELETE(progressDlg2);
		return;
	}


	char *lPart = strrchr(fullName, '\\');
	if (lPart) strcpy(shortName, lPart + 1);
	else strcpy(shortName, fullName);

	try {
		ClearFormula();
		ClearParameters();
		ClearAllSelections();
		ClearAllViews();		

		worker.LoadGeometry(fullName);

		Geometry *geom = worker.GetGeometry();

		// Default initialisation
		viewer[0]->SetWorker(&worker);
		viewer[1]->SetWorker(&worker);
		viewer[2]->SetWorker(&worker);
		viewer[3]->SetWorker(&worker);
		//UpdateModelParams();
		startSimu->SetEnabled(TRUE);
		compACBtn->SetEnabled(modeCombo->GetSelectedIndex() == 1);
		singleACBtn->SetEnabled(modeCombo->GetSelectedIndex() == 1);
		//resetSimu->SetEnabled(TRUE);
		ClearFacetParams();
		nbDesStart = worker.nbDesorption;
		nbHitStart = worker.nbHit;
		AddRecent(fullName);
		geom->viewStruct = -1;


		UpdateStructMenu();
		UpdateCurrentDir(fullName);

		// Check non simple polygon
		progressDlg2->SetMessage("Checking for non simple polygons...");

		geom->CheckCollinear();
		geom->CheckNonSimple();
		geom->CheckIsolatedVertex();
		// Set up view
		// Default
		viewer[0]->SetProjection(ORTHOGRAPHIC_PROJ);
		viewer[0]->ToFrontView();
		viewer[1]->SetProjection(ORTHOGRAPHIC_PROJ);
		viewer[1]->ToTopView();
		viewer[2]->SetProjection(ORTHOGRAPHIC_PROJ);
		viewer[2]->ToSideView();
		viewer[3]->SetProjection(PERSPECTIVE_PROJ);
		viewer[3]->ToFrontView();
		SelectViewer(0);

		ResetAutoSaveTimer();
		if (profilePlotter) profilePlotter->Refresh();
		if (pressureEvolution) pressureEvolution->Refresh();
		if (textureSettings) textureSettings->Update();
		if (texturePlotter) texturePlotter->Update(m_fTime, TRUE);
		if (outgassingMap) outgassingMap->Update(m_fTime, TRUE);
		if (facetDetails) facetDetails->Update();
		if (facetCoordinates) facetCoordinates->UpdateFromSelection();
		if (vertexCoordinates) vertexCoordinates->Update();
		if (movement) movement->Update();
		if (globalSettings && globalSettings->IsVisible()) globalSettings->Update();
		if (timewisePlotter) timewisePlotter->Refresh();
		if (timeSettings) mApp->timeSettings->RefreshMoments();
		if (momentsEditor) mApp->momentsEditor->Refresh();
		if (parameterEditor) mApp->parameterEditor->UpdateCombo();

	}
	catch (Error &e) {

		char errMsg[512];
		sprintf(errMsg, "%s\nFile:%s", e.GetMsg(), shortName);
		GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
		RemoveRecent(fName);

	}
	progressDlg2->SetVisible(FALSE);
	SAFE_DELETE(progressDlg2);
	changedSinceSave = FALSE;
}

//-----------------------------------------------------------------------------

void MolFlow::InsertGeometry(BOOL newStr, char *fName) {
	if (!AskToReset()) return;
	ResetSimulation(FALSE);

	char fullName[1024];
	char shortName[512];
	strcpy(fullName, "");

	if (fName == NULL) {
		FILENAME *fn = GLFileBox::OpenFile(currentDir, NULL, "Open File", fileInsFilters, 0);
		if (fn)
			strcpy(fullName, fn->fullName);
	}
	else {
		strcpy(fullName, fName);
	}

	GLProgress *progressDlg2 = new GLProgress("Loading file...", "Please wait");
	progressDlg2->SetVisible(TRUE);
	progressDlg2->SetProgress(0.0);
	//GLWindowManager::Repaint();

	if (strlen(fullName) == 0) {
		progressDlg2->SetVisible(FALSE);
		SAFE_DELETE(progressDlg2);
		return;
	}


	char *lPart = strrchr(fullName, '\\');
	if (lPart) strcpy(shortName, lPart + 1);
	else strcpy(shortName, fullName);

	try {

		//worker.InsertGeometry(newStr, fullName);
		worker.LoadGeometry(fullName, TRUE, newStr);

		Geometry *geom = worker.GetGeometry();
		geom->CalcTotalOutGassing();
		/*
		// Default initialisation
		viewer[0]->SetWorker(&worker);
		viewer[1]->SetWorker(&worker);
		viewer[2]->SetWorker(&worker);
		viewer[3]->SetWorker(&worker);*/
		//UpdateModelParams();
		startSimu->SetEnabled(TRUE);

		compACBtn->SetEnabled(modeCombo->GetSelectedIndex() == 1);
		singleACBtn->SetEnabled(modeCombo->GetSelectedIndex() == 1);
		//resetSimu->SetEnabled(TRUE);
		//ClearFacetParams();
		//nbDesStart = worker.nbDesorption;
		//nbHitStart = worker.nbHit;
		AddRecent(fullName);
		geom->viewStruct = -1;


		//worker.LoadTextures(fullName);
		UpdateStructMenu();
		if (profilePlotter) profilePlotter->Reset();

		if (pressureEvolution) pressureEvolution->Reset();
		if (timewisePlotter) timewisePlotter->Reset();
		//UpdateCurrentDir(fullName);

		geom->CheckCollinear();
		geom->CheckNonSimple();
		geom->CheckIsolatedVertex();

		/*
		// Set up view
		// Default
		viewer[0]->SetProjection(ORTHOGRAPHIC_PROJ);
		viewer[0]->ToFrontView();
		viewer[1]->SetProjection(ORTHOGRAPHIC_PROJ);
		viewer[1]->ToTopView();
		viewer[2]->SetProjection(ORTHOGRAPHIC_PROJ);
		viewer[2]->ToSideView();
		viewer[3]->SetProjection(PERSPECTIVE_PROJ);
		viewer[3]->ToFrontView();
		SelectViewer(0);
		*/
		if (profilePlotter) profilePlotter->Refresh();

		if (pressureEvolution) pressureEvolution->Refresh();
		if (timewisePlotter) timewisePlotter->Refresh();
		if (texturePlotter) texturePlotter->Update(m_fTime, TRUE);
		if (outgassingMap) outgassingMap->Update(m_fTime, TRUE);
		if (facetDetails) facetDetails->Update();
		if (facetCoordinates) facetCoordinates->UpdateFromSelection();
		if (vertexCoordinates) vertexCoordinates->Update();


	}
	catch (Error &e) {

		char errMsg[512];
		sprintf(errMsg, "%s\nFile:%s", e.GetMsg(), shortName);
		GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
		RemoveRecent(fName);

	}
	progressDlg2->SetVisible(FALSE);
	SAFE_DELETE(progressDlg2);
	changedSinceSave = TRUE;
}

void MolFlow::ClearParameters() {
	worker.parameters = std::vector<Parameter>();
	if (parameterEditor) parameterEditor->UpdateCombo();
}

void MolFlow::StartStopSimulation() {

	//if(nbSt<=10) BuildPipeStick((double)nbSt/10);
	//else         return;

	if (!(worker.nbHit>0) && !worker.calcConstantFlow && worker.moments.size() == 0) {
		BOOL ok = GLMessageBox::Display("Warning: in the Moments Editor, the option \"Calculate constant flow\" is disabled.\n"
			"This is useful for time-dependent simulations.\n"
			"However, you didn't define any moments, suggesting you're using steady-state mode.\n"
			"\nDo you want to continue?\n", "Strange time settings", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONWARNING) == GLDLG_OK;
		if (!ok) return;
	}

	worker.StartStop(m_fTime, modeCombo->GetSelectedIndex());
	if (profilePlotter) profilePlotter->Update(m_fTime, TRUE);
	if (pressureEvolution) pressureEvolution->Update(m_fTime, TRUE);
	if (timewisePlotter) timewisePlotter->Update(m_fTime, TRUE);
	if (texturePlotter) texturePlotter->Update(m_fTime, TRUE);
	if (autoUpdateFormulas) UpdateFormula();

	// Frame rate measurement
	lastMeasTime = m_fTime;
	dps = 0.0;
	hps = 0.0;
	lastHps = hps;
	lastDps = dps;
	lastNbHit = worker.nbHit;
	lastNbDes = worker.nbDesorption;
	lastUpdate = 0.0;

}

void MolFlow::DoEvents(BOOL forced)
{
	static float lastExec = 0;
	int time = SDL_GetTicks();
	if (forced || (time - lastExec>333)) { //Don't check for inputs more than 3 times a second
		SDL_Event sdlEvent;
		SDL_PollEvent(&sdlEvent);
		mApp->UpdateEventCount(&sdlEvent);
		/*if (GLWindowManager::ManageEvent(&sdlEvent)) {
			// Relay to GLApp EventProc
			mApp->EventProc(&sdlEvent);
		}*/
		GLWindowManager::ManageEvent(&sdlEvent);
		GLWindowManager::Repaint();
		GLToolkit::CheckGLErrors("GLApplication::Paint()");
		lastExec = time;
	}
}

//-----------------------------------------------------------------------------
// Name: EventProc()
// Desc: Message proc function to handle key and mouse input
//-----------------------------------------------------------------------------
void MolFlow::ProcessMessage(GLComponent *src, int message)
{
	Geometry *geom = worker.GetGeometry();
	char *input;
	char tmp[128];
	switch (message) {

		//MENU --------------------------------------------------------------------
	case MSG_MENU:
		switch (src->GetId()) {
		case MENU_FILE_LOAD:
			if (AskToSave()) {
				if (worker.running) worker.Stop_Public();
				LoadFile();
			}
			break;
		case MENU_FILE_IMPORTDES_SYN:

			if (geom->IsLoaded()) {
				Geometry *geom = worker.GetGeometry();
				if (!importDesorption) importDesorption = new ImportDesorption();
				importDesorption->SetGeometry(geom, &worker);
				importDesorption->SetVisible(TRUE);
			}
			else GLMessageBox::Display("No geometry loaded.", "No geometry", GLDLG_OK, GLDLG_ICONERROR);
			break;
		case MENU_FILE_IMPORTDES_DES:
			ImportDesorption_DES();
			break;
		case MENU_FILE_INSERTGEO:
			if (geom->IsLoaded()) {
				if (worker.running) worker.Stop_Public();
				InsertGeometry(FALSE);
			}
			else GLMessageBox::Display("No geometry loaded.", "No geometry", GLDLG_OK, GLDLG_ICONERROR);
			break;
		case MENU_FILE_INSERTGEO_NEWSTR:
			if (geom->IsLoaded()) {
				if (worker.running) worker.Stop_Public();
				InsertGeometry(TRUE);
			}
			else GLMessageBox::Display("No geometry loaded.", "No geometry", GLDLG_OK, GLDLG_ICONERROR);
			break;
		case MENU_FILE_SAVEAS:
			if (geom->IsLoaded()) {
				SaveFileAs();
			}
			else GLMessageBox::Display("No geometry loaded.", "No geometry", GLDLG_OK, GLDLG_ICONERROR);
			break;
		case MENU_FILE_EXPORT_SELECTION:
			ExportSelection();
			break;

		case MENU_FILE_EXPORTTEXTURE_AREA:
			ExportTextures(0, 0); break;
		case MENU_FILE_EXPORTTEXTURE_MCHITS:
			ExportTextures(0, 1); break;
		case MENU_FILE_EXPORTTEXTURE_IMPINGEMENT:
			ExportTextures(0, 2); break;
		case MENU_FILE_EXPORTTEXTURE_PART_DENSITY:
			ExportTextures(0, 3); break;
		case MENU_FILE_EXPORTTEXTURE_GAS_DENSITY:
			ExportTextures(0, 4); break;
		case MENU_FILE_EXPORTTEXTURE_PRESSURE:
			ExportTextures(0, 5); break;
		case MENU_FILE_EXPORTTEXTURE_AVG_V:
			ExportTextures(0, 6); break;
		case MENU_FILE_EXPORTTEXTURE_V_VECTOR:
			ExportTextures(1, 7); break;
		case MENU_FILE_EXPORTTEXTURE_N_VECTORS:
			ExportTextures(0, 8); break;
		case MENU_FILE_EXPORTTEXTURE_AREA_COORD:
			ExportTextures(1, 0); break;
		case MENU_FILE_EXPORTTEXTURE_MCHITS_COORD:
			ExportTextures(1, 1); break;
		case MENU_FILE_EXPORTTEXTURE_IMPINGEMENT_COORD:
			ExportTextures(1, 2); break;
		case MENU_FILE_EXPORTTEXTURE_PART_DENSITY_COORD:
			ExportTextures(1, 3); break;
		case MENU_FILE_EXPORTTEXTURE_GAS_DENSITY_COORD:
			ExportTextures(1, 4); break;
		case MENU_FILE_EXPORTTEXTURE_PRESSURE_COORD:
			ExportTextures(1, 5); break;
		case MENU_FILE_EXPORTTEXTURE_AVG_V_COORD:
			ExportTextures(1, 6); break;
		case MENU_FILE_EXPORTTEXTURE_V_VECTOR_COORD:
			ExportTextures(1, 7); break;
		case MENU_FILE_EXPORTTEXTURE_N_VECTORS_COORD:
			ExportTextures(1, 8); break;

		case MENU_FILE_EXPORTPROFILES:
			ExportProfiles();
			break;

		case MENU_FILE_SAVE:
			if (geom->IsLoaded()) SaveFile();
			else GLMessageBox::Display("No geometry loaded.", "No geometry", GLDLG_OK, GLDLG_ICONERROR);
			break;
		case MENU_FILE_EXIT:
			if (AskToSave()) Exit();
			break;
		case MENU_EDIT_TSCALING:
			if (!textureSettings || !textureSettings->IsVisible()) {
				SAFE_DELETE(textureSettings);
				textureSettings = new TextureSettings();
				textureSettings->Display(&worker, viewer);
			}
			break;
		case MENU_EDIT_ADDFORMULA:
			if (!formulaSettings) formulaSettings = new FormulaSettings();
			AddFormula(formulaSettings->NewFormula());
			break;
		case MENU_EDIT_UPDATEFORMULAS:
			UpdateFormula();
			break;
		case MENU_EDIT_MOVINGPARTS:
			if (!movement) movement = new Movement(geom, &worker);
			movement->Update();
			movement->SetVisible(TRUE);
			break;
		case MENU_EDIT_GLOBALSETTINGS:
			if (!globalSettings) globalSettings = new GlobalSettings(&worker);
			globalSettings->Update();
			globalSettings->SetVisible(TRUE);
			break;
		case MENU_FACET_COLLAPSE:
			if (geom->IsLoaded()) {
				DisplayCollapseDialog();
			}
			else GLMessageBox::Display("No geometry loaded.", "No geometry", GLDLG_OK, GLDLG_ICONERROR);
			break;
		case MENU_FACET_SWAPNORMAL:
			if (AskToReset()) {
				geom->SwapNormal();
				// Send to sub process
				try { worker.Reload(); }
				catch (Error &e) {
					GLMessageBox::Display((char *)e.GetMsg(), "Error", GLDLG_OK, GLDLG_ICONERROR);
				}
			}
			break;
		case MENU_FACET_EXTRUDE:
			if (!extrudeFacet || !extrudeFacet->IsVisible()) {
				SAFE_DELETE(extrudeFacet);
				extrudeFacet = new ExtrudeFacet(geom, &worker);
			}
			extrudeFacet->SetVisible(TRUE);
			break;
			
		case MENU_FACET_SHIFTVERTEX:
			if (AskToReset()) {
				geom->ShiftVertex();
				// Send to sub process
				worker.Reload();
			}
			break;
		case MENU_FACET_COORDINATES:

			if (!facetCoordinates) facetCoordinates = new FacetCoordinates();
			facetCoordinates->Display(&worker);
			break;
		case MENU_FACET_MOVE:
			if (!moveFacet || !moveFacet->IsVisible()) {
				SAFE_DELETE(moveFacet);
				moveFacet = new MoveFacet(geom, &worker);
			}
			moveFacet->SetVisible(TRUE);
			break;
		case MENU_FACET_SCALE:
			if (geom->IsLoaded()) {
				if (!scaleFacet) scaleFacet = new ScaleFacet(geom, &worker);

				scaleFacet->SetVisible(TRUE);

			}
			else GLMessageBox::Display("No geometry loaded.", "No geometry", GLDLG_OK, GLDLG_ICONERROR);
			break;
		case MENU_FACET_MIRROR:
			if (!mirrorFacet) mirrorFacet = new MirrorFacet(geom, &worker);
			mirrorFacet->SetVisible(TRUE);
			break;
		case MENU_FACET_SPLIT:			
			if (!splitFacet || !splitFacet->IsVisible()) {
				SAFE_DELETE(splitFacet);
				splitFacet = new SplitFacet(geom, &worker);
				splitFacet->SetVisible(TRUE);
			}
			break;
		case MENU_FACET_ROTATE:
			if (!rotateFacet) rotateFacet = new RotateFacet(geom, &worker);
			rotateFacet->SetVisible(TRUE);
			break;
		case MENU_FACET_ALIGN:
			if (!alignFacet) alignFacet = new AlignFacet(geom, &worker);
			alignFacet->MemorizeSelection();
			alignFacet->SetVisible(TRUE);
			break;
		case MENU_FACET_PROFPLOTTER:
			if (!profilePlotter) profilePlotter = new ProfilePlotter();
			profilePlotter->Display(&worker);
			break;


		case MENU_TIME_PRESSUREEVOLUTION:
			if (!pressureEvolution) pressureEvolution = new PressureEvolution();
			pressureEvolution->Display(&worker);
			break;
			/*case MENU_FACET_MESH:
				if (!facetMesh) facetMesh = new FacetMesh(&worker);
				facetMesh->SetVisible(!facetMesh->IsVisible());
				break;      */
		case MENU_FACET_TEXPLOTTER:
			if (!texturePlotter) texturePlotter = new TexturePlotter();
			texturePlotter->Display(&worker);
			break;
		case MENU_FACET_OUTGASSINGMAP:
			if (!outgassingMap) outgassingMap = new OutgassingMap();
			outgassingMap->Display(&worker);
			break;
		case MENU_FACET_REMOVESEL:
			if (GLMessageBox::Display("Remove selected facets?", "Question", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONINFO) == GLDLG_OK) {
				if (AskToReset()) {
					if (worker.running) worker.Stop_Public();
					geom->RemoveSelected();
					geom->CalcTotalOutGassing();
					//geom->CheckIsolatedVertex();
					UpdateModelParams();
					if (vertexCoordinates) vertexCoordinates->Update();
					if (facetCoordinates) facetCoordinates->UpdateFromSelection();
					if (profilePlotter) profilePlotter->Refresh();

					if (pressureEvolution) pressureEvolution->Refresh();
					if (timewisePlotter) timewisePlotter->Refresh();
					// Send to sub process
					try { worker.Reload(); }
					catch (Error &e) {
						GLMessageBox::Display((char *)e.GetMsg(), "Error reloading worker", GLDLG_OK, GLDLG_ICONERROR);
					}
				}
			}
			break;
		case MENU_FACET_EXPLODE:
			if (GLMessageBox::Display("Explode selected facet?", "Question", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONINFO) == GLDLG_OK) {
				if (AskToReset()) {
					int err;
					try {
						err = geom->ExplodeSelected();
					}
					catch (Error &e) {
						GLMessageBox::Display((char *)e.GetMsg(), "Error exploding", GLDLG_OK, GLDLG_ICONERROR);
					}
					if (err == -1) {
						GLMessageBox::Display("Empty selection", "Error", GLDLG_OK, GLDLG_ICONERROR);
					}
					else if (err == -2) {
						GLMessageBox::Display("All selected facets must have a mesh with boudary correction enabled", "Error", GLDLG_OK, GLDLG_ICONERROR);
					}
					else if (err == 0) {

						UpdateModelParams();
						UpdateFacetParams(TRUE);
						// Send to sub process
						try { worker.Reload(); }
						catch (Error &e) {
							GLMessageBox::Display((char *)e.GetMsg(), "Error reloading worker", GLDLG_OK, GLDLG_ICONERROR);
						}
					}
				}
			}
			break;
		case MENU_FACET_DETAILS:
			if (facetDetails == NULL) facetDetails = new FacetDetails();
			facetDetails->Display(&worker);
			break;

		case MENU_SELECTION_SMARTSELECTION:
			if (!smartSelection) smartSelection = new SmartSelection(worker.GetGeometry(),&worker);
			smartSelection->SetVisible(TRUE);
			break;
		case MENU_FACET_SELECTALL:
			geom->SelectAll();
			UpdateFacetParams(TRUE);
			break;
		case MENU_FACET_SELECTSTICK:
			geom->UnselectAll();
			for (int i = 0; i < geom->GetNbFacet(); i++)
				if (geom->GetFacet(i)->sh.sticking != 0.0 && !geom->GetFacet(i)->IsLinkFacet())
					geom->Select(i);
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;
		case MENU_FACET_SELECTTRANS:
			geom->UnselectAll();
			for (int i = 0; i < geom->GetNbFacet(); i++)
				if (geom->GetFacet(i)->sh.opacity != 1.0 && geom->GetFacet(i)->sh.opacity != 2.0)
					geom->Select(i);
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;
		case MENU_FACET_SELECTREFL:
			geom->UnselectAll();
			for (int i = 0; i < geom->GetNbFacet(); i++) {
				Facet *f = geom->GetFacet(i);
				if (f->sh.desorbType == DES_NONE && f->sh.sticking == 0.0 && f->sh.opacity > 0.0)
					geom->Select(i);
			}
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;
		case MENU_FACET_SELECT2SIDE:
			geom->UnselectAll();
			for (int i = 0; i < geom->GetNbFacet(); i++)
				if (geom->GetFacet(i)->sh.is2sided)
					geom->Select(i);
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;
		case MENU_FACET_SELECTVOL:
			geom->UnselectAll();
			for (int i = 0; i < geom->GetNbFacet(); i++)
				if (geom->GetFacet(i)->sh.isVolatile)
					geom->Select(i);
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;
		case MENU_FACET_SELECTTEXT:
			geom->UnselectAll();
			for (int i = 0; i < geom->GetNbFacet(); i++)
				if (geom->GetFacet(i)->sh.isTextured != NULL)
					geom->Select(i);
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;
		case MENU_FACET_SELECTPROF:
			geom->UnselectAll();
			for (int i = 0; i < geom->GetNbFacet(); i++)
				if (geom->GetFacet(i)->sh.isProfile != NULL)
					geom->Select(i);
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;

		case MENU_FACET_SELECTNONPLANAR:

			sprintf(tmp, "%g", planarityThreshold);
			//sprintf(title,"Pipe L/R = %g",L/R);
			input = GLInputBox::GetInput(tmp, "Planarity larger than:", "Select non planar facets");
			if (!input) return;
			if (!sscanf(input, "%lf", &planarityThreshold)) {
				GLMessageBox::Display("Invalid number", "Error", GLDLG_OK, GLDLG_ICONERROR);
				return;
			}
			geom->UnselectAll();
			for (int i = 0; i < geom->GetNbFacet(); i++)
				if (geom->GetFacet(i)->err >= planarityThreshold)
					geom->Select(i);
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;


		case MENU_FACET_SELECTERR:
			geom->UnselectAll();
			for (int i = 0; i < geom->GetNbFacet(); i++)

				if (geom->GetFacet(i)->sh.sign == 0.0)
					geom->Select(i);
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;

		case MENU_FACET_SELECTDEST:
			geom->UnselectAll();
			for (int i = 0; i < geom->GetNbFacet(); i++)

				if (geom->GetFacet(i)->sh.superDest != 0)
					geom->Select(i);
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;

		case MENU_FACET_SELECTTELEPORT:
			geom->UnselectAll();
			for (int i = 0; i < geom->GetNbFacet(); i++)

				if (geom->GetFacet(i)->sh.teleportDest != 0)
					geom->Select(i);
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;

		case MENU_FACET_SELECTABS:
			geom->UnselectAll();
			for (int i = 0; i < geom->GetNbFacet(); i++)
				if (geom->GetFacet(i)->counterCache.hit.nbAbsorbed > 0)
					geom->Select(i);
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;

		case MENU_FACET_SELECTHITS:
			geom->UnselectAll();
			for (int i = 0; i < geom->GetNbFacet(); i++)

				if (geom->GetFacet(i)->counterCache.hit.nbHit > 0)
					geom->Select(i);
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;

		case MENU_FACET_SELECTNOHITS_AREA:

			sprintf(tmp, "%g", largeArea);
			//sprintf(title,"Pipe L/R = %g",L/R);
			input = GLInputBox::GetInput(tmp, "Min.area (cm)", "Select large facets without hits");
			if (!input) return;
			if ((sscanf(input, "%lf", &largeArea) <= 0) || (largeArea <= 0.0)) {
				GLMessageBox::Display("Invalid number", "Error", GLDLG_OK, GLDLG_ICONERROR);
				return;
			}
			geom->UnselectAll();
			for (int i = 0; i < geom->GetNbFacet(); i++)
				if (geom->GetFacet(i)->counterCache.hit.nbHit == 0 && geom->GetFacet(i)->sh.area >= largeArea)
					geom->Select(i);
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;

		case MENU_FACET_SELECTDES:
			geom->UnselectAll();
			for (int i = 0; i < geom->GetNbFacet(); i++)
				if (geom->GetFacet(i)->sh.desorbType != DES_NONE)
					geom->Select(i);
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;
		case MENU_SELECT_HASDESFILE:
			geom->UnselectAll();
			for (int i = 0; i < geom->GetNbFacet(); i++)
				if (geom->GetFacet(i)->hasOutgassingFile)
					geom->Select(i);
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;
		case MENU_FACET_INVERTSEL:
			for (int i = 0; i < geom->GetNbFacet(); i++)
				geom->GetFacet(i)->selected = !geom->GetFacet(i)->selected;
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;
		case MENU_SELECTION_SELECTFACETNUMBER:
			if (!selectDialog) selectDialog = new SelectDialog(&worker);
			selectDialog->SetVisible(TRUE);
			break;
		case MENU_FACET_SAVESEL:
			SaveSelection();
			break;
		case MENU_FACET_LOADSEL:
			LoadSelection();
			break;
		case MENU_SELECTION_ADDNEW:
			AddSelection();
			break;
			break;
		case  MENU_SELECTION_CLEARALL:
			if (GLMessageBox::Display("Clear all selections ?", "Question", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONINFO) == GLDLG_OK) {
				ClearAllSelections();
			}
			break;
		case MENU_VERTEX_UNSELECTALL:
			geom->UnselectAllVertex();
			break;
		case MENU_VERTEX_SELECTALL:
			geom->SelectAllVertex();
			break;
		case MENU_SELECTION_ISOLATED_VERTEX:
			geom->SelectIsolatedVertices();
			break;
		case MENU_VERTEX_CLEAR_ISOLATED:
			geom->DeleteIsolatedVertices(FALSE);
			UpdateModelParams();
			break;
		case MENU_VERTEX_CREATE_POLY_CONVEX:
			if (AskToReset()) {
				try {
					geom->CreatePolyFromVertices_Convex();
				}
				catch (Error &e) {
					GLMessageBox::Display((char *)e.GetMsg(), "Error creating polygon", GLDLG_OK, GLDLG_ICONERROR);
				}
				//UpdateModelParams();
				try {
					worker.Reload();
				}
				catch (Error &e) {
					GLMessageBox::Display((char *)e.GetMsg(), "Error reloading worker", GLDLG_OK, GLDLG_ICONERROR);
				}
			}
			break;
		case MENU_VERTEX_CREATE_POLY_ORDER:
			if (AskToReset()) {
				try {
					geom->CreatePolyFromVertices_Order();
				}
				catch (Error &e) {
					GLMessageBox::Display((char *)e.GetMsg(), "Error creating polygon", GLDLG_OK, GLDLG_ICONERROR);
				}
				//UpdateModelParams();
				try {
					worker.Reload();
				}
				catch (Error &e) {

					GLMessageBox::Display((char *)e.GetMsg(), "Error reloading worker", GLDLG_OK, GLDLG_ICONERROR);
				}
			}
			break;
		case MENU_FACET_CREATE_DIFFERENCE:
			CreateOfTwoFacets(ClipperLib::ctDifference);
			break;
		case MENU_FACET_CREATE_DIFFERENCE2:
			CreateOfTwoFacets(ClipperLib::ctDifference,TRUE);
			break;
		case MENU_FACET_CREATE_UNION:
			CreateOfTwoFacets(ClipperLib::ctUnion);
			break;
		case MENU_FACET_CREATE_INTERSECTION:
			CreateOfTwoFacets(ClipperLib::ctIntersection);
			break;
		case MENU_FACET_CREATE_XOR:
			CreateOfTwoFacets(ClipperLib::ctXor);
			break;
		case MENU_FACET_LOFT:
			if (geom->GetNbSelected() != 2) {
				GLMessageBox::Display("Select exactly 2 facets", "Can't create loft", GLDLG_OK, GLDLG_ICONERROR);
				return;
			}
			if (AskToReset()) {
				geom->CreateLoft();
			}
			worker.Reload();
			mApp->UpdateModelParams();
			mApp->UpdateFacetlistSelected();
			mApp->UpdateViewers();
			break;
		case MENU_FACET_INTERSECT:
			if (!buildIntersection || !buildIntersection->IsVisible()) {
				SAFE_DELETE(buildIntersection);
				buildIntersection = new BuildIntersection(geom, &worker);
				buildIntersection->SetVisible(TRUE);
			}
			break;

		case MENU_VERTEX_SELECT_COPLANAR:
			char *input;


			if (geom->IsLoaded()) {
				if (geom->GetNbSelectedVertex() != 3) {
					GLMessageBox::Display("Select exactly 3 vertices", "Can't define plane", GLDLG_OK, GLDLG_ICONERROR);
					return;
				}
				sprintf(tmp, "%g", tolerance);
				//sprintf(title,"Pipe L/R = %g",L/R);
				input = GLInputBox::GetInput(tmp, "Tolerance (cm)", "Select coplanar vertices");
				if (!input) return;
				if ((sscanf(input, "%lf", &tolerance) <= 0) || (tolerance <= 0.0)) {
					GLMessageBox::Display("Invalid number", "Error", GLDLG_OK, GLDLG_ICONERROR);
					return;
				}
				try { viewer[curViewer]->SelectCoplanar(tolerance); }
				catch (Error &e) {
					GLMessageBox::Display((char *)e.GetMsg(), "Error selecting coplanar vertices", GLDLG_OK, GLDLG_ICONERROR);
				}
			}
			else GLMessageBox::Display("No geometry loaded.", "No geometry", GLDLG_OK, GLDLG_ICONERROR);
			break;
		case MENU_VERTEX_MOVE:
			if (geom->IsLoaded()) {
				if (!moveVertex) moveVertex = new MoveVertex(geom, &worker);

				//moveVertex->DoModal();
				moveVertex->SetVisible(TRUE);

				/*
				UpdateModelParams();
				try { worker.Reload(); } catch(Error &e) {
				GLMessageBox::Display((char *)e.GetMsg(),"Error reloading worker",GLDLG_OK,GLDLG_ICONERROR);
				*/

			}
			else GLMessageBox::Display("No geometry loaded.", "No geometry", GLDLG_OK, GLDLG_ICONERROR);
			break;
		case MENU_VERTEX_SCALE:
			if (geom->IsLoaded()) {
				if (!scaleVertex) scaleVertex = new ScaleVertex(geom, &worker);

				scaleVertex->SetVisible(TRUE);

			}
			else GLMessageBox::Display("No geometry loaded.", "No geometry", GLDLG_OK, GLDLG_ICONERROR);
			break;
		case MENU_VERTEX_COORDINATES:

			if (!vertexCoordinates) vertexCoordinates = new VertexCoordinates();
			vertexCoordinates->Display(&worker);
			break;

		case MENU_VERTEX_ADD:
			if (geom->IsLoaded()) {
				if (!addVertex) addVertex = new AddVertex(geom, &worker);
				addVertex->SetVisible(TRUE);
			}
			else GLMessageBox::Display("No geometry loaded.", "No geometry", GLDLG_OK, GLDLG_ICONERROR);
			break;

		case MENU_VERTEX_REMOVE:
			if (geom->IsLoaded()) {
				if (GLMessageBox::Display("Remove Selected vertices?\nNote: It will also affect facets that contain them!", "Question", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONINFO) == GLDLG_OK)  {
					if (AskToReset()) {
						if (worker.running) worker.Stop_Public();
						geom->RemoveSelectedVertex();
						geom->CalcTotalOutGassing();
						geom->Rebuild(); //Will recalculate facet parameters
						UpdateModelParams();
						if (vertexCoordinates) vertexCoordinates->Update();
						if (facetCoordinates) facetCoordinates->UpdateFromSelection();
						if (profilePlotter) profilePlotter->Refresh();

						if (pressureEvolution) pressureEvolution->Refresh();
						if (timewisePlotter) timewisePlotter->Refresh();
						// Send to sub process
						try { worker.Reload(); }
						catch (Error &e) {
							GLMessageBox::Display((char *)e.GetMsg(), "Error reloading worker", GLDLG_OK, GLDLG_ICONERROR);
						}
					}
				}

			}
			else GLMessageBox::Display("No geometry loaded.", "No geometry", GLDLG_OK, GLDLG_ICONERROR);
			break;

		case MENU_VIEW_FULLSCREEN:
			if (m_bWindowed) {
				ToggleFullscreen();
				PlaceComponents();
			}
			else {
				Resize(1024, 768, TRUE);
			}
			menu->GetSubMenu("View")->SetState(MENU_VIEW_FULLSCREEN, !m_bWindowed);
			break;

		case MENU_VIEW_ADDNEW:
			AddView();
			break;
			break;
		case  MENU_VIEW_CLEARALL:
			if (GLMessageBox::Display("Clear all views ?", "Question", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONINFO) == GLDLG_OK) {
				ClearAllViews();
			}
			break;

		case MENU_TIME_SETTINGS:
			if (!timeSettings) timeSettings = new TimeSettings(&worker);
			timeSettings->SetVisible(TRUE);
			break;
		case MENU_TIME_MOMENTS_EDITOR:
			if (momentsEditor == NULL || !momentsEditor->IsVisible()) {
				SAFE_DELETE(momentsEditor);
				momentsEditor = new MomentsEditor(&worker);
				momentsEditor->Refresh();
				momentsEditor->SetVisible(TRUE);
			}
			break;
		case MENU_TIME_PARAMETER_EDITOR:
			if (parameterEditor == NULL) parameterEditor = new ParameterEditor(&worker);
			parameterEditor->UpdateCombo();
			parameterEditor->SetVisible(TRUE);
			break;
		case MENU_TIMEWISE_PLOTTER:
			if (!timewisePlotter) timewisePlotter = new TimewisePlotter();
			timewisePlotter->Display(&worker);
			break;
		case MENU_TEST_PIPE0001:
			if (AskToSave()) BuildPipe(0.0001);
			break;
		case MENU_TEST_PIPE1:
			if (AskToSave()) BuildPipe(1.0);
			break;
		case MENU_TEST_PIPE10:
			if (AskToSave()) BuildPipe(10.0);
			break;
		case MENU_TEST_PIPE100:
			if (AskToSave()) BuildPipe(100.0);
			break;
		case MENU_TEST_PIPE1000:
			if (AskToSave()) BuildPipe(1000.0);
			break;
		case MENU_TEST_PIPE10000:
			if (AskToSave()) BuildPipe(10000.0);
			break;
		case MENU_QUICKPIPE:
			if (AskToSave()) QuickPipe();
			break;
		}
		// Load recent menu
		if (src->GetId() >= MENU_FILE_LOADRECENT && src->GetId() < MENU_FILE_LOADRECENT + nbRecent) {
			if (AskToSave()) {
				if (worker.running) worker.Stop_Public();
				LoadFile(recents[src->GetId() - MENU_FILE_LOADRECENT]);
			}
		}

		// Show structure menu
		if (src->GetId() >= MENU_VIEW_STRUCTURE && src->GetId() <= MENU_VIEW_STRUCTURE + geom->GetNbStructure()) {
			geom->viewStruct = src->GetId() - MENU_VIEW_STRUCTURE - 1;
			if (src->GetId() > MENU_VIEW_STRUCTURE) geom->UnselectAll();
			UpdateStructMenu();
		}
		if (src->GetId() == MENU_VIEW_NEWSTRUCT) {
			AddStruct();
			UpdateStructMenu();
		}
		if (src->GetId() == MENU_VIEW_DELSTRUCT) {
			DeleteStruct();
			UpdateStructMenu();
		}
		if (src->GetId() == MENU_VIEW_PREVSTRUCT) {
			geom->viewStruct = Remainder(geom->viewStruct - 1, geom->GetNbStructure());
			geom->UnselectAll();
			UpdateStructMenu();
		}
		if (src->GetId() == MENU_VIEW_NEXTSTRUCT) {
			geom->viewStruct = Remainder(geom->viewStruct + 1, geom->GetNbStructure());
			geom->UnselectAll();
			UpdateStructMenu();
		}

		// Select selection
		if (MENU_SELECTION_SELECTIONS + nbSelection > src->GetId() && src->GetId() >= MENU_SELECTION_SELECTIONS) { //Choose selection by number
			SelectSelection(src->GetId() - MENU_SELECTION_SELECTIONS);
		}
		else if (src->GetId() == (MENU_SELECTION_SELECTIONS + nbSelection)){ //Previous selection
			SelectSelection(Remainder(idSelection - 1, nbSelection));
		}
		else if (src->GetId() == (MENU_SELECTION_SELECTIONS + nbSelection + 1)){ //Next selection
			SelectSelection(Remainder(idSelection + 1, nbSelection));
		}

		// Clear selection
		if (src->GetId() >= MENU_SELECTION_CLEARSELECTIONS && src->GetId() < MENU_SELECTION_CLEARSELECTIONS + nbSelection) {
			char tmpname[256];
			sprintf(tmpname, "Clear %s?", selections[src->GetId() - MENU_SELECTION_CLEARSELECTIONS].name);
			if (GLMessageBox::Display(tmpname, "Confirmation", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONINFO) == GLDLG_OK) {
				ClearSelection(src->GetId() - MENU_SELECTION_CLEARSELECTIONS);
			}
		}

		// Memorize selection
		if (src->GetId() >= MENU_SELECTION_MEMORIZESELECTIONS && src->GetId() < MENU_SELECTION_MEMORIZESELECTIONS + nbSelection) {
			OverWriteSelection(src->GetId() - MENU_SELECTION_MEMORIZESELECTIONS);
		}

		// Select view
		if (src->GetId() >= MENU_VIEW_VIEWS && src->GetId() < MENU_VIEW_VIEWS + nbView) {
			SelectView(src->GetId() - MENU_VIEW_VIEWS);
		}
		// Clear view
		if (src->GetId() >= MENU_VIEW_CLEARVIEWS && src->GetId() < MENU_VIEW_CLEARVIEWS + nbView) {
			char tmpname[256];
			sprintf(tmpname, "Clear %s?", views[src->GetId() - MENU_VIEW_CLEARVIEWS].name);
			if (GLMessageBox::Display(tmpname, "Confirmation", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONINFO) == GLDLG_OK) {
				ClearView(src->GetId() - MENU_VIEW_CLEARVIEWS);
			}
		}
		// Memorize view
		if (src->GetId() >= MENU_VIEW_MEMORIZEVIEWS && src->GetId() < MENU_VIEW_MEMORIZEVIEWS + nbView) {
			OverWriteView(src->GetId() - MENU_VIEW_MEMORIZEVIEWS);
		}
		break;

		//TEXT --------------------------------------------------------------------
	case MSG_TEXT_UPD:
		if (src == facetSticking) {
			calcFlow();
			facetApplyBtn->SetEnabled(TRUE);
		}
		else if (src == facetOpacity) {
			facetApplyBtn->SetEnabled(TRUE);
		}
		else if (src == facetDesTypeN) {
			facetApplyBtn->SetEnabled(TRUE);
		}
		else if (src == facetTemperature) {
			calcFlow();
			facetApplyBtn->SetEnabled(TRUE);
		}
		else if (src == facetFlow) {
			double flow;
			double area;
			facetFlow->GetNumber(&flow);
			facetArea->GetNumber(&area);
			if (area == 0) facetFlowArea->SetText("#DIV0");
			else facetFlowArea->SetText(flow / area);
			facetApplyBtn->SetEnabled(TRUE);
			facetFILabel->SetState(TRUE);
			facetFIAreaLabel->SetState(FALSE);
			// //else if ( src == facetMass ) {
			//  facetApplyBtn->SetEnabled(TRUE);
		}
		else if (src == facetFlowArea) {
			double flowPerArea;
			double area;
			facetFlowArea->GetNumber(&flowPerArea);
			facetArea->GetNumber(&area);
			facetFlow->SetText(flowPerArea*area);
			facetApplyBtn->SetEnabled(TRUE);
			facetFIAreaLabel->SetState(TRUE);
			facetFILabel->SetState(FALSE);
			// //else if ( src == facetMass ) {
			//  facetApplyBtn->SetEnabled(TRUE);
		}
		else if (src == facetPumping) {
			calcSticking();
			facetApplyBtn->SetEnabled(TRUE);
			// //else if ( src == facetMass ) {
			//  facetApplyBtn->SetEnabled(TRUE);
		}
		break;

	case MSG_TEXT:
		if (src == facetSticking) {
			ApplyFacetParams();
		}
		else if (src == facetDesTypeN) {
			ApplyFacetParams();
		}
		else if (src == facetOpacity) {
			ApplyFacetParams();
		}
		else if (src == facetTemperature) {
			ApplyFacetParams();
		}
		else if (src == facetPumping) {
			ApplyFacetParams();
		}
		else if (src == facetFlow) {
			ApplyFacetParams();
		}
		else if (src == facetFlowArea) {
			ApplyFacetParams();
		}
		break;

		//COMBO -------------------------------------------------------------------
	case MSG_COMBO:

		if (src == facetDesType) {
			facetApplyBtn->SetEnabled(TRUE);
			BOOL hasDesorption = !facetDesType->GetSelectedIndex() == 0;

			facetFlow->SetEditable(hasDesorption);
			facetFlowArea->SetEditable(hasDesorption);

			int color = (hasDesorption) ? 0 : 110;
			facetFILabel->SetTextColor(color, color, color);
			facetFIAreaLabel->SetTextColor(color, color, color);
			facetFILabel->SetEnabled(hasDesorption);
			facetFIAreaLabel->SetEnabled(hasDesorption);
			facetDesTypeN->SetEditable(facetDesType->GetSelectedIndex() == 3);
		}
		else if (src == facetRecType) {
			facetApplyBtn->SetEnabled(TRUE);
		}
		else if (src == facetSideType) {
			facetApplyBtn->SetEnabled(TRUE);
		}
		else if (src == modeCombo) {

			compACBtn->SetEnabled(modeCombo->GetSelectedIndex() == 1);

			singleACBtn->SetEnabled(modeCombo->GetSelectedIndex() == 1);
			UpdateFacetHits();
		}
		break;

		//TOGGLE ------------------------------------------------------------------
	case MSG_TOGGLE:
		if (src == facetFILabel) {
			facetFILabel->SetState(TRUE);
			facetFIAreaLabel->SetState(FALSE);
		}
		else if (src == facetFIAreaLabel) {
			facetFILabel->SetState(FALSE);
			facetFIAreaLabel->SetState(TRUE);
		}
		else if (src == autoFrameMoveToggle) {
			autoFrameMove = autoFrameMoveToggle->GetState();
			forceFrameMoveButton->SetEnabled(!autoFrameMove);
		}

		else {
			// Update viewer flags
			UpdateViewerFlags();
		}
		break;

		//LIST --------------------------------------------------------------------
	case MSG_LIST:
		if (src == facetList && geom->IsLoaded()) {
			int *sels = (int *)malloc((geom->GetNbFacet())*sizeof(int));
			int nbSel;
			facetList->GetSelectedRows(&sels, &nbSel, TRUE);
			geom->UnselectAll();
			for (int i = 0; i < nbSel; i++)
				geom->Select(sels[i]);
			geom->UpdateSelection();
			UpdateFacetParams();
			SAFE_FREE(sels);
		}
		break;

		//GEOMVIEWER ------------------------------------------------------------------
	case MSG_GEOMVIEWER_MAXIMISE:
	{
		if (src == viewer[0]) {
			AnimateViewerChange(0);
		}
		else if (src == viewer[1]) {
			AnimateViewerChange(1);
		}
		else if (src == viewer[2]) {
			AnimateViewerChange(2);
		}
		else if (src == viewer[3]) {
			AnimateViewerChange(3);
		}
		Place3DViewer();

		BOOL neededTexture = needsTexture;

		BOOL neededMesh = needsMesh;
		CheckNeedsTexture();

		if (!needsTexture && neededTexture) { //We just disabled textures
			worker.GetGeometry()->ClearFacetTextures();
		}
		else if (needsTexture && !neededTexture) { //We just enabled textures
			worker.RebuildTextures();
		}

		if (!needsMesh && neededMesh) { //We just disabled mesh
			geom->ClearFacetMeshLists();
		}
		else if (needsMesh && !neededMesh) { //We just enabled mesh
			geom->BuildFacetMeshLists();
		}

		break;
	}
	case MSG_GEOMVIEWER_SELECT: {
		SelectViewer(src->GetId());
	}break;

		//BUTTON ------------------------------------------------------------------
	case MSG_BUTTON:
		if (src == startSimu) {
			changedSinceSave = TRUE;
			StartStopSimulation();
			resetSimu->SetEnabled(!worker.running);
		}
		else if (src == resetSimu) {
			changedSinceSave = TRUE;
			ResetSimulation();
		}
		else if (src == facetApplyBtn) {
			//changedSinceSave=TRUE;
			ApplyFacetParams();
		}
		else if (src == facetDetailsBtn) {
			if (facetDetails == NULL) facetDetails = new FacetDetails();
			facetDetails->Display(&worker);
		}
		else if (src == facetCoordBtn) {
			if (!facetCoordinates) facetCoordinates = new FacetCoordinates();
			facetCoordinates->Display(&worker);
		}
		else if (src == facetMoreBtn) {
			if (!facetMesh) {
				facetMesh = new FacetMesh(&worker);
				int *selection;
				int nbSel;
				worker.GetGeometry()->GetSelection(&selection, &nbSel);
				facetMesh->Refresh(nbSel, selection);
				SAFE_FREE(selection);
			}
			facetMesh->SetVisible(!facetMesh->IsVisible());
			facetMesh->Reposition();
		}
		else if (src == showMoreBtn) {
			if (!viewer3DSettings)	viewer3DSettings = new Viewer3DSettings();
			viewer3DSettings->SetVisible(!viewer3DSettings->IsVisible());
			viewer3DSettings->Reposition();
			viewer3DSettings->Refresh(geom, viewer[curViewer]);

		}
		else if (src == compACBtn) {
			try { lastUpdate = 0.0; worker.ComputeAC(m_fTime); }
			catch (Error &e) {
				GLMessageBox::Display((char *)e.GetMsg(), "Error", GLDLG_OK, GLDLG_ICONERROR);
				return;
			}
			break;
		}
		else if (src == singleACBtn) {
			try { lastUpdate = 0.0; worker.StepAC(m_fTime); }
			catch (Error &e) {
				GLMessageBox::Display((char *)e.GetMsg(), "Error", GLDLG_OK, GLDLG_ICONERROR);
				return;
			}
			break;
		}
		else if (src == textureScalingBtn) {
			if (!textureSettings || !textureSettings->IsVisible()) {
				SAFE_DELETE(textureSettings);
				textureSettings = new TextureSettings();
				textureSettings->Display(&worker, viewer);
			}
			else {
				textureSettings->SetVisible(FALSE);
			}
		}
		else if (src == profilePlotterBtn) {
			if (!profilePlotter) profilePlotter = new ProfilePlotter();
			if (!profilePlotter->IsVisible()) profilePlotter->Display(&worker);
			else profilePlotter->SetVisible(FALSE);
		}
		else if (src == texturePlotterBtn) {
			if (!texturePlotter) texturePlotter = new TexturePlotter();
			if (!texturePlotter->IsVisible()) texturePlotter->Display(&worker);
			else {
				texturePlotter->SetVisible(FALSE);
				SAFE_DELETE(texturePlotter);
			}
		}
		else if (src == globalSettingsBtn) {
			if (!globalSettings) globalSettings = new GlobalSettings(&worker);
			if (!globalSettings->IsVisible()) {
				globalSettings->Update();
				globalSettings->SetVisible(TRUE);
			}
			else globalSettings->SetVisible(FALSE);
		}
		else if (src == forceFrameMoveButton) {
			frameMoveRequested = TRUE;
			FrameMove();
		}
		else {
			ProcessFormulaButtons(src);
		}
		break;

		//Panel open/close ---------------------------------------------------------
	case MSG_PANELR:
		PlaceComponents();
		break;

	}

}

// ---------------------------------------------------------------------------
BOOL MolFlow::AskToReset(Worker *work) {
	if (work == NULL) work = &worker;
	if (work->nbHit > 0) {
		int rep = GLMessageBox::Display("This will reset simulation data.", "Geometry change", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONWARNING);
		if (rep == GLDLG_OK) {
			work->ResetStatsAndHits(m_fTime);
			nbDesStart = 0;
			nbHitStart = 0;

			//resetSimu->SetEnabled(FALSE);
			if (mApp->profilePlotter) mApp->profilePlotter->Update(m_fTime, TRUE);

			if (mApp->pressureEvolution) mApp->pressureEvolution->Update(m_fTime, TRUE);
			if (mApp->timewisePlotter) mApp->timewisePlotter->Update(m_fTime, TRUE);
			if (mApp->texturePlotter) mApp->texturePlotter->Update(m_fTime, TRUE);
			return TRUE;
		}
		else return FALSE;
	}
	else return TRUE;
}


//-----------------------------------------------------------------------------

void MolFlow::QuickPipe() {
	BuildPipe(5.0, 5);
}

void MolFlow::BuildPipe(double ratio, int steps) {

	char tmp[256];
	Geometry *geom = worker.GetGeometry();

	double R = 1.0;
	double L = ratio * R;
	int    step;
	if (steps) step = 5; //Quick Pipe
	else {
		sprintf(tmp, "100");
		//sprintf(title,"Pipe L/R = %g",L/R);
		char *nbF = GLInputBox::GetInput(tmp, "Number of facet", "Build Pipe");
		if (!nbF) return;
		if ((sscanf(nbF, "%d", &step) <= 0) || (step < 3)) {
			GLMessageBox::Display("Invalid number", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}
	}
	ResetSimulation(FALSE);

	try{
		geom->BuildPipe(L, R, 0, step);
		worker.CalcTotalOutgassing();
		//default values
		worker.enableDecay = FALSE;
		worker.halfLife = 1;
		worker.gasMass = 28;
		worker.ResetMoments();
	}
	catch (Error &e) {
		GLMessageBox::Display((char *)e.GetMsg(), "Error building pipe", GLDLG_OK, GLDLG_ICONERROR);
		geom->Clear();
		return;
	}
	//worker.nbDesorption = 0; //Already done by ResetWorkerStats
	//sprintf(tmp,"L|R %g",L/R);
	worker.SetFileName("");
	nbDesStart = 0;
	nbHitStart = 0;
	for (int i = 0; i < MAX_VIEWER; i++)
		viewer[i]->SetWorker(&worker);
	//UpdateModelParams();
	startSimu->SetEnabled(TRUE);
	compACBtn->SetEnabled(modeCombo->GetSelectedIndex() == 1);
	//resetSimu->SetEnabled(TRUE);
	ClearFacetParams();
	ClearFormula();
	ClearParameters();
	ClearAllSelections();
	ClearAllViews();
	worker.displayedMoment = 0;

	GLParser *f = new GLParser();
	f->SetExpression("A2/SUMDES");
	f->SetName("Trans. Prob.");
	f->Parse();
	AddFormula(f);

	UpdateStructMenu();
	// Send to sub process
	try { worker.Reload(); }
	catch (Error &e) {
		GLMessageBox::Display((char *)e.GetMsg(), "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}

	if (timeSettings) mApp->timeSettings->RefreshMoments();
	if (momentsEditor) mApp->momentsEditor->Refresh();
	if (parameterEditor) mApp->parameterEditor->UpdateCombo();
	if (timewisePlotter) mApp->timewisePlotter->refreshViews();
	if (profilePlotter) profilePlotter->Refresh();
	if (pressureEvolution) pressureEvolution->Refresh();
	if (textureSettings) textureSettings->Update();
	if (texturePlotter) texturePlotter->Update(m_fTime, TRUE);
	if (outgassingMap) outgassingMap->Update(m_fTime, TRUE);
	if (facetDetails) facetDetails->Update();
	if (facetCoordinates) facetCoordinates->UpdateFromSelection();
	if (vertexCoordinates) vertexCoordinates->Update();
	if (movement) movement->Update();
	if (globalSettings && globalSettings->IsVisible()) globalSettings->Update();
	UpdateTitle();
	changedSinceSave = FALSE;
	ResetAutoSaveTimer();
}


void MolFlow::LoadConfig() {

	FileReader *f = NULL;
	char *w;
	nbRecent = 0;

	try {

		f = new FileReader("molflow.cfg");
		Geometry *geom = worker.GetGeometry();

		f->ReadKeyword("showRules"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showRule = f->ReadInt();
		f->ReadKeyword("showNormals"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showNormal = f->ReadInt();
		f->ReadKeyword("showUV"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showUV = f->ReadInt();
		f->ReadKeyword("showLines"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showLine = f->ReadInt();
		f->ReadKeyword("showLeaks"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showLeak = f->ReadInt();
		f->ReadKeyword("showHits"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showHit = f->ReadInt();
		f->ReadKeyword("showVolume"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showVolume = f->ReadInt();
		f->ReadKeyword("showTexture"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showTexture = f->ReadInt();
		f->ReadKeyword("showFilter"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showFilter = f->ReadInt();
		f->ReadKeyword("showIndices"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showIndex = f->ReadInt();
		f->ReadKeyword("showVertices"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showVertex = f->ReadInt();
		f->ReadKeyword("showMode"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showBack = f->ReadInt();
		f->ReadKeyword("showMesh"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showMesh = f->ReadInt();
		f->ReadKeyword("showHidden"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showHidden = f->ReadInt();
		f->ReadKeyword("showHiddenVertex"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showHiddenVertex = f->ReadInt();
		f->ReadKeyword("showTimeOverlay"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showTime = f->ReadInt();
		f->ReadKeyword("texColormap"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showColormap = f->ReadInt();
		f->ReadKeyword("translation"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->transStep = f->ReadDouble();
		f->ReadKeyword("dispNumLines"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->dispNumHits = f->ReadInt();
		f->ReadKeyword("dispNumLeaks"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->dispNumLeaks = f->ReadInt();
		f->ReadKeyword("dirShow"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showDir = f->ReadInt();
		f->ReadKeyword("dirNorme"); f->ReadKeyword(":");
		geom->SetNormeRatio((float)f->ReadDouble());
		f->ReadKeyword("dirAutoNormalize"); f->ReadKeyword(":");
		geom->SetAutoNorme(f->ReadInt());
		f->ReadKeyword("dirCenter"); f->ReadKeyword(":");
		geom->SetCenterNorme(f->ReadInt());
		f->ReadKeyword("angle"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->angleStep = f->ReadDouble();
		f->ReadKeyword("autoScale"); f->ReadKeyword(":");
		geom->texAutoScale = f->ReadInt();
		f->ReadKeyword("autoScale_include_constant_flow"); f->ReadKeyword(":");
		geom->texAutoScaleIncludeConstantFlow = f->ReadInt();

		f->ReadKeyword("textures_min_pressure_all"); f->ReadKeyword(":");
		geom->texture_limits[0].autoscale.min.all = f->ReadDouble();
		f->ReadKeyword("textures_min_pressure_moments_only"); f->ReadKeyword(":");
		geom->texture_limits[0].autoscale.min.moments_only = f->ReadDouble();
		f->ReadKeyword("textures_max_pressure_all"); f->ReadKeyword(":");
		geom->texture_limits[0].autoscale.max.all = f->ReadDouble();
		f->ReadKeyword("textures_max_pressure_moments_only"); f->ReadKeyword(":");
		geom->texture_limits[0].autoscale.max.moments_only = f->ReadDouble();

		f->ReadKeyword("textures_min_impingement_all"); f->ReadKeyword(":");
		geom->texture_limits[1].autoscale.min.all = f->ReadDouble();
		f->ReadKeyword("textures_min_impingement_moments_only"); f->ReadKeyword(":");
		geom->texture_limits[1].autoscale.min.moments_only = f->ReadDouble();
		f->ReadKeyword("textures_max_impingement_all"); f->ReadKeyword(":");
		geom->texture_limits[1].autoscale.max.all = f->ReadDouble();
		f->ReadKeyword("textures_max_impingement_moments_only"); f->ReadKeyword(":");
		geom->texture_limits[1].autoscale.max.moments_only = f->ReadDouble();

		f->ReadKeyword("textures_min_density_all"); f->ReadKeyword(":");
		geom->texture_limits[2].autoscale.min.all = f->ReadDouble();
		f->ReadKeyword("textures_min_density_moments_only"); f->ReadKeyword(":");
		geom->texture_limits[2].autoscale.min.moments_only = f->ReadDouble();
		f->ReadKeyword("textures_max_density_all"); f->ReadKeyword(":");
		geom->texture_limits[2].autoscale.max.all = f->ReadDouble();
		f->ReadKeyword("textures_max_density_moments_only"); f->ReadKeyword(":");
		geom->texture_limits[2].autoscale.max.moments_only = f->ReadDouble();

		f->ReadKeyword("processNum"); f->ReadKeyword(":");
		nbProc = f->ReadInt();
#ifdef _DEBUG
		nbProc=1;
#endif
		if (nbProc <= 0) nbProc = 1;
		f->ReadKeyword("recents"); f->ReadKeyword(":"); f->ReadKeyword("{");
		w = f->ReadString();
		while (strcmp(w, "}") != 0 && nbRecent < MAX_RECENT)  {
			recents[nbRecent] = _strdup(w);
			nbRecent++;
			w = f->ReadString();
		}
		for (int i = nbRecent - 1; i >= 0; i--)
			menu->GetSubMenu("File")->GetSubMenu("Load recent")->Add(recents[i], MENU_FILE_LOADRECENT + i);

		f->ReadKeyword("cdir"); f->ReadKeyword(":");
		strcpy(currentDir, f->ReadString());
		f->ReadKeyword("cseldir"); f->ReadKeyword(":");
		strcpy(currentSelDir, f->ReadString());
		f->ReadKeyword("autonorme"); f->ReadKeyword(":");
		geom->SetAutoNorme(f->ReadInt());
		f->ReadKeyword("centernorme"); f->ReadKeyword(":");
		geom->SetCenterNorme(f->ReadInt());
		f->ReadKeyword("normeratio"); f->ReadKeyword(":");
		geom->SetNormeRatio((float)(f->ReadDouble()));
		f->ReadKeyword("showDirection"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showDir = f->ReadInt();

		f->ReadKeyword("autoSaveFrequency"); f->ReadKeyword(":");
		autoSaveFrequency = f->ReadDouble();
		f->ReadKeyword("autoSaveSimuOnly"); f->ReadKeyword(":");
		autoSaveSimuOnly = f->ReadInt();
		f->ReadKeyword("checkForUpdates"); f->ReadKeyword(":");
		checkForUpdates = f->ReadInt();
		f->ReadKeyword("autoUpdateFormulas"); f->ReadKeyword(":");
		autoUpdateFormulas = f->ReadInt();
		f->ReadKeyword("compressSavedFiles"); f->ReadKeyword(":");
		compressSavedFiles = f->ReadInt();
		f->ReadKeyword("gasMass"); f->ReadKeyword(":");
		worker.gasMass = f->ReadDouble();
		f->ReadKeyword("expandShortcutPanel"); f->ReadKeyword(":");
		BOOL isOpen = f->ReadInt();
		if (isOpen) shortcutPanel->Open();
		else shortcutPanel->Close();

	}
	catch (Error &err) {
		/*std::ostringstream tmp;
		tmp << err.GetMsg() << "\n\nThis is normal on the first launch and if you upgrade from an earlier version\n";
		tmp << "MolFlow will use default program settings.\nWhen you quit, a correct config file will be written\n";
		GLMessageBox::Display(tmp.str().c_str(), "Error loading config file", GLDLG_OK, GLDLG_ICONINFO);*/
		SAFE_DELETE(f);
		return;
	}

	SAFE_DELETE(f);


}

//-----------------------------------------------------------------------------

#define WRITEI(name,var) {             \
	f->Write(name);                      \
	f->Write(":");                       \
	for(int i=0;i<MAX_VIEWER;i++)        \
	f->WriteInt(viewer[i]->var," ");   \
	f->Write("\n");                      \
}                                      \

#define WRITED(name,var) {             \
	f->Write(name);                      \
	f->Write(":");                       \
	for(int i=0;i<MAX_VIEWER;i++)        \
	f->WriteDouble(viewer[i]->var," ");\
	f->Write("\n");                      \
}



void MolFlow::SaveConfig() {

	FileWriter *f = NULL;

	try {

		f = new FileWriter("molflow.cfg");
		Geometry *geom = worker.GetGeometry();

		// Save flags
		WRITEI("showRules", showRule);
		WRITEI("showNormals", showNormal);
		WRITEI("showUV", showUV);
		WRITEI("showLines", showLine);
		WRITEI("showLeaks", showLeak);
		WRITEI("showHits", showHit);
		WRITEI("showVolume", showVolume);
		WRITEI("showTexture", showTexture);
		WRITEI("showFilter", showFilter);
		WRITEI("showIndices", showIndex);
		WRITEI("showVertices", showVertex);
		WRITEI("showMode", showBack);
		WRITEI("showMesh", showMesh);
		WRITEI("showHidden", showHidden);
		WRITEI("showHiddenVertex", showHiddenVertex);
		WRITEI("showTimeOverlay", showTime);
		WRITEI("texColormap", showColormap);
		WRITED("translation", transStep);
		WRITEI("dispNumLines", dispNumHits);
		WRITEI("dispNumLeaks", dispNumLeaks);
		WRITEI("dirShow", showDir);
		f->Write("dirNorme:"); f->WriteDouble(geom->GetNormeRatio(), "\n");
		f->Write("dirAutoNormalize:"); f->WriteInt(geom->GetAutoNorme(), "\n");
		f->Write("dirCenter:"); f->WriteInt(geom->GetCenterNorme(), "\n");

		WRITED("angle", angleStep);
		f->Write("autoScale:"); f->WriteInt(geom->texAutoScale, "\n");
		f->Write("autoScale_include_constant_flow:"); f->WriteInt(geom->texAutoScaleIncludeConstantFlow, "\n");


		f->Write("textures_min_pressure_all:");
		f->WriteDouble(geom->texture_limits[0].autoscale.min.all, "\n");
		f->Write("textures_min_pressure_moments_only:");
		f->WriteDouble(geom->texture_limits[0].autoscale.min.moments_only, "\n");
		f->Write("textures_max_pressure_all:");
		f->WriteDouble(geom->texture_limits[0].autoscale.max.all, "\n");
		f->Write("textures_max_pressure_moments_only:");
		f->WriteDouble(geom->texture_limits[0].autoscale.max.moments_only, "\n");

		f->Write("textures_min_impingement_all:");
		f->WriteDouble(geom->texture_limits[1].autoscale.min.all, "\n");

		f->Write("textures_min_impingement_moments_only:");
		f->WriteDouble(geom->texture_limits[1].autoscale.min.moments_only, "\n");
		f->Write("textures_max_impingement_all:");
		f->WriteDouble(geom->texture_limits[1].autoscale.max.all, "\n");
		f->Write("textures_max_impingement_moments_only:");
		f->WriteDouble(geom->texture_limits[1].autoscale.max.moments_only, "\n");

		f->Write("textures_min_density_all:");
		f->WriteDouble(geom->texture_limits[2].autoscale.min.all, "\n");
		f->Write("textures_min_density_moments_only:");
		f->WriteDouble(geom->texture_limits[2].autoscale.min.moments_only, "\n");
		f->Write("textures_max_density_all:");
		f->WriteDouble(geom->texture_limits[2].autoscale.max.all, "\n");
		f->Write("textures_max_density_moments_only:");
		f->WriteDouble(geom->texture_limits[2].autoscale.max.moments_only, "\n");

#ifdef _DEBUG
		f->Write("processNum:");f->WriteInt(numCPU,"\n");
#else
		f->Write("processNum:"); f->WriteInt(worker.GetProcNumber(), "\n");
#endif
		f->Write("recents:{\n");
		for (int i = 0; i < nbRecent; i++) {
			f->Write("\"");
			f->Write(recents[i]);
			f->Write("\"\n");
		}
		f->Write("}\n");

		f->Write("cdir:\""); f->Write(currentDir); f->Write("\"\n");
		f->Write("cseldir:\""); f->Write(currentSelDir); f->Write("\"\n");
		f->Write("autonorme:"); f->WriteInt(geom->GetAutoNorme(), "\n");
		f->Write("centernorme:"); f->WriteInt(geom->GetCenterNorme(), "\n");
		f->Write("normeratio:"); f->WriteDouble((double)(geom->GetNormeRatio()), "\n");
		WRITEI("showDirection", showDir); f->Write("\n");

		f->Write("autoSaveFrequency:"); f->WriteDouble(autoSaveFrequency, "\n");
		f->Write("autoSaveSimuOnly:"); f->WriteInt(autoSaveSimuOnly, "\n");
		f->Write("checkForUpdates:"); f->WriteInt(checkForUpdates, "\n");
		f->Write("autoUpdateFormulas:"); f->WriteInt(autoUpdateFormulas, "\n");
		f->Write("compressSavedFiles:"); f->WriteInt(compressSavedFiles, "\n");
		f->Write("gasMass:"); f->WriteDouble(worker.gasMass, "\n");
		f->Write("expandShortcutPanel:"); f->WriteInt(!shortcutPanel->IsClosed(), "\n");

	}
	catch (Error &err) {
		GLMessageBox::Display(err.GetMsg(), "Error saving config file", GLDLG_OK, GLDLG_ICONWARNING);
	}

	SAFE_DELETE(f);

}

void MolFlow::calcFlow() {
	double sticking;
	double area;
	double flow;
	double temperature;
	//double mass;

	facetSticking->GetNumber(&sticking);
	facetArea->GetNumber(&area);
	facetTemperature->GetNumber(&temperature);
	//facetMass->GetNumber(&mass);

	flow = 1 * sticking*area / 10.0 / 4.0*sqrt(8.0*8.31*temperature / PI / (worker.gasMass*0.001));
	facetPumping->SetText(flow);
	return;
}

void MolFlow::calcSticking() {
	double sticking;
	double area;
	double flow;
	double temperature;
	//double mass;

	facetPumping->GetNumber(&flow);
	facetArea->GetNumber(&area);
	facetTemperature->GetNumber(&temperature);
	//facetMass->GetNumber(&mass);

	sticking = abs(flow / (area / 10.0)*4.0*sqrt(1.0 / 8.0 / 8.31 / (temperature)*PI*(worker.gasMass*0.001)));
	//if (sticking<=1.0) {
	facetSticking->SetText(sticking);
	//}
	//else { //Sticking: max. 1.0
	//	SetParam(facetSticking,1.0);
	//	calcFlow();
	//}
	return;
}

void MolFlow::CrashHandler(Error *e) {
	char tmp[1024];
	sprintf(tmp, "Well, that's emberassing. Molflow crashed and will exit now.\nBefore that, an autosave will be attempted.\nHere is the error info:\n\n%s", (char *)e->GetMsg());
	GLMessageBox::Display(tmp, "Main crash handler", GLDLG_OK, GLDGL_ICONDEAD);
	try {
		if (AutoSave(TRUE))
			GLMessageBox::Display("Good news, autosave worked!", "Main crash handler", GLDLG_OK, GLDGL_ICONDEAD);
		else
			GLMessageBox::Display("Sorry, I couldn't even autosave.", "Main crash handler", GLDLG_OK, GLDGL_ICONDEAD);
	}
	catch (Error &e) {
		e.GetMsg();
		GLMessageBox::Display("Sorry, I couldn't even autosave.", "Main crash handler", GLDLG_OK, GLDGL_ICONDEAD);
	}
}


BOOL MolFlow::EvaluateVariable(VLIST *v,Worker *w, Geometry *geom) {
	BOOL ok=TRUE;
	int nbFacet = geom->GetNbFacet();
	int idx;

	if ((idx = GetVariable(v->name, "A")) > 0) {
		ok = (idx>0 && idx <= nbFacet);
		if (ok) v->value = (double)geom->GetFacet(idx - 1)->counterCache.hit.nbAbsorbed;
	}
	else if ((idx = GetVariable(v->name, "D")) > 0) {
		ok = (idx>0 && idx <= nbFacet);
		if (ok) v->value = (double)geom->GetFacet(idx - 1)->counterCache.hit.nbDesorbed;
	}
	else if ((idx = GetVariable(v->name, "H")) > 0) {
		ok = (idx>0 && idx <= nbFacet);
		if (ok) v->value = (double)geom->GetFacet(idx - 1)->counterCache.hit.nbHit;
	}
	else if ((idx = GetVariable(v->name, "P")) > 0) {
		ok = (idx>0 && idx <= nbFacet);
		if (ok) v->value = geom->GetFacet(idx - 1)->counterCache.hit.sum_v_ort *
			worker.GetMoleculesPerTP()*1E4 / (geom->GetFacet(idx - 1)->sh.area*
			(geom->GetFacet(idx - 1)->sh.is2sided ? 2.0 : 1.0)) * (worker.gasMass / 1000 / 6E23)*0.0100;
	}
	else if ((idx = GetVariable(v->name, "DEN")) > 0) {
		ok = (idx>0 && idx <= nbFacet);
		if (ok) {
			double dCoef = 1.0;
			Facet *f = geom->GetFacet(idx - 1);
			if (f->counterCache.hit.nbHit>0 || f->counterCache.hit.nbDesorbed>0)
				if (f->counterCache.hit.nbAbsorbed >0 || f->counterCache.hit.nbDesorbed>0) //otherwise save calculation time
					dCoef *= 1.0 - ((double)f->counterCache.hit.nbAbsorbed + (double)f->counterCache.hit.nbDesorbed) / ((double)f->counterCache.hit.nbHit + (double)f->counterCache.hit.nbDesorbed) / 2.0;
			v->value = dCoef * f->counterCache.hit.sum_1_per_ort_velocity /
				f->GetArea() *
				worker.GetMoleculesPerTP()*1E4;
		}
	}
	else if ((idx = GetVariable(v->name, "Z")) > 0) {
		ok = (idx>0 && idx <= nbFacet);
		if (ok) v->value = geom->GetFacet(idx - 1)->counterCache.hit.nbHit /
			(geom->GetFacet(idx - 1)->sh.area*(geom->GetFacet(idx - 1)->sh.is2sided ? 2.0 : 1.0)) *
			worker.GetMoleculesPerTP()*1E4;
	}
	else if ((idx = GetVariable(v->name, "V")) > 0) {
		ok = (idx>0 && idx <= nbFacet);
		if (ok) v->value = 4.0*(double)(geom->GetFacet(idx - 1)->counterCache.hit.nbHit + geom->GetFacet(idx - 1)->counterCache.hit.nbDesorbed) /
			geom->GetFacet(idx - 1)->counterCache.hit.sum_1_per_ort_velocity;
	}
	else if ((idx = GetVariable(v->name, "T")) > 0) {
		ok = (idx>0 && idx <= nbFacet);
		if (ok) v->value = geom->GetFacet(idx - 1)->sh.temperature;
	}
	else if ((idx = GetVariable(v->name, "_A")) > 0) {
		ok = (idx>0 && idx <= nbFacet);
		if (ok) v->value = (double)geom->GetFacet(idx - 1)->counterCache.density.absorbed;
	}
	else if ((idx = GetVariable(v->name, "_D")) > 0) {
		ok = (idx>0 && idx <= nbFacet);
		if (ok) v->value = (double)geom->GetFacet(idx - 1)->counterCache.density.desorbed;
	}
	else if ((idx = GetVariable(v->name, "_H")) > 0) {
		ok = (idx>0 && idx <= nbFacet);
		if (ok) v->value = (double)geom->GetFacet(idx - 1)->counterCache.density.value;
	}
	else if ((idx = GetVariable(v->name, "AR")) > 0) {
		ok = (idx>0 && idx <= nbFacet);
		if (ok) v->value = geom->GetFacet(idx - 1)->sh.area;
	}
	else if (_stricmp(v->name, "SUMDES") == 0) {
		v->value = (double)worker.nbDesorption;
	}
	else if (_stricmp(v->name, "SUMABS") == 0) {
		v->value = (double)worker.nbAbsorption;
	}
	else if (_stricmp(v->name, "SUMHIT") == 0) {
		v->value = (double)worker.nbHit;
	}
	else if (_stricmp(v->name, "MPP") == 0) {
		v->value = worker.distTraveledTotal_total / (double)worker.nbDesorption;
	}
	else if (_stricmp(v->name, "MFP") == 0) {
		v->value = worker.distTraveledTotal_fullHitsOnly / (double)worker.nbHit;
	}
	else if (_stricmp(v->name, "DESAR") == 0) {
		double sumArea = 0.0;
		for (int i2 = 0; i2 < geom->GetNbFacet(); i2++) {
			Facet *f_tmp = geom->GetFacet(i2);
			if (f_tmp->sh.desorbType) sumArea += f_tmp->sh.area*(f_tmp->sh.is2sided ? 2.0 : 1.0);
		}
		v->value = sumArea;
	}
	else if (_stricmp(v->name, "ABSAR") == 0) {
		double sumArea = 0.0;

		for (int i2 = 0; i2<geom->GetNbFacet(); i2++) {
			Facet *f_tmp = geom->GetFacet(i2);
			if (f_tmp->sh.sticking>0.0) sumArea += f_tmp->sh.area*f_tmp->sh.opacity*(f_tmp->sh.is2sided ? 2.0 : 1.0);
		}
		v->value = sumArea;
	}
	else if (_stricmp(v->name, "QCONST") == 0) {
		v->value = worker.finalOutgassingRate_Pa_m3_sec*10.00; //10: Pa*m3/sec -> mbar*l/s
	}
	else if (_stricmp(v->name, "QCONST_N") == 0) {
		v->value = worker.finalOutgassingRate;
	}
	else if (_stricmp(v->name, "NTOT") == 0) {
		v->value = worker.totalDesorbedMolecules;
	}
	else if (_stricmp(v->name, "GASMASS") == 0) {
		v->value = worker.gasMass;
	}
	else if (_stricmp(v->name, "KB") == 0) {
		v->value = 1.3806504e-23;
	}
	else if (_stricmp(v->name, "R") == 0) {
		v->value = 8.314472;
	}
	else if (_stricmp(v->name, "Na") == 0) {
		v->value = 6.02214179e23;
	}
	else ok = FALSE;
	return ok;
}

void MolFlow::ResetSimulation(BOOL askConfirm){
	Interface::ResetSimulation(askConfirm);
	if (pressureEvolution) pressureEvolution->Update(m_fTime, TRUE);
	if (timewisePlotter) timewisePlotter->Update(m_fTime, TRUE);
}