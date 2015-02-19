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
#include "GLApp/GLSaveDialog.h"
#include "GLApp/GLInputBox.h"
#include "GLApp/GLFileBox.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLWindowManager.h"
#include "Utils.h" //for Remainder
#include "direct.h"
#include <vector>
#include <string>

#ifdef _DEBUG
#define APP_NAME "MolFlow+ development version 64-bit (Compiled "__DATE__" "__TIME__") DEBUG MODE"
#else
//#define APP_NAME "Molflow+ development version ("__DATE__")"
#define APP_NAME "Molflow+ 2.6.4 64-bit ("__DATE__")"
#endif

/*
//Leak detection
#ifdef _DEBUG
#define _CRTDBG_MAP_ALLOC
#define DEBUG_NEW new(_NORMAL_BLOCK, __FILE__, __LINE__)
#define new DEBUG_NEW
#endif
*/


// Some contants
/*
static const char *fileLFilters[] = { "All files" , "*.*" , "Geometry files" , "*.geo" , "Text files" , "*.txt" , "STR files" , "*.str" , "STL files" , "*.stl" , "ASE files" , "*.ase"};
static const int   nbLFilter = sizeof(fileLFilters) / (2*sizeof(char *));
static const char *fileSFilters[] = { "Geometry files" , "*.geo" , "Text files" , "*.txt" };
static const int   nbSFilter = sizeof(fileSFilters) / (2*sizeof(char *));
static const char *fileSelFilters[] = { "Selection files" , "*.sel" };
static const int   nbSelFilter = sizeof(fileSelFilters) / (2*sizeof(char *));
static const char *fileTexFilters[] = { "Text files" , "*.txt" , "Texture files" , "*.tex" };
static const int   nbTexFilter = sizeof(fileTexFilters) / (2*sizeof(char *));
*/

static const char *fileLFilters = "All MolFlow supported files\0*.txt;*.xml;*.zip;*.geo;*.geo7z;*.syn;*.syn7z;*.str;*.stl;*.ase\0XML files\0*.xml;*.zip;\0GEO files\0*.geo;*.geo7z;\0SYN files\0*.syn;*.syn7z;\0TXT files\0*.txt\0STR files\0*.str\0STL files\0*.stl\0ASE files\0*.ase\0";
static const int   nbLFilter = 9;
static const char *fileInsFilters = "\0GEO files\0*.geo;*.geo7z;SYN files\0*.syn;*.syn7z;\0Text files\0*.txt\0STL files\0*.stl\0";
static const int   nbInsFilter = 4;
static const char *fileSFilters = "All files\0*.*\0";
static const int   nbSFilter = 1;
static const char *fileSelFilters = "Selection files\0*.sel\0All files\0*.*\0";
static const int   nbSelFilter = 2;
static const char *fileTexFilters = "Text files\0*.txt\0Texture files\0*.tex\0All files\0*.*\0";
static const int   nbTexFilter = 3;
static const char *fileProfFilters = "CSV file\0*.csv\0Text files\0*.txt\0All files\0*.*\0";
static const int   nbProfFilter = 3;

static const char *fileDesFilters = "Desorption files\0*.des\0All files\0*.*\0";
static const int   nbDesFilter = 2;

static const int   cWidth[] = {30,56,50,50};
static const char *cName[] = {"#","Hits","Des","Abs"};

float m_fTime;
MolFlow *mApp;

#define MENU_FILE_LOAD       11
#define MENU_FILE_IMPORTDES_DES 120
#define MENU_FILE_IMPORTDES_SYN 121
#define MENU_FILE_SAVE       13
#define MENU_FILE_SAVEAS     14
#define MENU_FILE_INSERTGEO  140
#define MENU_FILE_INSERTGEO_NEWSTR  141
#define MENU_FILE_EXPORTMESH      16

#define MENU_FILE_EXPORTTEXTURES 150
#define MENU_FILE_EXPORTTEXTURE_AREA 151
#define MENU_FILE_EXPORTTEXTURE_MCHITS 152
#define MENU_FILE_EXPORTTEXTURE_IMPINGEMENT 153
#define MENU_FILE_EXPORTTEXTURE_PART_DENSITY 154
#define MENU_FILE_EXPORTTEXTURE_GAS_DENSITY 155
#define MENU_FILE_EXPORTTEXTURE_PRESSURE 156
#define MENU_FILE_EXPORTTEXTURE_AVG_V 157
#define MENU_FILE_EXPORTTEXTURE_V_VECTOR 158
#define MENU_FILE_EXPORTTEXTURE_N_VECTORS 159
#define MENU_FILE_EXPORTPROFILES 160

#define MENU_FILE_LOADRECENT 110
#define MENU_FILE_EXIT       17

#define MENU_EDIT_TSCALING     22
#define MENU_EDIT_ADDFORMULA   23
#define MENU_EDIT_UPDATEFORMULAS 24
#define MENU_EDIT_GLOBALSETTINGS 25
#define MENU_EDIT_MOVINGPARTS 26

#define MENU_FACET_COLLAPSE    300
#define MENU_FACET_SWAPNORMAL  301
#define MENU_FACET_SHIFTVERTEX 302
#define MENU_FACET_COORDINATES 303
#define MENU_FACET_PROFPLOTTER 304
#define MENU_FACET_DETAILS     305
//#define MENU_FACET_MESH        306
#define MENU_FACET_TEXPLOTTER  307
#define MENU_FACET_REMOVESEL   308
#define MENU_FACET_EXPLODE     309
#define MENU_FACET_SELECTALL   310
#define MENU_FACET_SELECTSTICK 311
#define MENU_FACET_SELECTDES   312
#define MENU_SELECT_HASDESFILE 313
#define MENU_FACET_SELECTABS   314
#define MENU_FACET_SELECTTRANS 315
#define MENU_FACET_SELECTREFL  316
#define MENU_FACET_SELECT2SIDE 317
#define MENU_FACET_SELECTTEXT  318
#define MENU_FACET_SELECTPROF  319
#define MENU_FACET_SELECTDEST  320

#define MENU_FACET_SELECTVOL   321
#define MENU_FACET_SELECTERR   322
#define MENU_FACET_SELECTNONPLANAR 323
#define MENU_FACET_SELECTHITS        3230
#define MENU_FACET_SELECTNOHITS_AREA 3231
#define MENU_FACET_SAVESEL     324
#define MENU_FACET_LOADSEL     325
#define MENU_FACET_INVERTSEL   326
#define MENU_FACET_MOVE		   327
#define MENU_FACET_SCALE       328
#define MENU_FACET_MIRROR	   329
#define MENU_FACET_ROTATE	   330
#define MENU_FACET_ALIGN       331
#define MENU_FACET_OUTGASSINGMAP 332
#define MENU_FACET_CREATE_DIFFERENCE 3321

#define MENU_SELECTION_ADDNEW             333
#define MENU_SELECTION_CLEARALL           334

#define MENU_SELECTION_MEMORIZESELECTIONS   3300
#define MENU_SELECTION_SELECTIONS           3400
#define MENU_SELECTION_CLEARSELECTIONS      3500

#define MENU_SELECTION_SELECTFACETNUMBER 360

#define MENU_VERTEX_SELECTALL   701
#define MENU_VERTEX_UNSELECTALL 702
#define MENU_VERTEX_SELECT_ISOLATED 703
#define MENU_VERTEX_CREATE_POLY_CONVEX   7040
#define MENU_VERTEX_CREATE_POLY_ORDER    7041
#define MENU_VERTEX_SELECT_COPLANAR   705
#define MENU_VERTEX_MOVE   706
#define MENU_VERTEX_ADD	   707
#define MENU_VERTEX_SCALE  708
#define MENU_VERTEX_REMOVE 709
#define MENU_VERTEX_COORDINATES 710

#define MENU_VIEW_STRUCTURE       4000
#define MENU_VIEW_STRUCTURE_P     40
#define MENU_VIEW_NEWSTRUCT       401
#define MENU_VIEW_DELSTRUCT       402
#define MENU_VIEW_PREVSTRUCT	  403		
#define MENU_VIEW_NEXTSTRUCT	  404
#define MENU_VIEW_FULLSCREEN      41

#define MENU_VIEW_ADDNEW          431
#define MENU_VIEW_CLEARALL        432

#define MENU_VIEW_MEMORIZEVIEWS   4300
#define MENU_VIEW_VIEWS           4400
#define MENU_VIEW_CLEARVIEWS      4500

#define MENU_TIME_SETTINGS          50
#define MENU_TIMEWISE_PLOTTER       51
#define MENU_TIME_PRESSUREEVOLUTION 52
#define MENU_TIME_MOMENTS_EDITOR    53
#define MENU_TIME_PARAMETER_EDITOR  54

#define MENU_TEST_PIPE0001        60
#define MENU_TEST_PIPE1           61
#define MENU_TEST_PIPE10          62
#define MENU_TEST_PIPE100         63
#define MENU_TEST_PIPE1000        64
#define MENU_TEST_PIPE10000       65

#define MENU_QUICKPIPE            66

//-----------------------------------------------------------------------------
// Name: WinMain()
// Desc: Entry point to the program. Initializes everything, and goes into a
//       message-processing loop. Idle time is used to render the scene.
//-----------------------------------------------------------------------------

INT WINAPI WinMain( HINSTANCE hInst, HINSTANCE, LPSTR, INT )
{
	MolFlow *mApp = new MolFlow();
	
	if( !mApp->Create( 1024 , 768 , FALSE ) ) {
		char *logs = GLToolkit::GetLogs();
#ifdef WIN64
		if(logs) MessageBox(NULL,logs,"Molflow [Fatal error]",MB_OK);
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
	} catch(Error &e) {
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

	//Get number of cores
	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);

	numCPU = (int)sysinfo.dwNumberOfProcessors;

	mApp = this; //to refer to the app as extern variable
	antiAliasing = TRUE;
	whiteBg = FALSE;
	checkForUpdates = FALSE;
	autoUpdateFormulas = FALSE;
	compressSavedFiles = TRUE;
	/*double gasMass=28;
	double totalOutgassing=0.0; //total outgassing in Pa*m3/sec (internally everything is in SI units)
	double totalInFlux = 0.0; //total incoming molecules per second. For anisothermal system, it is (totalOutgassing / Kb / T)*/
	autoSaveFrequency = 10.0; //in minutes
	autoSaveSimuOnly = FALSE;
	compressProcessHandle = NULL;

	lastSaveTime = 0.0f;
	lastSaveTimeSimu = 0.0f;
	changedSinceSave = FALSE;
	//lastHeartBeat=0.0f;
	nbDesStart = 0;
	nbHitStart = 0;

	lastUpdate=0.0;
	nbFormula = 0;
	nbRecent = 0;

	nbView = 0;
	nbSelection = 0;
	idView = 0;
	idSelection = 0;

#ifdef _DEBUG
	nbProc = 1;
#else
	nbProc = numCPU;
#endif

	curViewer = 0;
	strcpy(currentDir,".");
	strcpy(currentSelDir,".");
	memset(formulas,0,sizeof formulas);
	formulaSettings = NULL;
	collapseSettings = NULL;
	importDesorption = NULL;
	moveVertex = NULL;
	timeSettings = NULL;
	scaleVertex = NULL;
	scaleFacet = NULL;
	selectDialog = NULL;
	moveFacet = NULL;
	mirrorFacet = NULL;
	rotateFacet = NULL;
	movement = NULL;
	alignFacet = NULL;
	addVertex = NULL;
	facetMesh = NULL;
	facetDetails = NULL;
	viewer3DSettings = NULL;
	textureSettings = NULL;
	globalSettings = NULL;
	facetCoordinates = NULL;
	vertexCoordinates = NULL;
	profilePlotter = NULL;
	pressureEvolution = NULL;
	timewisePlotter = NULL;
	viewEditor = NULL;
	texturePlotter = NULL;
	outgassingMap = NULL;
	momentsEditor = NULL;
	parameterEditor = NULL;
	m_strWindowTitle = APP_NAME;
	wnd->SetBackgroundColor(212,208,200);
	m_bResizable = TRUE;
	m_minScreenWidth  = 800;
	m_minScreenHeight = 600;
	tolerance=1e-8;
	largeArea=1.0;
	planarityThreshold = 1e-5;
	/*totalOutgassing=0.0;
	totalInFlux = 0.0;*/
	//molflowHandle=GetCurrentProcess();
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

	GLToolkit::SetIcon32x32("images/app_icon.png");

	for(int i=0;i<MAX_VIEWER;i++) {
		viewer[i] = new GeometryViewer(i);
		Add(viewer[i]);
	}
	modeSolo = TRUE;
	//nbSt = 0;

	menu = new GLMenuBar(0);
	wnd->SetMenuBar(menu);
	menu->Add("File");
	menu->GetSubMenu("File")->Add("&Load",MENU_FILE_LOAD,SDLK_o,CTRL_MODIFIER);
	menu->GetSubMenu("File")->Add("Import desorption file");
	menu->GetSubMenu("File")->GetSubMenu("Import desorption file")->Add("SYN file", MENU_FILE_IMPORTDES_SYN);
	menu->GetSubMenu("File")->GetSubMenu("Import desorption file")->Add("DES file (deprecated)", MENU_FILE_IMPORTDES_DES);
	menu->GetSubMenu("File")->Add("&Save", MENU_FILE_SAVE, SDLK_s, CTRL_MODIFIER);
	menu->GetSubMenu("File")->Add("&Save as",MENU_FILE_SAVEAS);

	menu->GetSubMenu("File")->Add(NULL); //separator
	menu->GetSubMenu("File")->Add("&Insert geometry");
	menu->GetSubMenu("File")->GetSubMenu("Insert geometry")->Add("&To current structure",MENU_FILE_INSERTGEO);
	menu->GetSubMenu("File")->GetSubMenu("Insert geometry")->Add("&To new structure",MENU_FILE_INSERTGEO_NEWSTR);
	menu->GetSubMenu("File")->Add("Export selected facets",MENU_FILE_EXPORTMESH);

	menu->GetSubMenu("File")->Add("Export selected textures",MENU_FILE_EXPORTTEXTURES);

	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->Add("Cell Area (cm\262)",MENU_FILE_EXPORTTEXTURE_AREA);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->Add("# of MC Hits",MENU_FILE_EXPORTTEXTURE_MCHITS);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->Add("Impingement rate (1/s/m\262)", MENU_FILE_EXPORTTEXTURE_IMPINGEMENT);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->Add("Particle density (1/m\263)", MENU_FILE_EXPORTTEXTURE_PART_DENSITY);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->Add("Gas density (kg/m\263)", MENU_FILE_EXPORTTEXTURE_GAS_DENSITY);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->Add("Pressure (mbar)",MENU_FILE_EXPORTTEXTURE_PRESSURE);

	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->Add("Avg. Velocity (m/s)",MENU_FILE_EXPORTTEXTURE_AVG_V);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->Add("Velocity vector (m/s)",MENU_FILE_EXPORTTEXTURE_V_VECTOR);

	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->Add("# of velocity vectors",MENU_FILE_EXPORTTEXTURE_N_VECTORS);

	menu->GetSubMenu("File")->Add("Export selected profiles", MENU_FILE_EXPORTPROFILES);

	menu->GetSubMenu("File")->Add(NULL); // Separator
	menu->GetSubMenu("File")->Add("Load recent");
	
	for(int i=nbRecent-1;i>=0;i--)
		menu->GetSubMenu("File")->GetSubMenu("Load recent")->Add(recents[i],MENU_FILE_LOADRECENT+i);
	
	menu->GetSubMenu("File")->Add(NULL); // Separator
	menu->GetSubMenu("File")->Add("E&xit",MENU_FILE_EXIT);

	menu->GetSubMenu("File")->SetIcon(MENU_FILE_EXPORTMESH,126,77);//126,77
	menu->GetSubMenu("File")->SetIcon(MENU_FILE_EXPORTTEXTURES,126,77);//126,77
	menu->GetSubMenu("File")->SetIcon(MENU_FILE_SAVE,83,24);
	menu->GetSubMenu("File")->SetIcon(MENU_FILE_SAVEAS,101,24);
	menu->GetSubMenu("File")->SetIcon(MENU_FILE_LOAD,65,24);//65,24
	//menu->GetSubMenu("File")->SetIcon(MENU_FILE_LOADRECENT,83,24);//83,24

	menu->Add("Selection");
	menu->GetSubMenu("Selection")->Add("Select All Facets",MENU_FACET_SELECTALL,SDLK_a,CTRL_MODIFIER);
	menu->GetSubMenu("Selection")->Add("Select by Facet Number...",MENU_SELECTION_SELECTFACETNUMBER,SDLK_n,ALT_MODIFIER);
	menu->GetSubMenu("Selection")->Add("Select Sticking",MENU_FACET_SELECTSTICK);
	menu->GetSubMenu("Selection")->Add("Select Desoprtion",MENU_FACET_SELECTDES);
	menu->GetSubMenu("Selection")->Add("Select Outgassing Map",MENU_SELECT_HASDESFILE);
	menu->GetSubMenu("Selection")->Add("Select Transparent",MENU_FACET_SELECTTRANS);
	menu->GetSubMenu("Selection")->Add("Select Reflective",MENU_FACET_SELECTREFL);
	menu->GetSubMenu("Selection")->Add("Select 2 sided",MENU_FACET_SELECT2SIDE);
	menu->GetSubMenu("Selection")->Add("Select Texture",MENU_FACET_SELECTTEXT);
	menu->GetSubMenu("Selection")->Add("Select Profile",MENU_FACET_SELECTPROF);

	menu->GetSubMenu("Selection")->Add(NULL); // Separator
	menu->GetSubMenu("Selection")->Add("Select Abs > 0",MENU_FACET_SELECTABS);
	menu->GetSubMenu("Selection")->Add("Select Hit > 0",MENU_FACET_SELECTHITS);
	menu->GetSubMenu("Selection")->Add("Select large with no hits",MENU_FACET_SELECTNOHITS_AREA);
	menu->GetSubMenu("Selection")->Add(NULL); // Separator

	menu->GetSubMenu("Selection")->Add("Select link facets",MENU_FACET_SELECTDEST);
	menu->GetSubMenu("Selection")->Add("Select volatile facets", MENU_FACET_SELECTVOL);
	menu->GetSubMenu("Selection")->Add("Select non planar facets", MENU_FACET_SELECTNONPLANAR);
	menu->GetSubMenu("Selection")->Add("Select non simple facets",MENU_FACET_SELECTERR);
	//menu->GetSubMenu("Selection")->Add(NULL); // Separator
	//menu->GetSubMenu("Selection")->Add("Load selection",MENU_FACET_LOADSEL);
	//menu->GetSubMenu("Selection")->Add("Save selection",MENU_FACET_SAVESEL);
	menu->GetSubMenu("Selection")->Add("Invert selection",MENU_FACET_INVERTSEL,SDLK_i, CTRL_MODIFIER);
	menu->GetSubMenu("Selection")->Add(NULL); // Separator 

	menu->GetSubMenu("Selection")->Add("Memorize selection to");
	memorizeSelectionsMenu = menu->GetSubMenu("Selection")->GetSubMenu("Memorize selection to");
	memorizeSelectionsMenu->Add("Add new...",MENU_SELECTION_ADDNEW,SDLK_w,CTRL_MODIFIER);
	memorizeSelectionsMenu->Add(NULL); // Separator

	menu->GetSubMenu("Selection")->Add("Select memorized");
	selectionsMenu = menu->GetSubMenu("Selection")->GetSubMenu("Select memorized");

	menu->GetSubMenu("Selection")->Add("Clear memorized",MENU_SELECTION_CLEARSELECTIONS);
	clearSelectionsMenu = menu->GetSubMenu("Selection")->GetSubMenu("Clear memorized");
	clearSelectionsMenu->Add("Clear All",MENU_SELECTION_CLEARALL);
	clearSelectionsMenu->Add(NULL); // Separator

	menu->GetSubMenu("Selection")->Add(NULL); // Separator
	menu->GetSubMenu("Selection")->Add(NULL); // Separator
	menu->GetSubMenu("Selection")->Add("Select all vertex",MENU_VERTEX_SELECTALL);
	menu->GetSubMenu("Selection")->Add("Unselect all vertex",MENU_VERTEX_UNSELECTALL);
	menu->GetSubMenu("Selection")->Add("Select Isolated vertex",MENU_VERTEX_SELECT_ISOLATED);
	menu->GetSubMenu("Selection")->Add("Select Coplanar vertex (visible on screen)",MENU_VERTEX_SELECT_COPLANAR);


	menu->Add("Tools");
	
	menu->GetSubMenu("Tools")->Add("Add formula ..."   ,MENU_EDIT_ADDFORMULA);
	menu->GetSubMenu("Tools")->Add("Update formulas now!",MENU_EDIT_UPDATEFORMULAS,SDLK_f,ALT_MODIFIER);
	menu->GetSubMenu("Tools")->Add("Moving parts...", MENU_EDIT_MOVINGPARTS);
	menu->GetSubMenu("Tools")->Add(NULL); // Separator
	menu->GetSubMenu("Tools")->Add("Texture Plotter ...",MENU_FACET_TEXPLOTTER,SDLK_t,ALT_MODIFIER);
	menu->GetSubMenu("Tools")->Add("Profile Plotter ...",MENU_FACET_PROFPLOTTER,SDLK_p,ALT_MODIFIER);
	menu->GetSubMenu("Tools")->Add(NULL); // Separator
	menu->GetSubMenu("Tools")->Add("Texture scaling...",MENU_EDIT_TSCALING,SDLK_d,CTRL_MODIFIER);
	menu->GetSubMenu("Tools")->Add("Global Settings ..."   ,MENU_EDIT_GLOBALSETTINGS);

	menu->GetSubMenu("Tools")->SetIcon(MENU_EDIT_TSCALING,137,24);
	menu->GetSubMenu("Tools")->SetIcon(MENU_EDIT_ADDFORMULA,155,24);
	menu->GetSubMenu("Tools")->SetIcon(MENU_EDIT_GLOBALSETTINGS,0,77);

	menu->Add("Facet");
	menu->GetSubMenu("Facet")->Add("Collapse ...",MENU_FACET_COLLAPSE);
	menu->GetSubMenu("Facet")->Add("Swap normal",MENU_FACET_SWAPNORMAL,SDLK_n,CTRL_MODIFIER);
	menu->GetSubMenu("Facet")->Add("Shift vertex",MENU_FACET_SHIFTVERTEX,SDLK_h,CTRL_MODIFIER);
	menu->GetSubMenu("Facet")->Add("Edit coordinates ...",MENU_FACET_COORDINATES);
	menu->GetSubMenu("Facet")->Add("Move ...",MENU_FACET_MOVE);
	menu->GetSubMenu("Facet")->Add("Scale ...",MENU_FACET_SCALE);
	menu->GetSubMenu("Facet")->Add("Mirror ...",MENU_FACET_MIRROR);
	menu->GetSubMenu("Facet")->Add("Rotate ...",MENU_FACET_ROTATE);
	menu->GetSubMenu("Facet")->Add("Align ...",MENU_FACET_ALIGN);
	menu->GetSubMenu("Facet")->Add("Remove selected",MENU_FACET_REMOVESEL,SDLK_DELETE,CTRL_MODIFIER);
	menu->GetSubMenu("Facet")->Add("Explode selected",MENU_FACET_EXPLODE);
	menu->GetSubMenu("Facet")->Add("Create difference of 2",MENU_FACET_CREATE_DIFFERENCE);
	menu->GetSubMenu("Facet")->Add("Convert to outgassing map...",MENU_FACET_OUTGASSINGMAP);
	menu->GetSubMenu("Facet")->Add(NULL); // Separator

	menu->GetSubMenu("Facet")->Add("Facet Details ...",MENU_FACET_DETAILS);
	//menu->GetSubMenu("Facet")->Add("Facet Mesh ...",MENU_FACET_MESH);

	//facetMenu = menu->GetSubMenu("Facet");
	//facetMenu->SetEnabled(MENU_FACET_MESH,FALSE);

	menu->GetSubMenu("Facet")->SetIcon(MENU_FACET_COLLAPSE,173,24);
	menu->GetSubMenu("Facet")->SetIcon(MENU_FACET_SWAPNORMAL,191,24);
	menu->GetSubMenu("Facet")->SetIcon(MENU_FACET_SHIFTVERTEX,90,77);
	menu->GetSubMenu("Facet")->SetIcon(MENU_FACET_COORDINATES,209,24);
	menu->GetSubMenu("Facet")->SetIcon(MENU_FACET_PROFPLOTTER,227,24);
	menu->GetSubMenu("Facet")->SetIcon(MENU_FACET_DETAILS,54,77);
	//menu->GetSubMenu("Facet")->SetIcon(MENU_FACET_MESH,72,77);
	menu->GetSubMenu("Facet")->SetIcon(MENU_FACET_TEXPLOTTER,108,77);

	menu->Add("Vertex");
	menu->GetSubMenu("Vertex")->Add("Create Facet from Selected");
	menu->GetSubMenu("Vertex")->GetSubMenu("Create Facet from Selected")->Add("Convex Hull",MENU_VERTEX_CREATE_POLY_CONVEX,SDLK_v,ALT_MODIFIER);
	menu->GetSubMenu("Vertex")->GetSubMenu("Create Facet from Selected")->Add("Keep selection order",MENU_VERTEX_CREATE_POLY_ORDER);
	menu->GetSubMenu("Vertex")->Add("Remove selected",MENU_VERTEX_REMOVE);
	menu->GetSubMenu("Vertex")->Add("Vertex coordinates...",MENU_VERTEX_COORDINATES);
	menu->GetSubMenu("Vertex")->Add("Move selected...",MENU_VERTEX_MOVE);
	menu->GetSubMenu("Vertex")->Add("Scale selected...",MENU_VERTEX_SCALE);
	menu->GetSubMenu("Vertex")->Add("Add new...",MENU_VERTEX_ADD);

	menu->Add("View");

	menu->GetSubMenu("View")->Add("Structure",MENU_VIEW_STRUCTURE_P);
	structMenu = menu->GetSubMenu("View")->GetSubMenu("Structure");
	UpdateStructMenu();

	menu->GetSubMenu("View")->Add("Full Screen",MENU_VIEW_FULLSCREEN);

	menu->GetSubMenu("View")->Add(NULL); // Separator 

	menu->GetSubMenu("View")->Add("Memorize view to");
	memorizeViewsMenu = menu->GetSubMenu("View")->GetSubMenu("Memorize view to");
	memorizeViewsMenu->Add("Add new...",MENU_VIEW_ADDNEW,SDLK_q,CTRL_MODIFIER);
	memorizeViewsMenu->Add(NULL); // Separator

	menu->GetSubMenu("View")->Add("Select memorized");
	viewsMenu = menu->GetSubMenu("View")->GetSubMenu("Select memorized");

	menu->GetSubMenu("View")->Add("Clear memorized",MENU_VIEW_CLEARVIEWS);
	clearViewsMenu = menu->GetSubMenu("View")->GetSubMenu("Clear memorized");
	clearViewsMenu->Add("Clear All",MENU_VIEW_CLEARALL);

	//menu->GetSubMenu("View")->SetIcon(MENU_VIEW_STRUCTURE_P,0,77);
	menu->GetSubMenu("View")->SetIcon(MENU_VIEW_FULLSCREEN,18,77);
	//menu->GetSubMenu("View")->SetIcon(MENU_VIEW_ADD,36,77);

	menu->Add("Time");
	menu->GetSubMenu("Time")->Add("Time settings...",MENU_TIME_SETTINGS,SDLK_i,ALT_MODIFIER);
	menu->GetSubMenu("Time")->Add("Edit moments...",MENU_TIME_MOMENTS_EDITOR);
	menu->GetSubMenu("Time")->Add("Edit parameters...", MENU_TIME_PARAMETER_EDITOR);
	menu->GetSubMenu("Time")->Add(NULL);
	menu->GetSubMenu("Time")->Add("Timewise plotter",MENU_TIMEWISE_PLOTTER);
	menu->GetSubMenu("Time")->Add("Pressure evolution",MENU_TIME_PRESSUREEVOLUTION);

	menu->Add("Test");
	menu->GetSubMenu("Test")->Add("Pipe (L/R=0.0001)",MENU_TEST_PIPE0001);
	menu->GetSubMenu("Test")->Add("Pipe (L/R=1)",MENU_TEST_PIPE1);
	menu->GetSubMenu("Test")->Add("Pipe (L/R=10)",MENU_TEST_PIPE10);
	menu->GetSubMenu("Test")->Add("Pipe (L/R=100)",MENU_TEST_PIPE100);
	menu->GetSubMenu("Test")->Add("Pipe (L/R=1000)",MENU_TEST_PIPE1000);
	menu->GetSubMenu("Test")->Add("Pipe (L/R=10000)",MENU_TEST_PIPE10000);
	//Quick test pipe
	menu->GetSubMenu("Test")->Add(NULL);
	menu->GetSubMenu("Test")->Add("Quick Pipe",MENU_QUICKPIPE,SDLK_q,ALT_MODIFIER);

	/*profilePlotterShortcut=new GLButton(0,"Profile Plotter");
	Add(profilePlotterShortcut);

	texturePlotterShortcut=new GLButton(0,"Texture Plotter");
	Add(texturePlotterShortcut);*/


	geomNumber = new GLTextField(0,NULL);
	geomNumber->SetEditable(FALSE);
	Add(geomNumber);

	togglePanel = new GLTitledPanel("3D Viewer settings");
	togglePanel->SetClosable(TRUE);
	Add(togglePanel);

	showNormal = new GLToggle(0,"Normals");
	togglePanel->Add(showNormal);

	showRule = new GLToggle(0,"Rules");
	togglePanel->Add(showRule);

	showUV = new GLToggle(0,"\201,\202");
	togglePanel->Add(showUV);

	showLeak = new GLToggle(0,"Leaks");
	togglePanel->Add(showLeak);

	showHit = new GLToggle(0,"Hits");
	togglePanel->Add(showHit);

	showLine = new GLToggle(0,"Lines");
	togglePanel->Add(showLine);

	showVolume = new GLToggle(0,"Volume");
	togglePanel->Add(showVolume);

	showTexture = new GLToggle(0,"Texture");
	togglePanel->Add(showTexture);

	showFilter = new GLToggle(0,"Filtering");
	//togglePanel->Add(showFilter);

	showIndex = new GLToggle(0,"Indices");
	togglePanel->Add(showIndex);

	showVertex = new GLToggle(0,"Vertices");
	togglePanel->Add(showVertex);

	showMoreBtn = new GLButton(0,"<< More");
	togglePanel->Add(showMoreBtn);

	shortcutPanel = new GLTitledPanel("Shortcuts");
	shortcutPanel->SetClosable(TRUE);
	Add(shortcutPanel);

	profilePlotterBtn = new GLButton(0, "Profile pl.");
	shortcutPanel->Add(profilePlotterBtn);

	texturePlotterBtn = new GLButton(0, "Texture pl.");
	shortcutPanel->Add(texturePlotterBtn);

	textureScalingBtn = new GLButton(0, "Tex.scaling");
	shortcutPanel->Add(textureScalingBtn);

	simuPanel = new GLTitledPanel("Simulation");
	simuPanel->SetClosable(TRUE);
	Add(simuPanel);

	globalSettingsBtn = new GLButton(0,"<< More");
	simuPanel->Add(globalSettingsBtn);

	startSimu = new GLButton(0,"Start/Stop");
	//startSimu->SetEnabled(FALSE);
	simuPanel->Add(startSimu);

	resetSimu = new GLButton(0,"Reset");
	//resetSimu->SetEnabled(FALSE);
	simuPanel->Add(resetSimu);

	/*
	statusSimu = new GLButton(0,"...");
	simuPanel->Add(statusSimu);
	*/

	modeLabel = new GLLabel("Mode");
	simuPanel->Add(modeLabel);

	modeCombo = new GLCombo(0);
	modeCombo->SetEditable(TRUE);
	modeCombo->SetSize(2);
	modeCombo->SetValueAt(0,"Monte Carlo");
	modeCombo->SetValueAt(1,"Angular Coef");
	modeCombo->SetSelectedIndex(0);
	simuPanel->Add(modeCombo);

	compACBtn = new GLButton(0,"Calc AC");
	compACBtn->SetEnabled(FALSE);
	simuPanel->Add(compACBtn);

	singleACBtn = new GLButton(0,"1");
	singleACBtn->SetEnabled(FALSE);
	//simuPanel->Add(singleACBtn);

	hitLabel = new GLLabel("Hits");
	simuPanel->Add(hitLabel);

	hitNumber = new GLTextField(0,NULL);
	hitNumber->SetEditable(FALSE);
	simuPanel->Add(hitNumber);

	desLabel = new GLLabel("Des.");
	simuPanel->Add(desLabel);

	desNumber = new GLTextField(0,NULL);
	desNumber->SetEditable(FALSE);
	simuPanel->Add(desNumber);

	leakLabel = new GLLabel("Leaks");
	simuPanel->Add(leakLabel);

	leakNumber = new GLTextField(0,NULL);
	leakNumber->SetEditable(FALSE);
	simuPanel->Add(leakNumber);

	sTimeLabel = new GLLabel("Time");
	simuPanel->Add(sTimeLabel);

	sTime = new GLTextField(0,NULL);
	sTime->SetEditable(FALSE);
	simuPanel->Add(sTime);

	facetPanel = new GLTitledPanel("Selected Facet");
	facetPanel->SetClosable(TRUE);
	Add(facetPanel);

	inputPanel = new GLTitledPanel("Particles in");
	facetPanel->Add(inputPanel);

	facetDLabel = new GLLabel("Desorption");
	facetPanel->Add(facetDLabel);
	facetDesType = new GLCombo(0);
	facetDesType->SetSize(4);
	facetDesType->SetValueAt(0,"None");
	facetDesType->SetValueAt(1,"Uniform");
	facetDesType->SetValueAt(2,"Cosine");
	facetDesType->SetValueAt(3,"Cosine^N");
	inputPanel->Add(facetDesType);

	facetDesTypeN = new GLTextField(0,NULL);
	facetDesTypeN->SetEditable(FALSE);
	facetPanel->Add(facetDesTypeN);

	facetFILabel = new GLToggle(0,"Outgassing (mbar*l/s):");
	facetFILabel->SetEnabled(FALSE);
	facetFILabel->SetState(TRUE);
	inputPanel->Add(facetFILabel);
	facetFlow = new GLTextField(0,NULL);
	inputPanel->Add(facetFlow);

	facetFIAreaLabel = new GLToggle(1,"Outg/area(mbar*l/s/cm\262):");
	facetFIAreaLabel->SetEnabled(FALSE);
	inputPanel->Add(facetFIAreaLabel);
	facetFlowArea = new GLTextField(0,NULL);
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
	facetSticking = new GLTextField(0,NULL);
	outputPanel->Add(facetSticking);

	facetPumpingLabel = new GLLabel("Pumping Speed (l/s):");
	outputPanel->Add(facetPumpingLabel);
	facetPumping = new GLTextField(0,NULL);

	outputPanel->Add(facetPumping);

	facetSideLabel = new GLLabel("Sides:");
	facetPanel->Add(facetSideLabel);

	facetSideType = new GLCombo(0);
	facetSideType->SetSize(2);
	facetSideType->SetValueAt(0,"1 Sided");
	facetSideType->SetValueAt(1,"2 Sided");
	facetPanel->Add(facetSideType);

	facetTLabel = new GLLabel("Opacity:");
	facetPanel->Add(facetTLabel);
	facetOpacity = new GLTextField(0,NULL);
	facetPanel->Add(facetOpacity);

	facetTempLabel = new GLLabel("Temperature (\260K):");
	facetPanel->Add(facetTempLabel);
	facetTemperature = new GLTextField(0,NULL);
	facetPanel->Add(facetTemperature);

	facetAreaLabel = new GLLabel("Area (cm\262):");
	facetPanel->Add(facetAreaLabel);
	facetArea = new GLTextField(0,NULL);
	facetPanel->Add(facetArea);

	facetReLabel = new GLLabel("Profile:");
	facetPanel->Add(facetReLabel);
	facetRecType = new GLCombo(0);
	facetRecType->SetSize(6);
	facetRecType->SetValueAt(0,"None");
	facetRecType->SetValueAt(1,"Pressure (along \201)");
	facetRecType->SetValueAt(2,"Pressure (along \202)");
	facetRecType->SetValueAt(3,"Angular");
	facetRecType->SetValueAt(4,"Speed distribution");
	facetRecType->SetValueAt(5,"Orthogonal velocity");
	facetPanel->Add(facetRecType);

	facetMoreBtn = new GLButton(0,"<< More");
	facetPanel->Add(facetMoreBtn);

	facetDetailsBtn = new GLButton(0,"Details...");
	facetPanel->Add(facetDetailsBtn);

	facetCoordBtn = new GLButton(0,"Coord...");
	facetPanel->Add(facetCoordBtn);

	facetApplyBtn = new GLButton(0,"Apply");
	facetApplyBtn->SetEnabled(FALSE);
	facetPanel->Add(facetApplyBtn);

	facetList = new GLList(0);
	facetList->SetWorker(&worker);
	facetList->SetGrid(TRUE);
	facetList->SetSelectionMode(MULTIPLE_ROW);
	facetList->SetSize(4,1);
	facetList->SetColumnWidths((int*)cWidth);
	facetList->SetColumnLabels((char **)cName);
	facetList->SetColumnLabelVisible(TRUE);
	facetList->Sortable = TRUE;
	Add(facetList);

	ClearFacetParams();
	UpdateViewerParams();
	PlaceComponents();
	LoadConfig();
	//LoadFile();
	try {

		worker.SetProcNumber(nbProc);

	} catch (Error &e) {
		char errMsg[512];
		sprintf(errMsg,"Failed to start working sub-process(es), simulation not available\n%s",e.GetMsg());
		GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
	}

	//AnimateViewerChange(0,TRUE);

	PlaceComponents();

	//SelectViewer(0);

	//viewer[0]->Paint();

	if (checkForUpdates) {
		//Launch updater tool
		char command[1024];
		char CWD [MAX_PATH];
		_getcwd( CWD, MAX_PATH );
		//sprintf(tmp5,"%s\\molflow_updater_tmp.exe",CWD);
		if(FileUtils::Exist("molflow_updater_tmp.exe")) { //rename after new installation
			sprintf(command,"move \"%s\\molflow_updater_tmp.exe\" \"%s\\molflow_updater.exe\"",CWD,CWD);
			system(command);
		}

		if(FileUtils::Exist("molflow_updater.exe"))
			StartProc_background("molflow_updater.exe");
		else GLMessageBox::Display("molflow_updater.exe not found. You will not receive updates to Molflow."
			"\n(You can disable checking for updates in Tools/Global Settings)","Updater module missing.",GLDLG_OK,GLDLG_ICONINFO);
	}
	return GL_OK;
}

//-----------------------------------------------------------------------------
// Name: Resize()
// Desc: Called when the window is resized
//-----------------------------------------------------------------------------
int MolFlow::Resize(DWORD width, DWORD height, BOOL forceWindowed) {
	int r = GLApplication::Resize(width,height,forceWindowed);
	PlaceComponents();
	return r;
}

//-----------------------------------------------------------------------------
// Name: PlaceComponents()
// Desc: Place components on screen
//-----------------------------------------------------------------------------
void MolFlow::PlaceComponents() {

	int sx = m_screenWidth-205;  
	int sy = 30;

	Place3DViewer();

	geomNumber->SetBounds(sx,3,202,18);

	// Viewer settings ----------------------------------------
	togglePanel->SetBounds(sx,sy,202,112);

	togglePanel->SetCompBounds(showRule,5,20,60,18);
	togglePanel->SetCompBounds(showNormal,70,20,60,18);
	togglePanel->SetCompBounds(showUV,135,20,60,18);

	togglePanel->SetCompBounds(showLine,5,42,60,18);
	togglePanel->SetCompBounds(showLeak,70,42,60,18);
	togglePanel->SetCompBounds(showHit,135,42,60,18);

	togglePanel->SetCompBounds(showVolume,5,64,60,18);
	togglePanel->SetCompBounds(showTexture,70,64,60,18);
	togglePanel->SetCompBounds(showFilter,135,64,60,18);

	togglePanel->SetCompBounds(showMoreBtn,5,86,55,19);
	togglePanel->SetCompBounds(showVertex,70,86,60,18);
	togglePanel->SetCompBounds(showIndex,137,86,60,18);

	sy += (togglePanel->GetHeight() + 5);

	// Selected facet -----------------------------------------
	facetPanel->SetBounds(sx,sy,202,330);

	facetPanel->SetCompBounds(inputPanel,5,16,192,90);

	int cursorY = 15;
	inputPanel->SetCompBounds(facetDLabel, 5, cursorY, 60, 18);
	inputPanel->SetCompBounds(facetDesType, 65, cursorY, 80, 18);
	inputPanel->SetCompBounds(facetDesTypeN, 150, cursorY, 30, 18);

	inputPanel->SetCompBounds(facetFILabel, 5, cursorY+=25, 110, 18);
	inputPanel->SetCompBounds(facetFlow, 140, cursorY, 45, 18);

	inputPanel->SetCompBounds(facetFIAreaLabel, 5, cursorY+=25, 110, 18);
	inputPanel->SetCompBounds(facetFlowArea, 140, cursorY, 45, 18);

	//inputPanel->SetCompBounds(facetUseDesFileLabel,5,90,60,18);
	//inputPanel->SetCompBounds(facetUseDesFile,65,90,120,18);

	facetPanel->SetCompBounds(outputPanel, 5, cursorY+=45, 192, 65);

	outputPanel->SetCompBounds(facetSLabel, 7, cursorY=15, 100, 18);
	outputPanel->SetCompBounds(facetSticking, 140, cursorY, 45, 18);

	outputPanel->SetCompBounds(facetPumpingLabel, 7, cursorY+=25, 100, 18);
	outputPanel->SetCompBounds(facetPumping, 140, cursorY, 45, 18);

	facetPanel->SetCompBounds(facetSideLabel, 7, cursorY=180, 50, 18);
	facetPanel->SetCompBounds(facetSideType, 65, cursorY, 130, 18);

	facetPanel->SetCompBounds(facetTLabel, 7, cursorY+=25, 100, 18);
	facetPanel->SetCompBounds(facetOpacity, 110, cursorY, 82, 18);

	facetPanel->SetCompBounds(facetTempLabel, 7, cursorY+=25, 100, 18);
	facetPanel->SetCompBounds(facetTemperature, 110, cursorY, 82, 18);

	facetPanel->SetCompBounds(facetAreaLabel, 7, cursorY+=25, 100, 18);
	facetPanel->SetCompBounds(facetArea, 110, cursorY, 82, 18);

	facetPanel->SetCompBounds(facetReLabel, 7, cursorY+=25, 60, 18);
	facetPanel->SetCompBounds(facetRecType, 65, cursorY, 130, 18);

	facetPanel->SetCompBounds(facetMoreBtn, 5, cursorY+=25, 48, 18);
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
	simuPanel->SetBounds(sx,sy,202,169);

	simuPanel->SetCompBounds(globalSettingsBtn, 5, 20, 48, 19);
	simuPanel->SetCompBounds(startSimu,58,20,66,19);
	simuPanel->SetCompBounds(resetSimu,128,20,66,19);
	//simuPanel->SetCompBounds(statusSimu,175,20,20,19);
	simuPanel->SetCompBounds(modeLabel,5,45,30,18);
	simuPanel->SetCompBounds(modeCombo,40,45,85,18);
	simuPanel->SetCompBounds(compACBtn,130,45,65,19);
	//simuPanel->SetCompBounds(compACBtn,123,45,52,19);
	//simuPanel->SetCompBounds(singleACBtn,178,45,20,19);
	simuPanel->SetCompBounds(hitLabel,5,70,30,18);
	simuPanel->SetCompBounds(hitNumber,40,70,155,18);
	simuPanel->SetCompBounds(desLabel,5,95,30,18);
	simuPanel->SetCompBounds(desNumber,40,95,155,18);
	simuPanel->SetCompBounds(leakLabel,5,120,30,18);
	simuPanel->SetCompBounds(leakNumber,40,120,155,18);
	simuPanel->SetCompBounds(sTimeLabel,5,145,30,18);


	simuPanel->SetCompBounds(sTime,40,145,155,18);

	sy += (simuPanel->GetHeight() + 5);

	// ---------------------------------------------------------
	int lg = m_screenHeight - (nbFormula*25 + 23);

	facetList->SetBounds(sx,sy,202,lg-sy);

	// ---------------------------------------------------------

	for(int i=0;i<nbFormula;i++) {
		formulas[i].name->SetBounds(sx,lg+5,95,18);
		formulas[i].value->SetBounds(sx+90,lg+5,87,18);
		formulas[i].setBtn->SetBounds(sx+182,lg+5,20,18);
		lg+=25;
	}

}

void MolFlow::Place3DViewer() {

	int sx = m_screenWidth-205;

	// 3D Viewer ----------------------------------------------
	int fWidth = m_screenWidth-215;
	int fHeight = m_screenHeight-27;
	int Width2 = fWidth/2-1;
	int Height2 = fHeight/2-1;

	if( modeSolo ) {
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->SetVisible(FALSE);
		viewer[curViewer]->SetBounds(3,3,fWidth,fHeight);
		viewer[curViewer]->SetVisible(TRUE);
	} else {
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->SetVisible(TRUE);
		viewer[0]->SetBounds(3       ,3        ,Width2,Height2);
		viewer[1]->SetBounds(6+Width2,3        ,Width2,Height2);
		viewer[2]->SetBounds(3       ,6+Height2,Width2,Height2);
		viewer[3]->SetBounds(6+Width2,6+Height2,Width2,Height2);
	}

}

void MolFlow::UpdateViewers() {
	for(int i=0;i<MAX_VIEWER;i++)
		viewer[i]->UpdateMatrix();
}

void MolFlow::SetFacetSearchPrg(BOOL visible,char *text) {
	for(int i=0;i<MAX_VIEWER;i++) {
		viewer[i]->facetSearchState->SetVisible(visible);
		viewer[i]->facetSearchState->SetText(text);
	}
	GLWindowManager::Repaint();
}

void MolFlow::UpdateViewerParams() {

	showNormal->SetState(viewer[curViewer]->showNormal);
	showRule->SetState(viewer[curViewer]->showRule);
	showUV->SetState(viewer[curViewer]->showUV);
	showLeak->SetState(viewer[curViewer]->showLeak);
	showHit->SetState(viewer[curViewer]->showHit);
	showVolume->SetState(viewer[curViewer]->showVolume);
	showLine->SetState(viewer[curViewer]->showLine);
	showTexture->SetState(viewer[curViewer]->showTexture);
	showFilter->SetState(viewer[curViewer]->showFilter);
	showVertex->SetState(viewer[curViewer]->showVertex);
	showIndex->SetState(viewer[curViewer]->showIndex);

	// Force all views to have the same showColormap
	viewer[1]->showColormap = viewer[0]->showColormap;
	viewer[2]->showColormap = viewer[0]->showColormap;
	viewer[3]->showColormap = viewer[0]->showColormap;
	worker.GetGeometry()->texColormap = viewer[0]->showColormap;

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
	if( facetSticking->GetNumber(&sticking) ) {
		if( sticking<0.0 || sticking>1.0 ) {
			GLMessageBox::Display("Sticking must be in the range [0,1]","Error",GLDLG_OK,GLDLG_ICONERROR);
			UpdateFacetParams();
			return;
		}
		doSticking = TRUE;
		stickingNotNumber = FALSE;
	} else {
		if( strcmp(facetSticking->GetText(),"..." )==0 ) doSticking = FALSE;
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
	if( facetOpacity->GetNumber(&opacity) ) {
		if( opacity<0.0 || opacity>1.0 ) {
			GLMessageBox::Display("Opacity must be in the range [0,1]","Error",GLDLG_OK,GLDLG_ICONERROR);
			UpdateFacetParams();
			return;
		}
		doOpacity = TRUE;
		opacityNotNumber = FALSE;
	} else {
		if( strcmp(facetOpacity->GetText(),"..." )==0 ) doOpacity = FALSE;
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
	if( facetTemperature->GetNumber(&temperature) ) {
		if( temperature<0.0 ) {
			GLMessageBox::Display("Temperature must be positive","Error",GLDLG_OK,GLDLG_ICONERROR);
			UpdateFacetParams();
			return;
		}
		doTemperature = TRUE;
	} else {
		if( strcmp(facetTemperature->GetText(),"..." )==0 ) doTemperature = FALSE;
		else {
			GLMessageBox::Display("Invalid temperature number","Error",GLDLG_OK,GLDLG_ICONERROR);
			UpdateFacetParams();
			return;
		}
	}

	// Outgassing
	double flow=0.0;
	BOOL doFlow = FALSE;
	BOOL outgassingNotNumber;
	//Calculate flow
	if (facetFILabel->GetState() && strcmp(facetFlow->GetText(), "...") != 0 && facetDesType->GetSelectedIndex() != 0
		&& strcmp(facetDesType->GetSelectedValue(), "...") != 0) {  //We want outgassing
		if (facetFlow->GetNumber(&flow)) { //If we can parse the number
			if (!(flow > 0.0) && !(facetMesh && facetMesh->facetUseDesFile->GetSelectedIndex() == 1)) {
				GLMessageBox::Display("Outgassing must be positive", "Error", GLDLG_OK, GLDLG_ICONERROR);
				UpdateFacetParams();
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
	double flowA=0;
	BOOL doFlowA = FALSE;
	//Calculate flow
	
	if( facetFIAreaLabel->GetState() && strcmp(facetFlowArea->GetText(),"..." )!=0
		&& facetDesType->GetSelectedIndex()!=0 && strcmp(facetDesType->GetSelectedValue(),"..." )!=0) { //We want outgassing per area
		if (facetFlowArea->GetNumber(&flowA)) { //Can be parsed as number
			if( !(flowA>0.0) ) {
				GLMessageBox::Display("Outgassing per area must be positive","Error",GLDLG_OK,GLDLG_ICONERROR);
				UpdateFacetParams();
				return;
			}
			doFlowA = TRUE;
		} else {
			GLMessageBox::Display("Invalid outgassing per area number","Error",GLDLG_OK,GLDLG_ICONERROR);
			UpdateFacetParams();
			return;
		}
	}

	// Desorption type
	int desorbType = facetDesType->GetSelectedIndex();

	double desorbTypeN;
	BOOL doDesorbTypeN = FALSE;
	if (desorbType==3) {
		if( facetDesTypeN->GetNumber(&desorbTypeN) ) {
			if( !(desorbTypeN>0.0) ) {
				GLMessageBox::Display("Desorption type exponent must be greater than 0.0","Error",GLDLG_OK,GLDLG_ICONERROR);
				UpdateFacetParams();
				return;
			}
			doDesorbTypeN = TRUE;
		} else {
			if( strcmp(facetDesTypeN->GetText(),"..." )==0 ) doDesorbTypeN = FALSE;
			else {
				GLMessageBox::Display("Invalid desorption type exponent","Error",GLDLG_OK,GLDLG_ICONERROR);
				UpdateFacetParams();
				return;
			}
		}
	}

	// Record (profile) type
	int rType = facetRecType->GetSelectedIndex(); // -1 if "..."

	// 2sided
	int is2Sided = facetSideType->GetSelectedIndex();

	
	// Update facets (local)
	for(int i=0;i<nbFacet;i++) {
		Facet *f = geom->GetFacet(i);
		if( f->selected ) {
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
			if(doTemperature) f->sh.temperature = temperature;
			
			if(doFlowA/* && !useMapA*/) f->sh.flow = flowA*f->sh.area*(f->sh.is2sided?2.0:1.0)*0.100;
			if(desorbType>=0) {
				if(desorbType==0) f->sh.flow=0.0;
				if(desorbType!=3) f->sh.desorbTypeN=0.0;
				f->sh.desorbType = desorbType;
				if (doDesorbTypeN) f->sh.desorbTypeN = desorbTypeN;
			}

			
			if(rType>=0) {
				f->sh.profileType = rType;
				//f->sh.isProfile = (rType!=REC_NONE); //included below by f->UpdateFlags();
			}
			if(is2Sided>=0) f->sh.is2sided = is2Sided;
			
			f->sh.maxSpeed = 4.0 * sqrt(2.0*8.31*f->sh.temperature / 0.001 / worker.gasMass);
			f->UpdateFlags();
		}
	}
	if (facetMesh && facetMesh->IsVisible()) {
		if (!facetMesh->Apply()) {
			UpdateFacetParams();
			return;
		}
	}
	
	// Mark "needsReload"
	try { worker.Reload(); } catch(Error &e) {
		GLMessageBox::Display((char *)e.GetMsg(),"Error",GLDLG_OK,GLDLG_ICONERROR);
		return;
	}
	geom->CalcTotalOutGassing();
	UpdateFacetParams();
	if (profilePlotter) profilePlotter->Refresh();
	if (pressureEvolution) pressureEvolution->Refresh();
	if (timewisePlotter) timewisePlotter->Refresh();
	//if (facetMesh) facetMesh->Refresh();
}

void MolFlow::UpdateFacetlistSelected() {
	int nbSelected=0;
	Geometry *geom = worker.GetGeometry();
	int nbFacet=geom->GetNbFacet();
	int* selection=(int*)malloc(nbFacet*sizeof(int));
	for (int i=0;i<nbFacet;i++) {
		if (geom->GetFacet(i)->selected) {
			selection[nbSelected]=i;
			nbSelected++;
		}
	}

	//facetList->SetSelectedRows(selection,nbSelected,TRUE);
	if (nbSelected>1000) {
		facetList->ReOrder();
		facetList->SetSelectedRows(selection,nbSelected,FALSE);
	} else {
		facetList->SetSelectedRows(selection,nbSelected,TRUE);
	}
	SAFE_FREE(selection);
}

//-----------------------------------------------------------------------------
// Name: UpdateFacetParams()
// Desc: Update selected facet parameters.
//-----------------------------------------------------------------------------

void MolFlow::UpdateFacetParams(BOOL updateSelection) {

	char tmp[256];

	// Update params
	Geometry *geom = worker.GetGeometry();
	int nbSel = geom->GetNbSelected();
	if(nbSel>0) {

		Facet *f0;
		Facet *f;

		// Get list of selected facet
		int *selection = (int *)malloc(nbSel*sizeof(int));
		int count=0;
		for(int i=0;i<geom->GetNbFacet();i++)
			if( geom->GetFacet(i)->selected )
				selection[count++] = i;

		f0 = geom->GetFacet(selection[0]);

		double area = f0->sh.area*(f0->sh.is2sided?2.0:1.0);

		BOOL stickingE = TRUE;
		BOOL opacityE = TRUE;
		BOOL temperatureE = TRUE;
		BOOL flowE = TRUE;
		BOOL flowAreaE = TRUE;
		BOOL desorbTypeE = TRUE;
		BOOL desorbTypeNE = TRUE;
		BOOL recordE = TRUE;
		BOOL is2sidedE = TRUE;

		for(int i=1;i<count;i++) {
			f = geom->GetFacet(selection[i]);
			stickingE = stickingE && (f0->userSticking.compare(f->userSticking)==0 && abs(f0->sh.sticking - f->sh.sticking)<1e-7);
			opacityE = opacityE && (f0->userOpacity.compare(f->userOpacity) == 0 && abs(f0->sh.opacity - f->sh.opacity)<1e-7);
			temperatureE = temperatureE && (abs(f0->sh.temperature - f->sh.temperature)<1e-7);
			flowE = flowE && (f0->userOutgassing.compare(f->userOutgassing) == 0 && abs(f0->sh.flow - f->sh.flow)<1e-7);
			flowAreaE = flowAreaE && (abs(f0->sh.flow/f0->sh.area/(f0->sh.is2sided?2.0:1.0) - f->sh.flow/f->sh.area/(f->sh.is2sided?2.0:1.0))<1e-20);
			is2sidedE = is2sidedE && (f0->sh.is2sided == f->sh.is2sided);
			desorbTypeE = desorbTypeE && (f0->sh.desorbType == f->sh.desorbType);
			desorbTypeNE = desorbTypeNE && (abs(f0->sh.desorbTypeN - f->sh.desorbTypeN)<1e-7);
			recordE = recordE && (f0->sh.profileType == f->sh.profileType);  //profiles
			if (f->sh.area>0) area+=f->sh.area*(f->sh.is2sided?2.0:1.0);
		}

		if( nbSel==1 )
			sprintf(tmp,"Selected Facet (#%d)",selection[0]+1);
		else
			sprintf(tmp,"Selected Facet (%d selected)",count);

		// Old STR compatibility
		if( stickingE && f0->sh.superDest ) stickingE = FALSE;

		facetPanel->SetTitle(tmp);
		if (count>1) facetAreaLabel->SetText("Sum Area (cm\262):");
		else facetAreaLabel->SetText("Area (cm\262):");
		facetArea->SetText(area);
		if (stickingE) {
			if (f0->userSticking.length() == 0)
				facetSticking->SetText(f0->sh.sticking);
			else facetSticking->SetText(f0->userSticking.c_str());
		} else facetSticking->SetText("...");
		
		if (opacityE) {
			if (f0->userOpacity.length() == 0)
				facetOpacity->SetText(f0->sh.opacity);
			else facetOpacity->SetText(f0->userOpacity.c_str());
		}
		else facetOpacity->SetText("...");
		
		if(flowAreaE) facetFlowArea->SetText(f0->sh.flow/f0->sh.area/(f0->sh.is2sided?2.0:1.0)*10.00); else facetFlowArea->SetText("...");
		if(temperatureE) facetTemperature->SetText(f0->sh.temperature); else facetTemperature->SetText("...");
		if(is2sidedE) facetSideType->SetSelectedIndex(f0->sh.is2sided); else facetSideType->SetSelectedValue("...");
		if(desorbTypeE) facetDesType->SetSelectedIndex(f0->sh.desorbType); else facetDesType->SetSelectedValue("...");
		if (desorbTypeNE) facetDesTypeN->SetText(f0->sh.desorbTypeN); else facetDesTypeN->SetText("...");
		if(recordE) facetRecType->SetSelectedIndex(f0->sh.profileType); else facetRecType->SetSelectedValue("...");

		if (count==1) {
			facetPumping->SetEditable(TRUE);
			calcFlow();
		}
		else{
			facetPumping->SetEditable(FALSE);
			facetPumping->SetText("...");
		}

		if( desorbTypeE && f0->sh.desorbType>0 ) {

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
			if(flowAreaE) facetFlowArea->SetText(f0->sh.flow/f0->sh.area/(f0->sh.is2sided?2.0:1.0)*10.00); else facetFlowArea->SetText("...");
			if (f0->sh.desorbType==3) {
				facetDesTypeN->SetEditable(TRUE);
				if(desorbTypeNE) facetDesTypeN->SetText(f0->sh.desorbTypeN); else facetDesTypeN->SetText("...");
			} else {
				facetDesTypeN->SetText("");
				facetDesTypeN->SetEditable(FALSE);
			};

		} else {
			facetFILabel->SetEnabled(FALSE);
			facetFlow->SetEditable(FALSE);
			facetFIAreaLabel->SetEnabled(FALSE);
			facetFlowArea->SetEditable(FALSE);
			facetDesTypeN->SetText("");
			facetDesTypeN->SetEditable(FALSE);
			//facetFlow->SetText("");
		}

		if( updateSelection ) {
			if (facetMesh) facetMesh->Refresh(count, selection);
			if (nbSel>1000 || geom->GetNbFacet()>50000) { //If it would take too much time to look up every selected facet in the list
				facetList->ReOrder();
				facetList->SetSelectedRows(selection,nbSel,FALSE);
			} else {
				facetList->SetSelectedRows(selection,nbSel,TRUE);
			}
			facetList->lastRowSel=-1;
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
	} else {
		ClearFacetParams();
		if (facetMesh) facetMesh->Refresh(0,NULL); //Clear
		if( updateSelection ) facetList->ClearSelection();
	}
	
	if( facetDetails ) facetDetails->Update();
	if( facetCoordinates ) facetCoordinates->UpdateFromSelection();
	if( texturePlotter ) texturePlotter->Update(m_fTime,TRUE);
	if( outgassingMap ) outgassingMap->Update(m_fTime,TRUE);
}
//-----------------------------------------------------------------------------
// Name: Evaluate formula
// Desc: Update formula
//-----------------------------------------------------------------------------

int getVariable(char *name,char *preffix) {

	char tmp[256];
	int  idx;
	int lgthP = (int)strlen(preffix);
	int lgthN = (int)strlen(name);

	if( lgthP>=lgthN ) {
		return -1;
	} else {
		strcpy(tmp,name);
		tmp[lgthP]=0;
		if( _stricmp(tmp,preffix)==0 ) {
			strcpy(tmp,name+lgthP);
			int conv = sscanf(tmp,"%d",&idx);
			if(conv) {
				return idx;
			} else {
				return -1;
			}
		}
	}
	return -1;

}


void MolFlow::UpdateFormula() {

	char tmp[256];

	int idx;
	Geometry *geom = worker.GetGeometry();
	int nbFacet = geom->GetNbFacet();

	for(int i=0;i<nbFormula;i++) {

		GLParser *f = formulas[i].parser;

		// Variables
		int nbVar = f->GetNbVariable();
		BOOL ok = TRUE;
		BOOL iName = FALSE;
		for(int j=0;j<nbVar && ok;j++) {

			VLIST *v = f->GetVariableAt(j);
			if( (idx=getVariable(v->name,"A"))>0 ) {
				ok = (idx<=nbFacet);
				if( ok ) v->value = (double)geom->GetFacet(idx-1)->sh.counter.hit.nbAbsorbed;
			} else if( (idx=getVariable(v->name,"D"))>0 ) {
				ok = (idx<=nbFacet);
				if( ok ) v->value = (double)geom->GetFacet(idx-1)->sh.counter.hit.nbDesorbed;
			} else if( (idx=getVariable(v->name,"H"))>0 ) {
				ok = (idx<=nbFacet);
				if( ok ) v->value = (double)geom->GetFacet(idx-1)->sh.counter.hit.nbHit;
			} else if( (idx=getVariable(v->name,"_A"))>0 ) {
				ok = (idx<=nbFacet);
				if( ok ) v->value = (double)geom->GetFacet(idx-1)->sh.counter.density.absorbed;
			} else if( (idx=getVariable(v->name,"_D"))>0 ) {
				ok = (idx<=nbFacet);
				if( ok ) v->value = (double)geom->GetFacet(idx-1)->sh.counter.density.desorbed;
			} else if( (idx=getVariable(v->name,"_H"))>0 ) {
				ok = (idx<=nbFacet);
				if( ok ) v->value = (double)geom->GetFacet(idx-1)->sh.counter.density.value;
			} else if( (idx=getVariable(v->name,"AR"))>0 ) {
				ok = (idx<=nbFacet);
				if( ok ) v->value = geom->GetFacet(idx-1)->sh.area;
			} else if( _stricmp(v->name,"SUMDES")==0 ) {
				v->value = (double)worker.nbDesorption;
			} else if( _stricmp(v->name,"SUMABS")==0 ) {
				v->value = (double)worker.nbAbsorption;
			} else if( _stricmp(v->name,"SUMHIT")==0 ) {
				v->value = (double)worker.nbHit;
			} else if( _stricmp(v->name,"MPP")==0 ) {
				v->value = worker.distTraveledTotal/(double)worker.nbDesorption;
			} else if( _stricmp(v->name,"MFP")==0 ) {
				v->value = worker.distTraveledTotal/(double)worker.nbHit;
			}else if( _stricmp(v->name,"DESAR")==0 ) {
				double sumArea = 0.0;
				for(int i2=0;i2<geom->GetNbFacet();i2++) {
					Facet *f_tmp = geom->GetFacet(i2);
					if(f_tmp->sh.desorbType) sumArea += f_tmp->sh.area*(f_tmp->sh.is2sided?2.0:1.0);
				}
				v->value = sumArea;
			} else if( _stricmp(v->name,"ABSAR")==0 ) {
				double sumArea = 0.0;

				for(int i2=0;i2<geom->GetNbFacet();i2++) {
					Facet *f_tmp = geom->GetFacet(i2);
					if(f_tmp->sh.sticking>0.0) sumArea += f_tmp->sh.area*f_tmp->sh.opacity*(f_tmp->sh.is2sided?2.0:1.0);
				}
				v->value = sumArea;
			} else if( _stricmp(v->name,"QCONST")==0 ) {
				v->value = worker.finalOutgassingRate_Pa_m3_sec*10.00; //10: Pa*m3/sec -> mbar*l/s
			} else if (_stricmp(v->name, "QCONST_N") == 0) {
				v->value = worker.finalOutgassingRate;
			} else if (_stricmp(v->name, "NTOT") == 0) {
				v->value = worker.totalDesorbedMolecules;
			}else if( _stricmp(v->name,"KB")==0 ) {
				v->value = 1.3806504e-23;
			} else if( _stricmp(v->name,"R")==0 ) {
				v->value = 8.314472;
			} else if( _stricmp(v->name,"Na")==0 ) {
				v->value = 6.02214179e23;
			} else {
				formulas[i].value->SetText("Invalid var name");
				ok = FALSE;
				iName = TRUE;
			}

			if(!ok && !iName) formulas[i].value->SetText("Invalid var index");

		}

		// Evaluation
		if( ok ) {
			double r;
			if( f->Evaluate(&r) ) {
				sprintf(tmp,"%g",r);
				formulas[i].value->SetText(tmp);
			} else {
				formulas[i].value->SetText(f->GetErrorMsg());
			}
		}

	}

}

void MolFlow::OffsetFormula(char *expression,int offset,int filter) {
	//will increase or decrease facet numbers in a formula
	//only applies to facet numbers larger than "filter" parameter

	vector<string> prefixes;
	prefixes.push_back("A");
	prefixes.push_back("D");
	prefixes.push_back("H");
	prefixes.push_back("AR");
	prefixes.push_back("a");
	prefixes.push_back("d");
	prefixes.push_back("h");
	prefixes.push_back("ar");
	prefixes.push_back(","); //for sum formulas

	string expr=expression; //convert char array to string

	size_t pos=0; //analyzed until this position
	while (pos<expr.size()) { //while not end of expression

		vector<size_t> location; //for each prefix, we store where it was found

		for (int j=0;j<(int)prefixes.size();j++) { //try all expressions
			location.push_back(expr.find(prefixes[j],pos));
		}
		size_t minPos=string::npos;
		size_t maxLength=0;
		for (int j=0;j<(int)prefixes.size();j++)  //try all expressions, find first prefix location
			if (location[j]<minPos) minPos=location[j];
		for (int j=0;j<(int)prefixes.size();j++)  //try all expressions, find longest prefix at location
			if (location[j]==minPos && prefixes[j].size()>maxLength) maxLength=prefixes[j].size();
		int digitsLength=0;
		if (minPos!=string::npos) { //found expression, let's find tailing facet number digits
			while ((minPos+maxLength+digitsLength)<expr.length() && expr[minPos+maxLength+digitsLength]>='0' && expr[minPos+maxLength+digitsLength]<='9')
				digitsLength++;
			if (digitsLength>0) { //there was a digit after the prefix
				int facetNumber;
				if (sscanf(expr.substr(minPos+maxLength,digitsLength).c_str(),"%d",&facetNumber)){
					if ((facetNumber - 1)>filter) {
						char tmp[10];
						sprintf(tmp, "%d", facetNumber += offset);
						expr.replace(minPos + maxLength, digitsLength, tmp);
					}
					else if ((facetNumber - 1) == filter) {
						expr.replace(minPos + maxLength, digitsLength, "0");
					}
				}
			}
		}
		if (minPos!=string::npos) pos=minPos+maxLength+digitsLength;
		else pos=minPos;
	}
	strcpy(expression,expr.c_str());
}

void MolFlow::RenumberFormulas(int startId) {
	for (int i=0;i<nbFormula;i++) {
		char expression[1024];
		strcpy(expression,this->formulas[i].parser->GetExpression());
		OffsetFormula(expression,-1,startId);
		this->formulas[i].parser->SetExpression(expression);
		this->formulas[i].parser->Parse();
	}
}




void InitLog(char *tmp) {
	FILE *f = fopen("D:\\c++\\molflow\\tests\\data\\scan.txt","w");
	fprintf(f,tmp);
	fprintf(f,"\n");
	fclose(f);
}

void Log(char *tmp) {
	FILE *f = fopen("D:\\c++\\molflow\\tests\\data\\scan.txt","a");
	fprintf(f,tmp);
	fprintf(f,"\n");
	fclose(f);
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
	if (autoSaveSimuOnly) lastSaveTimeSimu=worker.simuTime+(m_fTime-worker.startTime);
	else lastSaveTime=m_fTime;
}

//-----------------------------------------------------------------------------
BOOL MolFlow::AutoSave(BOOL crashSave) {
	if (!changedSinceSave) return TRUE;
	GLProgress *progressDlg2 = new GLProgress("Peforming autosave...","Please wait");
	progressDlg2->SetProgress(0.0);
	progressDlg2->SetVisible(TRUE);
	//GLWindowManager::Repaint();
	char CWD [MAX_PATH];
	_getcwd( CWD, MAX_PATH );
	char filename[1024];
	sprintf(filename,"%s\\Molflow_AutoSave.zip",CWD);
	try {
		worker.SaveGeometry(filename,progressDlg2,FALSE,FALSE,TRUE,crashSave);
		ResetAutoSaveTimer(); //deduct saving time from interval
	} catch (Error &e) {
		char errMsg[512];
		sprintf(errMsg,"%s\nFile:%s",e.GetMsg(),worker.GetFileName());
		GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
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
	char CWD [MAX_PATH];
	_getcwd( CWD, MAX_PATH );
	char filename[1024];
	sprintf(filename,"%s\\Molflow_AutoSave.zip",CWD);
	if (FileUtils::Exist(filename)) {
		int rep = GLMessageBox::Display("Autosave file found. Load it now?\nIf you click CANCEL the file will be discarded.","Autosave recovery",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONWARNING);
		if( rep == GLDLG_OK ) {
			LoadFile(filename);
			RemoveRecent(filename);
		}
		return;
	}
	sprintf(filename,"%s\\Molflow_AutoSave.xml",CWD);
	if (FileUtils::Exist(filename)) {
		int rep = GLMessageBox::Display("Autosave file found. Load it now?\nIf you click CANCEL the file will be discarded.","Autosave recovery",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONWARNING);
		if( rep == GLDLG_OK ) {
			LoadFile(filename);
		}
	}
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
	if (geom->IsLoaded()) {
		if (autoSaveSimuOnly) {
			if (worker.running) {
				if(((worker.simuTime+(m_fTime-worker.startTime))-lastSaveTimeSimu)>=(float)autoSaveFrequency*60.0f) {
					AutoSave(); 
				}
			}
		} else {
			if((m_fTime-lastSaveTime)>=(float)autoSaveFrequency*60.0f) {
				AutoSave();
			}
		}
	}

	// Simulation monitoring
	if(globalSettings) globalSettings->SMPUpdate(m_fTime);
	if(profilePlotter) profilePlotter->Update(m_fTime);

	if(pressureEvolution) pressureEvolution->Update(m_fTime);
	if(timewisePlotter) timewisePlotter->Update(m_fTime);
	if(texturePlotter) texturePlotter->Update(m_fTime);  
	//if(facetDetails) facetDetails->Update();
	if(worker.running) {
		if( m_fTime-lastUpdate>=1.0f ) {

			// Update hits
			try {
				worker.Update(m_fTime);
			} catch(Error &e) {
				GLMessageBox::Display((char *)e.GetMsg(),"Error (Stop)",GLDLG_OK,GLDLG_ICONERROR);
			}
			if( textureSettings ) textureSettings->Update();
			lastUpdate = m_fTime;



			// Update timing measurements
			if( worker.nbHit!=lastNbHit || worker.nbDesorption!=lastNbDes ) {
				double dTime = (double)(m_fTime-lastMeasTime);
				hps = (double)(worker.nbHit-lastNbHit)/dTime;
				dps = (double)(worker.nbDesorption-lastNbDes)/dTime;
				if( lastHps!=0.0 ) {
					hps = 0.2*(hps) + 0.8*lastHps;
					dps = 0.2*(dps) + 0.8*lastDps;
				}
				lastHps = hps;
				lastDps = dps;
				lastNbHit = worker.nbHit;
				lastNbDes = worker.nbDesorption;
				lastMeasTime = m_fTime;
			}

		}
		if( worker.calcAC ) {
			sprintf(tmp,"Calc AC: %s (%d %%)",FormatTime(worker.simuTime+(m_fTime-worker.startTime)),
				worker.calcACprg);
		} else {
			sprintf(tmp,"Running: %s",FormatTime(worker.simuTime+(m_fTime-worker.startTime)));
		}
		sTime->SetText( tmp );
	} else {
		if( worker.simuTime>0.0 ) {
			hps = (double)(worker.nbHit-nbHitStart)/worker.simuTime;
			dps = (double)(worker.nbDesorption-nbDesStart)/worker.simuTime;
		} else {
			hps = 0.0;
			dps = 0.0;
		}
		sprintf(tmp,"Stopped: %s",FormatTime(worker.simuTime));
		sTime->SetText( tmp );
	}

	if( (m_fTime-worker.startTime <= 2.0f) && worker.running ) {
		hitNumber->SetText("Starting...");
		desNumber->SetText("Starting...");

	} else {

		if(worker.mode!=AC_MODE) {
			sprintf(tmp,"%s (%s)",FormatInt(worker.nbHit,"hit"),FormatPS(hps,"hit"));
			hitNumber->SetText(tmp);
		} else {
			hitNumber->SetText("");
		}
		sprintf(tmp,"%s (%s)",FormatInt(worker.nbDesorption,"des"),FormatPS(dps,"des"));
		desNumber->SetText(tmp);
	}


	if( worker.nbLeakTotal ) {
		sprintf(tmp,"%g (%.4f%%)",(double)worker.nbLeakTotal,(double)(worker.nbLeakTotal)*100.0/(double)worker.nbDesorption);
		leakNumber->SetText(tmp);
	} else {
		leakNumber->SetText("None");
	}

	resetSimu->SetEnabled(!worker.running&&worker.nbDesorption>0);

	if (worker.running) {
		startSimu->SetText("Pause");
	} else if (worker.nbHit>0) {
		startSimu->SetText("Resume");
	} else {
		startSimu->SetText("Begin");
	}

	// Facet parameters and hits
	if( viewer[0]->SelectionChanged() ||
		viewer[1]->SelectionChanged() ||
		viewer[2]->SelectionChanged() ||
		viewer[3]->SelectionChanged() ) {
			UpdateFacetParams(TRUE);
	}

	UpdateFacetHits();

	// Formulas
	if (autoUpdateFormulas) UpdateFormula();

	/*
	if(worker.running) {
	if( m_fTime - lastWrite > 1.0 ) {
	llong A = geom->GetFacet(73)->sh.counter.hit.nbAbsorbed;
	llong D = geom->GetFacet(72)->sh.counter.hit.nbDesorbed;
	double t = (double)A / (double)D;
	//double t = geom->GetFacet(73)->sh.counter.density.absorbed;
	sprintf(tmp,"%f %.16f",m_fTime-worker.startTime,t);
	Log(tmp);
	lastWrite=m_fTime;
	}
	}
	*/

	/*
	if(worker.running) {
	if( worker.nbDesorption != lastNbD ) {
	double D = geom->GetFacet(1)->sh.counter.density.absorbed;
	sprintf(tmp,"%I64d %.16f",worker.nbDesorption,D);
	Log(tmp);
	lastNbD=worker.nbDesorption;
	}
	}
	*/

	/*
	if(!worker.running) {
	if( worker.simuTime>0 ) {
	if(nbSt<=10) {
	LogProfile();
	nbSt++;
	StartStopSimulation();
	}
	}
	}
	*/

	// Sleep a bit to avoid unwanted CPU load
	if( viewer[0]->IsDragging() || 
		viewer[1]->IsDragging() || 
		viewer[2]->IsDragging() || 
		viewer[3]->IsDragging() || !worker.running )
		SDL_Delay(32);
	else
		SDL_Delay(60);

	return GL_OK;
}

// ----------------------------------------------------------------

void MolFlow::UpdateFacetHits() {
	char tmp[256];
	Geometry *geom = worker.GetGeometry();

	try{
		// Facet list
		if(geom->IsLoaded()) {

			int sR,eR;
			facetList->GetVisibleRows(&sR,&eR);


			for(int i=sR;i<=eR;i++) {
				int facetId=facetList->GetValueInt(i,0)-1;
				if (facetId==-2) facetId=i;
				if (i>=geom->GetNbFacet()) {
					char errMsg[512];
					sprintf(errMsg,"Molflow::UpdateFacetHits()\nError while updating facet hits. Was looking for facet #%d in list.\nMolflow will now autosave and crash.",i+1);
					GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
					AutoSave();
				}
				Facet *f = geom->GetFacet(facetId);
				sprintf(tmp,"%d",facetId+1);
				facetList->SetValueAt(0,i,tmp);
				switch(modeCombo->GetSelectedIndex()) {
				case MC_MODE:
					facetList->SetColumnLabel(1,"Hits");
					sprintf(tmp,"%I64d",f->sh.counter.hit.nbHit);
					facetList->SetValueAt(1,i,tmp);
					sprintf(tmp,"%I64d",f->sh.counter.hit.nbDesorbed);
					facetList->SetValueAt(2,i,tmp);
					sprintf(tmp,"%I64d",f->sh.counter.hit.nbAbsorbed);
					facetList->SetValueAt(3,i,tmp);
					break;
				case AC_MODE:
					facetList->SetColumnLabel(1,"Density");
					sprintf(tmp,"%g",f->sh.counter.density.value);
					facetList->SetValueAt(1,i,tmp);

					sprintf(tmp,"%g",f->sh.counter.density.desorbed);
					facetList->SetValueAt(2,i,tmp);
					sprintf(tmp,"%g",f->sh.counter.density.absorbed);
					facetList->SetValueAt(3,i,tmp);
					break;
				}
			}

		}
	} 	 catch (Error &e) {
		char errMsg[512];
		sprintf(errMsg,"%s\nError while updating facet hits",e.GetMsg());
		GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
	}

}

//-----------------------------------------------------------------------------
// Name: FormatInt()
// Desc: Format an integer in K,M,G,..
//-----------------------------------------------------------------------------
char *MolFlow::FormatInt(llong v,char *unit)
{

	double x = (double)v;

	static char ret[64];
	if(x<1E3) {
		sprintf(ret,"%g %s",(double)x,unit);
	} else if(x<1E6) {
		sprintf(ret,"%.1f K%s",x/1E3,unit);
	} else if(x<1E9) {
		sprintf(ret,"%.2f M%s",x/1E6,unit);
	} else if(x<1E12) {
		sprintf(ret,"%.2f G%s",x/1E9,unit);
	} else {
		sprintf(ret,"%.2f T%s",x/1E12,unit);
	}

	return ret;

}

//-----------------------------------------------------------------------------
// Name: FormatTime()
// Desc: Format time in HH:MM:SS
//-----------------------------------------------------------------------------
char *MolFlow::FormatTime(float t) {
	static char ret[64];
	int nbSec = (int)(t+0.5f);
	sprintf(ret,"%02d:%02d:%02d",nbSec/3600,(nbSec%3600)/60,nbSec%60);
	return ret;
}

//-----------------------------------------------------------------------------
// Name: FormatPS()
// Desc: Format a double in K,M,G,.. per sec
//-----------------------------------------------------------------------------
char *MolFlow::FormatPS(double v,char *unit)
{

	static char ret[64];
	if(v<1000.0) {
		sprintf(ret,"%.1f %s/s",v,unit);
	} else if(v<1000000.0) {
		sprintf(ret,"%.1f K%s/s",v/1000.0,unit);
	} else if(v<1000000000.0) {
		sprintf(ret,"%.1f M%s/s",v/1000000.0,unit);
	} else {
		sprintf(ret,"%.1f G%s/s",v/1000000000.0,unit);
	}

	return ret;

}

//-----------------------------------------------------------------------------
// Name: FormatSize()
// Desc: Format a double in K,M,G,.. per sec
//-----------------------------------------------------------------------------
char *MolFlow::FormatSize(DWORD size)
{

	static char ret[64];
	if(size<1024UL) {
		sprintf(ret,"%d Bytes",size);
	} else if(size<1048576UL) {
		sprintf(ret,"%.1f KB",(double)size/1024.0);
	} else if(size<1073741824UL) {
		sprintf(ret,"%.1f MB",(double)size/1048576.0);
	} else {
		sprintf(ret,"%.1f GB",(double)size/1073741824.0);
	}

	return ret;

}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Name: RestoreDeviceObjects()
// Desc: Initialize scene objects.
//-----------------------------------------------------------------------------
int MolFlow::RestoreDeviceObjects()
{
	Geometry *geom = worker.GetGeometry();
	geom->RestoreDeviceObjects();
	//worker.Update(0.0f);

	// Restore dialog which are not displayed
	// Those which are displayed are invalidated by the window manager
	RVALIDATE_DLG(formulaSettings);
	RVALIDATE_DLG(collapseSettings);
	RVALIDATE_DLG(importDesorption);
	RVALIDATE_DLG(moveVertex);
	RVALIDATE_DLG(timeSettings);
	RVALIDATE_DLG(scaleVertex);
	RVALIDATE_DLG(scaleFacet);
	RVALIDATE_DLG(selectDialog);
	RVALIDATE_DLG(moveFacet);

	RVALIDATE_DLG(mirrorFacet);
	RVALIDATE_DLG(rotateFacet);
	RVALIDATE_DLG(movement);
	RVALIDATE_DLG(alignFacet);
	RVALIDATE_DLG(addVertex);
	RVALIDATE_DLG(facetMesh);
	RVALIDATE_DLG(facetDetails);

	RVALIDATE_DLG(viewer3DSettings);
	RVALIDATE_DLG(textureSettings);
	RVALIDATE_DLG(globalSettings);
	RVALIDATE_DLG(facetCoordinates);
	RVALIDATE_DLG(vertexCoordinates);
	RVALIDATE_DLG(profilePlotter);
	RVALIDATE_DLG(pressureEvolution);
	RVALIDATE_DLG(timewisePlotter);
	RVALIDATE_DLG(viewEditor);
	RVALIDATE_DLG(texturePlotter);
	RVALIDATE_DLG(outgassingMap);
	RVALIDATE_DLG(parameterEditor);

	UpdateTitle();
	return GL_OK;
}

//-----------------------------------------------------------------------------
// Name: InvalidateDeviceObjects()
// Desc: Free all alocated resource
//-----------------------------------------------------------------------------

int MolFlow::InvalidateDeviceObjects()
{
	Geometry *geom = worker.GetGeometry();
	geom->InvalidateDeviceObjects();

	// Invalidate dialog which are not displayed
	// Those which are displayed are invalidated by the window manager
	IVALIDATE_DLG(formulaSettings);
	IVALIDATE_DLG(collapseSettings);
	IVALIDATE_DLG(importDesorption);
	IVALIDATE_DLG(moveVertex);
	IVALIDATE_DLG(timeSettings);
	IVALIDATE_DLG(scaleVertex);
	IVALIDATE_DLG(scaleFacet);
	IVALIDATE_DLG(selectDialog);
	IVALIDATE_DLG(moveFacet);

	IVALIDATE_DLG(mirrorFacet);
	IVALIDATE_DLG(rotateFacet);
	IVALIDATE_DLG(alignFacet);
	IVALIDATE_DLG(addVertex);
	IVALIDATE_DLG(facetDetails);

	IVALIDATE_DLG(viewer3DSettings);
	IVALIDATE_DLG(textureSettings);
	IVALIDATE_DLG(globalSettings);
	IVALIDATE_DLG(facetCoordinates);
	IVALIDATE_DLG(vertexCoordinates);
	IVALIDATE_DLG(profilePlotter);
	IVALIDATE_DLG(pressureEvolution);
	IVALIDATE_DLG(timewisePlotter);
	IVALIDATE_DLG(viewEditor);
	IVALIDATE_DLG(texturePlotter);
	IVALIDATE_DLG(outgassingMap);
	IVALIDATE_DLG(parameterEditor);

	return GL_OK;
}

//-----------------------------------------------------------------------------
// Name: InvalidateDeviceObjects()
// Desc: Called before exiting
//-----------------------------------------------------------------------------
int MolFlow::OnExit() {
	SaveConfig();
	worker.Exit();
	remove("Molflow_AutoSave.xml");
	remove("Molflow_AutoSave.zip");
	//empty TMP directory
	char tmp[1024];
	char CWD [MAX_PATH];
	_getcwd( CWD, MAX_PATH );
	sprintf(tmp,"del /Q \"%s\\tmp\\*.*\"",CWD);
	system(tmp);

	return GL_OK;
}

//-----------------------------------------------------------------------------

void MolFlow::UpdateCurrentDir(char *fileName) {

	strncpy(currentDir,fileName,1024);
	char *dp = strrchr(currentDir,'\\');
	if(!dp) dp = strrchr(currentDir,'/');
	if(dp) *dp=0;

}

//-----------------------------------------------------------------------------

void MolFlow::UpdateCurrentSelDir(char *fileName) {

	strncpy(currentSelDir,fileName,1024);
	char *dp = strrchr(currentSelDir,'\\');
	if(!dp) dp = strrchr(currentSelDir,'/');
	if(dp) *dp=0;

}

//-----------------------------------------------------------------------------

void MolFlow::UpdateTitle() {

	static char title[128];

	Geometry *geom = worker.GetGeometry();

	if( !geom->IsLoaded() ) {
		sprintf(title,"%s",APP_NAME);
	} else {
		if( geom->viewStruct<0 ) {
			sprintf(title,"%s [%s]",APP_NAME,worker.GetShortFileName());
		} else {
			sprintf(title,"%s [%s: Struct #%d %s]",APP_NAME,worker.GetShortFileName(),geom->viewStruct+1,geom->GetStructureName(geom->viewStruct));
		}
	}

	SetTitle(title);

}

//-----------------------------------------------------------------------------

void MolFlow::SaveFileAs() {

	FILENAME *fn = GLFileBox::SaveFile(currentDir,worker.GetShortFileName(),"Save File",fileSFilters,nbSFilter);

	GLProgress *progressDlg2 = new GLProgress("Saving file...","Please wait");
	progressDlg2->SetProgress(0.0);
	progressDlg2->SetVisible(TRUE);
	//GLWindowManager::Repaint();  
	if( fn ) {

		try {

			worker.SaveGeometry(fn->fullName,progressDlg2);
			ResetAutoSaveTimer();
			changedSinceSave=FALSE;
			UpdateCurrentDir(worker.fullFileName);
			UpdateTitle();
			AddRecent(worker.fullFileName);
		} catch (Error &e) {
			char errMsg[512];
			sprintf(errMsg,"%s\nFile:%s",e.GetMsg(),fn->fullName);
			GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
			RemoveRecent(fn->fullName);
		}

	}

	progressDlg2->SetVisible(FALSE);
	SAFE_DELETE(progressDlg2);
}

//-----------------------------------------------------------------------------

void MolFlow::ExportSelection() {

	Geometry *geom = worker.GetGeometry();
	if(geom->GetNbSelected()==0) {
		GLMessageBox::Display("Empty selection","Error",GLDLG_OK,GLDLG_ICONERROR);
		return;
	}

	FILENAME *fn = GLFileBox::SaveFile(currentDir,worker.GetShortFileName(),"Export selection",fileSFilters,nbSFilter);
	GLProgress *progressDlg2 = new GLProgress("Saving file...","Please wait");
	progressDlg2->SetProgress(0.0);
	progressDlg2->SetVisible(TRUE);
	//GLWindowManager::Repaint();
	if( fn ) {

		try {
			worker.SaveGeometry(fn->fullName,progressDlg2,TRUE,TRUE);
			AddRecent(fn->fullName);
			//UpdateCurrentDir(fn->fullName);
			//UpdateTitle();
		} catch (Error &e) {
			char errMsg[512];
			sprintf(errMsg,"%s\nFile:%s",e.GetMsg(),fn->fullName);
			GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
		}

	}

	progressDlg2->SetVisible(FALSE);
	SAFE_DELETE(progressDlg2);
}

//-----------------------------------------------------------------------------

void MolFlow::ExportTextures(int mode) {

	Geometry *geom = worker.GetGeometry();
	if(geom->GetNbSelected()==0) {
		GLMessageBox::Display("Empty selection","Error",GLDLG_OK,GLDLG_ICONERROR);
		return;
	}
	if (!worker.IsDpInitialized()) {
		GLMessageBox::Display("Worker Dataport not initialized yet", "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}

	FILENAME *fn = GLFileBox::SaveFile(currentDir,NULL,"Save File",fileTexFilters,nbTexFilter);

	if( fn ) {

		try {
			worker.ExportTextures(fn->fullName,mode,TRUE,TRUE);
			//UpdateCurrentDir(fn->fullName);
			//UpdateTitle();
		} catch (Error &e) {
			char errMsg[512];
			sprintf(errMsg,"%s\nFile:%s",e.GetMsg(),fn->fullName);
			GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
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

	FILENAME *fn = GLFileBox::SaveFile(currentDir, NULL, "Save File", fileProfFilters, nbProfFilter);

	if (fn) {

		try {
			worker.ExportProfiles(fn->fullName);
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

void MolFlow::ImportDesorption_DES() {

	FILENAME *fn = GLFileBox::OpenFile(currentDir, NULL, "Import desorption File", fileDesFilters, nbDesFilter);

	if (fn) {

		try {
			worker.ImportDesorption_DES(fn->fullName);
			//UpdateCurrentDir(fn->fullName);
			//UpdateTitle();
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
	if(strlen(worker.fullFileName)>0){

		GLProgress *progressDlg2 = new GLProgress("Saving...","Please wait");
		progressDlg2->SetProgress(0.5);
		progressDlg2->SetVisible(TRUE);
		//GLWindowManager::Repaint();

		try {
			worker.SaveGeometry(worker.fullFileName,progressDlg2,FALSE);
			ResetAutoSaveTimer();
		} catch (Error &e) {
			char errMsg[512];
			sprintf(errMsg,"%s\nFile:%s",e.GetMsg(),worker.GetFileName());
			GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
		}
		progressDlg2->SetVisible(FALSE);
		SAFE_DELETE(progressDlg2);
		changedSinceSave=FALSE;

	} else SaveFileAs();
}

//-----------------------------------------------------------------------------

void MolFlow::SaveSelection() {

	FileWriter *f = NULL;
	Geometry *geom = worker.GetGeometry();
	if(geom->GetNbSelected()==0) return;
	GLProgress *progressDlg2 = new GLProgress("Saving file","Please wait");
	progressDlg2->SetProgress(0.5);
	progressDlg2->SetVisible(TRUE);
	//GLWindowManager::Repaint();

	FILENAME *fn = GLFileBox::SaveFile(currentSelDir,worker.GetShortFileName(),"Save selection",fileSelFilters,nbSelFilter);

	if( fn ) {

		try {


			char *ext = fn->fullName+strlen(fn->fullName)-4;

			if(!(*ext=='.') ) {
				sprintf(fn->fullName,"%s.sel",fn->fullName); //set to default SEL format
				ext = strrchr(fn->fullName,'.');
			}
			ext++;

			f = new FileWriter(fn->fullName);
			int nbSelected = geom->GetNbSelected();
			int nbFacet = geom->GetNbFacet();
			for(int i=0;i<nbFacet;i++) {
				if(geom->GetFacet(i)->selected) f->WriteInt(i,"\n");
			}

		} catch (Error &e) {
			char errMsg[512];
			sprintf(errMsg,"%s\nFile:%s",e.GetMsg(),fn->fullName);
			GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
		}

		SAFE_DELETE(f);

	}
	progressDlg2->SetVisible(FALSE);
	SAFE_DELETE(progressDlg2);
	changedSinceSave=FALSE;
}

//-----------------------------------------------------------------------------

void MolFlow::LoadSelection(char *fName) {

	char fullName[1024];
	strcpy(fullName,"");
	FileReader *f=NULL;

	if(fName==NULL) {
		FILENAME *fn = GLFileBox::OpenFile(currentSelDir,NULL,"Load Selection",fileSelFilters,nbSelFilter);
		if( fn )
			strcpy(fullName,fn->fullName);
	} else {
		strcpy(fullName,fName);
	}

	if(strlen(fullName)==0) return;

	try {

		Geometry *geom = worker.GetGeometry();
		geom->Unselect();
		int nbFacet = geom->GetNbFacet();

		f = new FileReader(fullName);
		while(!f->IsEof()) {
			int s = f->ReadInt();
			if( s>=0 && s<nbFacet ) geom->Select(s);
		}
		geom->UpdateSelection();

		UpdateFacetParams(TRUE);
		UpdateCurrentSelDir(fullName);

	} catch (Error &e) {

		char errMsg[512];
		sprintf(errMsg,"%s\nFile:%s",e.GetMsg(),fullName);
		GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);

	}

	SAFE_DELETE(f);
	changedSinceSave=FALSE;

}

//-----------------------------------------------------------------------------

void MolFlow::LoadFile(char *fName) {


	char fullName[1024];
	char shortName[512];
	strcpy(fullName,"");

	if(fName==NULL) {
		FILENAME *fn = GLFileBox::OpenFile(currentDir,NULL,"Open File",fileLFilters,nbLFilter);
		if( fn )
			strcpy(fullName,fn->fullName);
	} else {
		strcpy(fullName,fName);
	}

	GLProgress *progressDlg2 = new GLProgress("Preparing to load file...","Please wait");
	progressDlg2->SetVisible(TRUE);
	progressDlg2->SetProgress(0.0);
	//GLWindowManager::Repaint();

	if(strlen(fullName)==0) {
		progressDlg2->SetVisible(FALSE);
		SAFE_DELETE(progressDlg2);
		return;
	}


	char *lPart = strrchr(fullName,'\\');
	if(lPart) strcpy(shortName,lPart+1);
	else strcpy(shortName,fullName);

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
		compACBtn->SetEnabled(modeCombo->GetSelectedIndex()==1);
		singleACBtn->SetEnabled(modeCombo->GetSelectedIndex()==1);
		//resetSimu->SetEnabled(TRUE);
		ClearFacetParams();
		nbDesStart = worker.nbDesorption;
		nbHitStart = worker.nbHit;
		AddRecent(fullName);
		geom->viewStruct = -1;


		UpdateStructMenu();
		if(profilePlotter) profilePlotter->Reset();
		if (timewisePlotter) timewisePlotter->Refresh();
		if(pressureEvolution) pressureEvolution->Reset();
		if(timewisePlotter) timewisePlotter->Reset();
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
		if(profilePlotter) profilePlotter->Refresh();

		if(pressureEvolution) pressureEvolution->Refresh();
		if (textureSettings) textureSettings->Update();
		if(texturePlotter) texturePlotter->Update(m_fTime,TRUE);
		if(outgassingMap) outgassingMap->Update(m_fTime,TRUE);
		if(facetDetails) facetDetails->Update();
		if(facetCoordinates) facetCoordinates->UpdateFromSelection();
		if(vertexCoordinates) vertexCoordinates->Update();
		if(movement) movement->Update();

	} catch (Error &e) {

		char errMsg[512];
		sprintf(errMsg,"%s\nFile:%s",e.GetMsg(),shortName);
		GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
		RemoveRecent(fName);

	}
	progressDlg2->SetVisible(FALSE);
	SAFE_DELETE(progressDlg2);
	changedSinceSave=FALSE;
}

//-----------------------------------------------------------------------------

void MolFlow::InsertGeometry(BOOL newStr,char *fName) {
	if (!AskToReset()) return;
	ResetSimulation(FALSE);

	char fullName[1024];
	char shortName[512];
	strcpy(fullName,"");

	if(fName==NULL) {
		FILENAME *fn = GLFileBox::OpenFile(currentDir,NULL,"Open File",fileInsFilters,nbInsFilter);
		if( fn )
			strcpy(fullName,fn->fullName);
	} else {
		strcpy(fullName,fName);
	}

	GLProgress *progressDlg2 = new GLProgress("Loading file...","Please wait");
	progressDlg2->SetVisible(TRUE);
	progressDlg2->SetProgress(0.0);
	//GLWindowManager::Repaint();

	if(strlen(fullName)==0) {
		progressDlg2->SetVisible(FALSE);
		SAFE_DELETE(progressDlg2);
		return;
	}


	char *lPart = strrchr(fullName,'\\');
	if(lPart) strcpy(shortName,lPart+1);
	else strcpy(shortName,fullName);

	try {

		worker.InsertGeometry(newStr,fullName);


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

		compACBtn->SetEnabled(modeCombo->GetSelectedIndex()==1);
		singleACBtn->SetEnabled(modeCombo->GetSelectedIndex()==1);
		//resetSimu->SetEnabled(TRUE);
		//ClearFacetParams();
		//nbDesStart = worker.nbDesorption;
		//nbHitStart = worker.nbHit;
		AddRecent(fullName);
		geom->viewStruct = -1;


		//worker.LoadTextures(fullName);
		UpdateStructMenu();
		if(profilePlotter) profilePlotter->Reset();

		if(pressureEvolution) pressureEvolution->Reset();
		if(timewisePlotter) timewisePlotter->Reset();
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
		if(profilePlotter) profilePlotter->Refresh();

		if(pressureEvolution) pressureEvolution->Refresh();
		if(timewisePlotter) timewisePlotter->Refresh();
		if(texturePlotter) texturePlotter->Update(m_fTime,TRUE);
		if(outgassingMap) outgassingMap->Update(m_fTime,TRUE);
		if(facetDetails) facetDetails->Update();
		if(facetCoordinates) facetCoordinates->UpdateFromSelection();
		if(vertexCoordinates) vertexCoordinates->Update();


	} catch (Error &e) {

		char errMsg[512];
		sprintf(errMsg,"%s\nFile:%s",e.GetMsg(),shortName);
		GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
		RemoveRecent(fName);

	}
	progressDlg2->SetVisible(FALSE);
	SAFE_DELETE(progressDlg2);
	changedSinceSave=TRUE;
}

// ----------------------------------------------------------------------------
void MolFlow::UpdateMeasurements() {
	char tmp[256];
	sprintf(tmp,"%s (%s)",FormatInt(worker.nbHit,"hit"),FormatPS(hps,"hit"));
	hitNumber->SetText(tmp);
	sprintf(tmp,"%s (%s)",FormatInt(worker.nbDesorption,"des"),FormatPS(dps,"des"));
	desNumber->SetText(tmp);
}

//-----------------------------------------------------------------------------
// Name: UpdateModelParams()
// Desc: Update displayed model parameter on geometry ghange
//-----------------------------------------------------------------------------
void MolFlow::UpdateModelParams() {

	Geometry *geom = worker.GetGeometry();
	char tmp[256];
	double sumArea=0;
	facetList->SetSize(4,geom->GetNbFacet(),TRUE);
	facetList->SetColumnWidths((int*)cWidth);
	facetList->SetColumnLabels((char **)cName);

	facetList->UpdateAllRows();
	AABB bb = geom->GetBB();

	for (int i=0;i<geom->GetNbFacet();i++) {
		Facet *f=geom->GetFacet(i);
		if (f->sh.area>0) sumArea+=f->sh.area*(f->sh.is2sided?2.0:1.0);
	}

	sprintf(tmp,"V:%d F:%d Dim:(%g,%g,%g) Area:%g",geom->GetNbVertex(),geom->GetNbFacet(),
		(bb.max.x-bb.min.x),(bb.max.y-bb.min.y),(bb.max.z-bb.min.z),sumArea);
	geomNumber->SetText(tmp);

}

//-----------------------------------------------------------------------------
// Name: AddFormula()
// Desc: Add a formula
//-----------------------------------------------------------------------------
void MolFlow::AddFormula(GLParser *f,BOOL doUpdate) {

	if( f ) {
		if( nbFormula<MAX_FORMULA ) {
			formulas[nbFormula].parser = f;
			formulas[nbFormula].name = new GLLabel(f->GetName());
			Add(formulas[nbFormula].name);
			formulas[nbFormula].value = new GLTextField(0,"");
			formulas[nbFormula].value->SetEditable(FALSE);
			Add(formulas[nbFormula].value);
			formulas[nbFormula].setBtn = new GLButton(0,"...");
			Add(formulas[nbFormula].setBtn);
			nbFormula++;
			PlaceComponents();
			if(doUpdate) UpdateFormula();
		} else {
			SAFE_DELETE(f);
		}
	}

}

void MolFlow::ClearFormula() {

	for(int i=0;i<nbFormula;i++) {
		wnd->PostDelete(formulas[i].name);
		wnd->PostDelete(formulas[i].value);
		wnd->PostDelete(formulas[i].setBtn);
		formulas[i].name=NULL;
		formulas[i].value=NULL;
		formulas[i].setBtn=NULL;
		SAFE_DELETE(formulas[i].parser);
	}
	nbFormula=0;
	PlaceComponents();

}

void MolFlow::AddFormula(const char *fName,const char *formula) {

	GLParser *f = new GLParser();
	f->SetExpression(formula);
	f->SetName(fName);
	f->Parse();
	AddFormula(f,FALSE);

}

void MolFlow::ClearParameters() {
	worker.parameters = std::vector<Parameter>();
	if (parameterEditor) parameterEditor->UpdateCombo();
}
//-----------------------------------------------------------------------------
// Name: ProcessFormulaButtons()
// Desc: Handle forumla button event
//-----------------------------------------------------------------------------
void MolFlow::ProcessFormulaButtons(GLComponent *src) {

	// Search formula buttons
	BOOL found = FALSE;
	int i=0;
	while(!found && i<nbFormula) {
		found = (src == formulas[i].setBtn);
		if(!found) i++;
	}
	if( found ) {
		if( !formulaSettings ) formulaSettings = new FormulaSettings();
		if( formulaSettings->EditFormula(formulas[i].parser) ) {
			// Apply change
			formulas[i].name->SetText(formulas[i].parser->GetName());
			UpdateFormula();
		} else {
			// Delete
			wnd->PostDelete(formulas[i].name);
			wnd->PostDelete(formulas[i].value);
			wnd->PostDelete(formulas[i].setBtn);
			formulas[i].name=NULL;
			formulas[i].value=NULL;
			formulas[i].setBtn=NULL;
			SAFE_DELETE(formulas[i].parser);
			for(int j=i;j<nbFormula-1;j++)
				formulas[j] = formulas[j+1];
			nbFormula--;
			PlaceComponents();
			UpdateFormula();
		}
	}

}

//-----------------------------------------------------------------------------

void MolFlow::StartStopSimulation() {

	//if(nbSt<=10) BuildPipeStick((double)nbSt/10);
	//else         return;

	worker.StartStop(m_fTime,modeCombo->GetSelectedIndex());
	if(profilePlotter) profilePlotter->Update(m_fTime,TRUE);
	if(pressureEvolution) pressureEvolution->Update(m_fTime,TRUE);
	if(timewisePlotter) timewisePlotter->Update(m_fTime,TRUE);
	if(texturePlotter) texturePlotter->Update(m_fTime,TRUE);

	// Frame rate measurement
	lastMeasTime = m_fTime;
	dps = 0.0;
	hps = 0.0;
	lastHps = hps;
	lastDps = dps;
	lastNbHit = worker.nbHit;
	lastNbDes = worker.nbDesorption;
	lastUpdate = 0.0;

	/*
	if( worker.running ) {
	InitLog("#convergence_speed");
	Log("0 0");
	lastWrite = 0.0f;
	}
	*/

	/*
	if( worker.running ) {
	InitLog("#Jacobi");
	Log("0 0");
	lastWrite = 0.0f;
	}
	*/

}

//-----------------------------------------------------------------------------

void MolFlow::ResetSimulation(BOOL askConfirm) {

	BOOL ok = TRUE;
	if( askConfirm ) 
		ok = GLMessageBox::Display("Reset simulation ?","Question",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONINFO)==GLDLG_OK;

	if(ok) {
		worker.Reset(m_fTime);
		nbDesStart = 0;
		nbHitStart = 0;
	}
	//resetSimu->SetEnabled(FALSE);
	if(profilePlotter) profilePlotter->Update(m_fTime,TRUE);
	if(pressureEvolution) pressureEvolution->Update(m_fTime,TRUE);
	if(timewisePlotter) timewisePlotter->Update(m_fTime,TRUE);
	if(texturePlotter) texturePlotter->Update(m_fTime,TRUE);

}

//-----------------------------------------------------------------------------

void MolFlow::SelectViewer(int s) {

	curViewer = s;
	for(int i=0;i<MAX_VIEWER;i++) viewer[i]->SetSelected(i==curViewer);
	UpdateViewerParams();

}

//-----------------------------------------------------------------------------
// Name: EventProc()
// Desc: Message proc function to handle key and mouse input
//-----------------------------------------------------------------------------
void MolFlow::ProcessMessage(GLComponent *src,int message)
{
	Geometry *geom = worker.GetGeometry();
	char *input;
	char tmp[128];
	switch(message) {

		//MENU --------------------------------------------------------------------
	case MSG_MENU:
		switch(src->GetId()) {
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
			if( geom->IsLoaded() ) {
				if (worker.running) worker.Stop_Public();
				InsertGeometry(FALSE);
			} else GLMessageBox::Display("No geometry loaded.","No geometry",GLDLG_OK,GLDLG_ICONERROR);
			break;
		case MENU_FILE_INSERTGEO_NEWSTR:
			if( geom->IsLoaded() ) {
				if (worker.running) worker.Stop_Public();
				InsertGeometry(TRUE);
			} else GLMessageBox::Display("No geometry loaded.","No geometry",GLDLG_OK,GLDLG_ICONERROR);
			break;
		case MENU_FILE_SAVEAS:
			if( geom->IsLoaded() ) {
				SaveFileAs();
			} else GLMessageBox::Display("No geometry loaded.","No geometry",GLDLG_OK,GLDLG_ICONERROR);
			break;
		case MENU_FILE_EXPORTMESH:
			ExportSelection();
			break;

		case MENU_FILE_EXPORTTEXTURE_AREA:
		case MENU_FILE_EXPORTTEXTURE_MCHITS:
		case MENU_FILE_EXPORTTEXTURE_IMPINGEMENT:
		case MENU_FILE_EXPORTTEXTURE_PART_DENSITY:
		case MENU_FILE_EXPORTTEXTURE_GAS_DENSITY:
		case MENU_FILE_EXPORTTEXTURE_PRESSURE:
		case MENU_FILE_EXPORTTEXTURE_AVG_V:
		case MENU_FILE_EXPORTTEXTURE_N_VECTORS:
			ExportTextures(src->GetId()-150); // 1..9
			break;

		case MENU_FILE_EXPORTPROFILES:
			ExportProfiles();
			break;

		case MENU_FILE_SAVE:
			if( geom->IsLoaded() ) SaveFile();
			else GLMessageBox::Display("No geometry loaded.","No geometry",GLDLG_OK,GLDLG_ICONERROR);
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
			if( !formulaSettings ) formulaSettings = new FormulaSettings();
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
			if( !globalSettings ) globalSettings = new GlobalSettings();
			globalSettings->Display(&worker);
			break;
		case MENU_FACET_COLLAPSE:
			if( geom->IsLoaded() ) {
				DisplayCollapseDialog();
			} else GLMessageBox::Display("No geometry loaded.","No geometry",GLDLG_OK,GLDLG_ICONERROR);
			break;
		case MENU_FACET_SWAPNORMAL:
			if (AskToReset()) {
				geom->SwapNormal();
				// Send to sub process
				try { worker.Reload(); } catch(Error &e) {
					GLMessageBox::Display((char *)e.GetMsg(),"Error",GLDLG_OK,GLDLG_ICONERROR);
				}
			}
			break;
		case MENU_FACET_SHIFTVERTEX:
			if (AskToReset()) {
				geom->ShiftVertex();
				// Send to sub process
				try { worker.Reload(); } catch(Error &e) {
					GLMessageBox::Display((char *)e.GetMsg(),"Error",GLDLG_OK,GLDLG_ICONERROR);
				}
			}
			break;
		case MENU_FACET_COORDINATES:

			if(!facetCoordinates) facetCoordinates = new FacetCoordinates();
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
			if( geom->IsLoaded() ) {
				if( !scaleFacet ) scaleFacet = new ScaleFacet(geom,&worker);

				scaleFacet->SetVisible(TRUE);

			} else GLMessageBox::Display("No geometry loaded.","No geometry",GLDLG_OK,GLDLG_ICONERROR);
			break;
		case MENU_FACET_MIRROR:
			if(!mirrorFacet) mirrorFacet = new MirrorFacet(geom,&worker);
			mirrorFacet->SetVisible(TRUE);
			break;
		case MENU_FACET_ROTATE:
			if(!rotateFacet) rotateFacet = new RotateFacet(geom,&worker);
			rotateFacet->SetVisible(TRUE);
			break;
		case MENU_FACET_ALIGN:
			if(!alignFacet) alignFacet = new AlignFacet(geom,&worker);
			alignFacet->MemorizeSelection();
			alignFacet->SetVisible(TRUE);
			break;
		case MENU_FACET_PROFPLOTTER:
			if(!profilePlotter) profilePlotter = new ProfilePlotter();
			profilePlotter->Display(&worker);
			break;


		case MENU_TIME_PRESSUREEVOLUTION:
			if(!pressureEvolution) pressureEvolution = new PressureEvolution();
			pressureEvolution->Display(&worker);
			break;
		/*case MENU_FACET_MESH:
			if (!facetMesh) facetMesh = new FacetMesh(&worker);
			facetMesh->SetVisible(!facetMesh->IsVisible());
			break;      */    
		case MENU_FACET_TEXPLOTTER:
			if( !texturePlotter ) texturePlotter = new TexturePlotter();
			texturePlotter->Display(&worker);
			break;
		case MENU_FACET_OUTGASSINGMAP:
			if( !outgassingMap ) outgassingMap = new OutgassingMap();
			outgassingMap->Display(&worker);
			break;
		case MENU_FACET_REMOVESEL:
			if( GLMessageBox::Display("Remove selected facets?","Question",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONINFO)==GLDLG_OK ) {
				if (AskToReset()) {
					if (worker.running) worker.Stop_Public();
					geom->RemoveSelected();
					//geom->CheckIsolatedVertex();
					UpdateModelParams();
					if (vertexCoordinates) vertexCoordinates->Update();
					if (facetCoordinates) facetCoordinates->UpdateFromSelection();
					if (profilePlotter) profilePlotter->Refresh();

					if (pressureEvolution) pressureEvolution->Refresh();
					if (timewisePlotter) timewisePlotter->Refresh();
					// Send to sub process
					try { worker.Reload(); } catch(Error &e) {
						GLMessageBox::Display((char *)e.GetMsg(),"Error reloading worker",GLDLG_OK,GLDLG_ICONERROR);
					}
				}
			}
			break;
		case MENU_FACET_EXPLODE:
			if( GLMessageBox::Display("Explode selected facet?","Question",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONINFO)==GLDLG_OK ) {
				if (AskToReset()) {
					int err;
					try {
						err = geom->ExplodeSelected();
					} catch(Error &e) {
						GLMessageBox::Display((char *)e.GetMsg(),"Error exploding",GLDLG_OK,GLDLG_ICONERROR);
					}
					if( err==-1 ) {
						GLMessageBox::Display("Empty selection","Error",GLDLG_OK,GLDLG_ICONERROR);
					} else if ( err==-2 ) {
						GLMessageBox::Display("All selected facets must have a mesh with boudary correction enabled","Error",GLDLG_OK,GLDLG_ICONERROR);
					} else if( err==0 ) {

						UpdateModelParams();
						UpdateFacetParams(TRUE);
						// Send to sub process
						try { worker.Reload(); } catch(Error &e) {
							GLMessageBox::Display((char *)e.GetMsg(),"Error reloading worker",GLDLG_OK,GLDLG_ICONERROR);
						}
					}
				}
			}
			break;
		case MENU_FACET_DETAILS:
			if( facetDetails==NULL ) facetDetails = new FacetDetails();
			facetDetails->Display(&worker);
			break;

		case MENU_FACET_SELECTALL:
			geom->SelectAll();
			UpdateFacetParams(TRUE);
			break;
		case MENU_FACET_SELECTSTICK:
			geom->Unselect();
			for(int i=0;i<geom->GetNbFacet();i++)
				if( geom->GetFacet(i)->sh.sticking != 0.0 && !geom->GetFacet(i)->IsLinkFacet() )
					geom->Select(i);
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;
		case MENU_FACET_SELECTTRANS:
			geom->Unselect();
			for(int i=0;i<geom->GetNbFacet();i++)
				if( geom->GetFacet(i)->sh.opacity != 1.0 && geom->GetFacet(i)->sh.opacity != 2.0 )
					geom->Select(i);
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;
		case MENU_FACET_SELECTREFL:
			geom->Unselect();
			for(int i=0;i<geom->GetNbFacet();i++) {
				Facet *f = geom->GetFacet(i);
				if( f->sh.desorbType==DES_NONE && f->sh.sticking==0.0 && f->sh.opacity>0.0 )
					geom->Select(i);
			}
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;
		case MENU_FACET_SELECT2SIDE:
			geom->Unselect();
			for(int i=0;i<geom->GetNbFacet();i++)
				if( geom->GetFacet(i)->sh.is2sided )
					geom->Select(i);
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;
		case MENU_FACET_SELECTVOL:
			geom->Unselect();
			for(int i=0;i<geom->GetNbFacet();i++)
				if( geom->GetFacet(i)->sh.isVolatile )
					geom->Select(i);
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;
		case MENU_FACET_SELECTTEXT:
			geom->Unselect();
			for(int i=0;i<geom->GetNbFacet();i++)
				if( geom->GetFacet(i)->sh.isTextured!=NULL )
					geom->Select(i);
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;
		case MENU_FACET_SELECTPROF:
			geom->Unselect();
			for(int i=0;i<geom->GetNbFacet();i++)
				if( geom->GetFacet(i)->sh.isProfile!=NULL )
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
		geom->Unselect();
		for (int i = 0; i < geom->GetNbFacet(); i++)
			if (geom->GetFacet(i)->err >= planarityThreshold)
				geom->Select(i);
		geom->UpdateSelection();
		UpdateFacetParams(TRUE); 
			break;


		case MENU_FACET_SELECTERR:
			geom->Unselect();
			for(int i=0;i<geom->GetNbFacet();i++)

				if( geom->GetFacet(i)->sh.sign==0.0 )
					geom->Select(i);
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;

		case MENU_FACET_SELECTDEST:
			geom->Unselect();
			for(int i=0;i<geom->GetNbFacet();i++)

				if( geom->GetFacet(i)->sh.superDest!=0 )
					geom->Select(i);
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;

		case MENU_FACET_SELECTABS:
			geom->Unselect();
			for(int i=0;i<geom->GetNbFacet();i++)
				if( geom->GetFacet(i)->sh.counter.hit.nbAbsorbed > 0 )
					geom->Select(i);
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;

		case MENU_FACET_SELECTHITS:
			geom->Unselect();
			for(int i=0;i<geom->GetNbFacet();i++)

				if( geom->GetFacet(i)->sh.counter.hit.nbHit > 0 )
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
		geom->Unselect();
		for (int i = 0; i < geom->GetNbFacet(); i++)
			if (geom->GetFacet(i)->sh.counter.hit.nbHit == 0 && geom->GetFacet(i)->sh.area >= largeArea)
				geom->Select(i);
		geom->UpdateSelection();
		UpdateFacetParams(TRUE); 
			break;

		case MENU_FACET_SELECTDES:
			geom->Unselect();
			for(int i=0;i<geom->GetNbFacet();i++)
				if( geom->GetFacet(i)->sh.desorbType!=DES_NONE )
					geom->Select(i);
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;
		case MENU_SELECT_HASDESFILE:
			geom->Unselect();
			for(int i=0;i<geom->GetNbFacet();i++)
				if( geom->GetFacet(i)->hasOutgassingMap )
					geom->Select(i);
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;
		case MENU_FACET_INVERTSEL:
			for(int i=0;i<geom->GetNbFacet();i++)
				geom->GetFacet(i)->selected = !geom->GetFacet(i)->selected;
			geom->UpdateSelection();
			UpdateFacetParams(TRUE);
			break;
		case MENU_SELECTION_SELECTFACETNUMBER:
			if( !selectDialog ) selectDialog = new SelectDialog(&worker);
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
			if( GLMessageBox::Display("Clear all selections ?","Question",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONINFO)==GLDLG_OK ) {
				ClearAllSelections();
			}
			break;
		case MENU_VERTEX_UNSELECTALL:
			geom->UnselectAllVertex();
			break;
		case MENU_VERTEX_SELECTALL:
			geom->SelectAllVertex();
			break;
		case MENU_VERTEX_SELECT_ISOLATED:
			geom->SelectIsolatedVertices();
			break;
		case MENU_VERTEX_CREATE_POLY_CONVEX:
			if (AskToReset()) {
				try {
					geom->CreatePolyFromVertices_Convex();
				} catch (Error &e) {
					GLMessageBox::Display((char *)e.GetMsg(),"Error creating polygon",GLDLG_OK,GLDLG_ICONERROR);
				}
				//UpdateModelParams();
				try {
					worker.Reload();
				} catch(Error &e) {
					GLMessageBox::Display((char *)e.GetMsg(),"Error reloading worker",GLDLG_OK,GLDLG_ICONERROR);
				}
			}
			break;
		case MENU_VERTEX_CREATE_POLY_ORDER:
			if (AskToReset()) {
				try {
					geom->CreatePolyFromVertices_Order();
				} catch (Error &e) {
					GLMessageBox::Display((char *)e.GetMsg(),"Error creating polygon",GLDLG_OK,GLDLG_ICONERROR);
				}
				//UpdateModelParams();
				try {
					worker.Reload();
				} catch(Error &e) {

					GLMessageBox::Display((char *)e.GetMsg(),"Error reloading worker",GLDLG_OK,GLDLG_ICONERROR);
				}
			}
			break;
		case MENU_FACET_CREATE_DIFFERENCE:
			if (geom->IsLoaded()){
				try { 
					if (AskToReset()) {
						geom->CreateDifference(); }} catch (Error &e) {
							GLMessageBox::Display((char *)e.GetMsg(),"Error creating polygon",GLDLG_OK,GLDLG_ICONERROR);
						}
						//UpdateModelParams();
						try { worker.Reload(); } catch(Error &e) {
							GLMessageBox::Display((char *)e.GetMsg(),"Error reloading worker",GLDLG_OK,GLDLG_ICONERROR);
						}
			}else GLMessageBox::Display("No geometry loaded.","No geometry",GLDLG_OK,GLDLG_ICONERROR);
			break;

		case MENU_VERTEX_SELECT_COPLANAR:
			char *input;


			if( geom->IsLoaded() ) {
				if (geom->GetNbSelectedVertex()!=3) {
					GLMessageBox::Display("Select exactly 3 vertices","Can't define plane",GLDLG_OK,GLDLG_ICONERROR);
					return;
				}
				sprintf(tmp,"%g",tolerance);
				//sprintf(title,"Pipe L/R = %g",L/R);
				input = GLInputBox::GetInput(tmp,"Tolerance (cm)","Select coplanar vertices");
				if( !input ) return;	
				if(( sscanf(input,"%lf",&tolerance)<=0 )||(tolerance<=0.0)) {
					GLMessageBox::Display("Invalid number","Error",GLDLG_OK,GLDLG_ICONERROR);
					return;
				}
				try { viewer[curViewer]->SelectCoplanar(tolerance);  } catch (Error &e) {
					GLMessageBox::Display((char *)e.GetMsg(),"Error selecting coplanar vertices",GLDLG_OK,GLDLG_ICONERROR);
				}
			} else GLMessageBox::Display("No geometry loaded.","No geometry",GLDLG_OK,GLDLG_ICONERROR);
			break;
		case MENU_VERTEX_MOVE:
			if( geom->IsLoaded() ) {
				if( !moveVertex ) moveVertex = new MoveVertex(geom,&worker);

				//moveVertex->DoModal();
				moveVertex->SetVisible(TRUE);

				/*
				UpdateModelParams();
				try { worker.Reload(); } catch(Error &e) {
				GLMessageBox::Display((char *)e.GetMsg(),"Error reloading worker",GLDLG_OK,GLDLG_ICONERROR);
				*/

			} else GLMessageBox::Display("No geometry loaded.","No geometry",GLDLG_OK,GLDLG_ICONERROR);
			break;
		case MENU_VERTEX_SCALE:
			if( geom->IsLoaded() ) {
				if( !scaleVertex ) scaleVertex = new ScaleVertex(geom,&worker);

				scaleVertex->SetVisible(TRUE);

			} else GLMessageBox::Display("No geometry loaded.","No geometry",GLDLG_OK,GLDLG_ICONERROR);
			break;
		case MENU_VERTEX_COORDINATES:

			if(!vertexCoordinates) vertexCoordinates = new VertexCoordinates();
			vertexCoordinates->Display(&worker);
			break;

		case MENU_VERTEX_ADD:
			if( geom->IsLoaded() ) {
				if( !addVertex ) addVertex = new AddVertex(geom,&worker);
				addVertex->SetVisible(TRUE);
			} else GLMessageBox::Display("No geometry loaded.","No geometry",GLDLG_OK,GLDLG_ICONERROR);
			break;

		case MENU_VERTEX_REMOVE:
			if( geom->IsLoaded() ) {
				if( GLMessageBox::Display("Remove Selected vertices?\nNote: It will also affect facets that contain them!","Question",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONINFO)==GLDLG_OK )  {        
					if (AskToReset()) {
						if (worker.running) worker.Stop_Public();
						geom->RemoveSelectedVertex();
						geom->Rebuild(); //Will recalculate facet parameters
						UpdateModelParams();
						if (vertexCoordinates) vertexCoordinates->Update();
						if (facetCoordinates) facetCoordinates->UpdateFromSelection();
						if (profilePlotter) profilePlotter->Refresh();

						if (pressureEvolution) pressureEvolution->Refresh();
						if (timewisePlotter) timewisePlotter->Refresh();
						// Send to sub process
						try { worker.Reload(); } catch(Error &e) {
							GLMessageBox::Display((char *)e.GetMsg(),"Error reloading worker",GLDLG_OK,GLDLG_ICONERROR);
						}
					}
				}

			} else GLMessageBox::Display("No geometry loaded.","No geometry",GLDLG_OK,GLDLG_ICONERROR);
			break;

		case MENU_VIEW_FULLSCREEN:
			if( m_bWindowed ) {
				ToggleFullscreen();
				PlaceComponents();
			} else {
				Resize(1024,768,TRUE);
			}
			menu->GetSubMenu("View")->SetState(MENU_VIEW_FULLSCREEN,!m_bWindowed);
			break;

		case MENU_VIEW_ADDNEW:
			AddView();
			break;
			break;
		case  MENU_VIEW_CLEARALL:
			if( GLMessageBox::Display("Clear all views ?","Question",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONINFO)==GLDLG_OK ) {
				ClearAllViews();
			}
			break;

		case MENU_TIME_SETTINGS:
			if( !timeSettings ) timeSettings = new TimeSettings(&worker);
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
			if(!timewisePlotter) timewisePlotter = new TimewisePlotter();
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
		if( src->GetId()>=MENU_FILE_LOADRECENT && src->GetId()<MENU_FILE_LOADRECENT+nbRecent ) {
			if (AskToSave()) {
				if (worker.running) worker.Stop_Public();
				LoadFile(recents[src->GetId()-MENU_FILE_LOADRECENT]);
			}
		}

		// Show structure menu
		if( src->GetId()>=MENU_VIEW_STRUCTURE && src->GetId()<=MENU_VIEW_STRUCTURE+geom->GetNbStructure() ) {
			geom->viewStruct = src->GetId()-MENU_VIEW_STRUCTURE-1;
			if(src->GetId()>MENU_VIEW_STRUCTURE) geom->Unselect();
			UpdateStructMenu();
		}
		if( src->GetId()==MENU_VIEW_NEWSTRUCT) {
			AddStruct();
			UpdateStructMenu();
		}
		if( src->GetId()==MENU_VIEW_DELSTRUCT) {
			DeleteStruct();
			UpdateStructMenu();
		}
		if( src->GetId()==MENU_VIEW_PREVSTRUCT) {
			geom->viewStruct=Remainder(geom->viewStruct-1,geom->GetNbStructure());
			geom->Unselect();
			UpdateStructMenu();
		}
		if( src->GetId()==MENU_VIEW_NEXTSTRUCT) {
			geom->viewStruct=Remainder(geom->viewStruct+1,geom->GetNbStructure());
			geom->Unselect();
			UpdateStructMenu();
		}

		// Select selection
		if( src->GetId()>=MENU_SELECTION_SELECTIONS && src->GetId()<MENU_SELECTION_SELECTIONS+nbSelection ) {
			SelectSelection( src->GetId() - MENU_SELECTION_SELECTIONS );
		}
		// Clear selection
		if( src->GetId()>=MENU_SELECTION_CLEARSELECTIONS && src->GetId()<MENU_SELECTION_CLEARSELECTIONS+nbSelection ) {
			char tmpname[256];
			sprintf(tmpname,"Clear %s?",selections[src->GetId() - MENU_SELECTION_CLEARSELECTIONS].name);
			if( GLMessageBox::Display(tmpname,"Confirmation",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONINFO)==GLDLG_OK ) {
				ClearSelection( src->GetId() - MENU_SELECTION_CLEARSELECTIONS );
			}
		}

		// Memorize selection
		if( src->GetId()>=MENU_SELECTION_MEMORIZESELECTIONS && src->GetId()<MENU_SELECTION_MEMORIZESELECTIONS+nbSelection ) {
			OverWriteSelection( src->GetId() - MENU_SELECTION_MEMORIZESELECTIONS );
		}

		// Select view
		if( src->GetId()>=MENU_VIEW_VIEWS && src->GetId()<MENU_VIEW_VIEWS+nbView ) {
			SelectView( src->GetId() - MENU_VIEW_VIEWS );
		}
		// Clear view
		if( src->GetId()>=MENU_VIEW_CLEARVIEWS && src->GetId()<MENU_VIEW_CLEARVIEWS+nbView ) {
			char tmpname[256];
			sprintf(tmpname,"Clear %s?",views[src->GetId() - MENU_VIEW_CLEARVIEWS].name);
			if( GLMessageBox::Display(tmpname,"Confirmation",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONINFO)==GLDLG_OK ) {
				ClearView( src->GetId() - MENU_VIEW_CLEARVIEWS );
			}
		}
		// Memorize view
		if( src->GetId()>=MENU_VIEW_MEMORIZEVIEWS && src->GetId()<MENU_VIEW_MEMORIZEVIEWS+nbView ) {
			OverWriteView( src->GetId() - MENU_VIEW_MEMORIZEVIEWS );
		}
		break;

		//TEXT --------------------------------------------------------------------
	case MSG_TEXT_UPD:
		if( src == facetSticking ) {
			calcFlow();
			facetApplyBtn->SetEnabled(TRUE);
		} else if ( src == facetOpacity ) {
			facetApplyBtn->SetEnabled(TRUE);
		} else if ( src == facetDesTypeN ) {
			facetApplyBtn->SetEnabled(TRUE);
		} else if ( src == facetTemperature ) {
			calcFlow();
			facetApplyBtn->SetEnabled(TRUE);
		} else if ( src == facetFlow ) {
			double flow;
			double area;
			facetFlow->GetNumber(&flow);
			facetArea->GetNumber(&area);
			if (area==0) facetFlowArea->SetText("#DIV0");
			else facetFlowArea->SetText(flow/area);
			facetApplyBtn->SetEnabled(TRUE);
			facetFILabel->SetState(TRUE);
			facetFIAreaLabel->SetState(FALSE);
			// //else if ( src == facetMass ) {
			//  facetApplyBtn->SetEnabled(TRUE);
		} else if ( src == facetFlowArea ) {
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
		} else if ( src == facetPumping ) {
			calcSticking();
			facetApplyBtn->SetEnabled(TRUE);
			// //else if ( src == facetMass ) {
			//  facetApplyBtn->SetEnabled(TRUE);
		}
		break;

	case MSG_TEXT:
		if( src == facetSticking ) {
			ApplyFacetParams();
		}  else if ( src == facetDesTypeN ) {
			ApplyFacetParams();
		} else if ( src == facetOpacity ) {
			ApplyFacetParams();
		} else if ( src == facetTemperature ) {
			ApplyFacetParams();
		} else if ( src == facetPumping ) {
			ApplyFacetParams();
		} else if ( src == facetFlow ) {
			ApplyFacetParams();
		} else if ( src == facetFlowArea ) {
			ApplyFacetParams();
		}
		break;

		//COMBO -------------------------------------------------------------------
	case MSG_COMBO:

		if ( src == facetDesType ) {
			facetApplyBtn->SetEnabled(TRUE);
			BOOL hasDesorption=!facetDesType->GetSelectedIndex()==0;

			facetFlow->SetEditable(hasDesorption);
			facetFlowArea->SetEditable(hasDesorption);

			int color=(hasDesorption)?0:110;
			facetFILabel->SetTextColor(color,color,color);
			facetFIAreaLabel->SetTextColor(color,color,color);
			facetFILabel->SetEnabled(hasDesorption);
			facetFIAreaLabel->SetEnabled(hasDesorption);
			facetDesTypeN->SetEditable(facetDesType->GetSelectedIndex()==3);
		} else if ( src == facetRecType ) {
			facetApplyBtn->SetEnabled(TRUE);
		} else if ( src == facetSideType ) {
			facetApplyBtn->SetEnabled(TRUE);
		} else if ( src == modeCombo ) {

			compACBtn->SetEnabled(modeCombo->GetSelectedIndex()==1);

			singleACBtn->SetEnabled(modeCombo->GetSelectedIndex()==1);
			UpdateFacetHits();
		} 
		break;

		//TOGGLE ------------------------------------------------------------------
	case MSG_TOGGLE:
		if (src==facetFILabel) {
			facetFILabel->SetState(TRUE);
			facetFIAreaLabel->SetState(FALSE);
		} else if (src==facetFIAreaLabel) {
			facetFILabel->SetState(FALSE);
			facetFIAreaLabel->SetState(TRUE);
		}
		
		else {
			// Update viewer flags
			UpdateViewerFlags();
		}
		break;

		//LIST --------------------------------------------------------------------
	case MSG_LIST:
		if( src == facetList && geom->IsLoaded()) {
			int *sels = (int *)malloc((geom->GetNbFacet())*sizeof(int));
			int nbSel;
			facetList->GetSelectedRows(&sels,&nbSel,TRUE);
			geom->Unselect();
			for(int i=0;i<nbSel;i++)
				geom->Select(sels[i]);
			geom->UpdateSelection();
			UpdateFacetParams();
			SAFE_FREE(sels);
		}
		break;

		//GEOMVIEWER ------------------------------------------------------------------
	case MSG_GEOMVIEWER_MAXIMISE:
		if( src==viewer[0] ) {
			AnimateViewerChange(0);
		} else if( src==viewer[1] ) {
			AnimateViewerChange(1);
		} else if( src==viewer[2] ) {
			AnimateViewerChange(2);
		} else if( src==viewer[3] ) {
			AnimateViewerChange(3);
		}
		Place3DViewer();
		break;

	case MSG_GEOMVIEWER_SELECT: {
		SelectViewer(src->GetId());
								}break;

		//BUTTON ------------------------------------------------------------------
	case MSG_BUTTON:
		if( src == startSimu ) {
			changedSinceSave=TRUE;
			StartStopSimulation();
			resetSimu->SetEnabled(!worker.running);
		} else if ( src == resetSimu ) {
			changedSinceSave=TRUE;
			ResetSimulation();
		} else if ( src == facetApplyBtn ) {
			//changedSinceSave=TRUE;
			ApplyFacetParams();
		} else if ( src == facetDetailsBtn ) {
			if( facetDetails==NULL ) facetDetails = new FacetDetails();
			facetDetails->Display(&worker);
		} else if ( src == facetCoordBtn ) {
			if(!facetCoordinates) facetCoordinates = new FacetCoordinates();
			facetCoordinates->Display(&worker);
		} else if ( src == facetMoreBtn ) {
			if (!facetMesh) {
				facetMesh = new FacetMesh(&worker);
				int *selection;
				int nbSel;
				worker.GetGeometry()->GetSelection(&selection, &nbSel);
				facetMesh->Refresh(nbSel,selection);
				SAFE_FREE(selection);
			}
			/*// Position dialog next to Facet parameters
			int facetX, facetY, facetW, facetH;
			facetPanel->GetBounds(&facetX, &facetY, &facetW, &facetH);
			facetMesh->SetPosition(facetX - facetMesh->GetWidth() - 10, facetY + 20);*/
			facetMesh->SetVisible(!facetMesh->IsVisible());
		} else if ( src == showMoreBtn ) {
			if (!viewer3DSettings)	viewer3DSettings = new Viewer3DSettings();
			viewer3DSettings->SetVisible(!viewer3DSettings->IsVisible());
			viewer3DSettings->Refresh(geom, viewer[curViewer]);
			
		} else if ( src==compACBtn ) {
			try { lastUpdate = 0.0;worker.ComputeAC(m_fTime); } catch(Error &e) {
				GLMessageBox::Display((char *)e.GetMsg(),"Error",GLDLG_OK,GLDLG_ICONERROR);
				return;
			}
			break;
		} else if ( src==singleACBtn ) {
			try { lastUpdate = 0.0;worker.StepAC(m_fTime); } catch(Error &e) {
				GLMessageBox::Display((char *)e.GetMsg(),"Error",GLDLG_OK,GLDLG_ICONERROR);
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
			if (!globalSettings) globalSettings = new GlobalSettings();
			if (!globalSettings->IsVisible()) globalSettings->Display(&worker);
			else globalSettings->SetVisible(FALSE);
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


/*
//----------------------------------------------------------------------------
void MolFlow::SendHeartBeat(BOOL forced) {

if (((m_fTime-lastHeartBeat)>1.0f) || (m_fTime==lastAppTime)) { //last heartbeat sent more than 1 seconds ago or app doesn't update
worker.SendHeartBeat();
m_fTime+=0.001f;
lastHeartBeat=m_fTime;
}
lastAppTime=m_fTime;
}
*/

//-----------------------------------------------------------------------------

void MolFlow::AnimateViewerChange(int next,BOOL init) {

	double xs1,ys1,xs2,ys2;
	double xe1,ye1,xe2,ye2;
	int sx = m_screenWidth-205;  
	int fWidth  = m_screenWidth-215;
	int fHeight = m_screenHeight-27;
	int Width2  = fWidth/2-1;
	int Height2 = fHeight/2-1;

	// Reset to layout and make all visible
	
	if (!init) {
		for(int i=0;i<MAX_VIEWER;i++)  viewer[i]->SetVisible(TRUE);
		viewer[0]->SetBounds(3       ,3        ,Width2,Height2);
		viewer[1]->SetBounds(6+Width2,3        ,Width2,Height2);
		viewer[2]->SetBounds(3       ,6+Height2,Width2,Height2);
		viewer[3]->SetBounds(6+Width2,6+Height2,Width2,Height2);


		if( modeSolo ) {

			// Go from single to layout
			xs1 = (double)3;
			ys1 = (double)3;
			xs2 = (double)fWidth+xs1;
			ys2 = (double)fHeight+ys1;

			switch( next ) {
			case 0:
				xe1 = (double)(3);
				ye1 = (double)(3);
				break;
			case 1:
				xe1 = (double)(5+Width2);
				ye1 = (double)(3);
				break;
			case 2:
				xe1 = (double)(3);
				ye1 = (double)(5+Height2);
				break;
			case 3:
				xe1 = (double)(5+Width2);
				ye1 = (double)(5+Height2);
				break;
			}

			xe2 = (double)(Width2)+xe1;
			ye2 = (double)(Height2)+ye1;

		} else {

			// Go from layout to single
			xe1 = (double)3;
			ye1 = (double)3;
			xe2 = (double)fWidth+xe1;
			ye2 = (double)fHeight+ye1;

			switch( next ) {
			case 0:
				xs1 = (double)(3);
				ys1 = (double)(3);
				break;
			case 1:
				xs1 = (double)(5+Width2);
				ys1 = (double)(3);
				break;
			case 2:
				xs1 = (double)(3);
				ys1 = (double)(5+Height2);
				break;
			case 3:
				xs1 = (double)(5+Width2);
				ys1 = (double)(5+Height2);
				break;
			}

			xs2 = (double)(Width2)+xs1;
			ys2 = (double)(Height2)+ys1;

		}

		double t0=(double)SDL_GetTicks()/1000.0;
		double t1 = t0;
		double T = 0.15;
		if (init) T=0;

		while((t1-t0)<T) {
			double t = (t1-t0) / T;
			int x1 = (int)(xs1 + t*(xe1-xs1) + 0.5);
			int y1 = (int)(ys1 + t*(ye1-ys1) + 0.5);
			int x2 = (int)(xs2 + t*(xe2-xs2) + 0.5);
			int y2 = (int)(ys2 + t*(ye2-ys2) + 0.5);
			viewer[next]->SetBounds(x1,y1,x2-x1,y2-y1);
			if (!init) wnd->Paint();
			// Overides moving component
			if (!init) viewer[next]->Paint();
			// Paint modeless
			int n;
			if (!init) n = GLWindowManager::GetNbWindow();
			if (!init) GLWindowManager::RepaintRange(1,n);
			t1 = (double)SDL_GetTicks()/1000.0;
		}

	} else { //init
		wnd->Paint();
		viewer[0]->Paint();
		int n = GLWindowManager::GetNbWindow();
		GLWindowManager::RepaintRange(1,n);
	}

	modeSolo = !modeSolo;
	SelectViewer(next);

}

//-----------------------------------------------------------------------------

/*
void MolFlow::BuildPipeStick(double s) {

char tmp[128];
Geometry *geom = worker.GetGeometry();

double R = 1.0;
double L = 1000.0;

ResetSimulation(FALSE);
geom->BuildPipe(L,R,s,100);
worker.nbDesorption = 0;
sprintf(tmp,"L|R %g S=%g",L/R,s);
worker.SetFileName(tmp);
nbDesStart = 0;
nbHitStart = 0;
for(int i=0;i<MAX_VIEWER;i++)
viewer[i]->SetWorker(&worker);
UpdateModelParams();
startSimu->SetEnabled(TRUE);
compACBtn->SetEnabled(modeCombo->GetSelectedIndex()==1);
resetSimu->SetEnabled(TRUE);
ClearFacetParams();
if( nbFormula==0 ) {
GLParser *f = new GLParser();
f->SetExpression("A29/D28");
f->SetName("Trans. Prob.");
f->Parse();
AddFormula(f);
}

// Send to sub process
worker.SetMaxDesorption(10000000*(nbSt+1));

UpdateTitle();

}
*/
// ---------------------------------------------------------------------------
BOOL MolFlow::AskToSave() {
	if (!changedSinceSave) return TRUE;
	int ret = GLSaveDialog::Display("Save current geometry first?","File not saved",GLDLG_SAVE | GLDLG_DISCARD | GLDLG_CANCEL_S,GLDLG_ICONINFO);
	if (ret==GLDLG_SAVE) {
		FILENAME *fn = GLFileBox::SaveFile(currentDir,worker.GetShortFileName(),"Save File",fileSFilters,nbSFilter);
		if( fn ) {
			GLProgress *progressDlg2 = new GLProgress("Saving file...","Please wait");
			progressDlg2->SetVisible(TRUE);
			progressDlg2->SetProgress(0.0);
			//GLWindowManager::Repaint();
			try {
				worker.SaveGeometry(fn->fullName,progressDlg2);
				changedSinceSave=FALSE;
				UpdateCurrentDir(fn->fullName);
				UpdateTitle();
				AddRecent(fn->fullName);
			} catch (Error &e) {
				char errMsg[512];
				sprintf(errMsg,"%s\nFile:%s",e.GetMsg(),fn->fullName);
				GLMessageBox::Display(errMsg,"Error",GLDLG_OK,GLDLG_ICONERROR);
				RemoveRecent(fn->fullName);
			}
			progressDlg2->SetVisible(FALSE);
			SAFE_DELETE(progressDlg2);
			return TRUE;
		} else return FALSE;
	} else if (ret==GLDLG_DISCARD) return TRUE;
	return FALSE;
}

// ---------------------------------------------------------------------------
BOOL MolFlow::AskToReset(Worker *work) {
	if (work==NULL) work=&worker;
	if (work->nbHit>0) {
		int rep = GLMessageBox::Display("This will reset simulation data.","Geometry change",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONWARNING);
		if( rep == GLDLG_OK ) {
			work->Reset(m_fTime);
			mApp->nbDesStart = 0;
			mApp->nbHitStart = 0;

			//resetSimu->SetEnabled(FALSE);
			if(mApp->profilePlotter) mApp->profilePlotter->Update(m_fTime,TRUE);

			if(mApp->pressureEvolution) mApp->pressureEvolution->Update(m_fTime,TRUE);
			if(mApp->timewisePlotter) mApp->timewisePlotter->Update(m_fTime,TRUE);
			if(mApp->texturePlotter) mApp->texturePlotter->Update(m_fTime,TRUE);
			return TRUE;
		} else return FALSE;
	} else return TRUE;
}


//-----------------------------------------------------------------------------

void MolFlow::QuickPipe() {

	BuildPipe(5.0, 5);
}

void MolFlow::BuildPipe(double ratio,int steps) {

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
		geom->BuildPipe(L,R,0,step);
	} catch(Error &e) {
		GLMessageBox::Display((char *)e.GetMsg(),"Error building pipe",GLDLG_OK,GLDLG_ICONERROR);
		geom->Clear();
		return;
	}
	worker.nbDesorption = 0;
	worker.needsReload = TRUE;
	//sprintf(tmp,"L|R %g",L/R);
	worker.SetFileName("");
	nbDesStart = 0;
	nbHitStart = 0;
	for(int i=0;i<MAX_VIEWER;i++)
		viewer[i]->SetWorker(&worker);
	//UpdateModelParams();
	startSimu->SetEnabled(TRUE);
	compACBtn->SetEnabled(modeCombo->GetSelectedIndex()==1);
	//resetSimu->SetEnabled(TRUE);
	ClearFacetParams();
	ClearFormula();
	ClearParameters();
	ClearAllSelections();
	ClearAllViews();

	GLParser *f = new GLParser();
	f->SetExpression("A2/SUMDES");
	f->SetName("Trans. Prob.");
	f->Parse();
	AddFormula(f);

	UpdateStructMenu();
	// Send to sub process
	try { worker.Reload(); } catch(Error &e) {
		GLMessageBox::Display((char *)e.GetMsg(),"Error",GLDLG_OK,GLDLG_ICONERROR);
		return;
	}
	if (pressureEvolution) pressureEvolution->Refresh();
	if (textureSettings) textureSettings->Update();
	if (texturePlotter) texturePlotter->Update(m_fTime, TRUE);
	if (outgassingMap) outgassingMap->Update(m_fTime, TRUE);
	if (facetDetails) facetDetails->Update();
	if (facetCoordinates) facetCoordinates->UpdateFromSelection();
	if (vertexCoordinates) vertexCoordinates->Update();
	if (movement) movement->Update();
	UpdateTitle();
	changedSinceSave=FALSE;
	ResetAutoSaveTimer();
}

//-----------------------------------------------------------------------------

void MolFlow::UpdateStructMenu() {

	char tmp[128];
	Geometry *geom = worker.GetGeometry();

	structMenu->Clear();
	structMenu->Add("New structure...",MENU_VIEW_NEWSTRUCT);
	structMenu->Add("Delete structure...",MENU_VIEW_DELSTRUCT);
	structMenu->Add(NULL); //Separator
	structMenu->Add("Show all",MENU_VIEW_STRUCTURE,SDLK_F1,CTRL_MODIFIER);
	structMenu->Add("Show previous",MENU_VIEW_PREVSTRUCT,SDLK_F11,CTRL_MODIFIER);
	structMenu->Add("Show next",MENU_VIEW_NEXTSTRUCT,SDLK_F12,CTRL_MODIFIER);
	structMenu->Add(NULL); //Separator

	for(int i=0;i<geom->GetNbStructure();i++) {
		sprintf(tmp,"Show #%d (%s)",i+1,geom->GetStructureName(i));
		if( i<10 )
			structMenu->Add(tmp,MENU_VIEW_STRUCTURE+(i+1),SDLK_F1+i+1,CTRL_MODIFIER);
		else
			structMenu->Add(tmp,MENU_VIEW_STRUCTURE+(i+1));
	}

	structMenu->SetState(MENU_VIEW_STRUCTURE+geom->viewStruct+1,TRUE);

	UpdateTitle();
}

//SELECTIONS
//-----------------------------------------------------------------------------

void MolFlow::SelectView(int v) {
	viewer[curViewer]->SetCurrentView(views[v]);
}

//-----------------------------------------------------------------------------

void MolFlow::SelectSelection(int v) {
	Geometry *geom = worker.GetGeometry();
	geom->SetSelection((&selections[v].selection),&(selections[v].nbSel));
}

//-----------------------------------------------------------------------------
void MolFlow::ClearSelectionMenus() {
	memorizeSelectionsMenu->Clear();
	memorizeSelectionsMenu->Add("Add new...",MENU_SELECTION_ADDNEW,SDLK_w,CTRL_MODIFIER);
	memorizeSelectionsMenu->Add(NULL); // Separator
	clearSelectionsMenu->Clear();
	clearSelectionsMenu->Add("Clear All",MENU_SELECTION_CLEARALL);
	clearSelectionsMenu->Add(NULL); // Separator
	selectionsMenu->Clear();
}

void MolFlow::RebuildSelectionMenus() {
	ClearSelectionMenus();
	for(int i=0;i<nbSelection;i++){
		if (i<=8) {
			selectionsMenu->Add(selections[i].name,MENU_SELECTION_SELECTIONS+i,SDLK_1+i,ALT_MODIFIER);
		} else {
			selectionsMenu->Add(selections[i].name,MENU_SELECTION_SELECTIONS+i); //no place for ALT+shortcut
		}
		clearSelectionsMenu->Add(selections[i].name,MENU_SELECTION_CLEARSELECTIONS+i);
		memorizeSelectionsMenu->Add(selections[i].name,MENU_SELECTION_MEMORIZESELECTIONS+i);
	}
}

void MolFlow::AddSelection(char *selectionName,ASELECTION s) {

	if(nbSelection<MAX_SELECTION) {
		selections[nbSelection] = s;
		selections[nbSelection].name = _strdup(selectionName);
		nbSelection++;
	} else {
		SAFE_FREE(selections[0].name);
		for(int i=0;i<MAX_SELECTION-1;i++) selections[i] = selections[i+1];
		selections[MAX_SELECTION-1] = s;
		selections[MAX_SELECTION-1].name = _strdup(selectionName);
	}
	RebuildSelectionMenus();
}

void MolFlow::ClearSelection(int idClr) {
	SAFE_FREE(selections[idClr].name);
	for(int i=idClr;i<nbSelection-1;i++) selections[i] = selections[i+1];
	nbSelection--;
	RebuildSelectionMenus();
}

void MolFlow::ClearAllSelections() {
	for(int i=0;i<nbSelection;i++) SAFE_FREE(selections[i].name);
	nbSelection=0;
	ClearSelectionMenus();
}

void MolFlow::OverWriteSelection(int idOvr) {
	Geometry *geom = worker.GetGeometry();
	char *selectionName = GLInputBox::GetInput(selections[idOvr].name,"Selection name","Enter selection name");
	if( !selectionName ) return;

	geom->GetSelection(&(selections[idOvr].selection),&(selections[idOvr].nbSel));
	selections[idOvr].name = _strdup(selectionName);
	RebuildSelectionMenus();
}

void MolFlow::AddSelection() {
	Geometry *geom = worker.GetGeometry();
	char tmp[32];
	sprintf(tmp,"Selection #%d",nbSelection+1);
	char *selectionName = GLInputBox::GetInput(tmp,"Selection name","Enter selection name");
	if( !selectionName ) return;

	if(nbSelection<MAX_SELECTION) {
		geom->GetSelection(&(selections[nbSelection].selection),&(selections[nbSelection].nbSel));
		selections[nbSelection].name = _strdup(selectionName);
		nbSelection++;
	} else {
		SAFE_FREE(selections[0].name);
		for(int i=0;i<MAX_SELECTION-1;i++) selections[i] = selections[i+1];
		geom->GetSelection(&(selections[MAX_SELECTION-1].selection),&(selections[MAX_SELECTION-1].nbSel));
		selections[MAX_SELECTION-1].name = _strdup(selectionName);
	}
	RebuildSelectionMenus();
}

//VIEWS
//-----------------------------------------------------------------------------
void MolFlow::ClearViewMenus() {
	memorizeViewsMenu->Clear();
	memorizeViewsMenu->Add("Add new...",MENU_VIEW_ADDNEW,SDLK_q,CTRL_MODIFIER);
	memorizeViewsMenu->Add(NULL); // Separator
	clearViewsMenu->Clear();
	clearViewsMenu->Add("Clear All",MENU_VIEW_CLEARALL);
	clearViewsMenu->Add(NULL); // Separator
	viewsMenu->Clear();
}

void MolFlow::RebuildViewMenus() {
	ClearViewMenus();
	for(int i=0;i<nbView;i++){
		int id=i;
		if (nbView>=10) id=i-nbView+8;
		if (id>=0 && id<=8) {
			viewsMenu->Add(views[i].name,MENU_VIEW_VIEWS+i,SDLK_F1+id,ALT_MODIFIER);
		} else {
			viewsMenu->Add(views[i].name,MENU_VIEW_VIEWS+i);
		}
		clearViewsMenu->Add(views[i].name,MENU_VIEW_CLEARVIEWS+i);
		memorizeViewsMenu->Add(views[i].name,MENU_VIEW_MEMORIZEVIEWS+i);
	}
}














void MolFlow::AddView(char *viewName,AVIEW v) {

	if(nbView<MAX_VIEW) {
		views[nbView] = v;
		views[nbView].name = _strdup(viewName);
		nbView++;
	} else {
		SAFE_FREE(views[0].name);
		for(int i=0;i<MAX_VIEW-1;i++) views[i] = views[i+1];
		views[MAX_VIEW-1] = v;
		views[MAX_VIEW-1].name = _strdup(viewName);
	}
	RebuildViewMenus();
}

void MolFlow::ClearView(int idClr) {
	SAFE_FREE(views[idClr].name);
	for(int i=idClr;i<nbView-1;i++) views[i] = views[i+1];
	nbView--;
	RebuildViewMenus();
}

void MolFlow::ClearAllViews() {
	for(int i=0;i<nbView;i++) SAFE_FREE(views[i].name);
	nbView=0;
	ClearViewMenus();
}

void MolFlow::OverWriteView(int idOvr) {
	Geometry *geom = worker.GetGeometry();
	char *viewName = GLInputBox::GetInput(views[idOvr].name,"View name","Enter view name");
	if( !viewName ) return;

	views[idOvr] = viewer[curViewer]->GetCurrentView();
	views[idOvr].name = _strdup(viewName);
	RebuildViewMenus();
}

void MolFlow::AddView() {
	Geometry *geom = worker.GetGeometry();
	char tmp[32];
	sprintf(tmp,"View #%d",nbView+1);
	char *viewName = GLInputBox::GetInput(tmp,"View name","Enter view name");
	if( !viewName ) return;

	if(nbView<MAX_VIEW) {
		views[nbView] = viewer[curViewer]->GetCurrentView();
		views[nbView].name = _strdup(viewName);
		nbView++;
	} else {
		SAFE_FREE(views[0].name);
		for(int i=0;i<MAX_VIEW-1;i++) views[i] = views[i+1];
		views[MAX_VIEW-1] = viewer[curViewer]->GetCurrentView();
		views[MAX_VIEW-1].name = _strdup(viewName);
	}
	RebuildViewMenus();
}

//-----------------------------------------------------------------------------

void MolFlow::RemoveRecent(char *fileName) {

	if( !fileName ) return;

	BOOL found = FALSE;
	int i = 0;
	while(!found && i<nbRecent) {
		found = strcmp(fileName,recents[i])==0;
		if(!found) i++;
	}
	if(!found) return;

	SAFE_FREE(recents[i]);
	for(int j=i;j<nbRecent-1;j++)
		recents[j] = recents[j+1];
	nbRecent--;

	// Update menu
	GLMenu *m = menu->GetSubMenu("File")->GetSubMenu("Load recent");
	m->Clear();
	for(i=nbRecent-1;i>=0;i--)
		m->Add(recents[i],MENU_FILE_LOADRECENT+i);
	SaveConfig();
}

//-----------------------------------------------------------------------------

void MolFlow::AddRecent(char *fileName) {

	// Check if already exists
	BOOL found = FALSE;
	int i = 0;
	while(!found && i<nbRecent) {
		found = strcmp(fileName,recents[i])==0;
		if(!found) i++;
	}
	if(found) {
		for (int j=i;j<nbRecent-1;j++) {
			recents[j]=recents[j+1];
		}
		recents[nbRecent-1]=_strdup(fileName);
		// Update menu
		GLMenu *m = menu->GetSubMenu("File")->GetSubMenu("Load recent");
		m->Clear();
		for(int i=nbRecent-1;i>=0;i--)
			m->Add(recents[i],MENU_FILE_LOADRECENT+i);
		SaveConfig();
		return;
	}

	// Add the new recent file
	if(nbRecent<MAX_RECENT) {
		recents[nbRecent] = _strdup(fileName);
		nbRecent++;
	} else {
		// Shift
		SAFE_FREE(recents[0]);
		for(int i=0;i<MAX_RECENT-1;i++)
			recents[i] = recents[i+1];
		recents[MAX_RECENT-1]=_strdup(fileName);
	}

	// Update menu
	GLMenu *m = menu->GetSubMenu("File")->GetSubMenu("Load recent");
	m->Clear();
	for(int i=nbRecent-1;i>=0;i--)
		m->Add(recents[i],MENU_FILE_LOADRECENT+i);
	SaveConfig();
}


void MolFlow::LoadConfig() {

	FileReader *f = NULL;
	char *w;
	nbRecent = 0;

	try {

		f = new FileReader("molflow.cfg");
		Geometry *geom = worker.GetGeometry();

		f->ReadKeyword("showRules");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->showRule = f->ReadInt();
		f->ReadKeyword("showNormals");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->showNormal = f->ReadInt();
		f->ReadKeyword("showUV");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->showUV = f->ReadInt();
		f->ReadKeyword("showLines");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->showLine = f->ReadInt();
		f->ReadKeyword("showLeaks");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->showLeak = f->ReadInt();
		f->ReadKeyword("showHits");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->showHit = f->ReadInt();
		f->ReadKeyword("showVolume");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->showVolume = f->ReadInt();
		f->ReadKeyword("showTexture");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->showTexture = f->ReadInt();
		f->ReadKeyword("showFilter");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->showFilter = f->ReadInt();
		f->ReadKeyword("showIndices");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->showIndex = f->ReadInt();
		f->ReadKeyword("showVertices");f->ReadKeyword(":");
		for (int i = 0; i<MAX_VIEWER; i++)
			viewer[i]->showVertex = f->ReadInt();
		f->ReadKeyword("showMode");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->showBack = f->ReadInt();
		f->ReadKeyword("showMesh");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->showMesh = f->ReadInt();
		f->ReadKeyword("showHidden");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->showHidden = f->ReadInt();
		f->ReadKeyword("showHiddenVertex");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->showHiddenVertex = f->ReadInt();
		f->ReadKeyword("showTimeOverlay");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->showTime = f->ReadInt();
		f->ReadKeyword("texColormap");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->showColormap = f->ReadInt();
		f->ReadKeyword("translation");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->transStep = f->ReadDouble();
		f->ReadKeyword("dispNumLines");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->dispNumHits = f->ReadInt();
		f->ReadKeyword("dispNumLeaks");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->dispNumLeaks = f->ReadInt();
		f->ReadKeyword("dirShow"); f->ReadKeyword(":");
		for (int i = 0; i<MAX_VIEWER; i++)
			viewer[i]->showDir = f->ReadInt();
		f->ReadKeyword("dirNorme"); f->ReadKeyword(":");
		geom->SetNormeRatio(f->ReadDouble());
		f->ReadKeyword("dirAutoNormalize"); f->ReadKeyword(":");
		geom->SetAutoNorme(f->ReadDouble());
		f->ReadKeyword("dirCenter"); f->ReadKeyword(":");
		geom->SetCenterNorme(f->ReadDouble());
		f->ReadKeyword("angle");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->angleStep = f->ReadDouble();
		f->ReadKeyword("autoScale");f->ReadKeyword(":");
		geom->texAutoScale = f->ReadInt();
		f->ReadKeyword("autoScale_include_constant_flow");f->ReadKeyword(":");
		geom->texAutoScaleIncludeConstantFlow = f->ReadInt();

		f->ReadKeyword("textures_min_pressure_all");f->ReadKeyword(":");
		geom->texture_limits[0].autoscale.min.all = f->ReadDouble();
		f->ReadKeyword("textures_min_pressure_moments_only");f->ReadKeyword(":");
		geom->texture_limits[0].autoscale.min.moments_only = f->ReadDouble();
		f->ReadKeyword("textures_max_pressure_all");f->ReadKeyword(":");
		geom->texture_limits[0].autoscale.max.all = f->ReadDouble();
		f->ReadKeyword("textures_max_pressure_moments_only");f->ReadKeyword(":");
		geom->texture_limits[0].autoscale.max.moments_only = f->ReadDouble();

		f->ReadKeyword("textures_min_impingement_all");f->ReadKeyword(":");
		geom->texture_limits[1].autoscale.min.all = f->ReadDouble();
		f->ReadKeyword("textures_min_impingement_moments_only");f->ReadKeyword(":");
		geom->texture_limits[1].autoscale.min.moments_only = f->ReadDouble();
		f->ReadKeyword("textures_max_impingement_all");f->ReadKeyword(":");
		geom->texture_limits[1].autoscale.max.all = f->ReadDouble();
		f->ReadKeyword("textures_max_impingement_moments_only");f->ReadKeyword(":");
		geom->texture_limits[1].autoscale.max.moments_only = f->ReadDouble();

		f->ReadKeyword("textures_min_density_all");f->ReadKeyword(":");
		geom->texture_limits[2].autoscale.min.all = f->ReadDouble();
		f->ReadKeyword("textures_min_density_moments_only");f->ReadKeyword(":");
		geom->texture_limits[2].autoscale.min.moments_only = f->ReadDouble();
		f->ReadKeyword("textures_max_density_all");f->ReadKeyword(":");
		geom->texture_limits[2].autoscale.max.all = f->ReadDouble();
		f->ReadKeyword("textures_max_density_moments_only");f->ReadKeyword(":");
		geom->texture_limits[2].autoscale.max.moments_only = f->ReadDouble();

		f->ReadKeyword("processNum");f->ReadKeyword(":");
		nbProc = f->ReadInt();
#ifdef _DEBUG
		nbProc=1;
#endif
		if(nbProc<=0) nbProc=1;
		f->ReadKeyword("recents");f->ReadKeyword(":");f->ReadKeyword("{");
		w = f->ReadString();
		while( strcmp(w,"}")!=0 && nbRecent<MAX_RECENT)  {
			recents[nbRecent] = _strdup(w);
			nbRecent++;
			w = f->ReadString();
		}

		f->ReadKeyword("cdir");f->ReadKeyword(":");
		strcpy(currentDir,f->ReadString());
		f->ReadKeyword("cseldir");f->ReadKeyword(":");
		strcpy(currentSelDir,f->ReadString());
		f->ReadKeyword("autonorme");f->ReadKeyword(":");
		geom->SetAutoNorme(f->ReadInt());
		f->ReadKeyword("centernorme");f->ReadKeyword(":");
		geom->SetCenterNorme(f->ReadInt());
		f->ReadKeyword("normeratio");f->ReadKeyword(":");
		geom->SetNormeRatio((float)(f->ReadDouble()));
		f->ReadKeyword("showDirection");f->ReadKeyword(":");
		for(int i=0;i<MAX_VIEWER;i++)
			viewer[i]->showDir = f->ReadInt();

		f->ReadKeyword("autoSaveFrequency");f->ReadKeyword(":");
		autoSaveFrequency = f->ReadDouble();
		f->ReadKeyword("autoSaveSimuOnly");f->ReadKeyword(":");
		autoSaveSimuOnly = f->ReadInt();
		f->ReadKeyword("checkForUpdates");f->ReadKeyword(":");
		checkForUpdates = f->ReadInt();
		f->ReadKeyword("autoUpdateFormulas");f->ReadKeyword(":");
		autoUpdateFormulas = f->ReadInt();
		f->ReadKeyword("compressSavedFiles");f->ReadKeyword(":");
		compressSavedFiles = f->ReadInt();
		f->ReadKeyword("gasMass");f->ReadKeyword(":");
		worker.gasMass = f->ReadDouble();
		f->ReadKeyword("expandShortcutPanel"); f->ReadKeyword(":");
		BOOL isOpen = f->ReadInt();
		if (isOpen) shortcutPanel->Open();
		else shortcutPanel->Close();

	} catch (Error &err) {
		std::ostringstream tmp;
		tmp << err.GetMsg() << "\n\nThis is normal on the first launch and if you upgrade from an earlier version\n";
		tmp << "MolFlow will use default program settings.\nWhen you quit, a correct config file will be written\n";
		GLMessageBox::Display(tmp.str().c_str(), "Error loading config file", GLDLG_OK, GLDLG_ICONINFO);
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

//----------------------------------------------------------------------------
void MolFlow::UpdateViewerFlags() {
	viewer[curViewer]->showNormal=showNormal->GetState();
	viewer[curViewer]->showRule=showRule->GetState();
	viewer[curViewer]->showUV=showUV->GetState();
	viewer[curViewer]->showLeak=showLeak->GetState();
	viewer[curViewer]->showHit=showHit->GetState();
	viewer[curViewer]->showLine=showLine->GetState();
	viewer[curViewer]->showVolume=showVolume->GetState();
	viewer[curViewer]->showTexture=showTexture->GetState();
	viewer[curViewer]->showFilter=showFilter->GetState();
	viewer[curViewer]->showVertex=showVertex->GetState();
	viewer[curViewer]->showIndex=showIndex->GetState();
	//worker.Update(0.0);
}

//-----------------------------------------------------------------------------

void MolFlow::SaveConfig() {

	FileWriter *f = NULL;

	try {

		f = new FileWriter("molflow.cfg");
		Geometry *geom = worker.GetGeometry();

		// Save flags
		WRITEI("showRules",showRule);
		WRITEI("showNormals",showNormal);
		WRITEI("showUV",showUV);
		WRITEI("showLines",showLine);
		WRITEI("showLeaks",showLeak);
		WRITEI("showHits",showHit);
		WRITEI("showVolume",showVolume);
		WRITEI("showTexture",showTexture);
		WRITEI("showFilter",showFilter);
		WRITEI("showIndices",showIndex);
		WRITEI("showVertices",showVertex);
		WRITEI("showMode",showBack);
		WRITEI("showMesh",showMesh);
		WRITEI("showHidden",showHidden);
		WRITEI("showHiddenVertex",showHiddenVertex);
		WRITEI("showTimeOverlay",showTime);
		WRITEI("texColormap",showColormap);
		WRITED("translation",transStep);
		WRITEI("dispNumLines",dispNumHits);
		WRITEI("dispNumLeaks",dispNumLeaks);
		WRITEI("dirShow", showDir);
		f->Write("dirNorme:");f->WriteDouble(geom->GetNormeRatio(), "\n");
		f->Write("dirAutoNormalize:"); f->WriteDouble(geom->GetAutoNorme(), "\n");
		f->Write("dirCenter:"); f->WriteDouble(geom->GetCenterNorme(), "\n");

		WRITED("angle",angleStep);
		f->Write("autoScale:");f->WriteInt(geom->texAutoScale,"\n");
		f->Write("autoScale_include_constant_flow:");f->WriteInt(geom->texAutoScaleIncludeConstantFlow,"\n");


		f->Write("textures_min_pressure_all:");
		f->WriteDouble(geom->texture_limits[0].autoscale.min.all,"\n");
		f->Write("textures_min_pressure_moments_only:");
		f->WriteDouble(geom->texture_limits[0].autoscale.min.moments_only,"\n");
		f->Write("textures_max_pressure_all:");
		f->WriteDouble(geom->texture_limits[0].autoscale.max.all,"\n");
		f->Write("textures_max_pressure_moments_only:");
		f->WriteDouble(geom->texture_limits[0].autoscale.max.moments_only,"\n");

		f->Write("textures_min_impingement_all:");
		f->WriteDouble(geom->texture_limits[1].autoscale.min.all,"\n");

		f->Write("textures_min_impingement_moments_only:");
		f->WriteDouble(geom->texture_limits[1].autoscale.min.moments_only,"\n");
		f->Write("textures_max_impingement_all:");
		f->WriteDouble(geom->texture_limits[1].autoscale.max.all ,"\n");
		f->Write("textures_max_impingement_moments_only:");
		f->WriteDouble(geom->texture_limits[1].autoscale.max.moments_only,"\n");

		f->Write("textures_min_density_all:");
		f->WriteDouble(geom->texture_limits[2].autoscale.min.all,"\n");
		f->Write("textures_min_density_moments_only:");
		f->WriteDouble(geom->texture_limits[2].autoscale.min.moments_only,"\n");
		f->Write("textures_max_density_all:");
		f->WriteDouble(geom->texture_limits[2].autoscale.max.all,"\n");
		f->Write("textures_max_density_moments_only:");
		f->WriteDouble(geom->texture_limits[2].autoscale.max.moments_only,"\n");

#ifdef _DEBUG
		f->Write("processNum:");f->WriteInt(numCPU,"\n");
#else
		f->Write("processNum:");f->WriteInt(worker.GetProcNumber(),"\n");
#endif
		f->Write("recents:{\n");
		for(int i=0;i<nbRecent;i++) {
			f->Write("\"");
			f->Write(recents[i]);
			f->Write("\"\n");
		}
		f->Write("}\n");

		f->Write("cdir:\"");f->Write(currentDir);f->Write("\"\n");
		f->Write("cseldir:\"");f->Write(currentSelDir);f->Write("\"\n");
		f->Write("autonorme:");f->WriteInt(geom->GetAutoNorme(),"\n");
		f->Write("centernorme:");f->WriteInt(geom->GetCenterNorme(),"\n");
		f->Write("normeratio:");f->WriteDouble((double)(geom->GetNormeRatio()),"\n");
		WRITEI("showDirection",showDir);f->Write("\n");

		f->Write("autoSaveFrequency:");f->WriteDouble(autoSaveFrequency,"\n");
		f->Write("autoSaveSimuOnly:");f->WriteInt(autoSaveSimuOnly,"\n");
		f->Write("checkForUpdates:");f->WriteInt(checkForUpdates,"\n");
		f->Write("autoUpdateFormulas:");f->WriteInt(autoUpdateFormulas,"\n");
		f->Write("compressSavedFiles:");f->WriteInt(compressSavedFiles,"\n");
		f->Write("gasMass:");f->WriteDouble(worker.gasMass,"\n");
		f->Write("expandShortcutPanel:");f->WriteInt(!shortcutPanel->IsClosed(), "\n");

	} catch (Error &err) {
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

	flow=1*sticking*area/10.0/4.0*sqrt(8.0*8.31*temperature/PI/(worker.gasMass*0.001));
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

	sticking=abs(flow/(area/10.0)*4.0*sqrt(1.0/8.0/8.31/(temperature)*PI*(worker.gasMass*0.001)));
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
	sprintf(tmp,"Well, that's emberassing. Molflow crashed and will exit now.\nBefore that, an autosave will be attempted.\nHere is the error info:\n\n%s",(char *)e->GetMsg());
	GLMessageBox::Display(tmp,"Main crash handler",GLDLG_OK,GLDGL_ICONDEAD);
	try {
		if (AutoSave(TRUE))
			GLMessageBox::Display("Good news, autosave worked!","Main crash handler",GLDLG_OK,GLDGL_ICONDEAD);
		else 
			GLMessageBox::Display("Sorry, I couldn't even autosave.","Main crash handler",GLDLG_OK,GLDGL_ICONDEAD);
	} catch(Error &e) {
		e.GetMsg();
		GLMessageBox::Display("Sorry, I couldn't even autosave.","Main crash handler",GLDLG_OK,GLDGL_ICONDEAD);
	}
}

void MolFlow::AddStruct() {
	Geometry *geom = worker.GetGeometry();
	char tmp[32];
	sprintf(tmp,"Structure #%d",geom->GetNbStructure()+1);
	char *structName = GLInputBox::GetInput(tmp,"Structure name","Enter name of new structure");
	if( !structName ) return;
	geom->AddStruct(structName);
	// Send to sub process
	try { worker.Reload(); } catch(Error &e) {
		GLMessageBox::Display((char *)e.GetMsg(),"Error",GLDLG_OK,GLDLG_ICONERROR);
		return;
	}
}

void MolFlow::DeleteStruct() {
	Geometry *geom = worker.GetGeometry();
	char *structNum = GLInputBox::GetInput("","Structure number","Number of structure to delete:");
	if( !structNum ) return;
	int structNumInt;
	if ( !sscanf(structNum,"%d",&structNumInt)) {
		GLMessageBox::Display("Invalid structure number");
		return;
	}
	if (structNumInt<1 || structNumInt>geom->GetNbStructure()) {
		GLMessageBox::Display("Invalid structure number");
		return;
	}
	BOOL hasFacets=FALSE;
	for (int i=0;i<geom->GetNbFacet() && !hasFacets;i++) {
		if (geom->GetFacet(i)->sh.superIdx==(structNumInt-1)) hasFacets=TRUE;
	}
	if (hasFacets) {
		int rep = GLMessageBox::Display("This structure has facets. They will be deleted with the structure.","Structure delete",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONWARNING);
		if( rep != GLDLG_OK ) return;
	}
	if (!AskToReset()) return;
	geom->DelStruct(structNumInt-1);
	// Send to sub process
	try { worker.Reload(); } catch(Error &e) {
		GLMessageBox::Display((char *)e.GetMsg(),"Error",GLDLG_OK,GLDLG_ICONERROR);
		return;
	}
}

void MolFlow::DisplayCollapseDialog() {
	Geometry *geom = worker.GetGeometry();
	if( !collapseSettings ) collapseSettings = new CollapseSettings();
	collapseSettings->SetGeometry(geom,&worker);
	collapseSettings->SetVisible(TRUE);
}




void MolFlow::RenumberSelections(int startFacetId){
	for (int i=0;i<nbSelection;i++) {
		BOOL found=FALSE;
		for (int j=0;j<selections[i].nbSel && !found;j++) { //remove from selection
			if (selections[i].selection[j]==startFacetId) {

				for (int k=j;k<(selections[i].nbSel-1);k++)
					selections[i].selection[k]=selections[i].selection[k+1];
				selections[i].nbSel--;
				if (selections[i].nbSel==0) ClearSelection(i); //last facet removed from selection
				found=TRUE;
			} 
		}
		for (int j=0;j<selections[i].nbSel;j++) { //renumber selection
			if (selections[i].selection[j]>startFacetId) {
				//decrease number by 1
				selections[i].selection[j]--; 
			}
		}
	}
}