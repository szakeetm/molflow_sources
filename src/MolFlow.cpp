#include "MolFlow.h"
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
#include <cmath>
#include "MolFlow.h"
#include "Facet_shared.h"
#include "MolflowGeometry.h"
#include "File.h"
#include "GLApp/GLMessageBox.h"
#include "GLApp/GLInputBox.h"
#include "NativeFileDialog/molflow_wrapper/nfd_wrapper.h"
#include "GLApp/GLWindowManager.h"
#include "Helper/MathTools.h"
#include <Helper/FormatHelper.h>

#include "GLApp/GLMenuBar.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLLabel.h"
#include "GLApp/GLCombo.h"
#include "GLApp/GLTextField.h"

#include "Interface/RecoveryDialog.h"
#include <vector>
#include <string>
#include <numeric> //std::iota
#include <filesystem>

#ifdef _WIN32
#include "direct.h"
#else
#include <unistd.h> //chdir
#endif

// Plotters
#include "Interface/ProfilePlotter.h"
#include "ProfileModes.h"
#include "Interface/PressureEvolution.h"
#include "Interface/TimewisePlotter.h"
#include "Interface/TexturePlotter.h"
#include "Interface/ConvergencePlotter.h"

#include "AppUpdater.h"
#include "Worker.h"
#include "Interface/ImportDesorption.h"
#include "Interface/TimeSettings.h"
#include "Interface/Movement.h"
#include "Interface/MeasureForce.h"
#include "Interface/FacetAdvParams.h"
#include "Interface/FacetDetails.h"
#include "Interface/Viewer3DSettings.h"
#include "Interface/TextureScaling.h"
#include "Interface/GlobalSettings.h"
#include "Interface/OutgassingMapWindow.h"
#include "Interface/MomentsEditor.h"
#include "Interface/FacetCoordinates.h"
#include "Interface/VertexCoordinates.h"
#include "Interface/ParameterEditor.h"
#include "Interface/SmartSelection.h"
#include "Interface/FormulaEditor.h"
#include "Interface/ParticleLogger.h"
#include "Interface/HistogramSettings.h"
#include "Interface/HistogramPlotter.h"
#include "FormulaEvaluator_MF.h"

/*
static const char *fileLFilters = "All MolFlow supported files\0*.txt;*.xml;*.zip;*.geo;*.geo7z;*.syn;*.syn7z;*.str;*.stl;*.ase\0"
"All files\0*.*\0";
static const char *fileInsFilters = "All insertable geometries\0*.txt;*.xml;*.zip;*.geo;*.geo7z;*.syn;*.syn7z;*.stl\0"
"All files\0*.*\0";
const char *fileSFilters = "MolFlow saveable files\0*.xml;*.zip;*.geo;*.geo7z;*.txt\0All files\0*.*\0";
static const char *fileDesFilters = "Desorption files\0*.des\0All files\0*.*\0";
*/

//NativeFileDialog compatible file filters
char fileLoadFilters[] = "txt,xml,zip,geo,syn,str,stl,geo7z,syn7z";
char fileInsertFilters[] = "txt,xml,zip,geo,syn,stl,geo7z,syn7z";
char fileSaveFilters[] = "zip,xml,txt,geo,stl,geo7z";
char fileSelFilters[] = "sel";
char fileTexFilters[] = "txt";
char fileProfFilters[] = "csv,txt";


int cSize = 4;
int   cWidth[] = { 30, 56, 50, 50 };
const char* cName[] = { "#", "Hits", "Des", "Abs" };

std::vector<std::string> formulaPrefixes = { "A","D","H","MCH","P","DEN","Z","V","T","AR","a","d","h","mch","p","den","z","v","t","ar","," };
char formulaSyntax[] =
R"(MC Variables: An (Absorption on facet n), Dn (Desorption on facet n), Hn (Hit on facet n)
Pn (Pressure [mbar] on facet n), DENn (Density [1/m3] on facet n)
Zn (Imp. rate on facet n), Vn (avg. speed [m/s] on facet n), Tn (temp[K] of facet n)

Forcen, ForceXn,ForceYn,ForceZn - the molecular force [N] on facet n (Norme or X,Y,Z component)
ForceSqrn, ForceSqrXn,ForceSqrYn,ForceSqrZn - square of mol. force [N2] on facet n
Torquen, TorqueXn,TorqueYn,TorqueZn - torque relative to ref. point [Nm] on facet n
(set reference point in Tools / Measure Forces...)

SUMABS (total absorbed), SUMDES (total desorbed), SUMHIT (total hit)

Sum over multiple facets:
SUM(H,3,8)    calculates the sum of hits on facets 3,4,... ...7,8.
SUM(H,S2)     calculates the sum of hits on selection group #2
SUM(H,SEL)    calculates the sum of hits on the current selection
SUM works with H,A,D,AR and Force, ForceSqr, Torque and their X,Y,Z components

Average over multiple facets:
same syntax as above, replace SUM with AVG in the formulas
AVG works (area-weighted averaging): P, DEN, Z
AVG works (equal weight per facet): Force, ForceSqr, Torque and their X,Y,Z components

Area variables: ARn (Area of facet n), DESAR (total desorption area), ABSAR (total absorption area)

Final (constant) outgassing rate [mbar*l/s]: QCONST
Final (constant) outgassing rate [molecules/s]: QCONST_N
Total desorbed molecules until last moment: [molecules]: NTOT
Gas mass [g/mol]: GASMASS

Mean Pumping Path: MPP (average path of molecules in the system before absorption)
Mean Free Path:      MFP (average path of molecules between two wall hits)

Math functions: sin(), cos(), tan(), sinh(), cosh(), tanh(), asin(), acos(),
                     atan(), exp(), ln(), pow(x,y), log2(), log10(), sqrt(), abs()

Constants:  Kb (Boltzmann's constant), R (Gas constant), Na (Avogadro's number), PI
Previous formula values: by index (ex. Formula3) or name (ex. "Trans.prob.")
)";
int formulaSyntaxHeight = 525;

MolFlow* mApp;

//Menu elements, Molflow specific:
//#define MENU_FILE_IMPORTDES_DES 140
#define MENU_FILE_IMPORTDES_SYN 141

#define MENU_FILE_EXPORTTEXTURE_AREA 151
#define MENU_FILE_EXPORTTEXTURE_MCHITS 152
#define MENU_FILE_EXPORTTEXTURE_IMPINGEMENT 153
#define MENU_FILE_EXPORTTEXTURE_PART_DENSITY 154
#define MENU_FILE_EXPORTTEXTURE_GAS_DENSITY 155
#define MENU_FILE_EXPORTTEXTURE_PRESSURE 156
#define MENU_FILE_EXPORTTEXTURE_AVG_V 157
#define MENU_FILE_EXPORTTEXTURE_V_VECTOR 158
#define MENU_FILE_EXPORTTEXTURE_N_VECTORS 159

#define MENU_FILE_EXPORTTEXTURE_AREA_COORD 171
#define MENU_FILE_EXPORTTEXTURE_MCHITS_COORD  172
#define MENU_FILE_EXPORTTEXTURE_IMPINGEMENT_COORD  173
#define MENU_FILE_EXPORTTEXTURE_PART_DENSITY_COORD  174
#define MENU_FILE_EXPORTTEXTURE_GAS_DENSITY_COORD  175
#define MENU_FILE_EXPORTTEXTURE_PRESSURE_COORD  176
#define MENU_FILE_EXPORTTEXTURE_AVG_V_COORD  177
#define MENU_FILE_EXPORTTEXTURE_V_VECTOR_COORD  178
#define MENU_FILE_EXPORTTEXTURE_N_VECTORS_COORD  179

#define MENU_TOOLS_MOVINGPARTS 410
#define MENU_TOOLS_MEASUREFORCE 420

#define MENU_SELECT_HASDESFILE 361
#define MENU_FACET_OUTGASSINGMAP 362

#define MENU_TIME_SETTINGS          900
#define MENU_TIMEWISE_PLOTTER       901
#define MENU_TIME_PRESSUREEVOLUTION 902
#define MENU_TIME_MOMENTS_EDITOR    903
#define MENU_TIME_PARAMETER_EDITOR  904

// Name: WinMain()
// Desc: Entry point to the program. Initializes everything, and goes into a
//       message-processing loop. Idle time is used to render the scene.

//INT WINAPI WinMain(HINSTANCE hInst, HINSTANCE, LPSTR, INT) //Can be replaced with main() if including SDL2Main.lib
int main(int argc, char* argv[])
{
	mApp = new MolFlow();

#ifndef _WIN32
	//Change working directory to executable path (if launched by dbl-click)
	std::string myPath = FileUtils::GetPath(argv[0]);
	if (!myPath.empty()) chdir(myPath.c_str());
#endif

	if (!mApp->Create(1024, 800, false)) {
		char* logs = GLToolkit::GetLogs();
#ifdef _WIN32
		if (logs) MessageBox(nullptr, logs, "Molflow [Fatal error]", MB_OK);
#else
		if (logs) {
			Log::console_error("Molflow [Fatal error]\n");
			Log::console_error("{}", logs);
		}
#endif
		SAFE_FREE(logs);
		delete mApp;
		return -1;
	}
	try {
		mApp->Run();
	}
	catch (const std::exception& e) {
		mApp->CrashHandler(e);
	}
	delete mApp;
	return 0;
}

// Name: MolFlow()
// Desc: Application constructor. Sets default attributes for the app.

MolFlow::MolFlow()
{
	mApp = this; //to refer to the app as extern variable

	auto eval = std::make_shared<FormulaEvaluator_MF>(&worker, (MolflowGeometry*)worker.GetGeometry(), &selections);
	appFormulas = std::make_shared<Formulas>(eval);
}

// Name: OneTimeSceneInit()
// Desc: Called during initial app startup, this function performs all the
//       permanent initialization.
int MolFlow::OneTimeSceneInit()
{
	/*
	//Enable memory check at EVERY malloc/free operation:
	int tmpFlag = _CrtSetDbgFlag( _CRTDBG_REPORT_FLAG );
	tmpFlag |= _CRTDBG_LEAK_CHECK_DF;
	_CrtSetDbgFlag( tmpFlag );
	*/

	OneTimeSceneInit_shared_pre();

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

	menu->GetSubMenu("File")->Add("Import desorption from SYN file", MENU_FILE_IMPORTDES_SYN);
	//menu->GetSubMenu("File")->GetSubMenu("Import desorption file")->Add("SYN file", );
	//menu->GetSubMenu("File")->GetSubMenu("Import desorption file")->Add("DES file (deprecated)", MENU_FILE_IMPORTDES_DES);

	menu->GetSubMenu("File")->Add(nullptr); // Separator
	menu->GetSubMenu("File")->Add("E&xit", MENU_FILE_EXIT); //Moved here from OnetimeSceneinit_shared to assert it's the last menu item

	menu->GetSubMenu("Selection")->Add(nullptr); // Separator
	menu->GetSubMenu("Selection")->Add("Select Desorption", MENU_FACET_SELECTDES);
	menu->GetSubMenu("Selection")->Add("Select Dyn. Desorption (Outg.Map)", MENU_SELECT_HASDESFILE);
	menu->GetSubMenu("Selection")->Add("Select Reflective", MENU_FACET_SELECTREFL);

	menu->GetSubMenu("Tools")->Add(nullptr);
	menu->GetSubMenu("Tools")->Add("Moving parts...", MENU_TOOLS_MOVINGPARTS);
	menu->GetSubMenu("Tools")->Add("Measure forces...", MENU_TOOLS_MEASUREFORCE);

	menu->GetSubMenu("Facet")->Add("Convert to outgassing map...", MENU_FACET_OUTGASSINGMAP);

	menu->Add("Time");
	menu->GetSubMenu("Time")->Add("Time settings...", MENU_TIME_SETTINGS, SDLK_i, ALT_MODIFIER);
	menu->GetSubMenu("Time")->Add("Edit moments...", MENU_TIME_MOMENTS_EDITOR);
	menu->GetSubMenu("Time")->Add("Edit parameters...", MENU_TIME_PARAMETER_EDITOR);
	menu->GetSubMenu("Time")->Add(nullptr);
	menu->GetSubMenu("Time")->Add("Timewise plotter", MENU_TIMEWISE_PLOTTER);
	menu->GetSubMenu("Time")->Add("Pressure evolution", MENU_TIME_PRESSUREEVOLUTION);

	viewerMoreButton = new GLButton(0, "<< View");
	togglePanel->Add(viewerMoreButton);

	shortcutPanel = new GLTitledPanel("Shortcuts");
	shortcutPanel->SetClosable(true);
	shortcutPanel->Close();
	Add(shortcutPanel);

	profilePlotterBtn = new GLButton(0, "Profile pl.");
	shortcutPanel->Add(profilePlotterBtn);

	texturePlotterBtn = new GLButton(0, "Texture pl.");
	shortcutPanel->Add(texturePlotterBtn);

	textureScalingBtn = new GLButton(0, "Tex.scaling");
	shortcutPanel->Add(textureScalingBtn);

	/*globalSettingsBtn = new GLButton(0, "<< Sim");
	simuPanel->Add(globalSettingsBtn);
*/
/*
statusSimu = new GLButton(0,"...");
simuPanel->Add(statusSimu);
*/

	inputPanel = new GLTitledPanel("Particles in");
	facetPanel->Add(inputPanel);

	facetDesTypeLabel = new GLLabel("Desorption");
	facetPanel->Add(facetDesTypeLabel);
	facetDesType = new GLCombo(0);
	facetDesType->SetSize(5);
	facetDesType->SetValueAt(0, "None");
	facetDesType->SetValueAt(1, "Uniform");
	facetDesType->SetValueAt(2, "Cosine");
	facetDesType->SetValueAt(3, "Cosine^N");
	facetDesType->SetValueAt(4, "Recorded");
	inputPanel->Add(facetDesType);

	facetDesTypeN = new GLTextField(0, nullptr);
	facetDesTypeN->SetEditable(false);
	facetPanel->Add(facetDesTypeN);

	facetOutgToggleLabel = new GLToggle(0, "Outgassing (mbar*l/s):");
	facetOutgToggleLabel->SetEnabled(false);
	facetOutgToggleLabel->SetState(true);
	inputPanel->Add(facetOutgToggleLabel);
	facetOutgassingText = new GLTextField(0, nullptr);
	inputPanel->Add(facetOutgassingText);

	facetOutgPerAreaToggleLabel = new GLToggle(1, "Outg/area(mbar*l/s/cm\262):");
	facetOutgPerAreaToggleLabel->SetEnabled(false);
	inputPanel->Add(facetOutgPerAreaToggleLabel);
	facetOutgPerAreaText = new GLTextField(0, nullptr);
	inputPanel->Add(facetOutgPerAreaText);

	outputPanel = new GLTitledPanel("Particles out");
	facetPanel->Add(outputPanel);

	facetStickingLabel = new GLLabel("Sticking factor:");
	outputPanel->Add(facetStickingLabel);
	facetStickingText = new GLTextField(0, nullptr);
	outputPanel->Add(facetStickingText);

	facetPumpingLabel = new GLLabel("Pumping Speed (l/s):");
	outputPanel->Add(facetPumpingLabel);
	facetPumpingSpeedText = new GLTextField(0, nullptr);
	outputPanel->Add(facetPumpingSpeedText);

	facetTempLabel = new GLLabel("Temperature (\260K):");
	facetPanel->Add(facetTempLabel);
	facetTemperatureText = new GLTextField(0, nullptr);
	facetPanel->Add(facetTemperatureText);



	facetProfileLabel = new GLLabel("Profile:");
	facetPanel->Add(facetProfileLabel);
	facetProfileCombo = new GLCombo(0);
	size_t nbRecModes = (size_t)ProfileRecordModes::NUMITEMS;
	facetProfileCombo->SetSize(nbRecModes);
	for (size_t i = 0; i < nbRecModes; i++) {
		facetProfileCombo->SetValueAt(i, profileRecordModeDescriptions[(ProfileRecordModes)i].first); //long description
	}
	facetPanel->Add(facetProfileCombo);

	facetAdvParamsBtn = new GLButton(0, "<< Adv");
	facetPanel->Add(facetAdvParamsBtn);

	facetList = new GLList(0);
	facetList->SetWorker(&worker);
	facetList->SetGrid(true);
	facetList->SetSelectionMode(MULTIPLE_ROW);
	facetList->SetSize(4, 0);
	facetList->SetColumnWidths((int*)cWidth);
	facetList->SetColumnLabels((const char**)cName);
	facetList->SetColumnLabelVisible(true);
	facetList->Sortable = true;
	Add(facetList);

	facetAdvParams = new FacetAdvParams(&worker); //To use its UpdatefacetParams() routines

	TimeDependentParameters::LoadParameterCatalog(worker.interfaceParameterCache);
	OneTimeSceneInit_shared_post();

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
	togglePanel->SetCompBounds(showFacetId, 135, 64, 60, 18);

	togglePanel->SetCompBounds(viewerMoreButton, 5, 86, 55, 18);
	togglePanel->SetCompBounds(showIndex, 70, 86, 60, 18);
	togglePanel->SetCompBounds(showVertexId, 135, 86, 60, 18);

	sy += (togglePanel->GetHeight() + 5);

	// Selected facet -----------------------------------------
	facetPanel->SetBounds(sx, sy, 202, 330);

	facetPanel->SetCompBounds(inputPanel, 5, 16, 192, 90);

	int cursorY = 15;
	inputPanel->SetCompBounds(facetDesTypeLabel, 5, cursorY, 60, 18);
	inputPanel->SetCompBounds(facetDesType, 65, cursorY, 80, 18);
	inputPanel->SetCompBounds(facetDesTypeN, 150, cursorY, 30, 18);

	inputPanel->SetCompBounds(facetOutgToggleLabel, 5, cursorY += 25, 110, 18);
	inputPanel->SetCompBounds(facetOutgassingText, 140, cursorY, 45, 18);

	inputPanel->SetCompBounds(facetOutgPerAreaToggleLabel, 5, cursorY += 25, 110, 18);
	inputPanel->SetCompBounds(facetOutgPerAreaText, 140, cursorY, 45, 18);

	facetPanel->SetCompBounds(outputPanel, 5, cursorY += 45, 192, 65);

	outputPanel->SetCompBounds(facetStickingLabel, 7, cursorY = 15, 100, 18);
	outputPanel->SetCompBounds(facetStickingText, 140, cursorY, 45, 18);

	outputPanel->SetCompBounds(facetPumpingLabel, 7, cursorY += 25, 100, 18);
	outputPanel->SetCompBounds(facetPumpingSpeedText, 140, cursorY, 45, 18);

	facetPanel->SetCompBounds(facetSideLabel, 7, cursorY = 180, 50, 18);
	facetPanel->SetCompBounds(facetSideType, 65, cursorY, 130, 18);

	facetPanel->SetCompBounds(facetTLabel, 7, cursorY += 25, 100, 18);
	facetPanel->SetCompBounds(facetOpacity, 110, cursorY, 82, 18);

	facetPanel->SetCompBounds(facetTempLabel, 7, cursorY += 25, 100, 18);
	facetPanel->SetCompBounds(facetTemperatureText, 110, cursorY, 82, 18);

	facetPanel->SetCompBounds(facetAreaLabel, 7, cursorY += 25, 100, 18);
	facetPanel->SetCompBounds(facetAreaText, 110, cursorY, 82, 18);

	facetPanel->SetCompBounds(facetProfileLabel, 7, cursorY += 25, 60, 18);
	facetPanel->SetCompBounds(facetProfileCombo, 65, cursorY, 130, 18);

	facetPanel->SetCompBounds(facetAdvParamsBtn, 5, cursorY += 25, 48, 18);
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
	//simuPanel->SetCompBounds(modeLabel, 5, 45, 30, 18);
	//simuPanel->SetCompBounds(compACBtn, 130, 45, 65, 19);
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


	int lg = m_screenHeight - 23 /*- (nbFormula * 25)*/;

	facetList->SetBounds(sx, sy, 202, lg - sy);

	/*
	for (int i = 0; i < nbFormula; i++) {
		formulas[i].name->SetBounds(sx, lg + 5, 95, 18);
		formulas[i].value->SetBounds(sx + 90, lg + 5, 87, 18);
		formulas[i].setBtn->SetBounds(sx + 182, lg + 5, 20, 18);
		lg += 25;
	}
	*/
}

// Name: ClearFacetParams()
// Desc: Reset selected facet parameters.

void MolFlow::ClearFacetParams() {
	facetPanel->SetTitle("Selected Facet (none)");
	facetStickingText->Clear();
	facetStickingText->SetEditable(false);
	facetOutgToggleLabel->SetEnabled(false);
	facetOutgPerAreaToggleLabel->SetEnabled(false);
	facetOutgassingText->Clear();
	facetOutgassingText->SetEditable(false);
	facetOutgPerAreaText->Clear();
	facetOutgPerAreaText->SetEditable(false);
	facetAreaText->SetEditable(false);
	facetAreaText->Clear();
	facetPumpingSpeedText->SetEditable(false);
	facetPumpingSpeedText->Clear();
	facetOpacity->Clear();
	facetOpacity->SetEditable(false);
	facetTemperatureText->Clear();
	facetTemperatureText->SetEditable(false);
	facetSideType->SetSelectedValue("");
	facetSideType->SetEditable(false);
	facetDesType->SetSelectedValue("");
	facetDesType->SetEditable(false);
	facetDesTypeN->SetText("");
	facetDesTypeN->SetEditable(false);
	facetProfileCombo->SetSelectedValue("");
	facetProfileCombo->SetEditable(false);
}

// Name: ApplyFacetParams()
// Desc: Apply facet parameters.

void MolFlow::ApplyFacetParams() {

	InterfaceGeometry* interfGeom = worker.GetGeometry();
	size_t nbFacet = interfGeom->GetNbFacet();

	// Sticking
	double sticking;
	bool stickingNotNumber;
	bool doSticking = false;
	if (facetStickingText->GetNumber(&sticking)) {
		if (sticking < 0.0 || sticking>1.0) {
			GLMessageBox::Display("Sticking must be in the range [0,1]", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}
		doSticking = true;
		stickingNotNumber = false;
	}
	else {
		if (facetStickingText->GetText() == "...") doSticking = false;
		else {/*
			GLMessageBox::Display("Invalid sticking number","Error",GLDLG_OK,GLDLG_ICONERROR);
			UpdateFacetParams();
			return;*/
			doSticking = true;
			stickingNotNumber = true;
		}
	}

	// opacity
	double opacity;
	bool doOpacity = false;
	bool opacityNotNumber;
	if (facetOpacity->GetNumber(&opacity)) {
		if (opacity < 0.0 || opacity>1.0) {
			GLMessageBox::Display("Opacity must be in the range [0,1]", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}
		doOpacity = true;
		opacityNotNumber = false;
	}
	else {
		if (facetOpacity->GetText() == "...") doOpacity = false;
		else {
			doOpacity = true;
			opacityNotNumber = true;
		}
	}

	// temperature
	double temperature;
	bool doTemperature = false;
	bool temperatureNotNumber;
	if (facetTemperatureText->GetNumber(&temperature)) {
		if (temperature <= 0.0) {
			GLMessageBox::Display("Temperature must be positive", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}
		doTemperature = true;
		temperatureNotNumber = false;
	}
	else {
		if (facetTemperatureText->GetText() == "...") doTemperature = false;
		else {
			doTemperature = true;
			temperatureNotNumber = true;
		}
	}

	// Outgassing
	double outgassing = 0.0;
	bool doFlow = false;
	bool outgassingNotNumber;
	//Calculate outgassing
	if (facetOutgToggleLabel->GetState() && facetOutgassingText->GetText() != "..." && facetDesType->GetSelectedIndex() != 0
		&& facetDesType->GetSelectedValue() != "..." && facetOutgassingText->IsEditable()) {  //We want outgassing
		if (facetOutgassingText->GetNumber(&outgassing)) { //If we can parse the number
			if (outgassing <= 0.0) {
				GLMessageBox::Display("Outgassing must be positive", "Error", GLDLG_OK, GLDLG_ICONERROR);
				return;
			}
			doFlow = true;
			outgassingNotNumber = false;
		}
		else { //could not parse as number
			doFlow = true;
			outgassingNotNumber = true;
		}
	}

	// Outgassing per area
	double flowA = 0;
	bool doFlowA = false;
	//Calculate outgassing

	if (facetOutgPerAreaToggleLabel->GetState() && facetOutgPerAreaText->GetText() != "..."
		&& facetDesType->GetSelectedIndex() != 0 && facetDesType->GetSelectedValue() != "..." && facetOutgPerAreaText->IsEditable()) { //We want outgassing per area
		if (facetOutgPerAreaText->GetNumber(&flowA)) { //Can be parsed as number
			if (flowA <= 0.0) {
				GLMessageBox::Display("Outgassing per area must be positive", "Error", GLDLG_OK, GLDLG_ICONERROR);
				return;
			}
			doFlowA = true;
		}
		else {
			GLMessageBox::Display("Invalid outgassing per area number", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}
	}

	// Desorption type
	int desorbType = facetDesType->GetSelectedIndex();

	double desorbTypeN;
	bool doDesorbTypeN = false;
	if (desorbType == 3) {
		if (facetDesTypeN->GetNumber(&desorbTypeN)) {
			if (desorbTypeN <= 0.0) {
				GLMessageBox::Display("Desorption type exponent must be greater than 0.0", "Error", GLDLG_OK, GLDLG_ICONERROR);
				UpdateFacetParams(false);
				return;
			}
			doDesorbTypeN = true;
		}
		else {
			if (facetDesTypeN->GetText() == "...") doDesorbTypeN = false;
			else {
				GLMessageBox::Display("Invalid desorption type exponent", "Error", GLDLG_OK, GLDLG_ICONERROR);
				return;
			}
		}
	}

	// Record (profile) type
	int rType = facetProfileCombo->GetSelectedIndex(); // -1 if "..."

	// 2sided
	int is2Sided = facetSideType->GetSelectedIndex();

	if (desorbType == DES_ANGLEMAP && facetAdvParams->IsAngleMapRecording()) {
		char errMsg[512];
		sprintf(errMsg, "Angle map can't be used for DESORPTION and RECORDING\n");
		GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}

	//Check complete, let's apply
	if (facetAdvParams && facetAdvParams->IsVisible()) {
		if (!facetAdvParams->Apply()) {
			return;
		}
	}

	if (!AskToReset()) return;

	changedSinceSave = true;

	// Update facets (local)
	for (int i = 0; i < nbFacet; i++) {
		InterfaceFacet* f = interfGeom->GetFacet(i);
		if (f->selected) {
			if (doSticking) {
				if (!stickingNotNumber) {
					f->sh.sticking = sticking;
					f->sh.stickingParam = "";
				}
				else {
					f->sh.stickingParam = facetStickingText->GetText();
				}
			}

			if (doOpacity) {
				if (!opacityNotNumber) {
					f->sh.opacity = opacity;
					f->sh.opacityParam = "";
				}
				else {
					f->sh.opacityParam = facetOpacity->GetText();
				}
			}
			if (doTemperature) {
				if (!temperatureNotNumber) {
					f->sh.temperature = temperature;
					f->sh.temperatureParam = "";
				}
				else {
					f->sh.temperatureParam = facetTemperatureText->GetText();
				}
			}
			if (doFlow) {
				if (!outgassingNotNumber) {
					f->sh.outgassing = outgassing * MBARLS_TO_PAM3S;
					f->sh.outgassingParam = "";
				}
				else {
					f->sh.outgassingParam = facetOutgassingText->GetText();
				}
			}
			if (doFlowA/* && !useMapA*/) f->sh.outgassing = flowA * f->GetArea() * MBARLS_TO_PAM3S;
			if (desorbType >= 0) {
				if (desorbType == 0) f->sh.outgassing = 0.0;
				if (desorbType != 3) f->sh.desorbTypeN = 0.0;
				f->sh.desorbType = desorbType;
				if (doDesorbTypeN) f->sh.desorbTypeN = desorbTypeN;
			}

			if (rType >= 0) {
				f->sh.profileType = rType;
				//f->sp.isProfile = (rType!=PROFILE_NONE); //included below by f->UpdateFlags();
			}
			if (is2Sided >= 0) f->sh.is2sided = is2Sided;

			f->sh.maxSpeed = 4.0 * sqrt(2.0 * 8.31 * f->sh.temperature / 0.001 / worker.model->sp.gasMass);
			f->UpdateFlags();
		}
	}

	// Mark "needsReload" to sync changes with workers on next simulation start
	worker.MarkToReload();

	//worker.CalcTotalOutgassing();
	UpdateFacetParams(false);
	if (profilePlotter) profilePlotter->Refresh();
	if (pressureEvolution) pressureEvolution->Refresh();
	if (timewisePlotter) timewisePlotter->Refresh();
	//if (facetAdvParams) facetAdvParams->Refresh();
}

// Name: UpdateFacetParams()
// Desc: Update selected facet parameters.

void MolFlow::UpdateFacetParams(bool updateSelection) { //Calls facetAdvParams->Refresh()

	char tmp[256];

	// Update params
	InterfaceGeometry* interfGeom = worker.GetGeometry();
	// Get list of selected facet
	auto selectedFacets = interfGeom->GetSelectedFacets();
	size_t nbSel = selectedFacets.size();

	if (nbSel > 0) {

		InterfaceFacet* f0;
		InterfaceFacet* f;



		f0 = interfGeom->GetFacet(selectedFacets[0]);

		double f0Area = f0->GetArea();
		double sumArea = f0Area; //sum facet area

		bool stickingE = true;
		bool opacityE = true;
		bool temperatureE = true;
		bool outgassingE = true;
		bool outgassingPerAreaE = true;
		bool desorbTypeE = true;
		bool desorbTypeNE = true;
		bool recordE = true;
		bool is2sidedE = true;

		for (size_t sel = 1; sel < selectedFacets.size(); sel++) {
			f = interfGeom->GetFacet(selectedFacets[sel]);
			double fArea = f->GetArea();
			stickingE = stickingE && (f0->sh.stickingParam == f->sh.stickingParam) && (!f0->sh.stickingParam.empty() || IsEqual(f0->sh.sticking, f->sh.sticking));
			opacityE = opacityE && (f0->sh.opacityParam == f->sh.opacityParam) && (!f0->sh.opacityParam.empty() || IsEqual(f0->sh.opacity, f->sh.opacity));
			temperatureE = temperatureE && (f0->sh.temperatureParam == f->sh.temperatureParam) && (!f0->sh.temperatureParam.empty() || IsEqual(f0->sh.temperature, f->sh.temperature));
			outgassingE = outgassingE && (f0->sh.outgassingParam == f->sh.outgassingParam) && (!f0->sh.outgassingParam.empty() || IsEqual(f0->sh.outgassing, f->sh.outgassing));
			outgassingPerAreaE = outgassingPerAreaE && IsEqual(f0->sh.outgassing / f0Area, f->sh.outgassing / fArea);
			is2sidedE = is2sidedE && (f0->sh.is2sided == f->sh.is2sided);
			desorbTypeE = desorbTypeE && (f0->sh.desorbType == f->sh.desorbType);
			desorbTypeNE = desorbTypeNE && IsEqual(f0->sh.desorbTypeN, f->sh.desorbTypeN);
			recordE = recordE && (f0->sh.profileType == f->sh.profileType);  //profiles
			sumArea += fArea;
		}

		if (nbSel == 1)
			sprintf(tmp, "Selected Facet (#%zd)", selectedFacets[0] + 1);
		else
			sprintf(tmp, "Selected Facets (%zd selected)", selectedFacets.size());

		facetPanel->SetTitle(tmp);
		if (selectedFacets.size() > 1) facetAreaLabel->SetText("Sum Area (cm\262):");
		else facetAreaLabel->SetText("Area (cm\262):");
		facetAreaText->SetText(sumArea);
		if (stickingE) {
			if (f0->sh.stickingParam.empty())
				facetStickingText->SetText(f0->sh.sticking);
			else facetStickingText->SetText(f0->sh.stickingParam);
		}
		else facetStickingText->SetText("...");

		if (opacityE) {
			if (f0->sh.opacityParam.empty())
				facetOpacity->SetText(f0->sh.opacity);
			else facetOpacity->SetText(f0->sh.opacityParam);
		}
		else facetOpacity->SetText("...");

		if (temperatureE) {
			if (f0->sh.temperatureParam.empty())
				facetTemperatureText->SetText(f0->sh.temperature);
			else facetTemperatureText->SetText(f0->sh.temperatureParam);
		}
		else facetTemperatureText->SetText("...");

		if (is2sidedE) facetSideType->SetSelectedIndex(f0->sh.is2sided); else facetSideType->SetSelectedValue("...");
		if (desorbTypeNE) facetDesTypeN->SetText(f0->sh.desorbTypeN); else facetDesTypeN->SetText("...");
		if (recordE) facetProfileCombo->SetSelectedIndex(f0->sh.profileType); else facetProfileCombo->SetSelectedValue("...");

		if (selectedFacets.size() == 1) {
			facetPumpingSpeedText->SetEditable(true);
			CalcPumpingSpeed();
		}
		else {
			facetPumpingSpeedText->SetEditable(false);
			facetPumpingSpeedText->SetText("...");
		}

		if (desorbTypeE) {
			facetDesType->SetSelectedIndex(f0->sh.desorbType);
			if (f0->sh.desorbType > DES_NONE) { //There is some desorption

				//facetOutgassingText->SetEnabled(true);
				facetOutgPerAreaToggleLabel->SetEnabled(true);
				facetOutgPerAreaText->SetEditable(true);
				facetOutgToggleLabel->SetEnabled(true);
				facetOutgassingText->SetEditable(true);
				if (outgassingE) {
					if (f0->sh.outgassingParam.empty())
						facetOutgassingText->SetText(f0->sh.outgassing * PAM3S_TO_MBARLS);
					else facetOutgassingText->SetText(f0->sh.outgassingParam);
				}
				else facetOutgassingText->SetText("...");
				if (outgassingPerAreaE) facetOutgPerAreaText->SetText(f0->sh.outgassing / f0Area * 10.00); else facetOutgPerAreaText->SetText("...");
				if (f0->sh.desorbType == 3) {
					facetDesTypeN->SetEditable(true);
					if (desorbTypeNE) facetDesTypeN->SetText(f0->sh.desorbTypeN); else facetDesTypeN->SetText("...");
				}
				else {
					facetDesTypeN->SetText("");
					facetDesTypeN->SetEditable(false);
				}

			}
			else { //No desorption
				facetOutgToggleLabel->SetEnabled(false);
				facetOutgassingText->SetEditable(false);
				facetOutgPerAreaToggleLabel->SetEnabled(false);
				facetOutgPerAreaText->SetEditable(false);
				facetDesTypeN->SetText("");
				facetDesTypeN->SetEditable(false);
				facetOutgassingText->SetText("");
				facetOutgPerAreaText->SetText("");
			}
		}
		else { //Mixed state
			facetDesType->SetSelectedValue("...");
			facetOutgToggleLabel->SetEnabled(false);
			facetOutgassingText->SetEditable(false);
			facetOutgPerAreaToggleLabel->SetEnabled(false);
			facetOutgPerAreaText->SetEditable(false);
			facetDesTypeN->SetText("");
			facetDesTypeN->SetEditable(false);
			facetOutgassingText->SetText("");
			facetOutgPerAreaText->SetText("");
		}

		if (facetAdvParams) facetAdvParams->Refresh(selectedFacets); //Refresh advanced facet parameters panel
		if (updateSelection) {

			if (nbSel > 1000 || interfGeom->GetNbFacet() > 50000) { //If it would take too much time to look up every selected facet in the list
				facetList->ReOrder();
				facetList->SetSelectedRows(selectedFacets, false);
			}
			else {
				facetList->SetSelectedRows(selectedFacets, true);
			}
			facetList->lastRowSel = -1;
		}

		//Enabled->Editable
		facetStickingText->SetEditable(true);
		facetOpacity->SetEditable(true);
		facetTemperatureText->SetEditable(true);
		facetSideType->SetEditable(true);
		facetDesType->SetEditable(true);
		facetProfileCombo->SetEditable(true);
		facetApplyBtn->SetEnabled(false);
	}
	else {
		ClearFacetParams();
		if (facetAdvParams) facetAdvParams->Refresh(std::vector<size_t>()); //Clear
		if (updateSelection) facetList->ClearSelection();
	}

	if (facetDetails) facetDetails->Update();
	if (facetCoordinates) facetCoordinates->UpdateFromSelection();
	if (texturePlotter) texturePlotter->Update(m_fTime, true); //Facet change
	if (outgassingMapWindow) outgassingMapWindow->Update(m_fTime, true);
	if (histogramSettings) histogramSettings->Refresh(selectedFacets);
}

// Name: FrameMove()
// Desc: Called once per frame, the call is the entry point for animating
//       the scene.

int MolFlow::FrameMove()
{
	bool runningState = worker.IsRunning();
	double elapsedTime = worker.simuTimer.Elapsed();
	if (runningState && ((m_fTime - lastUpdate) >= 1.0f)) {
		if (textureScaling) textureScaling->Update();
		//if (formulaEditor && formulaEditor->IsVisible()) formulaEditor->Refresh(); //Interface::Framemove does it already
	}
	Interface::FrameMove(); //might reset lastupdate
	char tmp[256];
	if (globalSettings) globalSettings->SMPUpdate();

	if ((elapsedTime <= 2.0f) && runningState) {
		hitNumber->SetText("Starting...");
		desNumber->SetText("Starting...");
	}
	else {
		double current_avg = hps.avg();
		if (!runningState)
			current_avg = hps_runtotal.avg();
		else
			current_avg = (current_avg != 0.0) ? current_avg : hps.last();

		sprintf(tmp, "%s (%s)", Util::formatInt(worker.globalStatCache.globalHits.nbMCHit, "hit"),
			Util::formatPs(current_avg, "hit"));
		hitNumber->SetText(tmp);

		current_avg = dps.avg();
		if (!runningState) current_avg = dps_runtotal.avg();
		else current_avg = (current_avg != 0.0) ? current_avg : dps.last();

		sprintf(tmp, "%s (%s)", Util::formatInt(worker.globalStatCache.globalHits.nbDesorbed, "des"),
			Util::formatPs(current_avg, "des"));
		desNumber->SetText(tmp);
	}

	//sprintf(tmp, "Running: %s", Util::formatTime(worker.simuTimer.Elapsed()));


	// Save previous state to react to changes
	prevRunningState = runningState;

	return GL_OK;
}

// Name: RestoreDeviceObjects()
// Desc: Initialize scene objects.

int MolFlow::RestoreDeviceObjects()
{
	RestoreDeviceObjects_shared();

	//Different Molflow implementations:
	RVALIDATE_DLG(facetAdvParams);
	RVALIDATE_DLG(facetDetails);
	RVALIDATE_DLG(smartSelection);
	RVALIDATE_DLG(viewer3DSettings);
	RVALIDATE_DLG(textureScaling);
	RVALIDATE_DLG(globalSettings);
	RVALIDATE_DLG(profilePlotter);
	RVALIDATE_DLG(texturePlotter);

	//Molflow only:
	RVALIDATE_DLG(importDesorption);
	RVALIDATE_DLG(timeSettings);
	RVALIDATE_DLG(movement);
	RVALIDATE_DLG(measureForces);
	RVALIDATE_DLG(outgassingMapWindow);
	RVALIDATE_DLG(parameterEditor);
	RVALIDATE_DLG(pressureEvolution);
	RVALIDATE_DLG(timewisePlotter);

	return GL_OK;
}

// Name: InvalidateDeviceObjects()
// Desc: Free all alocated resource

int MolFlow::InvalidateDeviceObjects()
{
	InvalidateDeviceObjects_shared();

	//Different Molflow implementations:
	IVALIDATE_DLG(facetAdvParams);
	IVALIDATE_DLG(facetDetails);
	IVALIDATE_DLG(smartSelection);
	IVALIDATE_DLG(viewer3DSettings);
	IVALIDATE_DLG(textureScaling);
	IVALIDATE_DLG(globalSettings);
	IVALIDATE_DLG(profilePlotter);
	IVALIDATE_DLG(texturePlotter);

	//Molflow only:
	IVALIDATE_DLG(importDesorption);
	IVALIDATE_DLG(timeSettings);
	IVALIDATE_DLG(movement);
	IVALIDATE_DLG(measureForces);
	IVALIDATE_DLG(outgassingMapWindow);
	IVALIDATE_DLG(parameterEditor);
	IVALIDATE_DLG(pressureEvolution);
	IVALIDATE_DLG(timewisePlotter);

	return GL_OK;
}

void MolFlow::ExportProfiles() {

	InterfaceGeometry* interfGeom = worker.GetGeometry();
	if (interfGeom->GetNbSelectedFacets() == 0) {
		GLMessageBox::Display("Empty selection", "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}

	std::string saveFile = NFD_SaveFile_Cpp(fileProfFilters, "");

	if (!saveFile.empty()) {
		try {
			worker.ExportProfiles(saveFile.c_str());
		}
		catch (const std::exception& e) {
			char errMsg[512];
			sprintf(errMsg, "%s\nFile:%s", e.what(), saveFile.c_str());
			GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
		}
	}
}

void MolFlow::ExportAngleMaps() {

	std::string profFile = NFD_SaveFile_Cpp(fileProfFilters, "");

	if (!profFile.empty()) {
		try {
			auto retVal = worker.ExportAngleMaps(profFile, false);
			if (!retVal) return; //user cancel or error
			std::vector<std::string> exportList = *retVal; //vector of fileNames, might be empty
			if (exportList.empty()) {
				GLMessageBox::Display("Select at least one facet with recorded angle map", "Error", GLDLG_OK, GLDLG_ICONERROR);
				return;
			}
		}
		catch (const std::exception& e) {
			char errMsg[512];
			sprintf(errMsg, "%s\nFile:%s", e.what(), profFile.c_str());
			GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}
	}
}

void MolFlow::ImportAngleMaps() {
	auto selFacets = worker.GetGeometry()->GetSelectedFacets();
	if (selFacets.empty()) {
		GLMessageBox::Display("Select at least one facet to import angle map to", "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}

	//std::vector<FILENAME> files = GLFileBox::OpenMultipleFiles("CSV files\0*.csv\0All files\0*.*\0", "Import angle map(s)");
	auto fileNames = NFD_OpenMultiple_Cpp("csv", "");
	if (fileNames.empty()) return;
	if (selFacets.size() != fileNames.size()) {
		GLMessageBox::Display("Select the same number of facets and files to import", "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}
	for (size_t i = 0; i < fileNames.size(); i++) {
		try {
			auto file = FileReader(fileNames[i]);
			std::vector<std::vector<std::string>> table = file.ImportCSV_string();
			AskToReset(&worker);
			worker.GetGeometry()->GetFacet(selFacets[i])->ImportAngleMap(table);
			worker.MarkToReload();
		}
		catch (const std::exception& e) {
			char errMsg[512];
			sprintf(errMsg, "%s\nFile:%s", e.what(), fileNames[i].c_str());
			GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}
	}
}

void MolFlow::CopyAngleMapToClipboard()
{
	InterfaceGeometry* interfGeom = worker.GetGeometry();
	size_t angleMapFacetIndex;
	bool found = false;
	for (size_t i = 0; i < interfGeom->GetNbFacet(); i++) {
		InterfaceFacet* f = interfGeom->GetFacet(i);
		if (f->selected && !f->angleMapCache.empty()) {
			if (found) {
				GLMessageBox::Display("More than one facet with recorded angle map selected", "Error", GLDLG_OK, GLDLG_ICONERROR);
				return;
			}
			angleMapFacetIndex = i;
			found = true;
		}
	}
	if (!found) {
		GLMessageBox::Display("No facet with recorded angle map selected", "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}

	/*
	if (!worker.IsDpInitialized()) {
		GLMessageBox::Display("Worker Dataport not initialized yet", "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}*/

	try {
		std::string map = interfGeom->GetFacet(angleMapFacetIndex)->GetAngleMap(2); //Clipboard format: tab-separated (format==2)
		if (map.length() > (10 * 1024 * 1024)) {
			if (GLMessageBox::Display("Angle map text over 10MB. Copy to clipboard?", "Warning", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONWARNING) != GLDLG_OK)
				return;
		}

		GLToolkit::CopyTextToClipboard(map);

	}
	catch (const std::exception& e) {
		GLMessageBox::Display(e.what(), "Error", GLDLG_OK, GLDLG_ICONERROR);
	}
}

void MolFlow::ClearAngleMapsOnSelection() {
	InterfaceGeometry* interfGeom = worker.GetGeometry();
	for (size_t i = 0; i < interfGeom->GetNbFacet(); i++) {
		InterfaceFacet* f = interfGeom->GetFacet(i);
		if (f->selected && !f->angleMapCache.empty()) {
			f->angleMapCache.clear();
		}
	}
}

/*
void MolFlow::ImportDesorption_DES() {

	FILENAME *fn = GLFileBox::OpenFile(currentDir, nullptr, "Import desorption File", fileDesFilters, 0);

	if (fn) {

		try {
			worker.ImportDesorption_DES(fn->fullName);
		}
		catch(const std::exception &e) {
			char errMsg[512];
			sprintf(errMsg, "%s\nFile:%s", e.what(), fn->fullName);
			GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
		}
		//worker.CalcTotalOutgassing();
		UpdateFacetParams();
	}

}
*/

void MolFlow::SaveFile() {
	if (!worker.fullFileName.empty()) {

		auto prg = GLProgress_GUI("Saving file...\nIn this beta version, you can see the progress in the console.", "Please wait");
		prg.SetVisible(true);
		prg.SetProgress(0.5);

		//GLWindowManager::Repaint();

		try {
			worker.SaveGeometry(worker.fullFileName, prg, false);
			ResetAutoSaveTimer();
		}
		catch (const std::exception& e) {
			char errMsg[512];
			sprintf(errMsg, "%s\nFile:%s", e.what(), worker.GetCurrentFileName().c_str());
			GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
		}
		changedSinceSave = false;

	}
	else SaveFileAs();
}

void MolFlow::LoadFile(const std::string& fileName) {

	std::string fileShortName, filePath;

	if (fileName.empty()) {
		filePath = NFD_OpenFile_Cpp(fileLoadFilters, "");
	}
	else {
		filePath = fileName;
	}

	if (filePath.empty()) return; //User closed Open... dialog

	if (!FileUtils::Exist(filePath)) {
		auto answer = GLMessageBox::Display(
			fmt::format("{}\nDoesn't exist. Remove from the Recent files menu?", filePath),
			"No such file",
			{ "Yes","No" },
			GLDLG_ICONERROR);
		if (answer == 0) { //"Yes"
			RemoveRecent(filePath.c_str());
		}
		return;
	}

	auto prg = GLProgress_GUI("Preparing to load file...", "Please wait");
	prg.SetVisible(true);
	//GLWindowManager::Repaint();

	fileShortName = FileUtils::GetFilename(filePath);

	try {
		ClearFacetParams();
		ClearFormulas();
		TimeDependentParameters::ClearParameters(worker.interfaceParameterCache);
		ClearAllSelections();
		ClearAllViews();
		ResetSimulation(false);
		
		SetDefaultViews();
		worker.LoadGeometry(filePath);

		InterfaceGeometry* interfGeom = worker.GetGeometry();


		// Default initialisation
		for (auto& view : viewers) {
			view->SetWorker(&worker);
		}

		//UpdateModelParams();
		startSimu->SetEnabled(true);
		//compACBtn->SetEnabled(modeCombo->GetSelectedIndex() == 1);
		//singleACBtn->SetEnabled(modeCombo->GetSelectedIndex() == 1);
		//resetSimu->SetEnabled(true);
		
	nbDesStart = worker.globalStatCache.globalHits.nbDesorbed;
		nbHitStart = worker.globalStatCache.globalHits.nbMCHit;
			AddRecent(filePath);
		interfGeom->viewStruct = -1;

		UpdateStructMenu();

		// Check non simple polygon
		prg.SetMessage("Checking for non simple polygons...");

		interfGeom->CheckCollinear();
		//interfGeom->CheckNonSimple();
		interfGeom->CheckIsolatedVertex();
		
		ResetAutoSaveTimer();
		//UpdatePlotters();

		if (timeSettings) timeSettings->RefreshMoments();
		if (momentsEditor) momentsEditor->Refresh();
		if (pressureEvolution) pressureEvolution->Reset();
		if (timewisePlotter) timewisePlotter->Refresh();
		if (histogramPlotter) histogramPlotter->Reset();
		if (profilePlotter) profilePlotter->Refresh();
		if (convergencePlotter) convergencePlotter->Refresh();
		if (texturePlotter) texturePlotter->Update(0.0, true);
		//if (parameterEditor) parameterEditor->UpdateCombo(); //Done by ClearParameters()
		if (textureScaling) textureScaling->Update();
		if (outgassingMapWindow) outgassingMapWindow->Update(m_fTime, true);
		if (facetDetails) facetDetails->Update();
		if (facetCoordinates) facetCoordinates->UpdateFromSelection();
		if (vertexCoordinates) vertexCoordinates->Update();
		if (movement) movement->Update();
		if (measureForces) measureForces->Update();
		if (globalSettings && globalSettings->IsVisible()) globalSettings->Update();
		if (formulaEditor) formulaEditor->Refresh();
		if (parameterEditor) parameterEditor->Refresh();
	}
	catch (const std::exception& e) {

		char errMsg[512];
		sprintf(errMsg, "%s\nFile:%s", e.what(), fileShortName.c_str());
		GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
		RemoveRecent(filePath.c_str());

	}
	changedSinceSave = false;
}

void MolFlow::InsertGeometry(bool newStr, const std::string& fileName) {


	std::string fileShortName, filePath;

	if (fileName.empty()) {
		filePath = NFD_OpenFile_Cpp(fileInsertFilters, "");
	}
	else {
		filePath = fileName;
	}


	auto prg = GLProgress_GUI("Preparing to load file...", "Please wait");
	prg.SetVisible(true);
	//GLWindowManager::Repaint();

	if (filePath.empty()) {
		return;
	}

	if (!AskToReset()) return;
	ResetSimulation(false);

	fileShortName = FileUtils::GetFilename(filePath);

	try {

		//worker.InsertGeometry(newStr, fullName);
		worker.LoadGeometry(filePath, true, newStr);

		InterfaceGeometry* interfGeom = worker.GetGeometry();

		//Increase BB
		for (auto& view : viewers)
			view->SetWorker(&worker);

		//UpdateModelParams();
		startSimu->SetEnabled(true);

		//compACBtn->SetEnabled(modeCombo->GetSelectedIndex() == 1);
		//singleACBtn->SetEnabled(modeCombo->GetSelectedIndex() == 1);
		//resetSimu->SetEnabled(true);
		//ClearFacetParams();
		//nbDesStart = worker.globalState->globalStats.globalStats.hit.nbDesorbed;
		//nbHitStart = worker.globalState->globalStats.globalStats.hit.nbMC;
		AddRecent(filePath);
		interfGeom->viewStruct = -1;

		//worker.LoadTexturesGEO(fullName);
		UpdateStructMenu();

		interfGeom->CheckCollinear();
		//interfGeom->CheckNonSimple();
		interfGeom->CheckIsolatedVertex();

		/*
		// Set up view
		// Default
		viewers[0]->SetProjection(ProjectionMode::Orthographic);
		viewers[0]->ToFrontView();
		viewers[1]->SetProjection(ProjectionMode::Orthographic);
		viewers[1]->ToTopView();
		viewers[2]->SetProjection(ProjectionMode::Orthographic);
		viewers[2]->ToSideView();
		viewers[3]->SetProjection(ProjectionMode::Perspective);
		viewers[3]->ToFrontView();
		SelectViewer(0);
		*/
		RefreshPlotterCombos();
		//UpdatePlotters();

		if (outgassingMapWindow) outgassingMapWindow->Update(m_fTime, true);
		if (facetDetails) facetDetails->Update();
		if (facetCoordinates) facetCoordinates->UpdateFromSelection();
		if (vertexCoordinates) vertexCoordinates->Update();
		if (formulaEditor) formulaEditor->Refresh();
		if (parameterEditor) parameterEditor->Refresh();
	}
	catch (const std::exception& e) {

		char errMsg[512];
		sprintf(errMsg, "%s\nFile:%s", e.what(), fileShortName.c_str());
		GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
		RemoveRecent(filePath.c_str());

	}
	changedSinceSave = true;
}

void MolFlow::StartStopSimulation() {

	if (worker.globalStatCache.globalHits.nbMCHit <= 0 && !worker.model->sp.calcConstantFlow && worker.interfaceMomentCache.empty()) {
		bool ok = GLMessageBox::Display("Warning: in the Moments Editor, the option \"Calculate constant flow\" is disabled.\n"
			"This is useful for time-dependent simulations.\n"
			"However, you didn't define any moments, suggesting you're using steady-state mode.\n"
			"\nDo you want to continue?\n", "Strange time settings", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONWARNING) == GLDLG_OK;
		if (!ok) return;
	}

	worker.StartStop(m_fTime);
	if (!worker.IsRunning()) { //Force update on simulation stop
		appFormulas->EvaluateFormulas(worker.globalStatCache.globalHits.nbDesorbed);
		UpdatePlotters();
		//if (autoUpdateFormulas) UpdateFormula();
		if (autoUpdateFormulas && formulaEditor && formulaEditor->IsVisible()) formulaEditor->UpdateValues();
		if (particleLogger && particleLogger->IsVisible()) particleLogger->UpdateStatus();
	}

	// Frame rate measurement
	// reset on start only
	lastMeasTime = m_fTime;
	if (worker.IsRunning()) {
		hps.clear();
		dps.clear();
		hps_runtotal.clear();
		dps_runtotal.clear();
	}
	lastNbHit = worker.globalStatCache.globalHits.nbMCHit;
	lastNbDes = worker.globalStatCache.globalHits.nbDesorbed;

	hps_runtotal.push(lastNbHit, lastMeasTime);
	dps_runtotal.push(lastNbDes, lastMeasTime);
	lastUpdate = 0.0;
}

// Name: EventProc()
// Desc: Message proc function to handle key and mouse input

void MolFlow::ProcessMessage(GLComponent* src, int message)
{

	if (ProcessMessage_shared(src, message)) return; //Already processed by common interface

	InterfaceGeometry* interfGeom = worker.GetGeometry();
	switch (message) {

		//MENU --------------------------------------------------------------------
	case MSG_MENU:
		switch (src->GetId()) {

		case MENU_FILE_IMPORTDES_SYN:

			if (interfGeom->IsLoaded()) {
				if (!importDesorption) importDesorption = new ImportDesorption();
				importDesorption->SetGeometry(interfGeom, &worker);
				importDesorption->SetVisible(true);
			}
			else GLMessageBox::Display("No geometry loaded.", "No geometry", GLDLG_OK, GLDLG_ICONERROR);
			break;
			/* //Deprecated
			case MENU_FILE_IMPORTDES_DES:
				ImportDesorption_DES();
				break;
			*/

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
			ExportTextures(0, 7); break;
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

		case MENU_TOOLS_MOVINGPARTS:
			if (!movement) movement = new Movement(interfGeom, &worker);
			movement->Update();
			movement->SetVisible(true);
			break;

		case MENU_TOOLS_MEASUREFORCE:
			if (!measureForces) measureForces = new MeasureForce(interfGeom, &worker);
			measureForces->Update();
			measureForces->SetVisible(true);
			break;

		case MENU_EDIT_TSCALING:
			if (!textureScaling || !textureScaling->IsVisible()) {
				SAFE_DELETE(textureScaling);
				textureScaling = new TextureScaling(&worker, viewers);
				textureScaling->Display();
			}
			break;

		case MENU_EDIT_GLOBALSETTINGS:
			if (!globalSettings) globalSettings = new GlobalSettings(&worker);
			globalSettings->Update();
			globalSettings->SetVisible(true);
			break;

		case MENU_TOOLS_PROFPLOTTER:
			if (!profilePlotter) profilePlotter = new ProfilePlotter(&worker);
			profilePlotter->Display(&worker);
			break;
		case MENU_TOOLS_TEXPLOTTER:
			if (!texturePlotter) texturePlotter = new TexturePlotter();
			texturePlotter->Display(&worker);
			break;
		case MENU_FACET_DETAILS:
			if (facetDetails == nullptr) facetDetails = new FacetDetails();
			facetDetails->Display(&worker);
			break;

		case MENU_TIME_PRESSUREEVOLUTION:
			if (!pressureEvolution) pressureEvolution = new PressureEvolution(&worker);
			pressureEvolution->SetVisible(true);
			break;

		case MENU_FACET_OUTGASSINGMAP:
			if (!outgassingMapWindow) outgassingMapWindow = new OutgassingMapWindow();
			outgassingMapWindow->Display(&worker);
			break;
		case MENU_FACET_REMOVESEL:
		{
			auto selectedFacets = interfGeom->GetSelectedFacets();
			if (selectedFacets.empty()) return; //Nothing selected
			if (GLMessageBox::Display("Remove selected facets?", "Question", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONINFO) == GLDLG_OK) {
				if (AskToReset()) {
					if (worker.IsRunning()) worker.Stop_Public();
					interfGeom->RemoveFacets(selectedFacets);
					//worker.CalcTotalOutgassing();
					//interfGeom->CheckIsolatedVertex();
					UpdateModelParams();
					RefreshPlotterCombos();
					//UpdatePlotters();
					if (vertexCoordinates) vertexCoordinates->Update();
					if (facetCoordinates) facetCoordinates->UpdateFromSelection();
					// Send to sub process
					worker.MarkToReload();
				}
			}
			break;
		}
		case MENU_FACET_SELECTSTICK:
			interfGeom->UnselectAll();
			for (int i = 0; i < interfGeom->GetNbFacet(); i++)
				if (!interfGeom->GetFacet(i)->sh.stickingParam.empty() || (interfGeom->GetFacet(i)->sh.sticking != 0.0 && !interfGeom->GetFacet(i)->IsTXTLinkFacet()))
					interfGeom->SelectFacet(i);
			interfGeom->UpdateSelection();
			UpdateFacetParams(true);
			break;

		case MENU_FACET_SELECTREFL:
			interfGeom->UnselectAll();
			for (int i = 0; i < interfGeom->GetNbFacet(); i++) {
				InterfaceFacet* f = interfGeom->GetFacet(i);
				if (f->sh.desorbType == DES_NONE && f->sh.sticking == 0.0 && f->sh.opacity > 0.0)
					interfGeom->SelectFacet(i);
			}
			interfGeom->UpdateSelection();
			UpdateFacetParams(true);
			break;

		case MENU_FACET_SELECTDES:
			interfGeom->UnselectAll();
			for (int i = 0; i < interfGeom->GetNbFacet(); i++)
				if (interfGeom->GetFacet(i)->sh.desorbType != DES_NONE)
					interfGeom->SelectFacet(i);
			interfGeom->UpdateSelection();
			UpdateFacetParams(true);
			break;
		case MENU_SELECT_HASDESFILE:
			interfGeom->UnselectAll();
			for (int i = 0; i < interfGeom->GetNbFacet(); i++)
				if (interfGeom->GetFacet(i)->hasOutgassingFile)
					interfGeom->SelectFacet(i);
			interfGeom->UpdateSelection();
			UpdateFacetParams(true);
			break;

		case MENU_TIME_SETTINGS:
			if (!timeSettings) timeSettings = new TimeSettings(&worker);
			timeSettings->SetVisible(true);
			break;
		case MENU_TIME_MOMENTS_EDITOR:
			if (!momentsEditor || !momentsEditor->IsVisible()) {
				SAFE_DELETE(momentsEditor);
				momentsEditor = new MomentsEditor(&worker);
				momentsEditor->Refresh();
				momentsEditor->SetVisible(true);
			}
			break;
		case MENU_TIME_PARAMETER_EDITOR:
			if (!parameterEditor) parameterEditor = new ParameterEditor(&worker);
			parameterEditor->SetVisible(true);
			break;
		case MENU_TIMEWISE_PLOTTER:
			if (!timewisePlotter) timewisePlotter = new TimewisePlotter();
			timewisePlotter->Display(&worker);
			break;
		case MENU_VERTEX_REMOVE:
			if (interfGeom->IsLoaded()) {
				if (GLMessageBox::Display("Remove Selected vertices?\nNote: It will also affect facets that contain them!", "Question", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONINFO) == GLDLG_OK) {
					if (AskToReset()) {
						if (worker.IsRunning()) worker.Stop_Public();
						interfGeom->RemoveSelectedVertex();
						//worker.CalcTotalOutgassing();
						interfGeom->Rebuild(); //Will recalculate facet parameters
						UpdateModelParams();
						if (vertexCoordinates) vertexCoordinates->Update();
						if (facetCoordinates) facetCoordinates->UpdateFromSelection();
						if (profilePlotter) profilePlotter->Refresh();
						if (pressureEvolution) pressureEvolution->Refresh();
						if (timewisePlotter) timewisePlotter->Refresh();
						if (facetCoordinates) facetCoordinates->UpdateFromSelection();
						if (vertexCoordinates) vertexCoordinates->Update();
						// Send to sub process
						worker.MarkToReload();
					}
				}

			}
			else GLMessageBox::Display("No geometry loaded.", "No geometry", GLDLG_OK, GLDLG_ICONERROR);
			break;
		}
		break;
		//TEXT --------------------------------------------------------------------
	case MSG_TEXT_UPD:
		if (src == facetStickingText || src == facetTemperatureText) {
			CalcPumpingSpeed();
			facetApplyBtn->SetEnabled(true);
		}
		else if (src == facetOpacity || src == facetDesTypeN) {
			facetApplyBtn->SetEnabled(true);
		}
		else if (src == facetOutgassingText) {
			double outgassing;
			double area;
			facetOutgassingText->GetNumber(&outgassing);
			facetAreaText->GetNumber(&area);
			if (area == 0) facetOutgPerAreaText->SetText("#DIV0");
			else facetOutgPerAreaText->SetText(outgassing / area);
			facetApplyBtn->SetEnabled(true);
			facetOutgToggleLabel->SetState(true);
			facetOutgPerAreaToggleLabel->SetState(false);
			// //else if ( src == facetMass ) {
			//  facetApplyBtn->SetEnabled(true);
		}
		else if (src == facetOutgPerAreaText) {
			double flowPerArea;
			double area;
			facetOutgPerAreaText->GetNumber(&flowPerArea);
			facetAreaText->GetNumber(&area);
			facetOutgassingText->SetText(flowPerArea * area);
			facetApplyBtn->SetEnabled(true);
			facetOutgPerAreaToggleLabel->SetState(true);
			facetOutgToggleLabel->SetState(false);
			// //else if ( src == facetMass ) {
			//  facetApplyBtn->SetEnabled(true);
		}
		else if (src == facetPumpingSpeedText) {
			CalcSticking();
			facetApplyBtn->SetEnabled(true);
			// //else if ( src == facetMass ) {
			//  facetApplyBtn->SetEnabled(true);
		}
		break;

	case MSG_TEXT:
		if (src == facetStickingText || src == facetDesTypeN || src == facetOpacity || src == facetTemperatureText
			|| src == facetPumpingSpeedText || src == facetOutgassingText || src == facetOutgPerAreaText) {
			ApplyFacetParams();
		}
		break;

		//COMBO -------------------------------------------------------------------
	case MSG_COMBO:

		if (src == facetDesType) {
			facetApplyBtn->SetEnabled(true);
			bool hasDesorption = facetDesType->GetSelectedIndex() != 0;

			facetOutgassingText->SetEditable(hasDesorption);
			facetOutgPerAreaText->SetEditable(hasDesorption);

			int color = (hasDesorption) ? 0 : 110;
			facetOutgToggleLabel->SetTextColor(color, color, color);
			facetOutgPerAreaToggleLabel->SetTextColor(color, color, color);
			facetOutgToggleLabel->SetEnabled(hasDesorption);
			facetOutgPerAreaToggleLabel->SetEnabled(hasDesorption);
			facetDesTypeN->SetEditable(facetDesType->GetSelectedIndex() == 3);
		}
		else if (src == facetProfileCombo || src == facetSideType) {
			facetApplyBtn->SetEnabled(true);
		}
		break;

		//TOGGLE ------------------------------------------------------------------
	case MSG_TOGGLE:
		if (src == facetOutgToggleLabel) {
			facetOutgToggleLabel->SetState(true);
			facetOutgPerAreaToggleLabel->SetState(false);
		}
		else if (src == facetOutgPerAreaToggleLabel) {
			facetOutgToggleLabel->SetState(false);
			facetOutgPerAreaToggleLabel->SetState(true);
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

		//BUTTON ------------------------------------------------------------------
	case MSG_BUTTON:
		if (src == startSimu) {
			changedSinceSave = true;
			StartStopSimulation();
			resetSimu->SetEnabled(!worker.IsRunning());
		}

		else if (src == facetApplyBtn) {
			//changedSinceSave=true;
			ApplyFacetParams();
		}
		else if (src == facetDetailsBtn) {
			if (facetDetails == nullptr) facetDetails = new FacetDetails();
			facetDetails->Display(&worker);
		}
		else if (src == viewerMoreButton) {
			if (!viewer3DSettings)	viewer3DSettings = new Viewer3DSettings();
			viewer3DSettings->SetVisible(!viewer3DSettings->IsVisible());
			viewer3DSettings->Reposition();
			viewer3DSettings->Refresh(interfGeom, viewers[curViewer]);
		}
		else if (src == textureScalingBtn) {
			if (!textureScaling || !textureScaling->IsVisible()) {
				SAFE_DELETE(textureScaling);
				textureScaling = new TextureScaling(&worker, viewers);
				textureScaling->Display();
			}
			else {
				textureScaling->SetVisible(false);
			}
		}
		else if (src == profilePlotterBtn) {
			if (!profilePlotter) profilePlotter = new ProfilePlotter(&worker);
			if (!profilePlotter->IsVisible()) profilePlotter->Display(&worker);
			else profilePlotter->SetVisible(false);
		}
		else if (src == texturePlotterBtn) {
			if (!texturePlotter) texturePlotter = new TexturePlotter();
			if (!texturePlotter->IsVisible()) texturePlotter->Display(&worker);
			else {
				texturePlotter->SetVisible(false);
				SAFE_DELETE(texturePlotter);
			}
		}
		else if (src == globalSettingsBtn) {
			if (!globalSettings) globalSettings = new GlobalSettings(&worker);
			if (!globalSettings->IsVisible()) {
				globalSettings->Update();
				globalSettings->SetVisible(true);
			}
			else globalSettings->SetVisible(false);
		}
		else if (src == facetAdvParamsBtn) {
			if (!facetAdvParams) {
				facetAdvParams = new FacetAdvParams(&worker);
				facetAdvParams->Refresh(interfGeom->GetSelectedFacets());
			}
			facetAdvParams->SetVisible(!facetAdvParams->IsVisible());
			facetAdvParams->Reposition();
		}
		/*else {
			ProcessFormulaButtons(src);
		}*/
		break;
	default:
		//GLMessageBox::Display("Corrupted menu item selected!", "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}
}

void MolFlow::SelectViewer(int s)
{
	Interface::SelectViewer(s);
	if (viewer3DSettings && viewer3DSettings->IsVisible()) viewer3DSettings->Refresh(worker.GetGeometry(),viewers[curViewer]);
}

void MolFlow::BuildPipe(double ratio, int steps) {

	char tmp[256];
	MolflowGeometry* interfGeom = worker.GetMolflowGeometry();

	double R = 1.0;
	double L = ratio * R;
	int    step;
	if (steps) step = 5; //Quick Pipe
	else {
		sprintf(tmp, "100");
		//sprintf(title,"Pipe L/R = %g",L/R);
		char* nbFacets = GLInputBox::GetInput(tmp, "Number of facet", "Build Pipe");
		if (!nbFacets) return;
		if ((sscanf(nbFacets, "%d", &step) <= 0) || (step < 3)) {
			GLMessageBox::Display("Invalid number", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}
	}
	std::ostringstream temp;
	temp << "PIPE" << L / R;
	interfGeom->UpdateName(temp.str().c_str());

	try {
		SetDefaultViews();
		interfGeom->BuildPipe(L, R, 0, step);
		worker.MarkToReload();
		//worker.CalcTotalOutgassing();
		//default values
		worker.model->sp = SimuParams(); //reset to default
		worker.model->otfParams = OntheflySimulationParams(); //reset to default
		worker.ResetMoments();
		ResetSimulation(false);
	}
	catch (const std::exception& e) {
		GLMessageBox::Display((char*)e.what(), "Error building pipe", GLDLG_OK, GLDLG_ICONERROR);
		interfGeom->Clear();
		ResetSimulation(false);
		return;
	}
	//worker.globalState->globalStats.globalStats.hit.nbDesorbed = 0; //Already done by ResetWorkerStats
	//sprintf(tmp,"L|R %g",L/R);
	worker.SetCurrentFileName("");
	nbDesStart = 0;
	nbHitStart = 0;

	for (auto& view : viewers) {
		view->SetWorker(&worker);
	}

	//UpdateModelParams();
	startSimu->SetEnabled(true);
	//compACBtn->SetEnabled(modeCombo->GetSelectedIndex() == 1);
	//resetSimu->SetEnabled(true);
	ClearFacetParams();
	ClearFormulas();
	TimeDependentParameters::ClearParameters(worker.interfaceParameterCache);
	ClearAllSelections();
	ClearAllViews();

	appFormulas->AddFormula("Trans.prob.", "A2/SUMDES");

	UpdateStructMenu();
	// Send to sub process
	worker.MarkToReload();

	//UpdatePlotters();

	if (timeSettings) timeSettings->RefreshMoments();
	if (momentsEditor) momentsEditor->Refresh();
	if (pressureEvolution) pressureEvolution->Reset();
	if (timewisePlotter) timewisePlotter->Refresh();
	if (profilePlotter) profilePlotter->Refresh();
	if (histogramSettings) histogramSettings->Refresh({});
	if (histogramPlotter) histogramPlotter->Reset();
	if (texturePlotter) texturePlotter->Update(0.0, true);
	//if (parameterEditor) parameterEditor->UpdateCombo(); //Done by ClearParameters()
	if (textureScaling) textureScaling->Update();
	if (outgassingMapWindow) outgassingMapWindow->Update(m_fTime, true);
	if (facetDetails) facetDetails->Update();
	if (facetCoordinates) facetCoordinates->UpdateFromSelection();
	if (vertexCoordinates) vertexCoordinates->Update();
	if (movement) movement->Update();
	if (measureForces) measureForces->Update();
	if (globalSettings && globalSettings->IsVisible()) globalSettings->Update();
	if (formulaEditor) formulaEditor->Refresh();
	if (parameterEditor) parameterEditor->Refresh();
	UpdateTitle();
	changedSinceSave = false;
	ResetAutoSaveTimer();
}

void MolFlow::EmptyGeometry() {

	InterfaceGeometry* interfGeom = worker.GetGeometry();
	ResetSimulation(false);

	try {
		SetDefaultViews();
		interfGeom->EmptyGeometry();
		//worker.CalcTotalOutgassing();
		//default values
		worker.model->sp = SimuParams(); //reset to default
		worker.model->otfParams = OntheflySimulationParams(); //reset to default
		worker.ResetMoments();
	}
	catch (const std::exception& e) {
		GLMessageBox::Display(e.what(), "Error resetting geometry", GLDLG_OK, GLDLG_ICONERROR);
		interfGeom->Clear();
		return;
	}
	worker.SetCurrentFileName("");
	nbDesStart = 0;
	nbHitStart = 0;

	for (auto& view : viewers) {
		view->SetWorker(&worker);
	}

	//UpdateModelParams();
	startSimu->SetEnabled(true);
	//compACBtn->SetEnabled(modeCombo->GetSelectedIndex() == 1);
	//resetSimu->SetEnabled(true);
	ClearFacetParams();
	ClearFormulas();
	TimeDependentParameters::ClearParameters(worker.interfaceParameterCache);
	ClearAllSelections();
	ClearAllViews();

	/*GLFormula *f = new GLFormula();
	f.SetExpression("A2/SUMDES");
	f.SetName("Trans. Prob.");
	f.Parse();*/
	//AddFormula("Trans.prob.", "A2/SUMDES");

	UpdateStructMenu();
	// Send to sub process
	worker.MarkToReload();

	//UpdatePlotters();

	if (timeSettings) timeSettings->RefreshMoments();
	if (momentsEditor) momentsEditor->Refresh();
	if (pressureEvolution) pressureEvolution->Refresh();
	if (timewisePlotter) timewisePlotter->Refresh();
	if (profilePlotter) profilePlotter->Refresh();
	if (histogramSettings) histogramSettings->Refresh({});
	if (histogramPlotter) histogramPlotter->Reset();
	if (texturePlotter) texturePlotter->Update(0.0, true);
	//if (parameterEditor) parameterEditor->UpdateCombo(); //Done by ClearParameters()
	if (outgassingMapWindow) outgassingMapWindow->Update(m_fTime, true);
	if (movement) movement->Update();
	if (measureForces) measureForces->Update();
	if (globalSettings && globalSettings->IsVisible()) globalSettings->Update();
	if (formulaEditor) formulaEditor->Refresh();

	if (textureScaling) textureScaling->Update();
	if (facetDetails) facetDetails->Update();
	if (facetCoordinates) facetCoordinates->UpdateFromSelection();
	if (vertexCoordinates) vertexCoordinates->Update();

	UpdateTitle();
	changedSinceSave = false;
	ResetAutoSaveTimer();
	RefreshPlotterCombos();
	//UpdatePlotters();
}

void MolFlow::LoadConfig() {

	try {

		auto file = FileReader("molflow.cfg");
		MolflowGeometry* interfGeom = worker.GetMolflowGeometry();

		file.ReadKeyword("showRules"); file.ReadKeyword(":");
		for (auto& view : viewers)
			view->showRule = file.ReadInt();
		file.ReadKeyword("showNormals"); file.ReadKeyword(":");
		for (auto& view : viewers)
			view->showNormal = file.ReadInt();
		file.ReadKeyword("showUV"); file.ReadKeyword(":");
		for (auto& view : viewers)
			view->showUV = file.ReadInt();
		file.ReadKeyword("showLines"); file.ReadKeyword(":");
		for (auto& view : viewers)
			view->showLine = file.ReadInt();
		file.ReadKeyword("showLeaks"); file.ReadKeyword(":");
		for (auto& view : viewers)
			view->showLeak = file.ReadInt();
		file.ReadKeyword("showHits"); file.ReadKeyword(":");
		for (auto& view : viewers)
			view->showHit = file.ReadInt();
		file.ReadKeyword("showVolume"); file.ReadKeyword(":");
		for (auto& view : viewers)
			view->showVolume = file.ReadInt();
		file.ReadKeyword("showTexture"); file.ReadKeyword(":");
		for (auto& view : viewers)
			view->showTexture = file.ReadInt();
		file.ReadKeyword("showFacetId"); file.ReadKeyword(":");
		for (auto& view : viewers)
			view->showFacetId = file.ReadInt();
		file.ReadKeyword("showFilter"); file.ReadKeyword(":");
		for (auto& view : viewers)
			view->showFilter = file.ReadInt();
		file.ReadKeyword("showIndices"); file.ReadKeyword(":");
		for (auto& view : viewers)
			view->showIndex = file.ReadInt();
		file.ReadKeyword("showVertices"); file.ReadKeyword(":");
		for (auto& view : viewers)
			view->showVertexId = file.ReadInt();
		file.ReadKeyword("showMode"); file.ReadKeyword(":");
		for (auto& view : viewers)
			view->volumeRenderMode = static_cast<VolumeRenderMode>(file.ReadInt());
		file.ReadKeyword("showMesh"); file.ReadKeyword(":");
		for (auto& view : viewers)
			view->showMesh = file.ReadInt();
		file.ReadKeyword("showHidden"); file.ReadKeyword(":");
		for (auto& view : viewers)
			view->showHiddenFacet = file.ReadInt();
		file.ReadKeyword("showHiddenVertex"); file.ReadKeyword(":");
		for (auto& view : viewers)
			view->showHiddenVertex = file.ReadInt();
		file.ReadKeyword("showTimeOverlay"); file.ReadKeyword(":");
		for (auto& view : viewers)
			view->showTime = file.ReadInt();
		file.ReadKeyword("texColormap"); file.ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			//viewers[i]->showColormap = 
			file.ReadInt();
		file.ReadKeyword("translation"); file.ReadKeyword(":");
		for (auto& view : viewers)
			view->transStep = file.ReadDouble();
		file.ReadKeyword("dispNumLines"); file.ReadKeyword(":");
		for (auto& view : viewers)
			view->dispNumHits = file.ReadSizeT();
		file.ReadKeyword("dispNumLeaks"); file.ReadKeyword(":");
		for (auto& view : viewers)
			view->dispNumLeaks = file.ReadSizeT();
		file.ReadKeyword("dirShow"); file.ReadKeyword(":");
		for (auto& view : viewers)
			view->showDir = file.ReadInt();
		file.ReadKeyword("dirNorme"); file.ReadKeyword(":");
		interfGeom->SetNormeRatio((float)file.ReadDouble());
		file.ReadKeyword("dirAutoNormalize"); file.ReadKeyword(":");
		interfGeom->SetAutoNorme(file.ReadInt());
		file.ReadKeyword("dirCenter"); file.ReadKeyword(":");
		interfGeom->SetCenterNorme(file.ReadInt());
		file.ReadKeyword("angle"); file.ReadKeyword(":");
		for (auto& view : viewers)
			view->angleStep = file.ReadDouble();
		file.ReadKeyword("autoScale"); file.ReadKeyword(":");
		interfGeom->texAutoScale = file.ReadInt();
		file.ReadKeyword("autoScale_include_constant_flow"); file.ReadKeyword(":");
		interfGeom->texAutoScaleMode = static_cast<AutoScaleMode>(file.ReadInt());

		file.ReadKeyword("textures_min_pressure_all"); file.ReadKeyword(":");
		interfGeom->texture_limits[0].autoscale.min.steady_state = file.ReadDouble();
		file.ReadKeyword("textures_min_pressure_moments_only"); file.ReadKeyword(":");
		interfGeom->texture_limits[0].autoscale.min.moments_only = file.ReadDouble();
		file.ReadKeyword("textures_max_pressure_all"); file.ReadKeyword(":");
		interfGeom->texture_limits[0].autoscale.max.steady_state = file.ReadDouble();
		file.ReadKeyword("textures_max_pressure_moments_only"); file.ReadKeyword(":");
		interfGeom->texture_limits[0].autoscale.max.moments_only = file.ReadDouble();

		file.ReadKeyword("textures_min_impingement_all"); file.ReadKeyword(":");
		interfGeom->texture_limits[1].autoscale.min.steady_state = file.ReadDouble();
		file.ReadKeyword("textures_min_impingement_moments_only"); file.ReadKeyword(":");
		interfGeom->texture_limits[1].autoscale.min.moments_only = file.ReadDouble();
		file.ReadKeyword("textures_max_impingement_all"); file.ReadKeyword(":");
		interfGeom->texture_limits[1].autoscale.max.steady_state = file.ReadDouble();
		file.ReadKeyword("textures_max_impingement_moments_only"); file.ReadKeyword(":");
		interfGeom->texture_limits[1].autoscale.max.moments_only = file.ReadDouble();

		file.ReadKeyword("textures_min_density_all"); file.ReadKeyword(":");
		interfGeom->texture_limits[2].autoscale.min.steady_state = file.ReadDouble();
		file.ReadKeyword("textures_min_density_moments_only"); file.ReadKeyword(":");
		interfGeom->texture_limits[2].autoscale.min.moments_only = file.ReadDouble();
		file.ReadKeyword("textures_max_density_all"); file.ReadKeyword(":");
		interfGeom->texture_limits[2].autoscale.max.steady_state = file.ReadDouble();
		file.ReadKeyword("textures_max_density_moments_only"); file.ReadKeyword(":");
		interfGeom->texture_limits[2].autoscale.max.moments_only = file.ReadDouble();

		file.ReadKeyword("processNum"); file.ReadKeyword(":");
		nbProc = file.ReadSizeT();
#if defined(_DEBUG)
		nbProc = 1;
#endif
		if (nbProc <= 0) nbProc = 1;
		file.ReadKeyword("recents"); file.ReadKeyword(":"); file.ReadKeyword("{");
		std::string path = file.ReadString();
		while (path != "}" && recentsList.size() < MAX_RECENT) {
			recentsList.push_back(path);
			path = file.ReadString();
		}
		file.ReadKeyword("autonorme"); file.ReadKeyword(":");
		interfGeom->SetAutoNorme(file.ReadInt());
		file.ReadKeyword("centernorme"); file.ReadKeyword(":");
		interfGeom->SetCenterNorme(file.ReadInt());
		file.ReadKeyword("normeratio"); file.ReadKeyword(":");
		interfGeom->SetNormeRatio((float)(file.ReadDouble()));
		file.ReadKeyword("autoSaveFrequency"); file.ReadKeyword(":");
		autoSaveFrequency = file.ReadDouble();
		file.ReadKeyword("autoSaveSimuOnly"); file.ReadKeyword(":");
		autoSaveSimuOnly = file.ReadInt();
		file.ReadKeyword("checkForUpdates"); file.ReadKeyword(":");
		/*checkForUpdates =*/ file.ReadInt(); //Old checkforupdates
		file.ReadKeyword("autoUpdateFormulas"); file.ReadKeyword(":");
		autoUpdateFormulas = file.ReadInt();
		file.ReadKeyword("compressSavedFiles"); file.ReadKeyword(":");
		compressSavedFiles = file.ReadInt();
		file.ReadKeyword("expandShortcutPanel"); file.ReadKeyword(":");
		bool isOpen = file.ReadInt();
		if (isOpen) shortcutPanel->Open();
		else shortcutPanel->Close();
		file.ReadKeyword("hideLot"); file.ReadKeyword(":");
		for (auto& view : viewers)
			view->hideLot = file.ReadInt();
		file.ReadKeyword("textureLogScale"); file.ReadKeyword(":");
		interfGeom->texLogScale = file.ReadInt();
		file.ReadKeyword("leftHandedView"); file.ReadKeyword(":");
		leftHandedView = file.ReadInt();
		file.ReadKeyword("highlightNonplanarFacets"); file.ReadKeyword(":");
		highlightNonplanarFacets = file.ReadInt();
		file.ReadKeyword("highlightSelection"); file.ReadKeyword(":");
		highlightSelection = file.ReadInt();
		file.ReadKeyword("useOldXMLFormat"); file.ReadKeyword(":");
		useOldXMLFormat = file.ReadInt();
		file.ReadKeyword("showTP"); file.ReadKeyword(":");
		for (auto& view : viewers)
			view->showTP = file.ReadInt();
	}
	catch (...) {
		/*std::ostringstream tmp;
		tmp << err.GetMsg() << "\n\nThis is normal on the first launch and if you upgrade from an earlier version\n";
		tmp << "MolFlow will use default program settings.\nWhen you quit, a correct config file will be written\n";
		GLMessageBox::Display(tmp.str().c_str(), "Error loading config file", GLDLG_OK, GLDLG_ICONINFO);*/
		return;
	}
}


#define WRITEI(name,var) {             \
	file.Write(name);                      \
	file.Write(":");                       \
	for(size_t i=0;i<MAX_VIEWER;i++)        \
	file.Write(viewers[i]->var," ");   \
	file.Write("\n");                      \
}


#define WRITED(name,var) {             \
	file.Write(name);                      \
	file.Write(":");                       \
	for(size_t i=0;i<MAX_VIEWER;i++)        \
	file.Write(viewers[i]->var," ");\
	file.Write("\n");                      \
}


void MolFlow::SaveConfig() {

	try {

		auto file = FileWriter("molflow.cfg");
		MolflowGeometry* interfGeom = worker.GetMolflowGeometry();

		// Save flags
		WRITEI("showRules", showRule);
		WRITEI("showNormals", showNormal);
		WRITEI("showUV", showUV);
		WRITEI("showLines", showLine);
		WRITEI("showLeaks", showLeak);
		WRITEI("showHits", showHit);
		WRITEI("showVolume", showVolume);
		WRITEI("showTexture", showTexture);
		WRITEI("showFacetId", showFacetId);
		WRITEI("showFilter", showFilter);
		WRITEI("showIndices", showIndex);
		WRITEI("showVertices", showVertexId);
		WRITEI("showMode", volumeRenderMode);
		WRITEI("showMesh", showMesh);
		WRITEI("showHidden", showHiddenFacet);
		WRITEI("showHiddenVertex", showHiddenVertex);
		WRITEI("showTimeOverlay", showTime);
		//WRITEI("texColormap", showColormap);
		file.Write("texColormap:1 1 1 1\n");
		WRITED("translation", transStep);
		WRITEI("dispNumLines", dispNumHits);
		WRITEI("dispNumLeaks", dispNumLeaks);
		WRITEI("dirShow", showDir);
		file.Write("dirNorme:"); file.Write(interfGeom->GetNormeRatio(), "\n");
		file.Write("dirAutoNormalize:"); file.Write(interfGeom->GetAutoNorme(), "\n");
		file.Write("dirCenter:"); file.Write(interfGeom->GetCenterNorme(), "\n");

		WRITED("angle", angleStep);
		file.Write("autoScale:"); file.Write(interfGeom->texAutoScale, "\n");
		file.Write("autoScale_include_constant_flow:"); file.Write(static_cast<int>(interfGeom->texAutoScaleMode), "\n");

		file.Write("textures_min_pressure_all:");
		file.Write(interfGeom->texture_limits[0].autoscale.min.steady_state, "\n");
		file.Write("textures_min_pressure_moments_only:");
		file.Write(interfGeom->texture_limits[0].autoscale.min.moments_only, "\n");
		file.Write("textures_max_pressure_all:");
		file.Write(interfGeom->texture_limits[0].autoscale.max.steady_state, "\n");
		file.Write("textures_max_pressure_moments_only:");
		file.Write(interfGeom->texture_limits[0].autoscale.max.moments_only, "\n");

		file.Write("textures_min_impingement_all:");
		file.Write(interfGeom->texture_limits[1].autoscale.min.steady_state, "\n");

		file.Write("textures_min_impingement_moments_only:");
		file.Write(interfGeom->texture_limits[1].autoscale.min.moments_only, "\n");
		file.Write("textures_max_impingement_all:");
		file.Write(interfGeom->texture_limits[1].autoscale.max.steady_state, "\n");
		file.Write("textures_max_impingement_moments_only:");
		file.Write(interfGeom->texture_limits[1].autoscale.max.moments_only, "\n");

		file.Write("textures_min_density_all:");
		file.Write(interfGeom->texture_limits[2].autoscale.min.steady_state, "\n");
		file.Write("textures_min_density_moments_only:");
		file.Write(interfGeom->texture_limits[2].autoscale.min.moments_only, "\n");
		file.Write("textures_max_density_all:");
		file.Write(interfGeom->texture_limits[2].autoscale.max.steady_state, "\n");
		file.Write("textures_max_density_moments_only:");
		file.Write(interfGeom->texture_limits[2].autoscale.max.moments_only, "\n");

#if defined(_DEBUG)
		file.Write("processNum:"); file.Write(numCPU, "\n");
#else
		file.Write("processNum:"); file.Write(worker.GetProcNumber(), "\n");
#endif
		file.Write("recents:{\n");
		for (auto& recent : recentsList) {
			file.Write("\"");
			file.Write(recent);
			file.Write("\"\n");
		}
		file.Write("}\n");

		file.Write("autonorme:"); file.Write(interfGeom->GetAutoNorme(), "\n");
		file.Write("centernorme:"); file.Write(interfGeom->GetCenterNorme(), "\n");
		file.Write("normeratio:"); file.Write((double)(interfGeom->GetNormeRatio()), "\n");

		file.Write("autoSaveFrequency:"); file.Write(autoSaveFrequency, "\n");
		file.Write("autoSaveSimuOnly:"); file.Write(autoSaveSimuOnly, "\n");
		file.Write("checkForUpdates:"); file.Write(/*checkForUpdates*/ 0, "\n"); //Deprecated
		file.Write("autoUpdateFormulas:"); file.Write(autoUpdateFormulas, "\n");
		file.Write("compressSavedFiles:"); file.Write(compressSavedFiles, "\n");
		file.Write("expandShortcutPanel:"); file.Write(!shortcutPanel->IsClosed(), "\n");

		WRITEI("hideLot", hideLot);
		file.Write("textureLogScale:"); file.Write(interfGeom->texLogScale, "\n");
		file.Write("leftHandedView:"); file.Write(leftHandedView, "\n");
		file.Write("highlightNonplanarFacets:"); file.Write(highlightNonplanarFacets, "\n");
		file.Write("highlightSelection:"); file.Write(highlightSelection, "\n");
		file.Write("useOldXMLFormat:"); file.Write(useOldXMLFormat, "\n");
		WRITEI("showTP", showTP);
	}
	catch (const std::exception& e) {
		GLMessageBox::Display(e.what(), "Error saving config file", GLDLG_OK, GLDLG_ICONWARNING);
	}
}

void MolFlow::CalcPumpingSpeed() {
	double sticking;
	double area;
	double pumpingSpeed;
	double temperature;
	//double mass;

	if (!facetStickingText->GetNumber(&sticking)) {
		facetPumpingSpeedText->SetText("sticking?");
		return;
	}
	facetAreaText->GetNumber(&area);
	if (!facetTemperatureText->GetNumber(&temperature)) {
		facetPumpingSpeedText->SetText("temp?");
		return;
	}

	pumpingSpeed = 1 * sticking * area / 10.0 / 4.0 * sqrt(8.0 * 8.31 * temperature / PI / (worker.model->sp.gasMass * 0.001));
	facetPumpingSpeedText->SetText(pumpingSpeed);
}

void MolFlow::CalcSticking() {
	double sticking;
	double area;
	double pumpingSpeed;
	double temperature;
	//double mass;

	if (!facetPumpingSpeedText->GetNumber(&pumpingSpeed)) {
		return;
	}
	facetAreaText->GetNumber(&area);
	if (!facetTemperatureText->GetNumber(&temperature)) {
		return;
	}

	sticking = std::abs(pumpingSpeed / (area / 10.0) * 4.0 * sqrt(1.0 / 8.0 / 8.31 / (temperature)*PI * (worker.model->sp.gasMass * 0.001)));
	facetStickingText->SetText(sticking);
}

void MolFlow::CrashHandler(const std::exception& e) {
	char tmp[1024];
	sprintf(tmp, "Well, that's embarrassing. Molflow crashed and will exit now.\nBefore that, an autosave will be attempted.\nHere is the error info:\n\n%s", e.what());
	GLMessageBox::Display(tmp, "Main crash handler", GLDLG_OK, GLDGL_ICONDEAD);
	try {
		if (AutoSave(true))
			GLMessageBox::Display("Good news, autosave worked!", "Main crash handler", GLDLG_OK, GLDGL_ICONDEAD);
		else
			GLMessageBox::Display("Sorry, I couldn't even autosave.", "Main crash handler", GLDLG_OK, GLDGL_ICONDEAD);
	}
	catch (const std::exception& e) {
		sprintf(tmp, "Sorry, I couldn't even autosave:\n\n%s", e.what());
		GLMessageBox::Display(tmp, "Main crash handler", GLDLG_OK, GLDGL_ICONDEAD);
	}
}

void MolFlow::UpdatePlotters() {
	const bool forceUpdate = true; //worker.IsRunning();
	if (pressureEvolution) pressureEvolution->Update(m_fTime, forceUpdate);
	if (timewisePlotter) timewisePlotter->Update(m_fTime, forceUpdate);
	if (profilePlotter) profilePlotter->Update(m_fTime, forceUpdate);
	if (texturePlotter) texturePlotter->Update(m_fTime, forceUpdate);
	if (histogramPlotter) histogramPlotter->Update(m_fTime, forceUpdate);
	if (convergencePlotter) convergencePlotter->Update(m_fTime);
}

void MolFlow::RefreshPlotterCombos() {
	//Removes non-present views, rebuilds combobox and refreshes plotted data
	if (pressureEvolution) pressureEvolution->Refresh();
	if (timewisePlotter) timewisePlotter->Refresh();
	if (profilePlotter) profilePlotter->Refresh();
	if (histogramPlotter) histogramPlotter->Refresh();
	if (convergencePlotter) convergencePlotter->Refresh();
}

void MolFlow::UpdateFacetHits(bool allRows) {
	char tmp[256];
	InterfaceGeometry* interfGeom = worker.GetGeometry();

	try {
		// Facet list
		if (interfGeom->IsLoaded()) {

			int sR, eR;
			if (allRows)
			{
				sR = 0;
				eR = (int)facetList->GetNbRow() - 1;
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
				if (facetId == -2) facetId = (int)i;
				if (i >= interfGeom->GetNbFacet()) {
					char errMsg[512];
					sprintf(errMsg, "Molflow::UpdateFacetHits()\nError while updating facet hits. Was looking for facet #%d (/%zu) in list.\nMolflow will now autosave and crash.", i + 1, interfGeom->GetNbFacet());
					GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
					AutoSave();
					throw Error(errMsg);
				}
				InterfaceFacet* f = interfGeom->GetFacet(facetId);
				sprintf(tmp, "%d", facetId + 1);
				facetList->SetValueAt(0, i, tmp);

				facetList->SetColumnLabel(1, "Hits");
				sprintf(tmp, "%zd", f->facetHitCache.nbMCHit);
				facetList->SetValueAt(1, i, tmp);
				sprintf(tmp, "%zd", f->facetHitCache.nbDesorbed);
				facetList->SetValueAt(2, i, tmp);
				sprintf(tmp, "%g", f->facetHitCache.nbAbsEquiv);
				facetList->SetValueAt(3, i, tmp);
			}

		}
	}
	catch (const std::exception& e) {
		char errMsg[512];
		sprintf(errMsg, "%s\nError while updating facet hits", e.what());
		GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
	}

}
