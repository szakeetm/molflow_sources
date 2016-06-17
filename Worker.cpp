/*
File:        Worker.cpp
Description: Sub processes handling
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

#include <Windows.h>
#include "Worker.h"
#include "GLApp/GLApp.h"
#include "GLApp/GLMessageBox.h"
#include <math.h>
#include <stdlib.h>
#include <Process.h>
#include "GLApp/GLUnitDialog.h"
#include "Molflow.h"
#include <direct.h>
#include "Distributions.h"
#include "ZipUtils/zip.h"
#include "ZipUtils/unzip.h"
#include "File.h" //File utils (Get extension, etc)

/*
//Leak detection
#ifdef _DEBUG
#define _CRTDBG_MAP_ALLOC
#define DEBUG_NEW new(_NORMAL_BLOCK, __FILE__, __LINE__)
#define new DEBUG_NEW
#endif
*/
using namespace pugi;
extern MolFlow *mApp;

Worker::Worker() {
	//test
	temperatures = std::vector<double>();
	moments = std::vector<double>();
	desorptionParameterIDs = std::vector<size_t>();
	userMoments = std::vector<std::string>(); //strings describing moments, to be parsed
	CDFs = std::vector<std::vector<std::pair<double, double>>>();
	IDs = std::vector<std::vector<std::pair<double, double>>>();
	parameters = std::vector<Parameter>();
	needsReload = TRUE;  //When main and subprocess have different geometries, needs to reload (synchronize)
	displayedMoment = 0; //By default, steady-state is displayed

	//desorptionStartTime=0.0;
	//desorptionStopTime=1;
	timeWindowSize = 0.1;
	useMaxwellDistribution = TRUE;
	calcConstantFlow = TRUE;
	//valveOpenMoment=99999.0;
	distTraveledTotal_total = 0.0;
	distTraveledTotal_fullHitsOnly = 0.0;
	gasMass = 28.0;
	enableDecay = FALSE;
	halfLife = 1;
	finalOutgassingRate = finalOutgassingRate_Pa_m3_sec = totalDesorbedMolecules = 0.0;

	motionType = 0;

	pid = _getpid();
	sprintf(ctrlDpName, "MFLWCTRL%d", pid);
	sprintf(loadDpName, "MFLWLOAD%d", pid);
	sprintf(hitsDpName, "MFLWHITS%d", pid);
	nbProcess = 0;
	maxDesorption = 0;
	ResetWorkerStats();
	geom = new Geometry();
	dpControl = NULL;
	dpHit = NULL;
	nbHHit = 0;
	nbHit = 0;
	nbLeakTotal = 0;
	nbLastLeaks = 0;
	startTime = 0.0f;
	stopTime = 0.0f;
	simuTime = 0.0f;
	running = FALSE;
	calcAC = FALSE;
	strcpy(fullFileName, "");
	//Debug memory check
	//_ASSERTE (!_CrtDumpMemoryLeaks());
	_ASSERTE(_CrtCheckMemory());
}

// -------------------------------------------------------------

Worker::~Worker() {
	CLOSEDP(dpHit);
	Exit();
	delete geom;
}

// -------------------------------------------------------------

Geometry *Worker::GetGeometry() {
	return geom;
}

BOOL Worker::IsDpInitialized(){

	return (dpHit != NULL);
}
// -------------------------------------------------------------

char *Worker::GetFileName() {
	return fullFileName;
}

char *Worker::GetShortFileName() {

	static char ret[512];
	char *r = strrchr(fullFileName, '/');
	if (!r) r = strrchr(fullFileName, '\\');
	if (!r) strcpy(ret, fullFileName);
	else   {
		r++;
		strcpy(ret, r);
	}

	return ret;

}

char *Worker::GetShortFileName(char* longFileName) {

	static char ret[512];
	char *r = strrchr(longFileName, '/');
	if (!r) r = strrchr(longFileName, '\\');
	if (!r) strcpy(ret, longFileName);
	else   {
		r++;
		strcpy(ret, r);
	}

	return ret;

}

// -------------------------------------------------------------

void Worker::SetFileName(char *fileName) {
	strncpy(fullFileName, fileName, 512);
}

void Worker::SaveGeometry(char *fileName, GLProgress *prg, BOOL askConfirm, BOOL saveSelected, BOOL autoSave, BOOL crashSave) {
	try {
		if (needsReload && (!crashSave && !saveSelected)) RealReload();
	}
	catch (Error &e) {
		char errMsg[512];
		sprintf(errMsg, "Error reloading worker. Trying safe save (geometry only):\n%s", e.GetMsg());
		GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
		crashSave = TRUE;
	}
	char tmp[10000]; //compress.exe command line
	char fileNameWithGeo[2048]; //file name with .geo extension (instead of .geo7z)
	char fileNameWithGeo7z[2048];
	char fileNameWithXML[2048];
	char fileNameWithZIP[2048];
	char fileNameWithoutExtension[2048]; //file name without extension
	//char *ext = fileName+strlen(fileName)-4;
	char *ext, *dir;

	dir = strrchr(fileName, '\\');
	ext = strrchr(fileName, '.');

	if (!(ext) || !(*ext == '.') || ((dir) && (dir > ext))) {
		sprintf(fileName, mApp->compressSavedFiles ? "%s.zip" : "%s.xml", fileName); //set to default XML/ZIP format
		ext = strrchr(fileName, '.');
	}

	ext++;

	// Read a file
	BOOL ok = TRUE;
	FileWriter *f = NULL;
	BOOL isTXT = _stricmp(ext, "txt") == 0;
	BOOL isSTR = _stricmp(ext, "str") == 0;
	BOOL isGEO = _stricmp(ext, "geo") == 0;
	BOOL isGEO7Z = _stricmp(ext, "geo7z") == 0;
	BOOL isXML = _stricmp(ext, "xml") == 0;
	BOOL isXMLzip = _stricmp(ext, "zip") == 0;

	if (isTXT || isGEO || isGEO7Z || isSTR || isXML || isXMLzip) {
		if ((isGEO7Z) && WAIT_TIMEOUT == WaitForSingleObject(mApp->compressProcessHandle, 0)) {
			GLMessageBox::Display("Compressing a previous save file is in progress. Wait until that finishes"
				"or close process \"compress.exe\"\nIf this was an autosave attempt,"
				"you have to lower the autosave frequency.", "Can't save right now.", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}
		if (isGEO) {
			memcpy(fileNameWithoutExtension, fileName, sizeof(char)*(strlen(fileName) - 4));
			fileNameWithoutExtension[strlen(fileName) - 4] = '\0';
			sprintf(fileNameWithGeo7z, "%s7z", fileName);
			memcpy(fileNameWithGeo, fileName, (strlen(fileName) + 1)*sizeof(char));
		}
		else if (isGEO7Z) {
			memcpy(fileNameWithoutExtension, fileName, sizeof(char)*(strlen(fileName) - 6));
			fileNameWithoutExtension[strlen(fileName) - 6] = '\0';
			memcpy(fileNameWithGeo, fileName, sizeof(char)*(strlen(fileName) - 2));
			fileNameWithGeo[strlen(fileName) - 2] = '\0';
			memcpy(fileNameWithGeo7z, fileName, (1 + strlen(fileName))*sizeof(char));
			sprintf(tmp, "A .geo file of the same name exists. Overwrite that file ?\n%s", fileNameWithGeo);
			if (!autoSave && FileUtils::Exist(fileNameWithGeo)) {
				ok = (GLMessageBox::Display(tmp, "Question", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONWARNING) == GLDLG_OK);
			}
		}

		if (isXML || isXMLzip) {
			memcpy(fileNameWithoutExtension, fileName, sizeof(char)*(strlen(fileName) - 4));
			fileNameWithoutExtension[strlen(fileName) - 4] = '\0';
			sprintf(fileNameWithXML, "%s.xml", fileNameWithoutExtension);
			sprintf(fileNameWithZIP, "%s.zip", fileNameWithoutExtension);
		}
		if (isXMLzip) {
			sprintf(tmp, "An .xml file of the same name exists. Overwrite that file ?\n%s", fileNameWithZIP);
			if (!autoSave && FileUtils::Exist(fileNameWithXML)) {
				ok = (GLMessageBox::Display(tmp, "Question", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONWARNING) == GLDLG_OK);
			}
		}

		if (!autoSave && ok && FileUtils::Exist(fileName)) {
			sprintf(tmp, "Overwrite existing file ?\n%s", fileName);
			if (askConfirm) ok = (GLMessageBox::Display(tmp, "Question", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONWARNING) == GLDLG_OK);
		}

		if (ok) {
			if (isSTR) {
				geom->SaveSTR(dpHit, saveSelected);
			}
			else {
				try {
					if (isGEO7Z) {
						f = new FileWriter(fileNameWithGeo); //We first write a GEO file, then compress it to GEO7Z later
					}
					else if (!(isXML || isXMLzip))
						f = new FileWriter(fileName);
				}
				catch (Error &e) {
					SAFE_DELETE(f);
					GLMessageBox::Display((char*)e.GetMsg(), "Error writing file.", GLDLG_OK, GLDLG_ICONERROR);
					return;
				}
				geom->tNbDesorptionMax = maxDesorption;
				if (isTXT) geom->SaveTXT(f, dpHit, saveSelected);
				else if (isGEO || isGEO7Z) {
					// Retrieve leak cache
					int nbLeakSave, nbHHitSave;
					LEAK pLeak[NBHLEAK];
					if (!crashSave && !saveSelected) GetLeak(pLeak, &nbLeakSave);
					// Retrieve hit cache (lines and dots)
					HIT pHits[NBHHIT];
					if (!crashSave && !saveSelected) GetHHit(pHits, &nbHHitSave);
					geom->SaveGEO(f, prg, dpHit, this->userMoments, this, saveSelected, pLeak, &nbLeakSave, pHits, &nbHHitSave, crashSave);
				}
				else if (isXML || isXMLzip) {
					xml_document saveDoc;
					geom->SaveXML_geometry(saveDoc, this, prg, saveSelected);
					xml_document geom_only; geom_only.reset(saveDoc);
					BOOL success = FALSE; //success: simulation state could be saved
					if (!crashSave && !saveSelected) {
						try {
							AccessDataport(dpHit);
							BYTE *buffer;
							buffer = (BYTE *)dpHit->buff;
							SHGHITS *gHits;
							gHits = (SHGHITS *)buffer;
							int nbLeakSave, nbHHitSave;
							LEAK pLeak[NBHLEAK];
							GetLeak(pLeak, &nbLeakSave);
							HIT pHits[NBHHIT];
							GetHHit(pHits, &nbHHitSave);

							success = geom->SaveXML_simustate(saveDoc, this, buffer, gHits, nbLeakSave, nbHHitSave, pLeak, pHits, prg, saveSelected);
							ReleaseDataport(dpHit);
						}
						catch (Error &e) {
							SAFE_DELETE(f);
							ReleaseDataport(dpHit);
							GLMessageBox::Display((char*)e.GetMsg(), "Error saving simulation state.", GLDLG_OK, GLDLG_ICONERROR);
							return;
						}
					}

					prg->SetMessage("Writing xml file...");
					if (success) {
						if (!saveDoc.save_file(fileNameWithXML)) throw Error("Error writing XML file."); //successful save
					}
					else {
						if (!geom_only.save_file(fileNameWithXML)) throw Error("Error writing XML file."); //simu state error
					}

					if (isXMLzip) {
						prg->SetProgress(0.75);
						prg->SetMessage("Compressing xml to zip...");
						//mApp->compressProcessHandle=CreateThread(0, 0, ZipThreadProc, 0, 0, 0);
						HZIP hz = CreateZip(fileNameWithZIP, 0);
						if (!hz) {
							throw Error("Error creating ZIP file");
						}
						if (!ZipAdd(hz, GetShortFileName(fileNameWithXML), fileNameWithXML)) remove(fileNameWithXML);
						else {
							CloseZip(hz);
							throw Error("Error compressing ZIP file.");
						}
						CloseZip(hz);
					}


				}
			}
			/*if (!autoSave && !saveSelected) {
				strcpy(fullFileName, fileName);
				remove("Molflow_AutoSave.zip");
				}*/
		}
	}
	else {
		SAFE_DELETE(f);
		throw Error("SaveGeometry(): Invalid file extension [only xml,zip,geo,geo7z,txt or str]");
	}

	SAFE_DELETE(f);

	//File written, compress it if the user wanted to
	if (ok && isGEO7Z) {
		if (FileUtils::Exist("compress.exe")) { //compress GEO file to GEO7Z using 7-zip launcher "compress.exe"
			sprintf(tmp, "compress.exe \"%s\"", fileNameWithGeo);
			int procId = StartProc_background(tmp);
			mApp->compressProcessHandle = OpenProcess(PROCESS_ALL_ACCESS, TRUE, procId);
			fileName = fileNameWithGeo7z;
		}
		else {
			GLMessageBox::Display("compress.exe (part of Molfow) not found.\n Will save as uncompressed GEO file.", "Compressor not found", GLDLG_OK, GLDLG_ICONERROR);
			fileName = fileNameWithGeo;
		}
	}
	else if (ok && isGEO) fileName = fileNameWithGeo;
	if (!autoSave && !saveSelected) {
		SetFileName(fileName);
		mApp->UpdateTitle();
	}
}



void Worker::ExportTextures(char *fileName, int grouping, int mode, BOOL askConfirm, BOOL saveSelected) {

	char tmp[512];

	// Read a file
	FILE *f = NULL;



	BOOL ok = TRUE;
	if (askConfirm) {
		if (FileUtils::Exist(fileName)) {
			sprintf(tmp, "Overwrite existing file ?\n%s", fileName);
			ok = (GLMessageBox::Display(tmp, "Question", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONWARNING) == GLDLG_OK);
		}
	}
	if (ok) {
		f = fopen(fileName, "w");
		if (!f) {
			char tmp[256];
			sprintf(tmp, "Cannot open file for writing %s", fileName);
			throw Error(tmp);
		}
		geom->ExportTextures(f, grouping, mode, dpHit, saveSelected);
		fclose(f);
	}

}

void Worker::ExportProfiles(char *fileName) {

	char tmp[512];

	// Read a file
	FILE *f = NULL;

	char *ext, *dir;

	dir = strrchr(fileName, '\\');
	ext = strrchr(fileName, '.');

	if (!(ext) || !(*ext == '.') || ((dir) && (dir > ext))) {
		sprintf(fileName, "%s.csv", fileName); //set to default CSV format
		ext = strrchr(fileName, '.');
	}
	ext++;
	BOOL isTXT = _stricmp(ext, "txt") == 0;

	BOOL ok = TRUE;

	if (FileUtils::Exist(fileName)) {
		sprintf(tmp, "Overwrite existing file ?\n%s", fileName);
		ok = (GLMessageBox::Display(tmp, "Question", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONWARNING) == GLDLG_OK);
	}

	if (ok) {
		f = fopen(fileName, "w");
		if (!f) {
			char tmp[256];
			sprintf(tmp, "Cannot open file for writing %s", fileName);
			throw Error(tmp);
		}
		geom->ExportProfiles(f, isTXT, dpHit, this);
		fclose(f);
	}
}

/*void Worker::ImportDesorption(char *fileName) {
	//if (needsReload) RealReload();

	// Read a file
	FileReader *f=new FileReader(fileName);
	geom->ImportDesorption(f,dpHit);
	SAFE_DELETE(f);
	changedSinceSave=TRUE;
	Reload();
	}*/

// -------------------------------------------------------------

void Worker::LoadGeometry(char *fileName,BOOL insert,BOOL newStr) {
	if (!insert) {
		needsReload = TRUE;
	}
	else {
		RealReload();
	}
	char CWD[MAX_PATH];
	_getcwd(CWD, MAX_PATH);

	std::string ext=FileUtils::GetExtension(fileName);

	if (ext == "")
		throw Error("LoadGeometry(): No file extension, can't determine type");

	// Read a file
	FileReader *f = NULL;
	GLProgress *progressDlg = new GLProgress("Reading file...", "Please wait");
	progressDlg->SetVisible(TRUE);
	progressDlg->SetProgress(0.0);

	if (!insert) {
		//Clear hits and leaks cache
		memset(hhitCache, 0, sizeof(HIT)*NBHHIT);
		memset(leakCache, 0, sizeof(LEAK)*NBHLEAK);
		ResetMoments();
		//default values
		enableDecay = FALSE;
		gasMass = 28;
	}

	/*
	BOOL isASE = (_stricmp(ext, "ase") == 0);
	BOOL isSTR = (_stricmp(ext, "str") == 0);
	BOOL isSTL = (_stricmp(ext, "stl") == 0);
	BOOL isTXT = (_stricmp(ext, "txt") == 0);
	BOOL isGEO = (_stricmp(ext, "geo") == 0);
	BOOL isGEO7Z = (_stricmp(ext, "geo7z") == 0);
	BOOL isSYN = (_stricmp(ext, "syn") == 0);
	BOOL isSYN7Z = (_stricmp(ext, "syn7z") == 0);
	BOOL isXML = (_stricmp(ext, "xml") == 0);
	BOOL isXMLzip = (_stricmp(ext, "zip") == 0);
	*/

	if (ext == "txt" || ext == "TXT") {

		try {
			if (!insert) ResetWorkerStats();
			else mApp->changedSinceSave = TRUE;

			f = new FileReader(fileName);
			
			if (!insert) {
				geom->LoadTXT(f, progressDlg);
				nbHit = geom->tNbHit;
				nbDesorption = geom->tNbDesorption;
				nbAbsorption = geom->tNbAbsorption;
				maxDesorption = geom->tNbDesorptionMax;
				nbLeakTotal = geom->tNbLeak;
				//RealReload();
				strcpy(fullFileName, fileName);
			}
			else { //insert
				geom->InsertTXT(f, progressDlg, newStr);
				nbHit = 0;
				nbDesorption = 0;
				maxDesorption = 0;
				nbLeakTotal = 0;
				Reload();
			}
		}
		catch (Error &e) {
			if (!insert) geom->Clear();
			SAFE_DELETE(f);
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			throw e;
		}

	}
	else if (ext == "stl" || ext == "STL") {
		try {
			int ret = GLUnitDialog::Display("", "Choose STL file units:", GLDLG_MM | GLDLG_CM | GLDLG_M | GLDLG_INCH | GLDLG_FOOT | GLDLG_CANCEL_U, GLDLG_ICONNONE);
			double scaleFactor = 1.0;
			switch (ret) {
			case GLDLG_MM:
				scaleFactor = 0.1;
				break;
			case GLDLG_CM:
				scaleFactor = 1.0;
				break;
			case GLDLG_M:
				scaleFactor = 100;
				break;
			case GLDLG_INCH:
				scaleFactor = 2.54;
				break;
			case GLDLG_FOOT:
				scaleFactor = 30.48;
				break;
			}
			if (ret != GLDLG_CANCEL_U) {
				progressDlg->SetMessage("Resetting worker...");
				progressDlg->SetVisible(TRUE);
				ResetWorkerStats();				
				progressDlg->SetMessage("Reading geometry...");
				f = new FileReader(fileName);
				if (!insert) {
					geom->LoadSTL(f, progressDlg, scaleFactor);
					strcpy(fullFileName, fileName);
					mApp->DisplayCollapseDialog();
				}
				else { //insert
					mApp->changedSinceSave = TRUE;
					geom->InsertSTL(f, progressDlg, scaleFactor, newStr);
					Reload();
				}
			}
		}
		catch (Error &e) {
			if (!insert) geom->Clear();
			SAFE_DELETE(f);
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			throw e;

		}

	}
	else if (ext == "str" || ext == "STR") {
		if (insert) throw Error("STR file inserting is not supported.");
		try {
			ResetWorkerStats();
			f = new FileReader(fileName);
			progressDlg->SetVisible(TRUE);
			geom->LoadSTR(f, progressDlg);
			//RealReload();
			strcpy(fullFileName, fileName);
		}
		catch (Error &e) {
			geom->Clear();
			SAFE_DELETE(f);
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			throw e;
		}

	}
	else if (ext=="syn" || ext=="syn7z") { //Synrad file
		int version;
		progressDlg->SetVisible(TRUE);
		try {
			if (ext=="syn7z") {
				//decompress file
				progressDlg->SetMessage("Decompressing file...");
				char tmp[1024];
				//char fileOnly[512];
				sprintf(tmp, "cmd /C \"pushd \"%s\"&&7za.exe x -t7z -aoa \"%s\" -otmp&&popd\"", CWD, fileName);
				//execCMD(tmp);
				system(tmp);

				/*filebegin = strrchr(fileName, '\\');
				if (filebegin) filebegin++;
				else filebegin = fileName;
				memcpy(fileOnly, filebegin, sizeof(char)*(strlen(filebegin) - 2));
				fileOnly[strlen(filebegin) - 2] = '\0';
				sprintf(tmp2, "%s\\tmp\\Geometry.syn", CWD);*/
				f = new FileReader((std::string)CWD + "\\tmp\\Geometry.syn"); //Open extracted file
			} else f = new FileReader(fileName); //syn file, open it directly
			if (!insert) {
				progressDlg->SetMessage("Resetting worker...");
				ResetWorkerStats();

				//leaks
				LEAK pLeak[NBHLEAK];
				//hits
				HIT pHits[NBHHIT];
			
				geom->LoadSYN(f, progressDlg, pLeak, &nbLastLeaks, pHits, &nbHHit, &version);
				maxDesorption = 0;
			}
			else { //insert
				geom->InsertSYN(f, progressDlg, newStr);
				nbLeakTotal = 0;
				nbHit = 0;
				nbDesorption = 0;
				maxDesorption = geom->tNbDesorptionMax;
			}

			progressDlg->SetMessage("Reloading worker with new geometry...");
			RealReload(); //for the loading of textures
			//SHGHITS *gHits = (SHGHITS *)dpHit->buff;
			SAFE_DELETE(f);
			if (!insert) strcpy(fullFileName, fileName);
		}
		catch (Error &e) {
			if (!insert) geom->Clear();
			SAFE_DELETE(f);
			//if (isSYN7Z) remove(tmp2);
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			throw e;
		}

	}
	else if (ext=="geo" || ext=="geo7z") {
		//char tmp2[1024];
		std::string toOpen;
		int version;
		progressDlg->SetVisible(TRUE);
		try {
			if (ext == "geo7z") {
				//decompress file
				progressDlg->SetMessage("Decompressing file...");
				char tmp[1024];
				//char fileOnly[512];
				sprintf(tmp, "cmd /C \"pushd \"%s\"&&7za.exe x -t7z -aoa \"%s\" -otmp&&popd\"", CWD, fileName);
				system(tmp);
				toOpen = (std::string)CWD + "\\tmp\\Geometry.geo"; //newer geo7z format: contain Geometry.geo
				if (!FileUtils::Exist(toOpen)) toOpen = ((std::string)fileName).substr(0, strlen(fileName) - 2); //Inside the zip, try original filename with extension changed from geo7z to geo
				f = new FileReader(toOpen);
			}
			else { //not geo7z
				toOpen = fileName;
				f = new FileReader(fileName); //geo file, open it directly
			}
			progressDlg->SetMessage("Resetting worker...");
			ResetWorkerStats();
			if (insert) mApp->changedSinceSave = TRUE;

			if (!insert) {
				//leaks
				LEAK pLeak[NBHLEAK];
				//hits
				HIT pHits[NBHHIT];

				geom->LoadGEO(f, progressDlg, pLeak, &nbLastLeaks, pHits, &nbHHit, &version, this);
				//copy temp values from geom to worker:
				nbLeakTotal = geom->tNbLeak;
				nbHit = geom->tNbHit;
				nbDesorption = geom->tNbDesorption;
				maxDesorption = geom->tNbDesorptionMax;
				nbAbsorption = geom->tNbAbsorption;
				distTraveledTotal_total = geom->distTraveledTotal_total;
				distTraveledTotal_fullHitsOnly = geom->distTraveledTotal_fullHitsOnly;
				progressDlg->SetMessage("Reloading worker with new geometry...");
				RealReload(); //for the loading of textures
				if (version >= 8) geom->LoadProfile(f, dpHit, version);
				SHGHITS *gHits = (SHGHITS *)dpHit->buff;
				SetLeak(pLeak, &nbLastLeaks, gHits);
				SetHHit(pHits, &nbHHit, gHits);
				SendHits(); //Global and facet hit counters
				SAFE_DELETE(f);
				progressDlg->SetMessage("Loading textures...");
				LoadTexturesGEO(toOpen, version);
				strcpy(fullFileName, fileName);
				//if (isGEO7Z) remove(tmp2);
			}
			else { //insert
				mApp->changedSinceSave = TRUE;
				geom->InsertGEO(f, progressDlg, newStr);
				Reload();
			}
		}
		catch (Error &e) {
			if (!insert) geom->Clear();
			SAFE_DELETE(f);
			//if (isGEO7Z) remove(tmp2);
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			throw e;
		}

	}
	else if (ext == "xml" || ext=="zip" ) { //XML file, optionally in ZIP container
		xml_document loadXML;
		xml_parse_result parseResult;
		progressDlg->SetVisible(TRUE);
		try {
			if (ext=="zip") { //compressed in ZIP container
				//decompress file
				progressDlg->SetMessage("Decompressing file...");

				HZIP hz = OpenZip(fileName, 0);
				if (!hz) {
					throw Error("Can't open ZIP file");
				}
				ZIPENTRY ze; GetZipItem(hz, -1, &ze); int numitems = ze.index;
				BOOL notFoundYet = TRUE;
				for (int i = 0; i < numitems && notFoundYet; i++) { //extract first xml file found in ZIP archive
					GetZipItem(hz, i, &ze);
					std::string zipFileName = ze.name;

					if (FileUtils::GetExtension(zipFileName) == "xml") { //if it's an .xml file
						notFoundYet = FALSE;
						std::string tmpFileName = "tmp/" + zipFileName;
						UnzipItem(hz, i, tmpFileName.c_str()); //unzip it to tmp directory
						CloseZip(hz);
						progressDlg->SetMessage("Reading and parsing XML file...");
						parseResult = loadXML.load_file(tmpFileName.c_str()); //load and parse it
					}
				}
				if (notFoundYet) {
					CloseZip(hz);
					throw Error("Didn't find any XML file in the ZIP file.");
				}
			} else parseResult = loadXML.load_file(fileName); //parse xml file directly
			ResetWorkerStats();
			if (!parseResult) {
				//Parse error
				std::stringstream err;
				err << "XML parsed with errors.\n";
				err << "Error description: " << parseResult.description() << "\n";
				err << "Error offset: " << parseResult.offset << "\n";
				throw Error(err.str().c_str());
			}

			progressDlg->SetMessage("Building geometry...");
			if (!insert) {
				geom->LoadXML_geom(loadXML, this, progressDlg);
				geom->UpdateName(fileName);

				progressDlg->SetMessage("Reloading worker with new geometry...");
				try {
					RealReload(); //for the loading of textures, profiles, etc...
					strcpy(fullFileName, fileName);

					if (ext == "xml" || ext == "zip")
						progressDlg->SetMessage("Restoring simulation state...");
					geom->LoadXML_simustate(loadXML, dpHit, this, progressDlg);
					SendHits(TRUE); //Send hits without sending facet counters, as they are directly written during the load process
					RebuildTextures();
				}
				catch (Error &e) {
					GLMessageBox::Display(e.GetMsg(), "Error while loading simulation state", GLDLG_CANCEL, GLDLG_ICONWARNING);
				}
			}
			else { //insert
				geom->InsertXML(loadXML, this, progressDlg, newStr);
				mApp->changedSinceSave = TRUE;
				ResetWorkerStats();
				Reload();
			}
		}
		catch (Error &e) {
			if (!insert) geom->Clear();
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			throw e;
		}

	}
	else if (ext=="ase" || ext=="ASE") {
		if (insert) throw Error("ASE file inserting is not supported.");
		try {
			ResetWorkerStats();
			f = new FileReader(fileName);
			progressDlg->SetVisible(TRUE);
			geom->LoadASE(f, progressDlg);
			//RealReload();
			strcpy(fullFileName, fileName);
		}
		catch (Error &e) {
			geom->Clear();
			SAFE_DELETE(f);
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			throw e;
		}

	}
	else {

		SAFE_DELETE(f);
		progressDlg->SetVisible(FALSE);
		SAFE_DELETE(progressDlg);
		throw Error("LoadGeometry(): Invalid file extension [Only xml,zip,geo,geo7z,syn.syn7z,txt,ase,stl or str]");
	}
	//geom->CalcTotalOutGassing();
	if (!insert)
	{
		CalcTotalOutgassing();
		if (mApp->momentsEditor) mApp->momentsEditor->Refresh();
		if (mApp->parameterEditor) mApp->parameterEditor->UpdateCombo();
		if (mApp->timeSettings) mApp->timeSettings->RefreshMoments();
		if (mApp->timewisePlotter) mApp->timewisePlotter->refreshViews();
	}
	progressDlg->SetVisible(FALSE);
	SAFE_DELETE(progressDlg);
	SAFE_DELETE(f);
	if (insert) {
		mApp->UpdateFacetlistSelected();
		mApp->UpdateViewers();
	}
}

void Worker::GetLeak(LEAK *buffer, int *nb) {

	*nb = 0;
	if (dpHit) {
		memcpy(buffer, leakCache, sizeof(LEAK)*NBHLEAK);
		*nb = (int)MIN(nbLastLeaks, NBHLEAK);
	}

}


void Worker::SetLeak(LEAK *buffer, int *nb, SHGHITS *gHits) { //When loading from file

	if (nb > 0) {
		memcpy(leakCache, buffer, sizeof(LEAK)*MIN(NBHLEAK, *nb));
		memcpy(gHits->pLeak, buffer, sizeof(LEAK)*MIN(NBHLEAK, *nb));
		gHits->nbLastLeaks = *nb;
	}

}


void Worker::LoadTexturesGEO(std::string fileName, int version) {

	if (FileUtils::GetExtension(fileName) == "geo") {
		GLProgress *progressDlg = new GLProgress("Loading textures", "Please wait");
		progressDlg->SetProgress(0.0);
		FileReader *f = NULL;
		try {
			f = new FileReader(fileName);
			progressDlg->SetVisible(TRUE);
			geom->LoadTextures(f, progressDlg, dpHit, version);
			RebuildTextures();
		}
		catch (Error &e) {
			char tmp[256];
			sprintf(tmp, "Couldn't load some textures. To avoid continuing a partially loaded state, it is recommended to reset the simulation.\n%s", e.GetMsg());
			GLMessageBox::Display(tmp, "Error while loading textures.", GLDLG_OK, GLDLG_ICONWARNING);
		}
		progressDlg->SetVisible(FALSE);
		SAFE_DELETE(progressDlg);
		SAFE_DELETE(f);
	}
}

void Worker::GetHHit(HIT *buffer, int *nb) {

	*nb = 0;
	if (dpHit) {
		memcpy(buffer, hhitCache, sizeof(HIT)*NBHHIT);
		*nb = nbHHit;
	}

}

// -------------------------------------------------------------

void Worker::SetHHit(HIT *buffer, int *nb, SHGHITS *gHits) {

	memcpy(hhitCache, buffer, sizeof(HIT)*MIN(*nb, NBHHIT));
	memcpy(gHits->pHits, buffer, sizeof(HIT)*MIN(*nb, NBHHIT));
	gHits->nbHHit = *nb;

}

// -------------------------------------------------------------

void Worker::InnerStop(float appTime) {

	stopTime = appTime;
	simuTime += appTime - startTime;
	running = FALSE;
	calcAC = FALSE;

}

// -------------------------------------------------------------

void Worker::OneStep() {

	if (nbProcess == 0)
		throw Error("No sub process found. (Simulation not available)");

	if (!running)  {
		if (!ExecuteAndWait(COMMAND_STEPAC, PROCESS_RUN, AC_MODE))
			ThrowSubProcError();
	}

}

// -------------------------------------------------------------

void Worker::StepAC(float appTime) {

	try {
		OneStep();
		Update(appTime);
	}
	catch (Error &e) {
		GLMessageBox::Display((char *)e.GetMsg(), "Error (Stop)", GLDLG_OK, GLDLG_ICONERROR);
	}

}

void Worker::Stop_Public() {
	// Stop
	InnerStop(m_fTime);
	try {
		Stop();
		Update(m_fTime);
	}
	catch (Error &e) {
		GLMessageBox::Display((char *)e.GetMsg(), "Error (Stop)", GLDLG_OK, GLDLG_ICONERROR);
	}
}

// -------------------------------------------------------------

void Worker::StartStop(float appTime, int mode) {

	if (running)  {

		// Stop
		InnerStop(appTime);
		try {
			Stop();
			Update(appTime);
		}
		catch (Error &e) {
			GLMessageBox::Display((char *)e.GetMsg(), "Error (Stop)", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}

	}
	else {

		// Start
		try {
			if (needsReload) RealReload(); //Synchronize subprocesses to main process
			startTime = appTime;
			running = TRUE;
			calcAC = FALSE;
			this->mode = mode;
			Start();

			// Particular case when simulation ends before getting RUN state
			if (allDone) {
				Update(appTime);
				GLMessageBox::Display("Max desorption reached", "Information (Start)", GLDLG_OK, GLDLG_ICONINFO);
			}
		}
		catch (Error &e) {
			running = FALSE;
			GLMessageBox::Display((char *)e.GetMsg(), "Error (Start)", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}

	}

}

// -------------------------------------------------------------

void Worker::Update(float appTime) {
	if (!needsReload) {
		// Check calculation ending
		BOOL done = TRUE;
		BOOL error = TRUE;
		if (dpControl) {
			if (AccessDataport(dpControl)) {
				int i = 0;
				SHMASTER *master = (SHMASTER *)dpControl->buff;
				for (int proc = 0; proc < nbProcess && done; proc++) {
					done = done && (master->states[proc] == PROCESS_DONE);
					error = error && (master->states[proc] == PROCESS_ERROR);
					if (master->states[proc] == PROCESS_RUNAC) calcACprg = master->cmdParam[proc];
				}
				ReleaseDataport(dpControl);
			}
		}

		// End of simulation reached (Stop GUI)
		if ((error || done) && running && appTime != 0.0f) {
			InnerStop(appTime);
			if (error) ThrowSubProcError();
		}

		// Retrieve hit count recording from the shared memory
		if (dpHit) {

			if (AccessDataport(dpHit)) {
				BYTE *buffer = (BYTE *)dpHit->buff;

				mApp->changedSinceSave = TRUE;
				// Globals
				SHGHITS *gHits = (SHGHITS *)buffer;

				// Copy Global hits and leaks
				nbHit = gHits->total.hit.nbHit;
				nbAbsorption = gHits->total.hit.nbAbsorbed;
				distTraveledTotal_total = gHits->distTraveledTotal_total;
				distTraveledTotal_fullHitsOnly = gHits->distTraveledTotal_fullHitsOnly;
				nbDesorption = gHits->total.hit.nbDesorbed;
				nbLeakTotal = gHits->nbLeakTotal;
				nbHHit = gHits->nbHHit;
				nbLastLeaks = gHits->nbLastLeaks;
				memcpy(hhitCache, gHits->pHits, sizeof(HIT)*NBHHIT);
				memcpy(leakCache, gHits->pLeak, sizeof(LEAK)*NBHHIT);

				// Refresh local facet hit cache for the displayed moment
				int nbFacet = geom->GetNbFacet();
				for (int i = 0; i < nbFacet; i++) {
					Facet *f = geom->GetFacet(i);
					f->counterCache=(*((SHHITS*)(buffer + f->sh.hitOffset+displayedMoment*sizeof(SHHITS))));
				}
				try {
					geom->BuildTexture(buffer);
				}
				catch (Error &e) {
					GLMessageBox::Display((char *)e.GetMsg(), "Error building texture", GLDLG_OK, GLDLG_ICONERROR);
					ReleaseDataport(dpHit);
					return;
				}
				ReleaseDataport(dpHit);
			}

		}
	}
}

// -------------------------------------------------------------

void Worker::SendHits(BOOL skipFacetHits ) {
	//if (!needsReload) {
	if (dpHit) {
		if (AccessDataport(dpHit)) {

			SHGHITS *gHits = (SHGHITS *)dpHit->buff;
			gHits->total.hit.nbHit = nbHit;
			gHits->nbLeakTotal = nbLeakTotal;
			gHits->total.hit.nbDesorbed = nbDesorption;
			gHits->total.hit.nbAbsorbed = nbAbsorption;
			gHits->distTraveledTotal_total = distTraveledTotal_total;
			gHits->distTraveledTotal_fullHitsOnly = distTraveledTotal_fullHitsOnly;

			if (!skipFacetHits) {
				int nbFacet = geom->GetNbFacet();
				for (int i = 0; i < nbFacet; i++) {
					Facet *f = geom->GetFacet(i);
					for (size_t m = 0;m <= moments.size();m++)
						*((SHHITS*)((BYTE*)dpHit->buff + f->sh.hitOffset + m * sizeof(SHHITS))) = f->counterCache;
				}
			}
			ReleaseDataport(dpHit);

		}
		else {
			throw Error("Failed to initialize 'hits' dataport");
		}
	}
}

// -------------------------------------------------------------

void Worker::ComputeAC(float appTime) {
	try {
		if (needsReload) RealReload();
	}
	catch (Error &e) {
		GLMessageBox::Display((char *)e.GetMsg(), "Error (Stop)", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}
	if (running)
		throw Error("Already running");

	// Send correction map to sub process
	// (correction map contains area a surface elements)
	size_t maxElem = geom->GetMaxElemNumber();
	if (!maxElem)
		throw Error("Mesh with boundary correction must be enabled on all polygons");
	int dpSize = maxElem*sizeof(SHELEM);

	Dataport *loader = CreateDataport(loadDpName, dpSize);
	if (!loader)
		throw Error("Failed to create 'loader' dataport");
	AccessDataport(loader);
	geom->CopyElemBuffer((BYTE *)loader->buff);
	ReleaseDataport(loader);

	// Load Elem area and send AC matrix calculation order
	// Send command
	if (!ExecuteAndWait(COMMAND_LOADAC, PROCESS_RUNAC, dpSize)) {
		CLOSEDP(loader);
		char errMsg[1024];
		sprintf(errMsg, "Failed to send geometry to sub process:\n%s", GetErrorDetails());
		GLMessageBox::Display(errMsg, "Warning (LoadAC)", GLDLG_OK, GLDLG_ICONWARNING);
		return;
	}

	CLOSEDP(loader);

	running = TRUE;
	calcAC = TRUE;
	startTime = appTime;

}

// -------------------------------------------------------------

void  Worker::ReleaseHits() {

	ReleaseDataport(dpHit);

}

// -------------------------------------------------------------

BYTE *Worker::GetHits() {
	try {
		if (needsReload) RealReload();
	}
	catch (Error &e) {
		GLMessageBox::Display((char *)e.GetMsg(), "Error (Stop)", GLDLG_OK, GLDLG_ICONERROR);
	}
	if (dpHit)
		if (AccessDataport(dpHit))
			return (BYTE *)dpHit->buff;

	return NULL;

}

// -------------------------------------------------------------

void Worker::ThrowSubProcError(char *message) {

	char errMsg[1024];
	if (!message)
		sprintf(errMsg, "Bad response from sub process(es):\n%s", GetErrorDetails());
	else
		sprintf(errMsg, "%s\n%s", message, GetErrorDetails());
	throw Error(errMsg);

}

void Worker::Reload(){
	needsReload = true;
}

// -------------------------------------------------------------

void Worker::RealReload() { //Sharing geometry with workers

	//Do preliminary calculations
	PrepareToRun();

	if (nbProcess == 0) return;

	GLProgress *progressDlg = new GLProgress("Asking subprocesses to clear geometry...", "Passing Geometry to workers");
	progressDlg->SetVisible(TRUE);
	progressDlg->SetProgress(0.0);

	// Clear geometry
	CLOSEDP(dpHit);
	if (!ExecuteAndWait(COMMAND_CLOSE, PROCESS_READY))
		ThrowSubProcError();

	if (!geom->IsLoaded()) {
		progressDlg->SetVisible(FALSE);
		SAFE_DELETE(progressDlg);
		return;
	}

	// Create the temporary geometry shared structure
	progressDlg->SetMessage("Creating dataport...");
	size_t loadSize = geom->GetGeometrySize();
	Dataport *loader = CreateDataport(loadDpName, loadSize);
	if( !loader )
		throw Error("Failed to create 'loader' dataport.\nMost probably out of memory.\nReduce number of subprocesses or texture size.");
	progressDlg->SetMessage("Accessing dataport...");
	AccessDataportTimed(loader, 3000 + nbProcess*(int)((double)loadSize / 10000.0));
	progressDlg->SetMessage("Assembling geometry to pass...");
	geom->CopyGeometryBuffer((BYTE *)loader->buff);
	progressDlg->SetMessage("Releasing dataport...");
	ReleaseDataport(loader);

	size_t hitSize = geom->GetHitsSize(&moments);
	dpHit = CreateDataport(hitsDpName, hitSize);ClearHits(TRUE);
	if (!dpHit) {
		CLOSEDP(loader);
		//GLMessageBox::Display("Failed to create 'hits' dataport: not enough memory.", "Warning (Load)", GLDLG_OK, GLDLG_ICONERROR);
		//return FALSE;
		progressDlg->SetVisible(FALSE);
		SAFE_DELETE(progressDlg);
		throw Error("Failed to create 'hits' dataport: out of memory.");
	}

	// Compute number of max desorption per process
	if (AccessDataportTimed(dpControl, 3000 + nbProcess*(int)((double)loadSize / 10000.0))) {
		SHMASTER *m = (SHMASTER *)dpControl->buff;
		llong common = maxDesorption / (llong)nbProcess;
		int remain = (int)(maxDesorption % (llong)nbProcess);
		for (int i = 0; i < nbProcess; i++) {
			m->cmdParam2[i] = common;
			if (i < remain) m->cmdParam2[i]++;
		}
		ReleaseDataport(dpControl);
	}

	// Load geometry
	progressDlg->SetMessage("Waiting for subprocesses to load geometry...");
	if (!ExecuteAndWait(COMMAND_LOAD, PROCESS_READY, loadSize, progressDlg)) {
		CLOSEDP(loader);
		char errMsg[1024];
		sprintf(errMsg, "Failed to send geometry to sub process:\n%s", GetErrorDetails());
		//GLMessageBox::Display(errMsg, "Warning (Load)", GLDLG_OK, GLDLG_ICONWARNING);
		//return FALSE;
		progressDlg->SetVisible(FALSE);
		SAFE_DELETE(progressDlg);
		throw Error(errMsg);
	}

	//Old send hits location

	progressDlg->SetMessage("Closing dataport...");
	CLOSEDP(loader);
	needsReload = false;
	progressDlg->SetVisible(FALSE);
	SAFE_DELETE(progressDlg);
}

// -------------------------------------------------------------

void Worker::SetMaxDesorption(llong max) {

	try {
		Reset(0.0);
		maxDesorption = max;
		Reload();
	}
	catch (Error &e) {
		GLMessageBox::Display((char *)e.GetMsg(), "Error", GLDLG_OK, GLDLG_ICONERROR);
	}

}

// -------------------------------------------------------------

char *Worker::GetErrorDetails() {

	static char err[1024];
	strcpy(err, "");

	AccessDataport(dpControl);
	SHMASTER *shMaster = (SHMASTER *)dpControl->buff;
	for (int i = 0; i < nbProcess; i++) {
		char tmp[256];
		if (pID[i] != 0) {
			int st = shMaster->states[i];
			if (st == PROCESS_ERROR) {
				sprintf(tmp, "[#%d] Process [PID %d] %s: %s\n", i, pID[i], prStates[st], shMaster->statusStr[i]);
			}
			else {
				sprintf(tmp, "[#%d] Process [PID %d] %s\n", i, pID[i], prStates[st]);
			}
		}
		else {
			sprintf(tmp, "[#%d] Process [PID ???] Not started\n", i);
		}
		strcat(err, tmp);
	}
	ReleaseDataport(dpControl);

	return err;
}

// -------------------------------------------------------------

void Worker::ClearHits(BOOL noReload) {
	try {
		if (!noReload && needsReload) RealReload();
	}
	catch (Error &e) {
		GLMessageBox::Display((char *)e.GetMsg(), "Error (Stop)", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}
	if (dpHit) {
		AccessDataport(dpHit);
		memset(dpHit->buff, 0, geom->GetHitsSize(&moments));
		ReleaseDataport(dpHit);
	}

}

// -------------------------------------------------------------

BOOL Worker::Wait(int waitState, int timeout,GLProgress *prg) {

	BOOL ok = FALSE;
	BOOL error = FALSE;
	int t = 0;
		int nbReady = 0;
	double initialProgress = 0.0;
	if (prg) initialProgress = prg->GetProgress();
	int nbError = 0;
	allDone = TRUE;

	// Wait for completion
	while (!ok && t < timeout) {

		ok = TRUE;
		AccessDataport(dpControl);
		SHMASTER *shMaster = (SHMASTER *)dpControl->buff;
		nbReady=nbError=0;
		for (int i = 0; i < nbProcess; i++) {
			if (shMaster->states[i]==waitState) nbReady++;
			ok = ok & (shMaster->states[i] == waitState || shMaster->states[i] == PROCESS_ERROR || shMaster->states[i] == PROCESS_DONE);
			if( shMaster->states[i]==PROCESS_ERROR ) {
				error = TRUE;
				nbError++;
			}
			allDone = allDone & (shMaster->states[i] == PROCESS_DONE);
		}
		ReleaseDataport(dpControl);

		if (!ok) {
			if (prg) prg->SetProgress(double(nbReady)/(double)nbProcess);
			Sleep(500);
			t += 500;
		}

	}

	if (t>=timeout) {
		if ((prg) && ((double)nbReady/(double)nbProcess)>initialProgress) //progress advanced, wait more
			return Wait(waitState,timeout,prg);
		char tmp[512];
		sprintf(tmp,"Total workers : %d\n"
			"%d are ready, %d reported errors\n"
			"Do you want to wait a bit more?\n"
			"(Loading continues while this dialog is visible)\n",nbProcess,nbReady,nbError);
		int waitmore = GLMessageBox::Display(tmp,"Info",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONINFO)==GLDLG_OK;
		if (waitmore) {
			t=0;
			return Wait(waitState,timeout,prg);
		}
	}

	return ok && !error;

}

BOOL Worker::ExecuteAndWait(int command, int waitState, size_t param,GLProgress *prg) {
	/*static BOOL alreadyreloading=FALSE;
	if (needsReload && !alreadyreloading) {
	alreadyreloading = TRUE;
	RealReload();
	alreadyreloading = FALSE;
	}*/
	if (!dpControl) return FALSE;

	// Send command
	AccessDataport(dpControl);
	SHMASTER *shMaster = (SHMASTER *)dpControl->buff;
	for (int i = 0; i < nbProcess; i++) {
		shMaster->states[i] = command;
		shMaster->cmdParam[i] = param;
	}
	ReleaseDataport(dpControl);

	Sleep(100);
	return Wait(waitState, 3000 + nbProcess * 500,prg);
}

// -------------------------------------------------------------

/*
void Worker::SendHeartBeat() {
if(!dpControl) return;

// Send command
AccessDataport(dpControl);
SHMASTER *shMaster = (SHMASTER *)dpControl->buff;
//int hb = shMaster->heartBeat;
//hb++;
//if (hb==1000) hb=0;
shMaster->heartBeat = m_fTime;
ReleaseDataport(dpControl);
}
*/

// -------------------------------------------------------------

void Worker::ResetWorkerStats() {

	nbAbsorption = 0;
	nbDesorption = 0;
	nbHit = 0;
	nbLeakTotal = 0;
	distTraveledTotal_total = 0.0;
	distTraveledTotal_fullHitsOnly = 0.0;

}

// -------------------------------------------------------------

void Worker::Reset(float appTime) {

	if (calcAC) {
		GLMessageBox::Display("Reset not allowed while calculating AC", "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}

	stopTime = 0.0f;
	startTime = 0.0f;
	simuTime = 0.0f;
	running = FALSE;
	if (nbProcess == 0)
		return;

	try {
		ResetWorkerStats();
		if (!ExecuteAndWait(COMMAND_RESET, PROCESS_READY))
			ThrowSubProcError();
		ClearHits(FALSE);
		Update(appTime);
	}
	catch (Error &e) {
		GLMessageBox::Display((char *)e.GetMsg(), "Error", GLDLG_OK, GLDLG_ICONERROR);
	}

	//Debug memory check
	//_ASSERTE (!_CrtDumpMemoryLeaks());;
	_ASSERTE(_CrtCheckMemory());

}

// -------------------------------------------------------------

void Worker::Start() {

	// Check that at least one desortion facet exists
	BOOL found = FALSE;
	int nbF = geom->GetNbFacet();
	int i = 0;
	while (i<nbF && !found) {
		found = (geom->GetFacet(i)->sh.desorbType != DES_NONE);
		if (!found) i++;
	}

	if (!found)
		throw Error("No desorption facet found");

	if (!(totalDesorbedMolecules>0.0))
		throw Error("Total outgassing is zero.");

	if (nbProcess == 0)
		throw Error("No sub process found. (Simulation not available)");

	if (!ExecuteAndWait(COMMAND_START, PROCESS_RUN, mode))
		ThrowSubProcError();

	//Debug memory check
	//_ASSERTE (!_CrtDumpMemoryLeaks());;
	_ASSERTE(_CrtCheckMemory());

}

// -------------------------------------------------------------

void Worker::Stop() {

	if (nbProcess == 0)
		throw Error("No sub process found. (Simulation not available)");

	if (!ExecuteAndWait(COMMAND_PAUSE, PROCESS_READY))
		ThrowSubProcError();

	//Debug memory check
	//_ASSERTE (!_CrtDumpMemoryLeaks());;
	_ASSERTE(_CrtCheckMemory());

}

// -------------------------------------------------------------

void Worker::KillAll() {

	if (dpControl && nbProcess > 0) {
		if (!ExecuteAndWait(COMMAND_EXIT, PROCESS_KILLED)) {
			AccessDataport(dpControl);
			SHMASTER *shMaster = (SHMASTER *)dpControl->buff;
			for (int i = 0; i < nbProcess; i++)
				if (shMaster->states[i] == PROCESS_KILLED) pID[i] = 0;
			ReleaseDataport(dpControl);
			// Force kill
			for (int i = 0; i < nbProcess; i++)
				if (pID[i]) KillProc(pID[i]);
		}
		CLOSEDP(dpHit);
	}
	nbProcess = 0;

}

// -------------------------------------------------------------

void Worker::Exit() {

	if (dpControl && nbProcess > 0) {
		KillAll();
		CLOSEDP(dpControl);
	}

}

// -------------------------------------------------------------

void Worker::GetProcStatus(int *states, char **status) {

	if (nbProcess == 0) return;

	AccessDataport(dpControl);
	SHMASTER *shMaster = (SHMASTER *)dpControl->buff;
	memcpy(states, shMaster->states, MAX_PROCESS*sizeof(int));
	memcpy(status, shMaster->statusStr, MAX_PROCESS * 64);
	ReleaseDataport(dpControl);

}



void Worker::SetProcNumber(int n) {

	char cmdLine[512];

	// Kill all sub process
	KillAll();

	// Create new control dataport
	if (!dpControl)
		dpControl = CreateDataport(ctrlDpName, sizeof(SHMASTER));
	if (!dpControl)
		throw Error("Failed to create 'control' dataport");
	AccessDataport(dpControl);
	memset(dpControl->buff, 0, sizeof(SHMASTER));
	ReleaseDataport(dpControl);

	// Launch n subprocess
	for (int i = 0; i < n; i++) {
		sprintf(cmdLine, "molflowSub.exe %d %d", pid, i);
		pID[i] = StartProc(cmdLine);
		Sleep(25); // Wait a bit
		if (pID[i] == 0) {
			nbProcess = 0;
			throw Error(cmdLine);
		}
	}

	nbProcess = n;

	if (!Wait(PROCESS_READY, 3000)){
		ThrowSubProcError("Sub process(es) starting failure");
	}
}

// -------------------------------------------------------------

int Worker::GetProcNumber() {
	return nbProcess;
}

// -------------------------------------------------------------

DWORD Worker::GetPID(int prIdx) {
	return pID[prIdx];
}

/*
std::string execCMD(char* cmd) {
FILE* pipe = _popen(cmd, "r");
if (!pipe) return "ERROR";
char buffer[128];
std::string result = "";
while(!feof(pipe)) {
if(fgets(buffer, 128, pipe) != NULL)
result += buffer;
}
_pclose(pipe);
return result;
}
*/

int Worker::AddMoment(std::vector<double> newMoments) {
	int nb = (int)newMoments.size();
	for (int i = 0; i < nb; i++)
		moments.push_back(newMoments[i]);
	return nb;
}

std::vector<double> Worker::ParseMoment(std::string userInput) {
	std::vector<double> parsedResult;
	double begin, interval, end;

	int nb = sscanf(userInput.c_str(), "%lf,%lf,%lf", &begin, &interval, &end);
	if (nb == 1 && (begin >= 0.0)) {
		//One moment
		parsedResult.push_back(begin);
		//} else if (nb==3 && (begin>0.0) && (end>begin) && (interval<(end-begin)) && ((end-begin)/interval<300.0)) {
	}
	else if (nb == 3 && (begin >= 0.0) && (end > begin) && (interval < (end - begin))) {
		//Range
		for (double time = begin; time <= end; time += interval)
			parsedResult.push_back(time);
	}
	return parsedResult;
}

void Worker::ResetMoments() {
	moments = std::vector<double>();
	userMoments = std::vector<std::string>();
}

void Worker::ImportDesorption_DES(char *fileName) {
	//if (needsReload) RealReload();
	// Read a file
	FileReader *f = new FileReader(fileName);
	geom->ImportDesorption_DES(f);
	SAFE_DELETE(f);
	mApp->changedSinceSave = TRUE;
	Reload();
}

void Worker::ImportDesorption_SYN(char *fileName, const size_t &source, const double &time,
	const size_t &mode, const double &eta0, const double &alpha, const double &cutoffdose,
	const std::vector<std::pair<double, double>> &convDistr,
	GLProgress *prg) {
	char *ext, *filebegin;
	char tmp2[1024];
	ext = strrchr(fileName, '.');
	char CWD[MAX_PATH];
	_getcwd(CWD, MAX_PATH);
	if (ext == NULL || (!(_stricmp(ext, ".syn7z") == 0)) && (!(_stricmp(ext, ".syn") == 0)))
		throw Error("ImportDesorption_SYN(): Invalid file extension [Only syn, syn7z]");
	ext++;

	// Read a file
	FileReader *f = NULL;

	GLProgress *progressDlg = new GLProgress("Analyzing SYN file...", "Please wait");
	progressDlg->SetProgress(0.0);
	progressDlg->SetVisible(TRUE);
	BOOL isSYN7Z = (_stricmp(ext, "syn7z") == 0);
	BOOL isSYN = (_stricmp(ext, "syn") == 0);

	if (isSYN || isSYN7Z) {
		progressDlg->SetVisible(TRUE);
		try {
			if (isSYN7Z) {
				//decompress file
				progressDlg->SetMessage("Decompressing file...");
				char tmp[1024];
				char fileOnly[512];
				sprintf(tmp, "cmd /C \"pushd \"%s\"&&7za.exe x -t7z -aoa \"%s\" -otmp&&popd\"", CWD, fileName);
				system(tmp);

				filebegin = strrchr(fileName, '\\');
				if (filebegin) filebegin++;
				else filebegin = fileName;
				memcpy(fileOnly, filebegin, sizeof(char)*(strlen(filebegin) - 2)); //remove ..7z from extension
				fileOnly[strlen(filebegin) - 2] = '\0';
				sprintf(tmp2, "%s\\tmp\\Geometry.syn", CWD);
				f = new FileReader(tmp2); //decompressed file opened
			}

			if (!isSYN7Z) f = new FileReader(fileName);  //original file opened

			geom->ImportDesorption_SYN(f, source, time, mode, eta0, alpha, cutoffdose, convDistr, prg);
			CalcTotalOutgassing();
			SAFE_DELETE(f);

		}
		catch (Error &e) {
			SAFE_DELETE(f);
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			throw e;
		}
		progressDlg->SetVisible(FALSE);
		SAFE_DELETE(progressDlg);
	}
}

void Worker::AnalyzeSYNfile(char *fileName, int *nbFacet, int *nbTextured, int *nbDifferent) {
	char *ext, *filebegin;
	char tmp2[1024];
	ext = strrchr(fileName, '.');
	char CWD[MAX_PATH];
	_getcwd(CWD, MAX_PATH);
	if (ext == NULL || (!(_stricmp(ext, ".syn7z") == 0)) && (!(_stricmp(ext, ".syn") == 0)))
		throw Error("AnalyzeSYNfile(): Invalid file extension [Only syn, syn7z]");
	ext++;

	// Read a file
	FileReader *f = NULL;

	GLProgress *progressDlg = new GLProgress("Analyzing SYN file...", "Please wait");
	progressDlg->SetProgress(0.0);
	progressDlg->SetVisible(TRUE);
	BOOL isSYN7Z = (_stricmp(ext, "syn7z") == 0);
	BOOL isSYN = (_stricmp(ext, "syn") == 0);

	if (isSYN || isSYN7Z) {
		progressDlg->SetVisible(TRUE);
		try {
			if (isSYN7Z) {
				//decompress file
				progressDlg->SetMessage("Decompressing file...");
				char tmp[1024];
				char fileOnly[512];
				sprintf(tmp, "cmd /C \"pushd \"%s\"&&7za.exe x -t7z -aoa \"%s\" -otmp&&popd\"", CWD, fileName);
				system(tmp);

				filebegin = strrchr(fileName, '\\');
				if (filebegin) filebegin++;
				else filebegin = fileName;
				memcpy(fileOnly, filebegin, sizeof(char)*(strlen(filebegin) - 2)); //remove ..7z from extension
				fileOnly[strlen(filebegin) - 2] = '\0';
				sprintf(tmp2, "%s\\tmp\\Geometry.syn", CWD);
				f = new FileReader(tmp2); //decompressed file opened
			}

			if (!isSYN7Z) f = new FileReader(fileName);  //original file opened

			geom->AnalyzeSYNfile(f, progressDlg, nbFacet, nbTextured, nbDifferent, progressDlg);

			SAFE_DELETE(f);

		}
		catch (Error &e) {
			SAFE_DELETE(f);
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			throw e;
		}
		progressDlg->SetVisible(FALSE);
		SAFE_DELETE(progressDlg);
	}
}

void Worker::PrepareToRun() {

	//determine latest moment
	latestMoment = 1E-10;
	for (size_t i = 0; i<moments.size(); i++)
		if (moments[i]>latestMoment) latestMoment = moments[i];
	latestMoment += timeWindowSize / 2.0;

	Geometry *g = GetGeometry();
	//Generate integrated desorption functions

	temperatures = std::vector<double>();
	desorptionParameterIDs = std::vector<size_t>();
	CDFs = std::vector<std::vector<std::pair<double, double>>>();
	IDs = std::vector<std::vector<std::pair<double, double>>>();

	for (int i = 0; i < g->GetNbFacet(); i++) {
		Facet *f = g->GetFacet(i);

		//match parameters
		if (f->userOutgassing.length() > 0) {
			int id = GetParamId(f->userOutgassing);
			if (id == -1) { //parameter not found
				char tmp[256];
				sprintf(tmp, "Facet #%d: Outgassing parameter \"%s\" isn't defined.", i + 1, f->userOutgassing.c_str());
				throw Error(tmp);
			}
			else f->sh.outgassing_paramId = id;
		}
		else f->sh.outgassing_paramId = -1;

		if (f->userOpacity.length() > 0) {
			int id = GetParamId(f->userOpacity);
			if (id == -1) { //parameter not found
				char tmp[256];
				sprintf(tmp, "Facet #%d: Opacity parameter \"%s\" isn't defined.", i + 1, f->userOpacity.c_str());
				throw Error(tmp);
			}
			else f->sh.opacity_paramId = id;
		}
		else f->sh.opacity_paramId = -1;

		if (f->userSticking.length() > 0) {
			int id = GetParamId(f->userSticking);
			if (id == -1) { //parameter not found
				char tmp[256];
				sprintf(tmp, "Facet #%d: Sticking parameter \"%s\" isn't defined.", i + 1, f->userSticking.c_str());
				throw Error(tmp);
			}
			else f->sh.sticking_paramId = id;
		}
		else f->sh.sticking_paramId = -1;

		if (f->sh.outgassing_paramId >= 0) { //if time-dependent desorption
			int id = GetIDId(f->sh.outgassing_paramId);
			if (id >= 0)
				f->sh.IDid = id; //we've already generated an ID for this temperature
			else
				f->sh.IDid = GenerateNewID(f->sh.outgassing_paramId);
		}

		//Generate speed distribution functions
		int id = GetCDFId(f->sh.temperature);
		if (id >= 0)
			f->sh.CDFid = id; //we've already generated a CDF for this temperature
		else
			f->sh.CDFid = GenerateNewCDF(f->sh.temperature);

	}

	CalcTotalOutgassing();
}


int Worker::GetCDFId(double temperature) {

	int i;
	for (i = 0; i<(int)temperatures.size() && (abs(temperature - temperatures[i])>1E-5); i++); //check if we already had this temperature
	if (i >= (int)temperatures.size()) i = -1; //not found
	return i;
}

int Worker::GenerateNewCDF(double temperature){
	size_t i = temperatures.size();
	temperatures.push_back(temperature);
	CDFs.push_back(Generate_CDF(temperature, gasMass, CDF_SIZE));
	return (int)i;
}

int Worker::GenerateNewID(int paramId){
	size_t i = desorptionParameterIDs.size();
	desorptionParameterIDs.push_back(paramId);
	IDs.push_back(Generate_ID(paramId));
	return (int)i;
}

int Worker::GetIDId(int paramId) {

	int i;
	for (i = 0; i < (int)desorptionParameterIDs.size() && (paramId != desorptionParameterIDs[i]); i++); //check if we already had this parameter Id
	if (i >= (int)desorptionParameterIDs.size()) i = -1; //not found
	return i;
}

void Worker::CalcTotalOutgassing() {
	// Compute the outgassing of all source facet
	totalDesorbedMolecules = finalOutgassingRate_Pa_m3_sec = finalOutgassingRate = 0.0;
	Geometry *g = GetGeometry();


	for (int i = 0; i < g->GetNbFacet(); i++) {
		Facet *f = g->GetFacet(i);
		if (f->sh.desorbType != DES_NONE) { //there is a kind of desorption
			if (f->sh.useOutgassingFile) { //outgassing file
				for (int l = 0; l < (f->sh.outgassingMapWidth*f->sh.outgassingMapHeight); l++) {
					totalDesorbedMolecules += latestMoment * f->outgassingMap[l] / (1.38E-23*f->sh.temperature);
					finalOutgassingRate += f->outgassingMap[l] / (1.38E-23*f->sh.temperature);
					finalOutgassingRate_Pa_m3_sec += f->outgassingMap[l];
				}
			}
			else { //regular outgassing
				if (f->sh.outgassing_paramId == -1) { //constant outgassing
					totalDesorbedMolecules += latestMoment * f->sh.flow / (1.38E-23*f->sh.temperature);
					finalOutgassingRate += f->sh.flow / (1.38E-23*f->sh.temperature);  //Outgassing molecules/sec
					finalOutgassingRate_Pa_m3_sec += f->sh.flow;
				}
				else { //time-dependent outgassing
					totalDesorbedMolecules += IDs[f->sh.IDid].back().second / (1.38E-23*f->sh.temperature);
					finalOutgassingRate += parameters[f->sh.outgassing_paramId].values.back().second *0.100 / (1.38E-23*f->sh.temperature); //0.1: mbar*l/s->Pa*m3/s
					finalOutgassingRate_Pa_m3_sec += parameters[f->sh.outgassing_paramId].values.back().second *0.100;
				}
			}
		}
	}
	if (mApp->globalSettings) mApp->globalSettings->UpdateOutgassing();
}


std::vector<std::pair<double, double>> Worker::Generate_CDF(double gasTempKelvins, double gasMassGramsPerMol, size_t size){
	std::vector<std::pair<double, double>> cdf; cdf.reserve(size);
	double Kb = 1.38E-23;
	double R = 8.3144621;
	double a = sqrt(Kb*gasTempKelvins / (gasMassGramsPerMol*1.67E-27)); //distribution a parameter. Converting molar mass to atomic mass

	//Generate cumulative distribution function
	double mostProbableSpeed = sqrt(2 * R*gasTempKelvins / (gasMassGramsPerMol / 1000.0));
	double binSize = 4.0*mostProbableSpeed / (double)size; //distribution generated between 0 and 4*V_prob
	/*double coeff1=1.0/sqrt(2.0)/a;
	double coeff2=sqrt(2.0/PI)/a;
	double coeff3=1.0/(2.0*pow(a,2));

	for (size_t i=0;i<size;i++) {
	double x=(double)i*binSize;
	cdf.push_back(std::make_pair(x,erf(x*coeff1)-coeff2*x*exp(-pow(x,2)*coeff3)));
	}*/
	for (size_t i = 0; i < size; i++) {
		double x = (double)i*binSize;
		double x_square_per_2_a_square = pow(x, 2) / (2 * pow(a, 2));
		cdf.push_back(std::make_pair(x, 1 - exp(-x_square_per_2_a_square)*(x_square_per_2_a_square + 1)));
	}

	/* //UPDATE: not generating inverse since it was introducing sampling problems at the large tail for high speeds
	//CDF created, let's generate its inverse
	std::vector<std::pair<double,double>> inverseCDF;inverseCDF.reserve(size);
	binSize=1.0/(double)size; //Divide probability to bins
	for (size_t i=0;i<size;i++) {
	double p=(double)i*binSize;
	//inverseCDF.push_back(std::make_pair(p,InterpolateX(p,cdf,TRUE)));
	inverseCDF.push_back(std::make_pair(p, InterpolateX(p, cdf, FALSE)));
	}
	return inverseCDF;
	*/
	return cdf;
}

std::vector<std::pair<double, double>> Worker::Generate_ID(int paramId){
	std::vector<std::pair<double, double>> ID;
	//First, let's check at which index is the latest moment
	size_t indexBeforeLastMoment;
	for (indexBeforeLastMoment = 0; indexBeforeLastMoment < parameters[paramId].values.size() &&
		(parameters[paramId].values[indexBeforeLastMoment].first < latestMoment); indexBeforeLastMoment++);
		if (indexBeforeLastMoment >= parameters[paramId].values.size()) indexBeforeLastMoment = parameters[paramId].values.size() - 1; //not found, set as last moment

	//Construct integral from 0 to latest moment
	//Zero
	ID.push_back(std::make_pair(0.0, 0.0));

	//First moment
	ID.push_back(std::make_pair(parameters[paramId].values[0].first,
		parameters[paramId].values[0].first*parameters[paramId].values[0].second*0.100)); //for the first moment (0.1: mbar*l/s -> Pa*m3/s)

	//Intermediate moments
	for (size_t pos = 1; pos <= indexBeforeLastMoment; pos++) {
		if (abs(parameters[paramId].values[pos].second - parameters[paramId].values[pos - 1].second) < 1E-10) //two equal values follow, simple integration by multiplying
			ID.push_back(std::make_pair(parameters[paramId].values[pos].first,
			ID.back().second +
			(parameters[paramId].values[pos].first - parameters[paramId].values[pos - 1].first)*parameters[paramId].values[pos].second*0.100));
		else { //difficult case, we'll integrate by dividing two 20equal sections
			for (double delta = 0.05; delta < 1.0001; delta += 0.05) {
				double delta_t = parameters[paramId].values[pos].first - parameters[paramId].values[pos - 1].first;
				double time = parameters[paramId].values[pos - 1].first + delta*delta_t;
				double avg_value = (InterpolateY(time - 0.05*delta_t, parameters[paramId].values)*0.100
					+ InterpolateY(time, parameters[paramId].values)*0.100) / 2.0;
				ID.push_back(std::make_pair(time,
					ID.back().second +
					0.05*delta_t*avg_value));
			}
		}
	}

	//latestMoment
	double valueAtLatestMoment = InterpolateY(latestMoment, parameters[paramId].values, TRUE);
	if ((valueAtLatestMoment - parameters[paramId].values[indexBeforeLastMoment].second) < 1E-10) //two equal values follow, simple integration by multiplying
		ID.push_back(std::make_pair(latestMoment,
		ID.back().second +
		(latestMoment - parameters[paramId].values[indexBeforeLastMoment].first)*parameters[paramId].values[indexBeforeLastMoment].second*0.100));
	else { //difficult case, we'll integrate by dividing two 5equal sections
		for (double delta = 0.0; delta < 1.0001; delta += 0.05) {
			double delta_t = latestMoment - parameters[paramId].values[indexBeforeLastMoment].first;
			double time = parameters[paramId].values[indexBeforeLastMoment].first + delta*delta_t;
			double avg_value = (parameters[paramId].values[indexBeforeLastMoment].second*0.100 + InterpolateY(time, parameters[paramId].values)*0.100) / 2.0;
			ID.push_back(std::make_pair(time,
				ID.back().second +
				0.05*delta_t*avg_value));
		}
	}

	return ID;
}

int Worker::GetParamId(const std::string name) {
	int foundId = -1;
	for (int i = 0; foundId == -1 && i < (int)parameters.size(); i++)
		if (name.compare(parameters[i].name) == 0) foundId = i;
	return foundId;
}

void Worker::RebuildTextures() {
	if (AccessDataport(dpHit)) {
		BYTE *buffer = (BYTE *)dpHit->buff;
		try{ geom->BuildTexture(buffer); }
		catch (Error &e) {
			ReleaseDataport(dpHit);
			throw e;
		}
		ReleaseDataport(dpHit);
	}
}


/*DWORD Worker::ZipThreadProc(char* fileNameWithXML, char* fileNameWithXMLzip)
{
HZIP hz = CreateZip(fileNameWithXMLzip, 0);
if (!hz) {
throw Error("Error creating ZIP file");
}
if (!ZipAdd(hz, GetShortFileName(fileNameWithXML), fileNameWithXML)) remove(fileNameWithXML);
else {
CloseZip(hz);
throw Error("Error compressing ZIP file.");
}
CloseZip(hz);
return 0;
}*/