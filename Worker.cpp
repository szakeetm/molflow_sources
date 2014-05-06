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

#include "Worker.h"
#include "GLApp/GLApp.h"
#include "GLApp/GLMessageBox.h"
#include <math.h>
#include <stdlib.h>
#include <Process.h>
#include "GLApp/GLUnitDialog.h"
#include "Molflow.h"
#include <direct.h>

/*
//Leak detection
#ifdef _DEBUG
#define _CRTDBG_MAP_ALLOC
#define DEBUG_NEW new(_NORMAL_BLOCK, __FILE__, __LINE__)
#define new DEBUG_NEW
#endif
*/

extern GLApplication *theApp;
extern HANDLE compressProcessHandle;
extern double totalOutgassing;
extern double totalInFlux;
extern int autoSaveSimuOnly;
extern int needsReload;
extern int compressSavedFiles;

// -------------------------------------------------------------

Worker::Worker() {
	//test
	moments=std::vector<double>();
	userMoments=std::vector<std::string>(); //strings describing moments, to be parsed

	displayedMoment=0; //By default, steady-state is displayed

	desorptionStartTime=0.0;
	desorptionStopTime=1;
	timeWindowSize=0.1;
	useMaxwellDistribution=TRUE;
	calcConstantFlow=TRUE;
	valveOpenMoment=99999.0;
	distTraveledTotal=0.0;

	pid = _getpid();
	sprintf(ctrlDpName,"MFLWCTRL%d",pid);
	sprintf(loadDpName,"MFLWLOAD%d",pid);
	sprintf(hitsDpName,"MFLWHITS%d",pid);
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
	strcpy(fullFileName,"");
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

// -------------------------------------------------------------

char *Worker::GetFileName() {
	return fullFileName;
}

char *Worker::GetShortFileName() {

	static char ret[512];
	char *r = strrchr(fullFileName,'/');
	if(!r) r = strrchr(fullFileName,'\\');
	if(!r) strcpy(ret,fullFileName);
	else   {
		r++;
		strcpy(ret,r);
	}

	return ret;

}

// -------------------------------------------------------------

void Worker::SetFileName(char *fileName) {
	strncpy(fullFileName,fileName,512);
}

void Worker::SaveGeometry(char *fileName,GLProgress *prg,BOOL askConfirm,BOOL saveSelected,BOOL autoSave,BOOL crashSave) {

	if (needsReload&&(!crashSave && !saveSelected)) RealReload();
	char tmp[1024];
	char fileNameWithGeo[1024]; //file name with .geo extension (instead of .geo7z)
	char fileNameWithGeo7z[1024]; //file name with .geo extension (instead of .geo7z)
	char fileNameWithoutExtension[1024]; //file name without extension
	//char *ext = fileName+strlen(fileName)-4;
	char *ext,*dir;
	MolFlow *mApp = (MolFlow *)theApp;
	
	dir = strrchr(fileName,'\\');
	ext = strrchr(fileName,'.');

	if(!(ext) || !(*ext=='.') || ((dir)&&(dir>ext)) ) { 
		sprintf(fileName, compressSavedFiles ? "%s.geo7z" : "%s.geo", fileName); //set to default GEO/GEO7Z format
		ext = strrchr(fileName,'.');
	}

	ext++;

	// Read a file
	BOOL ok = TRUE;
	FileWriter *f = NULL;
	BOOL isTXT = _stricmp(ext,"txt")==0;
	BOOL isSTR = _stricmp(ext,"str")==0;
	BOOL isGEO = _stricmp(ext,"geo")==0;
	BOOL isGEO7Z = _stricmp(ext,"geo7z")==0;

	if(isTXT || isGEO || isGEO7Z || isSTR) {
		if (WAIT_TIMEOUT==WaitForSingleObject(compressProcessHandle,0)) {
			GLMessageBox::Display("Compressing a previous save file is in progress. Wait until that finishes"
				"or close process \"compress.exe\"\nIf this was an autosave attempt,"
				"you have to lower the autosave frequency.","Can't save right now.",GLDLG_OK,GLDLG_ICONERROR);
			return;
		}
		if (isGEO) {
			memcpy(fileNameWithoutExtension,fileName,sizeof(char)*(strlen(fileName)-4));
			fileNameWithoutExtension[strlen(fileName)-4]='\0';
		} else if (isGEO7Z) { //geo7z
			memcpy(fileNameWithoutExtension,fileName,sizeof(char)*(strlen(fileName)-6));
			fileNameWithoutExtension[strlen(fileName)-6]='\0';
		}		
		
			if (isGEO) {
				sprintf(fileNameWithGeo7z,"%s7z",fileName);
				memcpy(fileNameWithGeo,fileName,(strlen(fileName)+1)*sizeof(char));

				if(!autoSave && FileUtils::Exist(fileNameWithGeo7z) && isGEO7Z) {
					sprintf(tmp,"A .geo7z file of the same name exists. Overwrite that file ?\n%s",fileNameWithGeo7z);
					ok = ( GLMessageBox::Display(tmp,"Question",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONWARNING)==GLDLG_OK );
				}
			} 
			if (isGEO7Z) {
				memcpy(fileNameWithGeo,fileName,sizeof(char)*(strlen(fileName)-2));
				fileNameWithGeo[strlen(fileName)-2]='\0';
				memcpy(fileNameWithGeo7z,fileName,(1+strlen(fileName))*sizeof(char));
				sprintf(tmp,"A .geo file of the same name exists. Overwrite that file ?\n%s",fileNameWithGeo);
				if(!autoSave && FileUtils::Exist(fileNameWithGeo) ) {
					ok = ( GLMessageBox::Display(tmp,"Question",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONWARNING)==GLDLG_OK );
				}
			} 
			if(!autoSave && ok && FileUtils::Exist(fileName) ) {
				sprintf(tmp,"Overwrite existing file ?\n%s",fileName);
				if (askConfirm) ok = ( GLMessageBox::Display(tmp,"Question",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONWARNING)==GLDLG_OK );
			}
		
		if( ok ) {
			if( isSTR ) {
				geom->SaveSTR(dpHit,saveSelected);
			} else {
				try {
				if (isGEO7Z) {
					memcpy(fileNameWithGeo,fileName,sizeof(char)*(strlen(fileName)-2));
					fileNameWithGeo[strlen(fileName)-2]='\0';
					f = new FileWriter(fileNameWithGeo);
				} else
					f = new FileWriter(fileName);
				} catch(Error &e) {
					SAFE_DELETE(f);
					GLMessageBox::Display((char*)e.GetMsg(),"Error in Worker::SaveGeometry()",GLDLG_OK,GLDLG_ICONERROR);
					return;
				}
				geom->tNbDesorptionMax = maxDesorption;
				if( isTXT ) geom->SaveTXT(f,dpHit,saveSelected);
				else if( isGEO || isGEO7Z ) {
					// Retrieve leak cache
					int nbLeakSave,nbHHitSave;
					LEAK pLeak[NBHLEAK];
					if (!crashSave&& !saveSelected) GetLeak(pLeak, &nbLeakSave);
					// Retrieve hit cache (lines and dots)
					HIT pHits[NBHHIT];
					if (!crashSave&& !saveSelected) GetHHit(pHits, &nbHHitSave);
					geom->SaveGEO(f,prg,dpHit,this->userMoments,this,saveSelected,pLeak,&nbLeakSave,pHits,&nbHHitSave,crashSave);
				}
			}
			if (!autoSave && !saveSelected) strcpy(fullFileName, fileName);
			if (!autoSave && !saveSelected) {
				remove("Molflow_AutoSave.geo");
				remove("Molflow_AutoSave.geo7z");
			}
		}
	} else {
		SAFE_DELETE(f);
		throw Error("SaveGeometry(): Invalid file extension [only geo,txt,str]");
	}

	SAFE_DELETE(f);
	

		if (isGEO7Z) {
			if (FileUtils::Exist("compress.exe")) { //compress GEO file to GEO7Z using 7-zip launcher "compress.exe"
				sprintf(tmp,"compress.exe \"%s\"",fileNameWithGeo);
				int procId = StartProc_background(tmp);
				compressProcessHandle=OpenProcess(PROCESS_ALL_ACCESS, TRUE, procId);
				fileName=fileNameWithGeo7z;
			} else {
				GLMessageBox::Display("compress.exe (part of Molfow) not found.\n Will save as uncompressed GEO file.","Compressor not found",GLDLG_OK,GLDLG_ICONERROR);
				fileName=fileNameWithGeo;
			}
		} else fileName=fileNameWithGeo;
		if (!autoSave && !saveSelected) {
			SetFileName(fileName);
			mApp->UpdateTitle();
		}
	
	//Debug memory check
	 //_ASSERTE (!_CrtDumpMemoryLeaks());;
	 _ASSERTE(_CrtCheckMemory());
}

void Worker::ExportTextures(char *fileName,int mode,BOOL askConfirm,BOOL saveSelected) {

	char tmp[512];

	// Read a file
	FILE *f = NULL;



	BOOL ok = TRUE;
	if( askConfirm ) {
		if( FileUtils::Exist(fileName) ) {
			sprintf(tmp,"Overwrite existing file ?\n%s",fileName);
			ok = ( GLMessageBox::Display(tmp,"Question",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONWARNING)==GLDLG_OK );
		}
	}
	if( ok ) {
		f=fopen(fileName,"w");
		geom->ExportTextures(f,mode,dpHit,saveSelected);
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

void Worker::LoadGeometry(char *fileName) {
	needsReload=TRUE;
	MolFlow *mApp = (MolFlow *)theApp;
	char *ext,*filebegin;
	BOOL isGEO7Z,isSYN7Z;
	char CWD [MAX_PATH];
	_getcwd( CWD, MAX_PATH );

	ext = strrchr(fileName,'.');

	if( ext==NULL )
		throw Error("LoadGeometry(): Invalid file extension [Only txt,stl,str,geo,geo7z or ase]");
	ext++;

	// Read a file
	FileReader *f = NULL;
	GLProgress *progressDlg = new GLProgress("Reading file...","Please wait");
	progressDlg->SetVisible(TRUE);
	progressDlg->SetProgress(0.0);

	//Clear hits and leaks cache
	memset(hhitCache,0,sizeof(HIT)*NBHHIT);
	memset(leakCache,0,sizeof(LEAK)*NBHLEAK);
	ResetMoments();

	isGEO7Z=(_stricmp(ext,"geo7z")==0);
	isSYN7Z=(_stricmp(ext,"syn7z")==0);
	if(_stricmp(ext,"txt")==0) {

		try {
			ResetWorkerStats();
			f = new FileReader(fileName);
			geom->LoadTXT(f,progressDlg);
			nbHit = geom->tNbHit;
			nbDesorption = geom->tNbDesorption;
			nbAbsorption = geom->tNbAbsorption;
			maxDesorption = geom->tNbDesorptionMax;
			nbLeakTotal = geom->tNbLeak;
			//RealReload();
			strcpy(fullFileName,fileName);
		} catch(Error &e) {
			geom->Clear();
			SAFE_DELETE(f);
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			throw e;
		}

	} else if(_stricmp(ext,"stl")==0) {

		try {
			int ret = GLUnitDialog::Display("","Choose STL file units:",GLDLG_MM | GLDLG_CM| GLDLG_M| GLDLG_INCH| GLDLG_FOOT | GLDLG_CANCEL_U,GLDLG_ICONNONE);
			double scaleFactor=1.0;
			switch(ret) {
			case GLDLG_MM:
				scaleFactor=0.1;
				break;
			case GLDLG_CM:
				scaleFactor=1.0;
				break;
			case GLDLG_M:
				scaleFactor=100;
				break;
			case GLDLG_INCH:
				scaleFactor=2.54;
				break;
			case GLDLG_FOOT:
				scaleFactor=30.48;
				break;
			}
			if (ret!=GLDLG_CANCEL_U) {
				progressDlg->SetMessage("Resetting worker...");
				ResetWorkerStats();
				f = new FileReader(fileName);
				progressDlg->SetVisible(TRUE);
				progressDlg->SetMessage("Reading geometry...");
				geom->LoadSTL(f,progressDlg,scaleFactor);
				//progressDlg->SetMessage("Reloading worker...");
				//RealReload();
				strcpy(fullFileName,fileName);
				mApp->DisplayCollapseDialog();
			}
		} catch(Error &e) {
			geom->Clear();
			SAFE_DELETE(f);
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			throw e;

		}

	} else if(_stricmp(ext,"str")==0) {

		try {
			ResetWorkerStats();
			f = new FileReader(fileName);
			progressDlg->SetVisible(TRUE);
			geom->LoadSTR(f,progressDlg);
			//RealReload();
			strcpy(fullFileName,fileName);
		} catch(Error &e) {
			geom->Clear();
			SAFE_DELETE(f);
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			throw e;
		}

	} else if(_stricmp(ext,"syn")==0 || isSYN7Z) {
		char tmp2[1024];
		int version;
		progressDlg->SetVisible(TRUE);
		try {
			if (isSYN7Z) {
				//decompress file
				progressDlg->SetMessage("Decompressing file...");
				char tmp[1024];
				char fileOnly[512];
				sprintf(tmp,"cmd /C \"pushd \"%s\"&&7za.exe x -t7z -aoa \"%s\" -otmp&&popd\"",CWD,fileName);
				//execCMD(tmp);
				system(tmp);

				filebegin= strrchr(fileName,'\\');
				if (filebegin) filebegin++;
				else filebegin=fileName;
				memcpy(fileOnly,filebegin,sizeof(char)*(strlen(filebegin)-2));
				fileOnly[strlen(filebegin)-2]='\0';
				sprintf(tmp2,"%s\\tmp\\Geometry.syn",CWD);
				f = new FileReader(tmp2);
			}
			progressDlg->SetMessage("Resetting worker...");
			ResetWorkerStats();
			if (!isSYN7Z) f = new FileReader(fileName);

			//leaks
			LEAK pLeak[NBHLEAK];
			//hits
			HIT pHits[NBHHIT];

			geom->LoadSYN(f,progressDlg,pLeak,&nbLastLeaks,pHits,&nbHHit,&version);
			maxDesorption = 0;
			
			progressDlg->SetMessage("Reloading worker with new geometry...");
			RealReload(); //for the loading of textures
			SHGHITS *gHits = (SHGHITS *)dpHit->buff;
			SAFE_DELETE(f);
			strcpy(fullFileName,fileName);
		} catch(Error &e) {
			geom->Clear();
			SAFE_DELETE(f);
			//if (isSYN7Z) remove(tmp2);
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			throw e;
		}

	} else if(_stricmp(ext,"geo")==0 || isGEO7Z) {
		char tmp2[1024];
		int version;
		progressDlg->SetVisible(TRUE);
		try {
			if (isGEO7Z) {
				//decompress file
				progressDlg->SetMessage("Decompressing file...");
				char tmp[1024];
				char fileOnly[512];
				sprintf(tmp,"cmd /C \"pushd \"%s\"&&7za.exe x -t7z -aoa \"%s\" -otmp&&popd\"",CWD,fileName);
				system(tmp);
				
				filebegin= strrchr(fileName,'\\');
				if (filebegin) filebegin++;
				else filebegin=fileName;
				memcpy(fileOnly,filebegin,sizeof(char)*(strlen(filebegin)-2));
				fileOnly[strlen(filebegin)-2]='\0';
				sprintf(tmp2,"%s\\tmp\\Geometry.geo",CWD);
				if (!FileUtils::Exist(tmp2)) //fall back to old geo7z format
					sprintf(tmp2,"%s\\tmp\\%s",CWD,fileOnly);
				f = new FileReader(tmp2);
			}
			progressDlg->SetMessage("Resetting worker...");
			ResetWorkerStats();
			if (!isGEO7Z) f = new FileReader(fileName);
			
			//leaks
			LEAK pLeak[NBHLEAK];
			//hits
			HIT pHits[NBHHIT];
			
			geom->LoadGEO(f,progressDlg,pLeak,&nbLastLeaks,pHits,&nbHHit,&version,this);
			//copy temp values from geom to worker:
			nbLeakTotal = geom->tNbLeak;
			nbHit = geom->tNbHit;
			nbDesorption = geom->tNbDesorption;
			maxDesorption = geom->tNbDesorptionMax;
			nbAbsorption = geom->tNbAbsorption;
			distTraveledTotal = geom->distTraveledTotal;
			
			progressDlg->SetMessage("Reloading worker with new geometry...");
			RealReload(); //for the loading of textures
			if (version>=8) geom->LoadProfile(f,dpHit,version);
			SHGHITS *gHits = (SHGHITS *)dpHit->buff;
			SetLeak(pLeak,&nbLastLeaks,gHits);
			SetHHit(pHits,&nbHHit,gHits);
			SAFE_DELETE(f);
			progressDlg->SetMessage("Loading textures...");
			LoadTextures((isGEO7Z)?tmp2:fileName,version);
			strcpy(fullFileName,fileName);
			//if (isGEO7Z) remove(tmp2);
		} catch(Error &e) {
			geom->Clear();
			SAFE_DELETE(f);
			//if (isGEO7Z) remove(tmp2);
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			throw e;
		}

	} else if(_stricmp(ext,"ase")==0) {

		try {
			ResetWorkerStats();
			f = new FileReader(fileName);
			progressDlg->SetVisible(TRUE);
			geom->LoadASE(f,progressDlg);
			//RealReload();
			strcpy(fullFileName,fileName);
		} catch(Error &e) {
			geom->Clear();
			SAFE_DELETE(f);
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			throw e;
		}

	} else {

		SAFE_DELETE(f);
		progressDlg->SetVisible(FALSE);
		SAFE_DELETE(progressDlg);
		throw Error("LoadGeometry(): Invalid file extension [Only txt,geo,geo7z,ase,stl or str]");

	}
	geom->CalcTotalOutGassing();
	if (mApp->momentsEditor) mApp->momentsEditor->Refresh();
	if (mApp->timeSettings) mApp->timeSettings->RefreshMoments();
	if (mApp->timewisePlotter) mApp->timewisePlotter->refreshViews();
	progressDlg->SetVisible(FALSE);
	SAFE_DELETE(progressDlg);
	SAFE_DELETE(f);

	//Debug memory check
	 //_ASSERTE (!_CrtDumpMemoryLeaks());;
	 _ASSERTE(_CrtCheckMemory());
}

// -------------------------------------------------------------
// -------------------------------------------------------------

void Worker::InsertGeometry(BOOL newStr,char *fileName) {
		MolFlow *mApp = (MolFlow *)theApp;
		if (needsReload) RealReload();
	char *ext,*filebegin;
	char tmp2[1024];
	ext = strrchr(fileName,'.');
	char CWD [MAX_PATH];
	_getcwd( CWD, MAX_PATH );

	if( ext==NULL )
		throw Error("InsertGeometry(): Invalid file extension [Only txt,stl,geo]");
	ext++;

	// Read a file
	FileReader *f = NULL;

	GLProgress *progressDlg = new GLProgress("Loading file...","Please wait");
	progressDlg->SetProgress(0.0);
	progressDlg->SetVisible(TRUE);
	BOOL isGEO7Z=(_stricmp(ext,"geo7z")==0);
	BOOL isSYN7Z=(_stricmp(ext,"syn7z")==0);

	if(_stricmp(ext,"txt")==0) {

		try {
			f=new FileReader(fileName);
			changedSinceSave = TRUE;
			geom->InsertTXT(f,progressDlg,newStr);
			//geom->sh.nbSuper++;
			//nbHit = geom->tNbHit;
			//nbDesorption = geom->tNbDesorption;
			//maxDesorption = geom->tNbDesorptionMax;
			//nbLeak = geom->tNbLeak;
			nbHit = 0;
			nbDesorption = 0;
			maxDesorption = 0;
			nbLeakTotal = 0;
			Reload();
			//strcpy(fullFileName,fileName);
		} catch(Error &e) {
			//geom->Clear();
			SAFE_DELETE(f);
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			throw e;
		}

	} else if(_stricmp(ext,"stl")==0) {
		try {
			int ret = GLUnitDialog::Display("","Choose STL file units:",GLDLG_MM | GLDLG_CM| GLDLG_M| GLDLG_INCH| GLDLG_FOOT | GLDLG_CANCEL_U,GLDLG_ICONNONE);
			double scaleFactor=1.0;
			switch(ret) {
			case GLDLG_MM:
				scaleFactor=0.1;
				break;
			case GLDLG_CM:
				scaleFactor=1.0;
				break;
			case GLDLG_M:
				scaleFactor=100;
				break;
			case GLDLG_INCH:
				scaleFactor=2.54;
				break;
			case GLDLG_FOOT:
				scaleFactor=30.48;
				break;
			}
			if (ret!=GLDLG_CANCEL_U) {
				changedSinceSave = TRUE;
				ResetWorkerStats();
				f = new FileReader(fileName);
				geom->InsertSTL(f,progressDlg,scaleFactor,newStr);
				nbHit = 0;
				nbLeakTotal = 0;
				nbDesorption = 0;
				maxDesorption = 0;
				nbLeakTotal = 0;
				Reload();}
			//strcpy(fullFileName,fileName);
		} catch(Error &e) {
			//geom->Clear();
			SAFE_DELETE(f);
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			throw e;
		}
	} else if(_stricmp(ext,"str")==0) {
		throw Error("STR file inserting is not supported.");
	} else if(_stricmp(ext,"syn")==0 || isSYN7Z) {
		//char tmp2[1024];
		
		progressDlg->SetVisible(TRUE);
		try {
			if (isSYN7Z) {
				//decompress file
				progressDlg->SetMessage("Decompressing file...");
				char tmp[1024];
				char fileOnly[512];
				sprintf(tmp,"cmd /C \"pushd \"%s\"&&7za.exe x -t7z -aoa \"%s\" -otmp&&popd\"",CWD,fileName);
				//execCMD(tmp);
				system(tmp);

				filebegin= strrchr(fileName,'\\');
				if (filebegin) filebegin++;
				else filebegin=fileName;
				memcpy(fileOnly,filebegin,sizeof(char)*(strlen(filebegin)-2));
				fileOnly[strlen(filebegin)-2]='\0';
				sprintf(tmp2,"%s\\tmp\\Geometry.syn",CWD);
				f = new FileReader(tmp2);
			}
			//progressDlg->SetMessage("Resetting worker...");
			//ResetWorkerStats();
			if (!isSYN7Z) f = new FileReader(fileName);


			geom->InsertSYN(f,progressDlg,newStr);
			nbLeakTotal = 0;
			nbHit = 0;
			nbDesorption = 0;
			maxDesorption = geom->tNbDesorptionMax;

			progressDlg->SetMessage("Reloading worker with new geometry...");
			RealReload();
			SAFE_DELETE(f);
			strcpy(fullFileName,fileName);
		} catch(Error &e) {
			//geom->Clear();
			SAFE_DELETE(f);
			//if (isSYN7Z) remove(tmp2);
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			throw e;
		}
	} else if(isGEO7Z || _stricmp(ext,"geo")==0) {
		try {
			if (isGEO7Z) {
				//decompress file
				progressDlg->SetMessage("Decompressing file...");
				char tmp[1024];
				char fileOnly[512];
				sprintf(tmp,"cmd /C \"pushd \"%s\"&&7za.exe x -t7z -aoa \"%s\" -otmp&&popd\"",CWD,fileName);
				//execCMD(tmp);
				system(tmp);
				
				filebegin= strrchr(fileName,'\\');
				if (filebegin) filebegin++;
				else filebegin=fileName;
				memcpy(fileOnly,filebegin,sizeof(char)*(strlen(filebegin)-2));
				fileOnly[strlen(filebegin)-2]='\0';
				sprintf(tmp2,"%s\\tmp\\Geometry.geo",CWD);
				if (!FileUtils::Exist(tmp2)) //fall back to old geo7z format
					sprintf(tmp2,"%s\\tmp\\%s",CWD,fileOnly);
				f = new FileReader(tmp2);
			}
			progressDlg->SetMessage("Resetting worker...");
			ResetWorkerStats();
			if (!isGEO7Z) f = new FileReader(fileName);
			changedSinceSave = TRUE;
			geom->InsertGEO(f,progressDlg,newStr);
			nbHit = 0;
			nbDesorption = 0;
			maxDesorption = 0;
			nbLeakTotal = 0;
			Reload();
			//strcpy(fullFileName,fileName);
		} catch(Error &e) {
			//geom->Clear();
			SAFE_DELETE(f);
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);
			throw e;
		}
	} else if(_stricmp(ext,"ase")==0) {
		throw Error("ASE file inserting is not supported.");
	} else {

		SAFE_DELETE(f);
		progressDlg->SetVisible(FALSE);
		SAFE_DELETE(progressDlg);
		throw Error("InsertGeometry(): Invalid file extension [Only txt,stl,geo,ase or str]");

	}

	progressDlg->SetVisible(FALSE);
	SAFE_DELETE(progressDlg);
	SAFE_DELETE(f);
	mApp->UpdateFacetlistSelected();	
	mApp->UpdateViewers();

	//Debug memory check
	 //_ASSERTE (!_CrtDumpMemoryLeaks());;
	 _ASSERTE(_CrtCheckMemory());
}

// -------------------------------------------------------------

void Worker::GetLeak(LEAK *buffer,int *nb) {

	*nb=0;
	if( dpHit ) {
		memcpy(buffer,leakCache,sizeof(LEAK)*NBHLEAK);
		*nb = (int)MIN(nbLastLeaks,NBHLEAK);
	}

}

// -------------------------------------------------------------

void Worker::SetLeak(LEAK *buffer,int *nb,SHGHITS *gHits) { //When loading from file

	if( nb>0 ) {
		memcpy(leakCache,buffer,sizeof(LEAK)*MIN(NBHLEAK,*nb));
		memcpy(gHits->pLeak,buffer,sizeof(LEAK)*MIN(NBHLEAK,*nb));
		gHits->nbLastLeaks = *nb;
	}

}

// -------------------------------------------------------------
void Worker::LoadTextures(char *fileName,int version) {
	// Read a file
	char *ext;

	ext = strrchr(fileName,'.');
	FileReader *f = NULL;

	GLProgress *progressDlg = new GLProgress("Loading textures","Please wait");
	progressDlg->SetProgress(0.0);


	if(_stricmp(ext,".geo")==0) {

		try {

			f = new FileReader(fileName);
			progressDlg->SetVisible(TRUE);
			geom->LoadTextures(f,progressDlg,dpHit,version);
			if( AccessDataport(dpHit) ) {
				BYTE *buffer = (BYTE *)dpHit->buff;
				try{geom->BuildTexture(buffer);} catch (Error &e) {
					ReleaseDataport(dpHit);
					throw e;
				}
				ReleaseDataport(dpHit);
			}
			progressDlg->SetVisible(FALSE);
			SAFE_DELETE(progressDlg);

		} catch(Error &e) {
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
	}

	SAFE_DELETE(f);

	//Debug memory check
	 //_ASSERTE (!_CrtDumpMemoryLeaks());;
	 _ASSERTE(_CrtCheckMemory());
}



// -------------------------------------------------------------

void Worker::GetHHit(HIT *buffer,int *nb) {

	*nb=0;
	if( dpHit ) {
		memcpy(buffer,hhitCache,sizeof(HIT)*NBHHIT);
		*nb = nbHHit;
	}

}

// -------------------------------------------------------------

void Worker::SetHHit(HIT *buffer,int *nb,SHGHITS *gHits) {
	
	/*if( nb>0 ) {
		memcpy(hhitCache,buffer,sizeof(HIT)*MIN(*nb,NBHHIT));
		memcpy(gHits->pHits,buffer,sizeof(LEAK)*MIN(NBHLEAK,*nb));
		gHits->nbHHit = *nb;
	}*/
	memcpy(hhitCache,buffer,sizeof(HIT)*MIN(*nb,NBHHIT));
	memcpy(gHits->pHits,buffer,sizeof(HIT)*MIN(*nb,NBHHIT));
	gHits->nbHHit = *nb;

}

// -------------------------------------------------------------

void Worker::InnerStop(float appTime) {

	stopTime =appTime;
	simuTime+=appTime-startTime;
	running  = FALSE;
	calcAC = FALSE;

}

// -------------------------------------------------------------

void Worker::OneStep() {

	if( nbProcess==0 )
		throw Error("No sub process found. (Simulation not available)");

	if( !running )  {
		if( !ExecuteAndWait(COMMAND_STEPAC,PROCESS_RUN,AC_MODE) )
			ThrowSubProcError();
	}

}

// -------------------------------------------------------------

void Worker::StepAC(float appTime) {

	try {
		OneStep();
		Update(appTime);
	} catch(Error &e) {
		GLMessageBox::Display((char *)e.GetMsg(),"Error (Stop)",GLDLG_OK,GLDLG_ICONERROR);
	}

}

void Worker::Stop_Public() {
	// Stop
	InnerStop(m_fTime);
	try {
		Stop();
		Update(m_fTime);
	} catch(Error &e) {
		GLMessageBox::Display((char *)e.GetMsg(),"Error (Stop)",GLDLG_OK,GLDLG_ICONERROR);
	}
}

// -------------------------------------------------------------

void Worker::StartStop(float appTime,int mode) {

	if( running )  {

		// Stop
		InnerStop(appTime);
		try {
			Stop();
			Update(appTime);
		} catch(Error &e) {
			GLMessageBox::Display((char *)e.GetMsg(),"Error (Stop)",GLDLG_OK,GLDLG_ICONERROR);
		}

	} else {

		// Start
		try {
			if (needsReload) RealReload();
			startTime = appTime;
			running = TRUE;
			calcAC = FALSE;
			this->mode = mode;
			Start();
			
			// Particular case when simulation ends before getting RUN state
			if(allDone) {
				Update(appTime);
				GLMessageBox::Display("Max desorption reached","Information (Start)",GLDLG_OK,GLDLG_ICONINFO);
			}
		} catch(Error &e) {
			running = FALSE;
			GLMessageBox::Display((char *)e.GetMsg(),"Error (Start)",GLDLG_OK,GLDLG_ICONERROR);
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
		if( dpControl ) {
			if( AccessDataport(dpControl) ) {
				int i = 0;
				SHMASTER *master =(SHMASTER *)dpControl->buff;
				for(int proc=0;proc<nbProcess && done;proc++) {
					done = done && (master->states[proc]==PROCESS_DONE);
					error = error && (master->states[proc]==PROCESS_ERROR);
					if(master->states[proc]==PROCESS_RUNAC) calcACprg = master->cmdParam[proc];
				}
				ReleaseDataport(dpControl);
			}
		}

		// End of simulation reached (Stop GUI)
		if( (error || done) && running && appTime!=0.0f ) {
			InnerStop(appTime);
			if (error) ThrowSubProcError();
		}

		// Retrieve hit count recording from the shared memory
		if( dpHit ) {

			if( AccessDataport(dpHit) ) {
				BYTE *buffer = (BYTE *)dpHit->buff;

				changedSinceSave = TRUE;
				// Globals
				SHGHITS *gHits = (SHGHITS *)buffer;

				// Global hits and leaks
				nbHit = gHits->total.hit.nbHit;
				nbAbsorption = gHits->total.hit.nbAbsorbed;
				distTraveledTotal = gHits->distTraveledTotal;
				nbDesorption = gHits->total.hit.nbDesorbed;
				nbLeakTotal = gHits->nbLeakTotal;
				nbHHit = gHits->nbHHit;
				nbLastLeaks = gHits->nbLastLeaks;
				memcpy(hhitCache,gHits->pHits,sizeof(HIT)*NBHHIT);
				memcpy(leakCache,gHits->pLeak,sizeof(LEAK)*NBHHIT);

				// Facets
				int nbFacet = geom->GetNbFacet();
				for(int i=0;i<nbFacet;i++) {    
					Facet *f = geom->GetFacet(i);
					memcpy(&(f->sh.counter),buffer+f->sh.hitOffset,sizeof(SHHITS));
					//f = f;
				}
				try {
					geom->BuildTexture(buffer);
				} catch (Error &e) {
					GLMessageBox::Display((char *)e.GetMsg(),"Error",GLDLG_OK,GLDLG_ICONERROR);
					ReleaseDataport(dpHit);
					return;
				}
				ReleaseDataport(dpHit);
			}

		}
	}
}

// -------------------------------------------------------------

void Worker::SendHits() {
	//if (!needsReload) {
	if(dpHit) {
		if( AccessDataport(dpHit) ) {

			// Store initial hit counts in the shared memory
			BYTE *pBuff = (BYTE *)dpHit->buff;
			memset(pBuff,0,geom->GetHitsSize(&moments));

			SHGHITS *gHits = (SHGHITS *)pBuff;
			gHits->total.hit.nbHit = nbHit;
			gHits->nbLeakTotal = nbLeakTotal;
			gHits->total.hit.nbDesorbed = nbDesorption;
			gHits->total.hit.nbAbsorbed = nbAbsorption;
			gHits->distTraveledTotal = distTraveledTotal;

			int nbFacet = geom->GetNbFacet();
			for(int i=0;i<nbFacet;i++) {
				Facet *f = geom->GetFacet(i);
				memcpy(pBuff+f->sh.hitOffset,&(f->sh.counter),sizeof(SHHITS));
			}
			ReleaseDataport(dpHit);

		} else {
			throw Error("Failed to initialize 'hits' dataport");
		}
	}
}

// -------------------------------------------------------------

void Worker::ComputeAC(float appTime) {
	if (needsReload) RealReload();
	if( running )
		throw Error("Already running");

	// Send correction map to sub process
	// (correction map contains area a surface elements)
	int maxElem = geom->GetMaxElemNumber();
	if( !maxElem )
		throw Error("Mesh with boundary correction must be enabled on all polygons");
	int dpSize = maxElem*sizeof(SHELEM);

	Dataport *loader = CreateDataport(loadDpName,dpSize);
	if( !loader )
		throw Error("Failed to create 'loader' dataport");
	AccessDataport(loader);
	geom->CopyElemBuffer((BYTE *)loader->buff);
	ReleaseDataport(loader);

	// Load Elem area and send AC matrix calculation order
	// Send command
	if( !ExecuteAndWait(COMMAND_LOADAC,PROCESS_RUNAC,dpSize) ) {
		CLOSEDP(loader);
		char errMsg[1024];
		sprintf(errMsg,"Failed to send geometry to sub process:\n%s",GetErrorDetails());
		GLMessageBox::Display(errMsg,"Warning (LoadAC)",GLDLG_OK,GLDLG_ICONWARNING);
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
	if (needsReload) RealReload();
	if( dpHit )
		if( AccessDataport(dpHit) ) 
			return (BYTE *)dpHit->buff;

	return NULL;

}

// -------------------------------------------------------------

void Worker::ThrowSubProcError(char *message) {

	char errMsg[1024];
	if( !message )
		sprintf(errMsg,"Bad response from sub process(es):\n%s",GetErrorDetails());
	else
		sprintf(errMsg,"%s\n%s",message,GetErrorDetails());
	throw Error(errMsg);

}

void Worker::Reload (){
	needsReload=true;
}

// -------------------------------------------------------------

void Worker::RealReload() { //Sharing geometry with workers

	if(nbProcess==0) return;

	// Clear geometry
	CLOSEDP(dpHit);
	if( !ExecuteAndWait(COMMAND_CLOSE,PROCESS_READY) )
		ThrowSubProcError();

	if(!geom->IsLoaded()) return;

	// Create the temporary geometry shared structure
	int loadSize = geom->GetGeometrySize(&moments);
	Dataport *loader = CreateDataport(loadDpName,loadSize);
	if( !loader )
		throw Error("Failed to create 'loader' dataport");
	AccessDataportTimed(loader,3000+nbProcess*(int)((double)loadSize/10000.0));
	geom->CopyGeometryBuffer((BYTE *)loader->buff,&moments);
	ReleaseDataport(loader);

	int hitSize = geom->GetHitsSize(&moments);
	dpHit = CreateDataport(hitsDpName,hitSize);
	if( !dpHit ) {
		CLOSEDP(loader);
		GLMessageBox::Display("Failed to create 'hits' dataport: not enough memory.","Warning (Load)",GLDLG_OK,GLDLG_ICONERROR);
		return;
		//throw Error("Failed to create 'hits' dataport. Probably out of memory.");
	}

	// Compute number of max desorption per process
	if(AccessDataportTimed(dpControl,3000+nbProcess*(int)((double)loadSize/10000.0))) {
		SHMASTER *m = (SHMASTER *)dpControl->buff;
		llong common = maxDesorption / (llong)nbProcess;
		int remain = (int)(maxDesorption % (llong)nbProcess);
		for(int i=0;i<nbProcess;i++) {
			m->cmdParam2[i] = common;
			if(i<remain) m->cmdParam2[i]++;
		}
		ReleaseDataport(dpControl);
	}

	// Load geometry
	if( !ExecuteAndWait(COMMAND_LOAD,PROCESS_READY,loadSize) ) {
		CLOSEDP(loader);
		char errMsg[1024];
		sprintf(errMsg,"Failed to send geometry to sub process:\n%s",GetErrorDetails());
		GLMessageBox::Display(errMsg,"Warning (Load)",GLDLG_OK,GLDLG_ICONWARNING);
		return;
	}

	//Send hit counts
	try {
		SendHits();
	} catch(Error &e) {
		// Geometry not loaded !
		CLOSEDP(dpHit);
		CLOSEDP(loader);
		throw e;
	}

	CLOSEDP(loader);
	needsReload=false;

	//Debug memory check
	 //_ASSERTE (!_CrtDumpMemoryLeaks());;
	 _ASSERTE(_CrtCheckMemory());
}

// -------------------------------------------------------------

void Worker::SetMaxDesorption(llong max) {

	try {
		Reset(0.0);
		maxDesorption = max;
		Reload();
	} catch(Error &e) {
		GLMessageBox::Display((char *)e.GetMsg(),"Error",GLDLG_OK,GLDLG_ICONERROR);
	}

}

// -------------------------------------------------------------

char *Worker::GetErrorDetails() {

	static char err[1024];
	strcpy(err,"");

	AccessDataport(dpControl);
	SHMASTER *shMaster = (SHMASTER *)dpControl->buff;
	for(int i=0;i<nbProcess;i++) {
		char tmp[256];
		if( pID[i]!=0 ) {
			int st = shMaster->states[i];
			if(st==PROCESS_ERROR) {
				sprintf(tmp,"[#%d] Process [PID %d] %s: %s\n",i,pID[i],prStates[st],shMaster->statusStr[i]);
			} else {
				sprintf(tmp,"[#%d] Process [PID %d] %s\n",i,pID[i],prStates[st]);
			}
		} else {
			sprintf(tmp,"[#%d] Process [PID ???] Not started\n",i);
		}
		strcat(err,tmp);
	}
	ReleaseDataport(dpControl);

	return err;
}

// -------------------------------------------------------------

void Worker::ClearHits() {
	if (needsReload) RealReload();
	if(dpHit) {
		AccessDataport(dpHit);
		memset(dpHit->buff,0,geom->GetHitsSize(&moments));
		ReleaseDataport(dpHit);
	}

}

// -------------------------------------------------------------

BOOL Worker::Wait(int waitState,int timeout) {

	BOOL ok = FALSE;
	BOOL error = FALSE;
	int t = 0;
	allDone = TRUE;

	// Wait for completion
	while(!ok && t<timeout) {

		ok = TRUE;
		AccessDataport(dpControl);
		SHMASTER *shMaster = (SHMASTER *)dpControl->buff;
		for(int i=0;i<nbProcess;i++) {
			ok = ok & (shMaster->states[i]==waitState || shMaster->states[i]==PROCESS_ERROR || shMaster->states[i]==PROCESS_DONE);
			if( shMaster->states[i]==PROCESS_ERROR ) error = TRUE;
			allDone = allDone & (shMaster->states[i]==PROCESS_DONE);
		}
		ReleaseDataport(dpControl);

		if(!ok) {
			Sleep(500);
			t+=500;
		}

	}

	if (t>=timeout) {
			int waitmore = GLMessageBox::Display("A subprocess is not responding. This might be normal for large geometries,\n"
				"otherwise it indicates an error. Click OK to wait more, or Cancel to stop trying.\n"
				"You can always restart subprocesses in the Tools / Global Settings... menu.\n"
				"Do you want to wait a bit more?","No reply from molflowSub.exe",GLDLG_OK|GLDLG_CANCEL,GLDLG_ICONWARNING)==GLDLG_OK;
			if (waitmore) {
				t=0;
				return Wait(waitState,timeout);
			}
	}

	return ok && !error;

}

// -------------------------------------------------------------

BOOL Worker::ExecuteAndWait(int command,int waitState,int param) {
	/*static BOOL alreadyreloading=FALSE;
	if (needsReload && !alreadyreloading) {
		alreadyreloading = TRUE;
		RealReload();
		alreadyreloading = FALSE;
	}*/
	if(!dpControl) return FALSE;

	// Send command
	AccessDataport(dpControl);
	SHMASTER *shMaster = (SHMASTER *)dpControl->buff;
	for(int i=0;i<nbProcess;i++) {
		shMaster->states[i]=command;
		shMaster->cmdParam[i]=param;
	}
	ReleaseDataport(dpControl);

	Sleep(100);
	return Wait(waitState,3000+nbProcess*500);
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
	distTraveledTotal = 0.0;

}

// -------------------------------------------------------------

void Worker::Reset(float appTime) {

	if( calcAC ) {
		GLMessageBox::Display("Reset not allowed while calculating AC","Error",GLDLG_OK,GLDLG_ICONERROR);
		return;
	}

	stopTime = 0.0f;
	startTime = 0.0f;
	simuTime = 0.0f;
	running = FALSE;
	if( nbProcess==0 ) 
		return;

	try {
		ResetWorkerStats();
		if( !ExecuteAndWait(COMMAND_RESET,PROCESS_READY) )
			ThrowSubProcError();
		ClearHits();
		Update(appTime);
	} catch(Error &e) {
		GLMessageBox::Display((char *)e.GetMsg(),"Error",GLDLG_OK,GLDLG_ICONERROR);
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
	int i=0;
	while(i<nbF && !found) {
		found = (geom->GetFacet(i)->sh.desorbType!=DES_NONE);
		if(!found) i++;
	}
	
	if( !found )
		throw Error("No desorption facet found");

	if(!(totalOutgassing>0.0))
		throw Error("Total outgassing is zero.");

	if( nbProcess==0 )
		throw Error("No sub process found. (Simulation not available)");

	if( !ExecuteAndWait(COMMAND_START,PROCESS_RUN,mode) )
		ThrowSubProcError();

	//Debug memory check
	 //_ASSERTE (!_CrtDumpMemoryLeaks());;
	 _ASSERTE(_CrtCheckMemory());

}

// -------------------------------------------------------------

void Worker::Stop() {

	if( nbProcess==0 )
		throw Error("No sub process found. (Simulation not available)");

	if( !ExecuteAndWait(COMMAND_PAUSE,PROCESS_READY) )
		ThrowSubProcError();

	//Debug memory check
	 //_ASSERTE (!_CrtDumpMemoryLeaks());;
	 _ASSERTE(_CrtCheckMemory());

}

// -------------------------------------------------------------

void Worker::KillAll() {

	if( dpControl && nbProcess>0 ) {
		if( !ExecuteAndWait(COMMAND_EXIT,PROCESS_KILLED) ) {
			AccessDataport(dpControl);
			SHMASTER *shMaster = (SHMASTER *)dpControl->buff;
			for(int i=0;i<nbProcess;i++) 
				if(shMaster->states[i]==PROCESS_KILLED) pID[i]=0;
			ReleaseDataport(dpControl);
			// Force kill
			for(int i=0;i<nbProcess;i++)
				if(pID[i]) KillProc(pID[i]);
		}
		CLOSEDP(dpHit);
	}
	nbProcess = 0;

}

// -------------------------------------------------------------

void Worker::Exit() {

	if( dpControl && nbProcess>0 ) {
		KillAll();
		CLOSEDP(dpControl);
	}

}

// -------------------------------------------------------------

void Worker::GetProcStatus(int *states,char **status) {

	if(nbProcess==0) return;

	AccessDataport(dpControl);
	SHMASTER *shMaster = (SHMASTER *)dpControl->buff;
	memcpy(states,shMaster->states,MAX_PROCESS*sizeof(int));
	memcpy(status,shMaster->statusStr,MAX_PROCESS*64);
	ReleaseDataport(dpControl);

}

// -------------------------------------------------------------

void Worker::SetProcNumber(int n) {

	char cmdLine[512];

	// Kill all sub process
	KillAll();

	// Create new control dataport
	if( !dpControl ) 
		dpControl = CreateDataport(ctrlDpName,sizeof(SHMASTER));
	if( !dpControl )
		throw Error("Failed to create 'control' dataport");
	AccessDataport(dpControl);
	memset(dpControl->buff,0,sizeof(SHMASTER));
	ReleaseDataport(dpControl);

	// Launch n subprocess
	for(int i=0;i<n;i++) {
		sprintf(cmdLine,"molflowSub.exe %d %d",pid,i);
		pID[i] = StartProc(cmdLine);
		Sleep(25); // Wait a bit
		if( pID[i]==0 ) {
			nbProcess = 0;
			throw Error(cmdLine);
		}
	}

	nbProcess = n;

	if( !Wait(PROCESS_READY,3000) ){
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
	int nb=(int)newMoments.size();
	for (int i=0;i<nb;i++)
		moments.push_back(newMoments[i]);
	return nb;
}

std::vector<double> Worker::ParseMoment(std::string userInput) {
	std::vector<double> parsedResult;
	double begin,interval,end;

	int nb=sscanf(userInput.c_str(),"%lf,%lf,%lf",&begin,&interval,&end);
	if (nb==1 && (begin>0.0)) {
		//One moment
		parsedResult.push_back(begin);
	//} else if (nb==3 && (begin>0.0) && (end>begin) && (interval<(end-begin)) && ((end-begin)/interval<300.0)) {
	} else if (nb==3 && (begin>0.0) && (end>begin) && (interval<(end-begin))) {
		//Range
		for (double time=begin;time<=end;time+=interval)
			parsedResult.push_back(time);
	}
	return parsedResult;
}

void Worker::ResetMoments() {
	moments=std::vector<double>();
	userMoments=std::vector<std::string>();
}

void Worker::ImportDesorption_DES(char *fileName) {
	//if (needsReload) RealReload();

	// Read a file
	FileReader *f = new FileReader(fileName);
	geom->ImportDesorption_DES(f);
	SAFE_DELETE(f);
	changedSinceSave = TRUE;
	Reload();
}

void Worker::ImportDesorption_SYN(char *fileName, const size_t &source, const double &time,
	const size_t &mode, const double &eta0, const double &alpha,
	const std::vector<std::pair<double, double>> &convDistr,
	GLProgress *prg) {
	MolFlow *mApp = (MolFlow *)theApp;
	char *ext, *filebegin;
	char tmp2[1024];
	ext = strrchr(fileName, '.');
	char CWD[MAX_PATH];
	_getcwd(CWD, MAX_PATH);
	if (ext == NULL || (!(_stricmp(ext, ".syn7z") == 0)) && (!(_stricmp(ext, ".syn") == 0)))
		throw Error("InsertGeometry(): Invalid file extension [Only syn, syn7z]");
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

			geom->ImportDesorption_SYN(f, source, time, mode, eta0, alpha, convDistr, prg);

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
	MolFlow *mApp = (MolFlow *)theApp;
	char *ext, *filebegin;
	char tmp2[1024];
	ext = strrchr(fileName, '.');
	char CWD[MAX_PATH];
	_getcwd(CWD, MAX_PATH);
	if (ext == NULL || (!(_stricmp(ext, ".syn7z") == 0)) && (!(_stricmp(ext, ".syn") == 0)))
		throw Error("InsertGeometry(): Invalid file extension [Only syn, syn7z]");
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