/*
  File:        molflowSub.c
  Description: Main function on the working sub process
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

#include <windows.h>

#include <tlhelp32.h>

//#include <iostream>

#include <stdio.h>
#include <math.h>
#include <time.h>

#include "Simulation.h"
#ifdef WIN
#include <Process.h> // For _getpid()
#endif

// -------------------------------------------------
// Global process variables
// -------------------------------------------------
#define WAITTIME    100  // Answer in STOP mode
//#define TIMEOUT     300  // Process kills itself after no heartbeat (seconds)

static Dataport *dpControl=NULL;
static Dataport *dpHit=NULL;
//static int       noHeartBeatSince;
static int       prIdx;
static int       prState;
static size_t       prParam;
static llong     prParam2;
static DWORD     hostProcessId;
//static float       heartBeat;
//static HANDLE    masterHandle;
static char      ctrlDpName[32];
static char      loadDpName[32];
static char      hitsDpName[32];

BOOL end = FALSE;
BOOL IsProcessRunning(DWORD pid);
//int FIND_PROC_BY_NAME(const char *szToFind);
// -------------------------------------------------

//void GetState(int sleepTime) {
void GetState() {
  prState = PROCESS_READY;
  prParam = 0;

  if( AccessDataport(dpControl) ) {
    SHCONTROL *master = (SHCONTROL *)dpControl->buff;
    prState = master->states[prIdx];
    prParam = master->cmdParam[prIdx];
    prParam2 = master->cmdParam2[prIdx];
    master->cmdParam[prIdx] = 0;
    master->cmdParam2[prIdx] = 0;
	/*
	//debug process states
	for (i=0;i<8;i++) {
		printf(" %d",master->states[i]);
		if (i==prIdx) printf("*");
	}
	printf("\n");
	*/
    ReleaseDataport(dpControl);
	/*if (FIND_PROC_BY_NAME("synrad.exe")!=1) { //if there is no running Windows process called "synrad.exe"
		printf("synrad.exe not running. Closing.");
		SetErrorSub("synrad.exe not running. Closing subprocess.");
		end = TRUE;
	}*/
	if (!IsProcessRunning(hostProcessId)) {
		printf("Host synrad.exe (process id %d) not running. Closing.",hostProcessId);
		SetErrorSub("Host synrad.exe not running. Closing subprocess.");
		end = TRUE;
	}
	/*
	if (master->heartBeat>heartBeat) {
		heartBeat = master->heartBeat;
		noHeartBeatSince=0;
		printf("Heartbeat: %g sec\n",heartBeat);
	} else {
		noHeartBeatSince+=sleepTime;
		if (noHeartBeatSince%1000==0 && noHeartBeatSince>=2000)printf("No heartbeat since: %d ms.\n",noHeartBeatSince);
		if (noHeartBeatSince >= TIMEOUT*1000) {
			printf("No heartbeat from main program. Closing.");
			SetErrorSub("No heartbeat from main program. Closing subprocess.");
			Sleep(5000);
			end = TRUE;
		}
	}*/
  } else {
	  printf("Subprocess couldn't connect to Molflow.\n");
	  SetErrorSub("No connection to main program. Closing subprocess.");
	  Sleep(5000);
	  end = TRUE;
  }
}

// -------------------------------------------------

int GetLocalState() {
  return prState;
}

// -------------------------------------------------

void SetState(int state,const char *status) {

	prState = state;
	printf("\n setstate %d \n",state);
	if( AccessDataport(dpControl) ) {
		SHCONTROL *master = (SHCONTROL *)dpControl->buff;
		master->states[prIdx] = state;
		strncpy(master->statusStr[prIdx],status,63);
		master->statusStr[prIdx][63]=0;
		if( state==PROCESS_RUNAC ) {
			master->cmdParam[prIdx] = sHandle->prgAC;
		}
		ReleaseDataport(dpControl);
	}

}

// -------------------------------------------------

void SetErrorSub(const char *message) {

  printf("Error: %s\n",message);
  SetState(PROCESS_ERROR,message);

}

// -------------------------------------------------

char *GetSimuStatus() {

  int mode;
  static char ret[128];
  llong count = sHandle->nbDesorbed;
  llong max   = sHandle->maxDesorption;

  mode = sHandle->sMode;
  if( GetLocalState()==PROCESS_RUNAC ) mode = AC_MODE;

  switch(mode) {

    case MC_MODE:
      if( max!=0 ) {
        double percent = (double)(count)*100.0 / (double)(max);
        sprintf(ret,"(%s) MC %I64d/%I64d (%.1f%%)",sHandle->name,count,max,percent);
      } else {
        sprintf(ret,"(%s) MC %I64d",sHandle->name,count);
      }
      break;

    case AC_MODE:
      if( sHandle->prgAC<100 ) {
          sprintf(ret,"(%s) AC (%dx%d) (%zd%%)",sHandle->name,
                      sHandle->nbAC,sHandle->nbAC,sHandle->prgAC);
      } else {
        if( max!=0 ) {
          double percent = (double)(count)*100.0 / (double)(max);
          sprintf(ret,"(%s) AC (%dx%d) %I64d/%I64d (%.1f%%)",sHandle->name,
                      sHandle->nbAC,sHandle->nbAC,count,max,percent);
        } else {
          sprintf(ret,"(%s) AC (%dx%d) %I64d",sHandle->name,sHandle->nbAC,
                      sHandle->nbAC,count);
        }
      }
      break;

  }

  return ret;

}

void SetReady() {

  if(sHandle->loadOK)
    SetState(PROCESS_READY,GetSimuStatus());
  else
    SetState(PROCESS_READY,"(No geometry)");

}

// -------------------------------------------------

void SetStatus(char *message) {

  if( AccessDataport(dpControl) ) {
    SHCONTROL *master = (SHCONTROL *)dpControl->buff;
    strcpy(master->statusStr[prIdx],message);
    ReleaseDataport(dpControl);
  }

}

// -------------------------------------------------

void LoadAC() {

  Dataport *loader;
  SHELEM *map;

  if( !sHandle->loadOK ) {
    SetErrorSub("No geometry loaded");
    return;
  }

  ClearACMatrix();

  // Load mesh
  loader = OpenDataport(loadDpName,prParam);
  if( !loader ) {
    char err[512];
    sprintf(err,"Failed to open 'loader' dataport %s (%zd Bytes)",loadDpName, prParam);
    SetErrorSub(err);
    return;
  }

  // Connect the dataport
  if( !AccessDataport(loader) ) {
    SetErrorSub("Failed to connect to DP");
    return;
  }
  printf("Connected to %s\n",loadDpName);

  map = (SHELEM *)malloc( prParam );
  memcpy(map,loader->buff,prParam);
  ReleaseDataport(loader);
  CLOSEDP(loader);

  SetState(PROCESS_RUNAC,GetSimuStatus());
  ComputeACMatrix(map);
  if( GetLocalState()!=PROCESS_ERROR ) {
    char s[128];
    if( sHandle->prgAC==100 ) {
      sprintf(s,"AC matrix calcultion ok (%.3f s)",sHandle->calcACTime);
    } else {
      sprintf(s,"AC matrix calcultion interrupted");
    }
    SetState(PROCESS_DONE,s);
  }
  free(map);

}

// -------------------------------------------------

void Load() {

  Dataport *loader;
  size_t hSize;

  // Load geometry
  loader = OpenDataport(loadDpName,prParam);
  if( !loader ) {
    char err[512];
    sprintf(err,"Failed to connect to 'loader' dataport %s (%zd Bytes)",loadDpName, prParam);
    SetErrorSub(err);
    return;
  }
  
  printf("Connected to %s\n",loadDpName);

  if( !LoadSimulation(loader) ) {
    CLOSEDP(loader);
    return;
  }
  CLOSEDP(loader);

  // Connect to hit dataport
  hSize = GetHitsSize();
  dpHit = OpenDataport(hitsDpName,hSize);
  if( !dpHit ) {
	  char err[512];
	  sprintf(err, "Failed to connect to 'hits' dataport (%zd Bytes)", hSize);
	  SetErrorSub(err);
	sHandle->loadOK = FALSE;
    return;
  }

  printf("Connected to %s (%zd bytes)\n",hitsDpName,hSize);

}

// -------------------------------------------------

int main(int argc,char* argv[])
{
  BOOL eos = FALSE;

  if(argc!=3) {
    printf("Usage: molflowSub peerId index\n");
    return 1;
  }

  hostProcessId=atoi(argv[1]);
  prIdx = atoi(argv[2]);

  sprintf(ctrlDpName,"MFLWCTRL%s",argv[1]);
  sprintf(loadDpName,"MFLWLOAD%s",argv[1]);
  sprintf(hitsDpName,"MFLWHITS%s",argv[1]);

  dpControl = OpenDataport(ctrlDpName,sizeof(SHCONTROL));
  if( !dpControl ) {
    printf("Usage: Cannot connect to MFLWCTRL%s\n",argv[1]);
    return 1;
  }

  printf("Connected to %s (%zd bytes), molflowSub.exe #%d\n",ctrlDpName,sizeof(SHCONTROL),prIdx);

  InitSimulation();

  // Sub process ready
  SetReady();

  // Main loop
  while( !end ) {

    //GetState((prState==PROCESS_READY)?100:1000);
	GetState();
    switch(prState) {

      case COMMAND_LOAD:
        printf("COMMAND: LOAD (%zd,%llu)\n",prParam,prParam2);
        Load();
        if( sHandle->loadOK ) {
          sHandle->maxDesorption = prParam2; // 0 for endless
          SetReady();
        }
        break;

      case COMMAND_LOADAC:
        printf("COMMAND: LOADAC (%zd)\n",prParam);
        LoadAC();
        break;

      case COMMAND_START:
        printf("COMMAND: START (%zd,%llu)\n",prParam,prParam2);
        if( sHandle->loadOK ) {
          if( StartSimulation(prParam) )
            SetState(PROCESS_RUN,GetSimuStatus());
          else {
            if( GetLocalState()!=PROCESS_ERROR )
              SetState(PROCESS_DONE,GetSimuStatus());
          }
        } else
          SetErrorSub("No geometry loaded");
        break;

      case COMMAND_PAUSE:
        printf("COMMAND: PAUSE (%zd,%llu)\n",prParam,prParam2);
        if( !sHandle->lastUpdateOK ) {
          // Last update not successful, retry with a longer tomeout
          if(dpHit) UpdateHits(dpHit,prIdx,8000);
        }
        SetReady();
        break;

      case COMMAND_RESET:
        printf("COMMAND: RESET (%zd,%llu)\n",prParam,prParam2);
        ResetSimulation();
        SetReady();
        break;

      case COMMAND_EXIT:
        printf("COMMAND: EXIT (%zd,%llu)\n",prParam,prParam2);
        end = TRUE;
        break;

      case COMMAND_CLOSE:
        printf("COMMAND: CLOSE (%zd,%llu)\n",prParam,prParam2);
        ClearSimulation();
        CLOSEDP(dpHit);
        SetReady();
        break;

      case COMMAND_STEPAC:
        // Debug command
        printf("COMMAND: STEPAC (%zd,%llu)\n",prParam,prParam2);
        if( sHandle->loadOK ) {
          if( StartSimulation(prParam) ) {
            SetState(PROCESS_RUN,GetSimuStatus());
            SimulationACStep(1);
            if(dpHit) UpdateHits(dpHit,prIdx,20);
            SetReady();
          } else {
            if( GetLocalState()!=PROCESS_ERROR )
              SetState(PROCESS_DONE,GetSimuStatus());
          }
        } else
          SetErrorSub("No geometry loaded");
        break;

      case PROCESS_RUN:
        SetStatus(GetSimuStatus()); //update hits only
        eos = SimulationRun();                // Run during 1 sec
        if(dpHit) UpdateHits(dpHit,prIdx,20); // Update hit with 20ms timeout
        if(eos) {
          if( GetLocalState()!=PROCESS_ERROR ) {
            // Max desorption reached
            SetState(PROCESS_DONE,GetSimuStatus());
            printf("COMMAND: PROCESS_DONE (Max reached)\n");
          }
        }
        break;

      default:
        Sleep(WAITTIME);
        break;
    }

  }

  // Release
  SetState(PROCESS_KILLED,"");
  CLOSEDP(dpControl);
  CLOSEDP(dpHit);
  return 0;

}

/*
int FIND_PROC_BY_NAME(const char *szToFind)

// Created: 12/29/2000  (RK)

// Last modified: 6/16/2003  (RK)

// Please report any problems or bugs to kochhar@physiology.wisc.edu

// The latest version of this routine can be found at:

//     http://www.neurophys.wisc.edu/ravi/software/killproc/

// Check whether the process "szToFind" is currently running in memory

// This works for Win/95/98/ME and also Win/NT/2000/XP

// The process name is case-insensitive, i.e. "notepad.exe" and "NOTEPAD.EXE"

// will both work (for szToFind)

// Return codes are as follows:

//   0   = Process was not found

//   1   = Process was found

//   605 = Unable to search for process

//   606 = Unable to identify system type

//   607 = Unsupported OS

//   632 = Process name is invalid

// Change history:

//  3/10/2002   - Fixed memory leak in some cases (hSnapShot and

//                and hSnapShotm were not being closed sometimes)

//  6/13/2003   - Removed iFound (was not being used, as pointed out

//                by John Emmas)

{

    BOOL bResult,bResultm;
    DWORD aiPID[1000],iCb=1000,iNumProc,iV2000=0;
    DWORD iCbneeded,i;
    char szName[MAX_PATH],szToFindUpper[MAX_PATH];
    HANDLE hProc,hSnapShot,hSnapShotm;
    OSVERSIONINFO osvi;
    HINSTANCE hInstLib;
    int iLen,iLenP,indx;
    HMODULE hMod;
    PROCESSENTRY32 procentry;      
    MODULEENTRY32 modentry;

    // PSAPI Function Pointers.
     BOOL (WINAPI *lpfEnumProcesses)( DWORD *, DWORD cb, DWORD * );
     BOOL (WINAPI *lpfEnumProcessModules)( HANDLE, HMODULE *,
        DWORD, LPDWORD );
     DWORD (WINAPI *lpfGetModuleBaseName)( HANDLE, HMODULE,
        LPTSTR, DWORD );

      // ToolHelp Function Pointers.
      HANDLE (WINAPI *lpfCreateToolhelp32Snapshot)(DWORD,DWORD) ;
      BOOL (WINAPI *lpfProcess32First)(HANDLE,LPPROCESSENTRY32) ;
      BOOL (WINAPI *lpfProcess32Next)(HANDLE,LPPROCESSENTRY32) ;
      BOOL (WINAPI *lpfModule32First)(HANDLE,LPMODULEENTRY32) ;
      BOOL (WINAPI *lpfModule32Next)(HANDLE,LPMODULEENTRY32) ;

    // Transfer Process name into "szToFindUpper" and
    // convert it to upper case
    iLenP=strlen(szToFind);
    if(iLenP<1 || iLenP>MAX_PATH) return 632;
    for(indx=0;indx<iLenP;indx++)
        szToFindUpper[indx]=toupper(szToFind[indx]);
    szToFindUpper[iLenP]=0;

    // First check what version of Windows we're in
    osvi.dwOSVersionInfoSize = sizeof(OSVERSIONINFO);
    bResult=GetVersionEx(&osvi);
    if(!bResult)     // Unable to identify system version
        return 606;

    // At Present we only support Win/NT/2000 or Win/9x/ME
    if((osvi.dwPlatformId != VER_PLATFORM_WIN32_NT) &&
        (osvi.dwPlatformId != VER_PLATFORM_WIN32_WINDOWS))
        return 607;

    if(osvi.dwPlatformId==VER_PLATFORM_WIN32_NT)
    {
        // Win/NT or 2000 or XP

         // Load library and get the procedures explicitly. We do
         // this so that we don't have to worry about modules using
         // this code failing to load under Windows 95, because
         // it can't resolve references to the PSAPI.DLL.
         hInstLib = LoadLibraryA("PSAPI.DLL");
         if(hInstLib == NULL)
            return 605;

         // Get procedure addresses.
         lpfEnumProcesses = (BOOL(WINAPI *)(DWORD *,DWORD,DWORD*))
            GetProcAddress( hInstLib, "EnumProcesses" ) ;
         lpfEnumProcessModules = (BOOL(WINAPI *)(HANDLE, HMODULE *,
            DWORD, LPDWORD)) GetProcAddress( hInstLib,
            "EnumProcessModules" ) ;
         lpfGetModuleBaseName =(DWORD (WINAPI *)(HANDLE, HMODULE,
            LPTSTR, DWORD )) GetProcAddress( hInstLib,
            "GetModuleBaseNameA" ) ;

         if( lpfEnumProcesses == NULL ||
            lpfEnumProcessModules == NULL ||
            lpfGetModuleBaseName == NULL)
            {
               FreeLibrary(hInstLib);
               return 605;
            }

        bResult=lpfEnumProcesses(aiPID,iCb,&iCbneeded);
        if(!bResult)
        {
            // Unable to get process list, EnumProcesses failed
            FreeLibrary(hInstLib);
            return 605;
        }

        // How many processes are there?
        iNumProc=iCbneeded/sizeof(DWORD);

        // Get and match the name of each process
        for(i=0;i<iNumProc;i++)
        {
            // Get the (module) name for this process

            strcpy(szName,"Unknown");
            // First, get a handle to the process
            hProc=OpenProcess(PROCESS_QUERY_INFORMATION|PROCESS_VM_READ,FALSE,
                aiPID[i]);
            // Now, get the process name
            if(hProc)
            {
               if(lpfEnumProcessModules(hProc,&hMod,sizeof(hMod),&iCbneeded) )
               {
                  iLen=lpfGetModuleBaseName(hProc,hMod,szName,MAX_PATH);
               }
            }
            CloseHandle(hProc);
            // Match regardless of lower or upper case
            if(strcmp(_strupr(szName),szToFindUpper)==0)
            {
                // Process found
                FreeLibrary(hInstLib);
                return 1;
            }
        }
    }

    if(osvi.dwPlatformId==VER_PLATFORM_WIN32_WINDOWS)
    {
        // Win/95 or 98 or ME

        hInstLib = LoadLibraryA("Kernel32.DLL");
        if( hInstLib == NULL )
            return FALSE ;

        // Get procedure addresses.
        // We are linking to these functions of Kernel32
        // explicitly, because otherwise a module using
        // this code would fail to load under Windows NT,
        // which does not have the Toolhelp32
        // functions in the Kernel 32.
        lpfCreateToolhelp32Snapshot=
            (HANDLE(WINAPI *)(DWORD,DWORD))
            GetProcAddress( hInstLib,
            "CreateToolhelp32Snapshot" ) ;
        lpfProcess32First=
            (BOOL(WINAPI *)(HANDLE,LPPROCESSENTRY32))
            GetProcAddress( hInstLib, "Process32First" ) ;
        lpfProcess32Next=
            (BOOL(WINAPI *)(HANDLE,LPPROCESSENTRY32))
            GetProcAddress( hInstLib, "Process32Next" ) ;
        lpfModule32First=
            (BOOL(WINAPI *)(HANDLE,LPMODULEENTRY32))
            GetProcAddress( hInstLib, "Module32First" ) ;
        lpfModule32Next=
            (BOOL(WINAPI *)(HANDLE,LPMODULEENTRY32))
            GetProcAddress( hInstLib, "Module32Next" ) ;
        if( lpfProcess32Next == NULL ||
            lpfProcess32First == NULL ||
            lpfModule32Next == NULL ||
            lpfModule32First == NULL ||
            lpfCreateToolhelp32Snapshot == NULL )
        {
            FreeLibrary(hInstLib);
            return 605;
        }

        // The Process32.. and Module32.. routines return names in all uppercase

        // Get a handle to a Toolhelp snapshot of all the systems processes.

        hSnapShot = lpfCreateToolhelp32Snapshot(
            TH32CS_SNAPPROCESS, 0 ) ;
        if( hSnapShot == INVALID_HANDLE_VALUE )
        {
            FreeLibrary(hInstLib);
            return 605;
        }

        // Get the first process' information.
        procentry.dwSize = sizeof(PROCESSENTRY32);
        bResult=lpfProcess32First(hSnapShot,&procentry);

        // While there are processes, keep looping and checking.
        while(bResult)
        {
            // Get a handle to a Toolhelp snapshot of this process.
            hSnapShotm = lpfCreateToolhelp32Snapshot(
                TH32CS_SNAPMODULE, procentry.th32ProcessID) ;
            if( hSnapShotm == INVALID_HANDLE_VALUE )
            {
                CloseHandle(hSnapShot);
                FreeLibrary(hInstLib);
                return 605;
            }
            // Get the module list for this process
            modentry.dwSize=sizeof(MODULEENTRY32);
            bResultm=lpfModule32First(hSnapShotm,&modentry);

            // While there are modules, keep looping and checking
            while(bResultm)
            {
                if(strcmp(modentry.szModule,szToFindUpper)==0)
                {
                    // Process found
                    CloseHandle(hSnapShotm);
                    CloseHandle(hSnapShot);
                    FreeLibrary(hInstLib);
                    return 1;
                }
                else
                {  // Look for next modules for this process
                    modentry.dwSize=sizeof(MODULEENTRY32);
                    bResultm=lpfModule32Next(hSnapShotm,&modentry);
                }
            }

            //Keep looking
            CloseHandle(hSnapShotm);
            procentry.dwSize = sizeof(PROCESSENTRY32);
            bResult = lpfProcess32Next(hSnapShot,&procentry);
        }
        CloseHandle(hSnapShot);
    }
    FreeLibrary(hInstLib);
    return 0;

}
*/

BOOL IsProcessRunning(DWORD pid)
{
	HANDLE process = OpenProcess(SYNCHRONIZE, FALSE, pid);
	DWORD ret = WaitForSingleObject(process, 0);
	CloseHandle(process);
	return ret == WAIT_TIMEOUT;
}