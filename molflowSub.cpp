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

#define NOMINMAX
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

// Global process variables

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

bool end = false;
bool IsProcessRunning(DWORD pid);

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

    ReleaseDataport(dpControl);
	
	if (!IsProcessRunning(hostProcessId)) {
		printf("Host synrad.exe (process id %d) not running. Closing.",hostProcessId);
		SetErrorSub("Host synrad.exe not running. Closing subprocess.");
		end = true;
	}
  } else {
	  printf("Subprocess couldn't connect to Molflow.\n");
	  SetErrorSub("No connection to main program. Closing subprocess.");
	  Sleep(5000);
	  end = true;
  }
}

int GetLocalState() {
  return prState;
}

void SetState(int state,const char *status,bool changeState, bool changeStatus) {

	prState = state;
	if (changeState) printf("\n setstate %d \n",state);
	if( AccessDataport(dpControl) ) {
		SHCONTROL *master = (SHCONTROL *)dpControl->buff;
		if (changeState) master->states[prIdx] = state;
		if (changeStatus) {
			strncpy(master->statusStr[prIdx], status, 127);
			master->statusStr[prIdx][127] = 0;
		}
		if( state==PROCESS_RUNAC ) {
			master->cmdParam[prIdx] = sHandle->prgAC;
		}
		ReleaseDataport(dpControl);
	}

}

void SetErrorSub(const char *message) {

  printf("Error: %s\n",message);
  SetState(PROCESS_ERROR,message);

}

char *GetSimuStatus() {

  size_t mode;
  static char ret[128];
  llong count = sHandle->totalDesorbed;
  llong max   = sHandle->desorptionLimit;

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
          sprintf(ret,"(%s) AC (%zdx%zd) (%zd%%)",sHandle->name,
                      sHandle->nbAC,sHandle->nbAC,sHandle->prgAC);
      } else {
        if( max!=0 ) {
          double percent = (double)(count)*100.0 / (double)(max);
          sprintf(ret,"(%s) AC (%zdx%zd) %I64d/%I64d (%.1f%%)",sHandle->name,
                      sHandle->nbAC,sHandle->nbAC,count,max,percent);
        } else {
          sprintf(ret,"(%s) AC (%zdx%zd) %I64d",sHandle->name,sHandle->nbAC,
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

void SetStatus(char *status) {

  if( AccessDataport(dpControl) ) {
    SHCONTROL *master = (SHCONTROL *)dpControl->buff;
	strncpy(master->statusStr[prIdx], status, 127);
	master->statusStr[prIdx][127] = 0;
    ReleaseDataport(dpControl);
  }

}

void LoadAC() {

  Dataport *loader;
  SHELEM_OLD *map;

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

  map = (SHELEM_OLD *)malloc( prParam );
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
	sHandle->loadOK = false;
    return;
  }

  printf("Connected to %s (%zd bytes)\n",hitsDpName,hSize);

}

int main(int argc,char* argv[])
{
  bool eos = false;

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
	GetState();
    switch(prState) {

      case COMMAND_LOAD:
        printf("COMMAND: LOAD (%zd,%llu)\n",prParam,prParam2);
        Load();
        if( sHandle->loadOK ) {
          sHandle->desorptionLimit = prParam2; // 0 for endless
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
          // Last update not successful, retry with a longer timeout
			if (dpHit && (GetLocalState() != PROCESS_ERROR)) UpdateHits(dpHit,prIdx,60000);
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
        end = true;
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
		if (dpHit && (GetLocalState() != PROCESS_ERROR)) UpdateHits(dpHit,prIdx,20); // Update hit with 20ms timeout. If fails, probably an other subprocess is updating, so we'll keep calculating and try it later (latest when the simulation is stopped).
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
  //CLOSEDP(dpControl);
  //CLOSEDP(dpHit);
  //Why bother closing dataports? Windows will release handles automatically.
  return 0;

}

bool IsProcessRunning(DWORD pid)
{
	HANDLE process = OpenProcess(SYNCHRONIZE, false, pid);
	DWORD ret = WaitForSingleObject(process, 0);
	CloseHandle(process);
	return ret == WAIT_TIMEOUT;
}