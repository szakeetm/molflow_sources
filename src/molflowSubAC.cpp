/*
Program:     MolFlow+ / Synrad+
Description: Monte Carlo simulator for ultra-high vacuum and synchrotron radiation
Authors:     Jean-Luc PONS / Roberto KERSEVAN / Marton ADY
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
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
#define NOMINMAX
//#include <windows.h>
//#include <tlhelp32.h>
//#include <Process.h> // For _getpid()
#else
#include <cstring>
#endif


//#include <iostream>

#include <cstdio>
#include <cmath>

#include "Simulation.h"
#include "Parameter.h"


// Global process variables
Simulation* sHandle; //Global handle to simulation, one per subprocess

#define WAITTIME    100  // Answer in STOP mode
//#define TIMEOUT     300  // Process kills itself after no heartbeat (seconds)

static Dataport *dpControl=NULL;
static Dataport *dpHit=NULL;
static Dataport *dpLog = NULL;
//static int       noHeartBeatSince;
static int       prIdx;
static size_t       prState;
static size_t       prParam;
static size_t     prParam2;
static DWORD     hostProcessId;
//static float       heartBeat;
//static HANDLE    masterHandle;
static char      ctrlDpName[32];
static char      loadDpName[32];
static char		 logDpName[32];
static char      hitsDpName[32];

bool endState = false;
bool IsProcessRunning(DWORD pid);

void GetState() {
    prState = PROCESS_READY;
    prParam = 0;

    if( AccessDataport(dpControl) ) {
        SHCONTROL *master = (SHCONTROL *)dpControl->buff;
        prState = master->states[prIdx];
        prParam = master->cmdParam[prIdx];
        prParam2 = master->oldStates[prIdx];
        master->cmdParam[prIdx] = 0;
        master->oldStates[prIdx] = 0;

        ReleaseDataport(dpControl);

        if (!IsProcessRunning(hostProcessId)) {
            printf("Host molflow.exe (process id %d) not running. Closing.",hostProcessId);
            SetErrorSub("Host molflow.exe not running. Closing subprocess.");
            endState = true;
        }
    } else {
	    printf("Subprocess %d couldn't connect to Molflow.\n", prIdx);
	    SetErrorSub("No connection to main program. Closing subprocess.");
        ProcessSleep(5000);
	    endState = true;
  }
}

size_t GetLocalState() {
  return prState;
}

void SetState(size_t state,const char *status,bool changeState, bool changeStatus) {

	prState = state;
	if (changeState) printf("\n setstate %zd \n",state);
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

  static char ret[128];
  size_t count = sHandle->totalDesorbed;
  size_t max   = 0;
  if(sHandle->ontheflyParams.nbProcess)
      max = sHandle->ontheflyParams.desorptionLimit / sHandle->ontheflyParams.nbProcess;

  
  if( GetLocalState()==PROCESS_RUNAC ) sHandle->wp.sMode = AC_MODE;

  switch(sHandle->wp.sMode) {

    case MC_MODE:
      if( max!=0 ) {
        double percent = (double)(count)*100.0 / (double)(max);
        sprintf(ret, "(%s) MC %zd/%zd (%.1f%%)", sHandle->sh.name.c_str(), count, max, percent);
      } else {
        sprintf(ret,"(%s) MC %zd",sHandle->sh.name.c_str(),count);
      }
      break;

    case AC_MODE:
      if( sHandle->prgAC<100 ) {
          sprintf(ret,"(%s) AC (%zdx%zd) (%zd%%)",sHandle->sh.name.c_str(),
                      sHandle->nbAC,sHandle->nbAC,sHandle->prgAC);
      } else {
        if( max!=0 ) {
          double percent = (double)(count)*100.0 / (double)(max);
          sprintf(ret,"(%s) AC (%zdx%zd) %zd/%zd (%.1f%%)",sHandle->sh.name.c_str(),
                      sHandle->nbAC,sHandle->nbAC,count,max,percent);
        } else {
          sprintf(ret,"(%s) AC (%zdx%zd) %zd",sHandle->sh.name.c_str(),sHandle->nbAC,
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
	if (AccessDataport(dpControl)) {
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
  CLOSEDPSUB(loader);

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
    CLOSEDPSUB(loader);
    return;
  }
  CLOSEDPSUB(loader);

  //Connect to log dataport
  if (sHandle->ontheflyParams.enableLogging) {
	  dpLog = OpenDataport(logDpName, sizeof(size_t) + sHandle->ontheflyParams.logLimit * sizeof(ParticleLoggerItem));
	  if (!dpLog) {
		  char err[512];
		  sprintf(err, "Failed to connect to 'dpLog' dataport %s (%zd Bytes)", logDpName, sizeof(size_t) + sHandle->ontheflyParams.logLimit * sizeof(ParticleLoggerItem));
		  SetErrorSub(err);
		  sHandle->loadOK = false;
		  return;
	  }
	  //*((size_t*)dpLog->buff) = 0; //Autofill with 0. Besides, we don't write without access!
  }

  // Connect to hit dataport
  hSize = GetHitsSize();
  dpHit = OpenDataport(hitsDpName,hSize);
  if(!dpHit){
      dpHit = CreateDataport(hitsDpName,hSize);
  }
  if( !dpHit ) {
	  char err[512];
	  sprintf(err, "Failed to connect to 'hits' dataport (%zd Bytes)", hSize);
	  SetErrorSub(err);
	sHandle->loadOK = false;
    return;
  }

  printf("Connected to %s (%zd bytes)\n",hitsDpName,hSize);

}

bool UpdateParams() {

	// Load geometry
	Dataport *loader = OpenDataport(loadDpName, prParam);
	if (!loader) {
		char err[512];
		sprintf(err, "Failed to connect to 'loader' dataport %s (%zd Bytes)", loadDpName, prParam);
		SetErrorSub(err);
		return false;
	}
	printf("Connected to %s\n", loadDpName);

	bool result = UpdateOntheflySimuParams(loader);
	CLOSEDPSUB(loader);

	if (sHandle->ontheflyParams.enableLogging) {
		dpLog = OpenDataport(logDpName, sizeof(size_t) + sHandle->ontheflyParams.logLimit * sizeof(ParticleLoggerItem));
		if (!dpLog) {
			char err[512];
			sprintf(err, "Failed to connect to 'dpLog' dataport %s (%zd Bytes)", logDpName, sizeof(size_t) + sHandle->ontheflyParams.logLimit * sizeof(ParticleLoggerItem));
			SetErrorSub(err);
			return false;
		}
		//*((size_t*)dpLog->buff) = 0; //Autofill with 0, besides we would need access first
	}
	sHandle->tmpParticleLog.clear();
	sHandle->tmpParticleLog.shrink_to_fit();
	if (sHandle->ontheflyParams.enableLogging) sHandle->tmpParticleLog.reserve(sHandle->ontheflyParams.logLimit / sHandle->ontheflyParams.nbProcess);

	return result;
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

    {
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
        const char* dpPrefix = "MFLW";
#else
        const char* dpPrefix = "/MFLW"; // creates semaphore as /dev/sem/%s_sema
#endif
        sprintf(ctrlDpName,"%sCTRL%s",dpPrefix,argv[1]);
        sprintf(loadDpName,"%sLOAD%s",dpPrefix,argv[1]);
        sprintf(hitsDpName,"%sHITS%s",dpPrefix,argv[1]);
        sprintf(logDpName, "%sLOG%s",dpPrefix,argv[1]);
    }
  dpControl = OpenDataport(ctrlDpName,sizeof(SHCONTROL));
  if( !dpControl ) {
    printf("Cannot connect to %s\n",ctrlDpName);
    return 1;
  }

  printf("Connected to %s (%zd bytes), molflowSub.exe #%d\n",ctrlDpName,sizeof(SHCONTROL),prIdx);

  InitSimulation(); //Creates sHandle instance

  // Sub process ready
  SetReady();

  // Main loop
  while( !endState ) {
	GetState();
    switch(prState) {

      case COMMAND_LOAD:
        printf("[%d] COMMAND: LOAD (%zd,%zu)\n", prIdx, prParam,prParam2);
        Load();
        if( sHandle->loadOK ) {
          //sHandle->desorptionLimit = prParam2; // 0 for endless
          SetReady();
        }
        break;

      case COMMAND_LOADAC:
        printf("[%d] COMMAND: LOADAC (%zd)\n", prIdx,prParam);
        LoadAC();
        break;

	  case COMMAND_UPDATEPARAMS:
		  printf("[%d] COMMAND: UPDATEPARAMS (%zd,%zd)\n", prIdx, prParam, prParam2);
		  if (UpdateParams()) {
			  SetState(prParam2, GetSimuStatus());
		  }
		  break;

	  case COMMAND_RELEASEDPLOG:
		  printf("[%d] COMMAND: RELEASEDPLOG (%zd,%zd)\n", prIdx, prParam, prParam2);
		  CLOSEDPSUB(dpLog);
		  SetState(prParam, GetSimuStatus());
		  break;

      case COMMAND_START:
        printf("[%d] COMMAND: START (%zd,%zu)\n", prIdx,prParam,prParam2);
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
        printf("[%d] COMMAND: PAUSE (%zd,%zu)\n", prIdx,prParam,prParam2);
        if( !sHandle->lastHitUpdateOK ) {
          // Last update not successful, retry with a longer timeout
			if (dpHit && (GetLocalState() != PROCESS_ERROR)) UpdateHits(dpHit,dpLog,prIdx,60000);
        }
        SetReady();
        break;

      case COMMAND_RESET:
        printf("[%d] COMMAND: RESET (%zd,%zu)\n", prIdx,prParam,prParam2);
        ResetSimulation();
        SetReady();
        break;

      case COMMAND_EXIT:
        printf("[%d] COMMAND: EXIT (%zd,%zu)\n", prIdx,prParam,prParam2);
        endState = true;
        break;

      case COMMAND_CLOSE:
        printf("[%d] COMMAND: CLOSE (%zd,%zu)\n", prIdx,prParam,prParam2);
        ClearSimulation();
        CLOSEDPSUB(dpHit);
		CLOSEDPSUB(dpLog);
        SetReady();
        break;

      case COMMAND_STEPAC:
        // Debug command
        printf("[%d] COMMAND: STEPAC (%zd,%zu)\n", prIdx,prParam,prParam2);
        if( sHandle->loadOK ) {
          if( StartSimulation(prParam) ) {
            SetState(PROCESS_RUN,GetSimuStatus());
            SimulationACStep(1);
            if(dpHit) UpdateHits(dpHit,dpLog,prIdx,20);
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
        eos = SimulationRun();      // Run during 1 sec
		if (dpHit && (GetLocalState() != PROCESS_ERROR)) {
            UpdateHits(dpHit, dpLog, prIdx,20); // Update hit with 20ms timeout. If fails, probably an other subprocess is updating, so we'll keep calculating and try it later (latest when the simulation is stopped).
		}
		if(eos) {
          if( GetLocalState()!=PROCESS_ERROR ) {
            // Max desorption reached
            SetState(PROCESS_DONE,GetSimuStatus());
            printf("[%d] COMMAND: PROCESS_DONE (Max reached)\n",prIdx);
          }
        }
        break;

      default:

        ProcessSleep(WAITTIME);
        break;
    }
  }
  // Release and clear memory
  SetState(PROCESS_KILLED,"");
  delete sHandle;
  //CLOSEDP(dpControl);
  //CLOSEDP(dpHit);
  //Why bother closing dataports? Windows will release handles automatically.
  return 0;

}

