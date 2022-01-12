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
#include <thread>
#include <time.h>
#include "SimulationGPU.h"


// Global process variables
//Simulation* sHandle; //Global handle to simulation, one per subprocess

#define WAITTIME    100  // Answer in STOP mode
//#define TIMEOUT     300  // Process kills itself after no heartbeat (seconds)

/*static Dataport *sHandle->dpControl=NULL;
static Dataport *sHandle->dpHit=NULL;
static Dataport *sHandle->dpLog = NULL;*/
//static int       noHeartBeatSince;
/*static int prIdx;
static size_t prState;
static size_t prParam;
static size_t prParam2;
static DWORD hostProcessId;*/
//static float       heartBeat;
//static HANDLE    masterHandle;
static char ctrlDpName[32];
static char loadDpName[32];
static char logDpName[32];
static char hitsDpName[32];

/*bool SimulationGPU::UpdateParams() {

    // Load geometry
    Dataport *loader = OpenDataport(loadDpName, procInfo.cmdParam);
    if (!loader) {
        char err[512];
        sprintf(err, "Failed to connect to 'loader' dataport %s (%zd Bytes)", loadDpName, procInfo.cmdParam);
        SetErrorSub(err);
        return false;
    }
    printf("Connected to %s\n", loadDpName);

    bool result = UpdateOntheflySimuParams(loader);
    CLOSEDPSUB(loader);

    if (this->ontheflyParams.enableLogging) {
        this->dpLog = OpenDataport(logDpName,
                                   sizeof(size_t) + this->ontheflyParams.logLimit * sizeof(ParticleLoggerItem));
        if (!this->dpLog) {
            char err[512];
            sprintf(err, "Failed to connect to 'this->dpLog' dataport %s (%zd Bytes)", logDpName,
                    sizeof(size_t) + this->ontheflyParams.logLimit * sizeof(ParticleLoggerItem));
            SetErrorSub(err);
            return false;
        }
        //*((size_t*)sHandle->dpLog->buff) = 0; //Autofill with 0, besides we would need access first
    }
    return result;
}*/



int main(int argc, char *argv[]) {
    ProcessSleep(10000);

    if (argc != 3) {
        printf("Usage: gpuSub peerId index\n");
        return 1;
    }

    size_t hostProcessId = atoi(argv[1]);
    size_t prIdx = atoi(argv[2]);

    /*SimulationController simController = {"molflow", "MFLW", hostProcessId, prIdx, new SimulationGPU()};
    {hostProcessId, prIdx, 0, new SimulationGPU(), &procInformation};
    //Simulation *sHandle = new Simulation();

    {
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
        const char* dpPrefix = "MFLW";
#else
        const char *dpPrefix = "/MFLW"; // creates semaphore as /dev/sem/%s_sema
#endif
        sprintf(ctrlDpName, "%sCTRL%s", dpPrefix, argv[1]);
        sprintf(loadDpName, "%sLOAD%s", dpPrefix, argv[1]);
        sprintf(hitsDpName, "%sHITS%s", dpPrefix, argv[1]);
        sprintf(logDpName, "%sLOG%s", dpPrefix, argv[1]);
    }

    InitTick();

    simController.controlledLoop();*/

    return 0;

}

