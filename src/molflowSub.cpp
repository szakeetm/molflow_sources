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
#include "Simulation.h"
#include "Parameter.h"


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



/*void Simulation::GetState() {
    prState = PROCESS_READY;
    prParam = 0;

    if (AccessDataport(this->dpControl)) {
        SHCONTROL *master = (SHCONTROL *) this->dpControl->buff;
        prState = master->states[prIdx];
        prParam = master->cmdParam[prIdx];
        prParam2 = master->oldStates[prIdx];
        //master->cmdParam[prIdx] = 0;
        //master->oldStates[prIdx] = 0;

        ReleaseDataport(this->dpControl);

        if (!IsProcessRunning(hostProcessId)) {
            printf("Host molflow.exe (process id %d) not running. Closing.", hostProcessId);
            SetErrorSub("Host molflow.exe not running. Closing subprocess.");
            endState = true;
        }
    } else {
        printf("Subprocess %d couldn't connect to Molflow.\n", prIdx);
        SetErrorSub("No connection to main program. Closing subprocess.");
        ProcessSleep(5000);
        endState = true;
    }
}*/

/*size_t Simulation::GetLocalState() {
    return prState;
}*/

/*void Simulation::SetState(size_t state, const char *status, bool changeState, bool changeStatus) {

    prState = state;
    if (changeState) printf("\n setstate %zd \n", state);
    if (AccessDataport(this->dpControl)) {
        SHCONTROL *master = (SHCONTROL *) this->dpControl->buff;
        if (changeState) master->states[prIdx] = state;
        if (changeStatus) {
            strncpy(master->statusStr[prIdx], status, 127);
            master->statusStr[prIdx][127] = 0;
        }
        ReleaseDataport(this->dpControl);
    }

}

void Simulation::SetErrorSub(const char *message) {
    printf("Error: %s\n", message);
    SetState(PROCESS_ERROR, message);
}

char *Simulation::GetSimuStatus() {
    static char ret[128];
    size_t count = this->totalDesorbed;
    size_t max = 0;
    if (this->ontheflyParams.nbProcess)
        max = this->ontheflyParams.desorptionLimit / this->ontheflyParams.nbProcess;

    if (max != 0) {
        double percent = (double) (count) * 100.0 / (double) (max);
        sprintf(ret, "(%s) MC %zd/%zd (%.1f%%)", this->sh.name.c_str(), count, max, percent);
    } else {
        sprintf(ret, "(%s) MC %zd", this->sh.name.c_str(), count);
    }

    return ret;
}
*/
/*void Simulation::SetReady() {

    if (this->loadOK)
        SetState(PROCESS_READY, this->GetSimuStatus());
    else
        SetState(PROCESS_READY, "(No geometry)");

}*/

/*void Simulation::SetStatus(char *status) {
    if (AccessDataport(this->dpControl)) {
        SHCONTROL *master = (SHCONTROL *) this->dpControl->buff;
        strncpy(master->statusStr[prIdx], status, 127);
        master->statusStr[prIdx][127] = 0;
        ReleaseDataport(this->dpControl);
    }
}*/

bool Simulation::Load() {

    Dataport *loader;
    size_t hSize;

    // Load geometry
    loader = OpenDataport(loadDpName, procInfo.cmdParam);
    if (!loader) {
        char err[512];
        sprintf(err, "Failed to connect to 'loader' dataport %s (%zd Bytes)", loadDpName, procInfo.cmdParam);
        SetErrorSub(err);
        return this->loadOK;
    }

    printf("Connected to %s\n", loadDpName);

    if (!LoadSimulation(loader)) {
        CLOSEDPSUB(loader);
        return this->loadOK;
    }
    CLOSEDPSUB(loader);

    //Connect to log dataport
    if (this->ontheflyParams.enableLogging) {
        this->dpLog = OpenDataport(logDpName,
                                   sizeof(size_t) + this->ontheflyParams.logLimit * sizeof(ParticleLoggerItem));
        if (!this->dpLog) {
            char err[512];
            sprintf(err, "Failed to connect to 'this->dpLog' dataport %s (%zd Bytes)", logDpName,
                    sizeof(size_t) + this->ontheflyParams.logLimit * sizeof(ParticleLoggerItem));
            SetErrorSub(err);
            this->loadOK = false;
            return this->loadOK;
        }
        //*((size_t*)this->dpLog->buff) = 0; //Autofill with 0. Besides, we don't write without access!
    }

    // Connect to hit dataport
    hSize = GetHitsSize();
    this->dpHit = OpenDataport(hitsDpName, hSize);

    if (!this->dpHit) { // in case of unknown size, create the DP itself
        this->dpHit = CreateDataport(hitsDpName, hSize);
    }
    if (!this->dpHit) {
        char err[512];
        sprintf(err, "Failed to connect to 'hits' dataport (%zd Bytes)", hSize);
        SetErrorSub(err);
        this->loadOK = false;
        return this->loadOK;
    }

    printf("Connected to %s (%zd bytes)\n", hitsDpName, hSize);

    return this->loadOK;
}

bool Simulation::UpdateParams() {

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
    this->tmpParticleLog.clear();
    this->tmpParticleLog.shrink_to_fit();
    if (this->ontheflyParams.enableLogging)
        this->tmpParticleLog.reserve(this->ontheflyParams.logLimit / this->ontheflyParams.nbProcess);

    return result;
}



int main(int argc, char *argv[]) {

    if (argc != 3) {
        printf("Usage: molflowSub peerId index\n");
        return 1;
    }

    size_t hostProcessId = atoi(argv[1]);
    size_t prIdx = atoi(argv[2]);

    Simulation sHandles = {"molflow", "MFLW", hostProcessId, prIdx};
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

    /*{
        Dataport* dpControl = OpenDataport(ctrlDpName, sizeof(SHCONTROL));
        if (!dpControl) {
            printf("Cannot connect to %s\n", ctrlDpName);
            return 1;
        }

        printf("Connected to %s (%zd bytes), molflowSub.exe #%d\n", ctrlDpName, sizeof(SHCONTROL), prIdx);

        //InitSimulation(); //Creates sHandle instance

        // Sub process ready
        for(int i=0;i<1;i++) {
            sHandles.dpControl = dpControl;
            //sHandles[i].SetReady();
        }
    }*/
    InitTick();

    sHandles.controlledLoop();

//Set priority to idle
    /*std::vector<std::thread> threads;
    threads.resize(1);
    for(auto & thread : threads) {
        thread =
                std::thread(&Simulation::controlledLoop, &sHandles, NULL, nullptr); //Launch main loop
        auto myHandle = thread.native_handle();
#ifdef _WIN32
        SetThreadPriority(myHandle, THREAD_PRIORITY_IDLE);
#else
        int policy;
        struct sched_param param;

        pthread_getschedparam(myHandle, &policy, &param);
        param.sched_priority = sched_get_priority_min(policy);
        pthread_setschedparam(myHandle, policy, &param);
        //Check! Some documentation says it's always 0
#endif

        //sHandle-mainL>oop();
    }

    // Release and clear memory
    for(int i=0;i<1;i++) {
        //sHandles[i].SetState(PROCESS_KILLED, "");
        //threads[i].detach();
        ProcessSleep(1000*1);

    }
*/

    //sHandle->SetState(PROCESS_KILLED, "");
    //delete sHandle;
    //CLOSEDP(sHandle->dpControl);
    //CLOSEDP(sHandle->dpHit);
    //Why bother closing dataports? Windows will release handles automatically.
    return 0;

}

