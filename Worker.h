/*
  File:        Worker.h
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

#ifndef _WORKERH_
#define _WORKERH_

#include "Geometry.h"
#include <string>

extern float m_fTime;
class Worker
{

public:

  // Constructor
  Worker();
  ~Worker();

  // Return a handle to the currently loaded geometry
  Geometry *GetGeometry();

  // Load a geometry (throws Error)
  void LoadGeometry(char *fileName);

  // Inserts a new geometry (throws Error)
  void InsertGeometry(BOOL newStr,char *fileName);

  // Load a textures(throws Error)
  void LoadTextures(char *fileName,int version);

  // Save a geometry (throws Error)
  void SaveGeometry(char *fileName,GLProgress *prg,BOOL askConfirm=TRUE,BOOL saveSelected=FALSE,BOOL autoSave=FALSE,BOOL crashSave=FALSE);

  // Save textures (throws Error)
  void ExportTextures(char *fileName,int mode,BOOL askConfirm=TRUE,BOOL saveSelected=FALSE);

  //Import desorption map
  void ImportDesorption(char *fileName);

  // Save a geometry using the current file name (throws Error)
  void SaveGeometry(GLProgress *prg);

  // Return/Set the current filename
  char *GetFileName();
  char *GetShortFileName();
  void  SetFileName(char *fileName);

  // Set number of processes [1..16] (throws Error)
  void SetProcNumber(int n);

  // Get number of processes
  int GetProcNumber();

  // Set the number of maximum desorption
  void SetMaxDesorption(llong max);

  // Get PID
  DWORD GetPID(int prIdx);

  // Reset simulation
  void Reset(float appTime);

  // Reload simulation (throws Error)
  void Reload();
  void RealReload();

  // Switch running/stopped
  void StartStop(float appTime,int mode);

    // Switch running/stopped
  void Stop_Public();

  // AC iteration single step
  void StepAC(float appTime);

  // Free all allocated resource
  void Exit();

  // Kill all sub processes
  void KillAll();

  // Get hit counts for sub process
  void Update(float appTime);

  // Send total and facet hit counts to subprocesses
  void SendHits();

  // Send heartbeat to subprocesses, otherwise they close
  //void SendHeartBeat();
 
  // Get Leak
  void GetLeak(LEAK *buffer,int *nb);

    // Set Leak
  void SetLeak(LEAK *buffer,int *nb,SHGHITS *gHits);

  // Get HHit
  void GetHHit(HIT *buffer,int *nb);

  // Set HHit
  void SetHHit(HIT *buffer,int *nb,SHGHITS *gHits);

  // Get process status
  void GetProcStatus(int *states,char **status);

  // Access to dataport (HIT)
  BYTE *GetHits();
  void  ReleaseHits();

  // Send Compute AC matrix order
  void ComputeAC(float appTime);

  int AddMoment(std::vector<double> newMoments); //Adds a time serie to moments and returns the number of elements
  std::vector<double> ParseMoment(std::string userInput); //Parses a user input and returns a vector of time moments
  void ResetMoments();

  // Global simulation parameters
  llong  nbAbsorption;      // Total number of molecules absorbed (64 bit integer)
  llong  nbDesorption;      // Total number of molecules generated (64 bit integer)
  llong  nbHit;             // Total number of hit (64 bit integer)
  llong  maxDesorption;     // Number of desoprtion before halting
  llong  nbLeakTotal;       // Total number of leak
  double distTraveledTotal; // Total distance traveled by particles (for mean free path calc.)
  int    nbHHit;            // Last hits
  int    nbLastLeaks;       // Last leaks
  BOOL   running;           // Started/Stopped state
  float  startTime;         // Start time
  float  stopTime;          // Stop time
  float  simuTime;          // Total simulation time
  int    mode;              // Simulation mode
  BOOL   calcAC;            // Calculating AC matrix
  int    calcACprg;         // AC matrix progress
  
  std::vector<double> moments;             //moments when a time-dependent simulation state is recorded
  std::vector<std::string> userMoments;    //user-defined text values for defining time moments (can be time or time series)
  std::vector<unsigned int> testStructures; //structures which are test-cubes
  int displayedMoment;

  //gas pulse parameters
  double desorptionStartTime;
  double desorptionStopTime;
  double timeWindowSize;
  BOOL useMaxwellDistribution;
  BOOL calcConstantFlow;
  double valveOpenMoment;

	// Current loaded file
  char fullFileName[512];

private:

  // Process management
  int    nbProcess;
  DWORD  pID[MAX_PROCESS];
  DWORD  pid;
  BOOL   allDone;

  // Geometry handle
  Geometry *geom;

  

  // Dataport handles and names
  Dataport *dpControl;
  Dataport *dpHit;
  char      ctrlDpName[32];
  char      loadDpName[32];
  char      hitsDpName[32];

  // Caches
  HIT  hhitCache[NBHHIT];
  LEAK leakCache[NBHHIT];

  // Methods
  BOOL ExecuteAndWait(int command,int waitState,int param=0);
  BOOL Wait(int waitState,int timeout);
  void ResetWorkerStats();
  void ClearHits();
  char *GetErrorDetails();
  void ThrowSubProcError(char *message=NULL);
  void Start();
  void Stop();
  void OneStep();
  void InnerStop(float appTime);

};

#endif /* _WORKERH_ */
