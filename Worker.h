#pragma once

#include "Parameter.h"
#include "PugiXML/pugixml.hpp"
#include <string>
#include "GLApp/GLTypes.h"
#include "Vector.h"
#include "Smp.h"
#include "MolflowTypes.h"
#include "Buffer_shared.h" //LEAK, HIT

#define CDF_SIZE 100 //points in a cumulative distribution function

class Geometry;
class MolflowGeometry;
class GLProgress;
class LoadStatus;
class FileReader;

class Worker
{

public:

  // Constructor
  Worker();
  ~Worker();

  // Return a handle to the currently loaded geometry
  Geometry *GetGeometry();
  MolflowGeometry* GetMolflowGeometry();

  // Load a geometry (throws Error)

  void LoadGeometry(char *fileName, bool insert=false,bool newStr=false); // Load or insert a geometry (throws Error)
  bool IsDpInitialized();

  // Inserts a new geometry (throws Error)
  //void InsertGeometry(bool newStr,char *fileName);

  // Load a textures(throws Error)
  void LoadTexturesGEO(FileReader *f,int version);
  void RebuildTextures();

  // Save a geometry (throws Error)
  void SaveGeometry(char *fileName,GLProgress *prg,bool askConfirm=true,bool saveSelected=false,bool autoSave=false,bool crashSave=false);

  // Save textures (throws Error)

  void ExportTextures(char *fileName,int grouping,int mode,bool askConfirm=true,bool saveSelected=false);
  void ExportProfiles(char *fileName);
  void ExportAngleMaps(std::vector<size_t> faceList, std::string fileName);

  //Import desorption map
  void ImportDesorption_DES(char *fileName);
  void ImportDesorption_SYN(char *fileName, const size_t &source, const double &time,
	  const size_t &mode, const double &eta0, const double &alpha, const double &cutoffdose,
	  const std::vector<std::pair<double, double>> &convDistr,
	  GLProgress *prg);
  void AnalyzeSYNfile(char *fileName, size_t *nbFacet, size_t *nbTextured, size_t *nbDifferent);
  void ImportCSV(FileReader *file, std::vector<std::vector<std::string>>& table);

  // Return/Set the current filename
  char *GetFileName();
  char *GetShortFileName();
  char *GetShortFileName(char* longFileName);
  void  SetFileName(char *fileName);

  // Set number of processes [1..32] (throws Error)
  void SetProcNumber(size_t n);

  // Get number of processes
  size_t GetProcNumber();

  // Set the number of maximum desorption
  void SetMaxDesorption(llong max);

  // Get PID
  DWORD GetPID(size_t prIdx);

  // Reset simulation
  void ResetStatsAndHits(float appTime);
  void Reload(); //Mark geometry as out of sync with subprocess
  void RealReload(); // Send geometry to subprocess (throws Error)

  // Switch running/stopped
  void StartStop(float appTime,int mode);

    // Switch running/stopped

  void Stop_Public();

  // AC iteration single step
  void StepAC(float appTime);

  // Kill all sub processes
  void KillAll();

  // Get hit counts for sub process
  void Update(float appTime);

  // Send total and facet hit counts to subprocesses
  void SendHits(bool skipFacetHits=false);
  void SetLeakCache(LEAK *buffer,size_t *nb,Dataport* dpHit);
  void SetHitCache(HIT *buffer, size_t *nb,Dataport* dpHit);

  // Get process status
  void GetProcStatus(int *states,std::vector<std::string>& statusStrings);

  //Do calculations necessary before launching simulation
  void PrepareToRun();

  //Get ID of parameter name
  int GetParamId(const std::string); 

  // Access to dataport (HIT)
  BYTE *GetHits();
  void  ReleaseHits();

  // Send Compute AC matrix order
  void ComputeAC(float appTime);

  int AddMoment(std::vector<double> newMoments); //Adds a time serie to moments and returns the number of elements
  std::vector<double> ParseMoment(std::string userInput); //Parses a user input and returns a vector of time moments
  void ResetMoments();
  double GetMoleculesPerTP(size_t moment);

  std::vector<std::pair<double,double>> Generate_ID(int paramId);
  int GenerateNewID(int paramId);
  std::vector<std::pair<double,double>> Generate_CDF(double gasTempKelvins,double gasMassGramsPerMol,size_t size);
  int GenerateNewCDF(double temperature);
  void CalcTotalOutgassing();
  int GetCDFId(double temperature);
  int GetIDId(int paramId);

  // Global simulation parameters
  llong  nbAbsorption;      // Total number of molecules absorbed (64 bit integer)
  llong  nbDesorption;      // Total number of molecules generated (64 bit integer)
  llong  nbHit;             // Total number of hit (64 bit integer)

  llong  desorptionLimit;     // Number of desoprtion before halting

  llong  nbLeakTotal;       // Total number of leak
  double distTraveledTotal_total; // Total distance traveled by particles (for mean pumping path calc.)
  double distTraveledTotal_fullHitsOnly; // Total distance traveled by particles between full hits (for mean free path calc.)

  bool   running;           // Started/Stopped state
  float  startTime;         // Start time
  float  stopTime;          // Stop time
  float  simuTime;          // Total simulation time
  int    mode;              // Simulation mode
  bool   calcAC;            // Calculating AC matrix
  size_t    calcACprg;         // AC matrix progress
  
  std::vector<std::vector<std::pair<double,double>>> CDFs; //cumulative distribution function for each temperature
  std::vector<std::vector<std::pair<double,double>>> IDs; //integrated distribution function for each time-dependent desorption type
  std::vector<double> temperatures; //keeping track of all temperatures that have a CDF already generated
  std::vector<double> moments;             //moments when a time-dependent simulation state is recorded
  std::vector<size_t> desorptionParameterIDs; //time-dependent parameters which are used as desorptions, therefore need to be integrated
  double latestMoment;
  std::vector<std::string> userMoments;    //user-defined text values for defining time moments (can be time or time series)
  //std::vector<unsigned int> testStructures; //structures which are test-cubes
  
  double totalDesorbedMolecules; //Number of molecules desorbed between t=0 and latest_moment
  double finalOutgassingRate; //Number of outgassing molecules / second at latest_moment (constant flow)
  double finalOutgassingRate_Pa_m3_sec; //For the user to see on Global Seetings and in formulas. Not shared with workers
  double gasMass;
  bool   enableDecay;
  double halfLife;
  double timeWindowSize;
  bool useMaxwellDistribution; //true: Maxwell-Boltzmann distribution, false: All molecules have the same (V_avg) speed
  bool calcConstantFlow;

  int motionType;
  Vector3d motionVector1; //base point for rotation
  Vector3d motionVector2; //rotation vector or velocity vector

  bool needsReload;

  std::vector<Parameter> parameters;
  bool abortRequested; //Signal to stop current operation (Collapse, Analyze, etc.)
  int displayedMoment;
  
	// Current loaded file
  char fullFileName[512];

  // Caches
  HIT  hitCache[HITCACHESIZE];
  LEAK leakCache[HITCACHESIZE];
  size_t    hitCacheSize;
  size_t    leakCacheSize;

private:

  // Process management
  size_t    nbProcess;
  DWORD  pID[MAX_PROCESS];
  DWORD  pid;
  bool   allDone;

  // Geometry handle
  MolflowGeometry *geom;

  // Dataport handles and names
  Dataport *dpControl;
  Dataport *dpHit;
  char      ctrlDpName[32];
  char      loadDpName[32];
  char      hitsDpName[32];

  // Methods
  bool ExecuteAndWait(int command, int readyState, size_t param = 0);
  bool Wait(int waitState, LoadStatus *statusWindow = NULL);
  void ResetWorkerStats();
  void ClearHits(bool noReload);
  char *GetErrorDetails();
  void ThrowSubProcError(std::string message);
  void ThrowSubProcError(char *message=NULL);
  void Start();
  void Stop();
  void OneStep();
  void InnerStop(float appTime);

};