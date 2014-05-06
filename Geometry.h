/*
  File:        Geometry.h
  Description: Main geometry class (Handles sets of facets)
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

#ifndef _GEOMETRYH_
#define _GEOMETRYH_

#include "Facet.h"
#include "File.h"
#include "Types.h"
#include "GLApp/GLToolkit.h"
//#include "GLApp/GLGradient.h"
#include "GLApp/GLProgress.h"
#include "smp/SMP.h"
#include "Utils.h"
#include "GrahamScan.h"
#include <vector>

extern int changedSinceSave;
class Worker;

#define SEL_HISTORY  100
#define MAX_SUPERSTR 256
#define GEOVERSION   13
#define SYNVERSION 7

#define TEXTURE_MODE_PRESSURE 0
#define TEXTURE_MODE_IMPINGEMENT 1
#define TEXTURE_MODE_DENSITY 2

class Geometry {

public:

  // Constructor/Destructor
  Geometry();
  ~Geometry();

  // Clear this geometry
  void Clear();

  // Load
  void LoadTXT(FileReader *file,GLProgress *prg);
  void LoadSTR(FileReader *file,GLProgress *prg);
  void LoadSTL(FileReader *file,GLProgress *prg,double scaleFactor);
  void LoadASE(FileReader *file,GLProgress *prg);
  void LoadGEO(FileReader *file,GLProgress *prg,LEAK *pleak,int *nbleakLoad,HIT *pHits,int *nbHHitLoad,int *version,Worker *worker);
  void LoadSYN(FileReader *file,GLProgress *prg,LEAK *pleak,int *nbleakLoad,HIT *pHits,int *nbHHitLoad,int *version);
  bool LoadTextures(FileReader *file,GLProgress *prg,Dataport *dpHit,int version);
  void ImportDesorption_DES(FileReader *file);
  void ImportDesorption_SYN(FileReader *synFile, const size_t &source, const double &time,
	  const size_t &mode, const double &eta0, const double &alpha,
	  const std::vector<std::pair<double, double>> &convDistr,
	  GLProgress *prg);
  void AnalyzeSYNfile(FileReader *f, GLProgress *progressDlg, int *nbFacet,
	  int *nbTextured, int *nbDifferent, GLProgress *prg);
  BOOL IsLoaded();

  // Insert
  void InsertTXT(FileReader *file,GLProgress *prg,BOOL newStr);
  void InsertGEO(FileReader *file,GLProgress *prg,BOOL newStr);
  void InsertSYN(FileReader *file,GLProgress *prg,BOOL newStr);
  void InsertSTL(FileReader *file,GLProgress *prg,double scaleFactor,BOOL newStr);

  // Save
  void SaveTXT(FileWriter *file,Dataport *dhHit,BOOL saveSelected);
  void ExportTextures(FILE *file,int mode,Dataport *dhHit,BOOL saveSelected);
  void SaveGEO(FileWriter *file,GLProgress *prg,Dataport *dpHit,std::vector<std::string> userMoments,Worker *worker,
	  BOOL saveSelected,LEAK *pleak,int *nbleakSave,HIT *pHits,int *nbHHitSave,BOOL crashSave=FALSE);
  void SaveSTR(Dataport *dhHit,BOOL saveSelected);


  // Selection (drawing stuff)
  void SelectAll();
  void UnSelectAll();
  void SelectArea(int x1,int y1,int x2,int y2,BOOL clear,BOOL unselect,BOOL vertexBound,BOOL circularSelection);
  void Select(int x,int y,BOOL clear,BOOL unselect,BOOL vertexBound,int width,int height);
  void Select(int facet);
  void Select(Facet *f);
  void Unselect();
  void CheckIsolatedVertex();
  void CheckNonSimple();
  void CheckCollinear();
  int  GetNbSelected();
  void UpdateSelection();
  void SwapNormal();
  void RemoveSelected();
  void RemoveSelectedVertex();
  void RemoveFromStruct(int numToDel);
  void RemoveCollinear();
  int  ExplodeSelected(BOOL toMap=FALSE,int desType=1,double exponent=0.0,double *values=NULL);
  void SelectCoplanar(int width,int height,double tolerance);
  void MoveSelectedVertex(double dX,double dY,double dZ,BOOL copy,Worker *worker);
  void ScaleSelectedVertices(VERTEX3D invariant,double factor,BOOL copy,Worker *worker);
  void ScaleSelectedFacets(VERTEX3D invariant,double factorX,double factorY,double factorZ,BOOL copy,Worker *worker);
  void MoveSelectedFacets(double dX,double dY,double dZ,BOOL copy,Worker *worker);
  void MirrorSelectedFacets(VERTEX3D P0,VERTEX3D N,BOOL copy,Worker *worker);
  void RotateSelectedFacets(VERTEX3D AXIS_P0,VERTEX3D AXIS_DIR,double theta,BOOL copy,Worker *worker);
  void AlignFacets(int* selection,int nbSelected,int Facet_source,int Facet_dest,int Anchor_source,int Anchor_dest,
	  int Aligner_source,int Aligner_dest,BOOL invertNormal,BOOL invertDir1,BOOL invertDir2,BOOL copy,Worker *worker);
  void CloneSelectedFacets();
  void AddVertex(double X,double Y,double Z);
  void CorrectNonSimple(int *nonSimpleList,int nbNonSimple);

  void AddStruct(char *name);
  void DelStruct(int numToDel);

    // Vertex Selection (drawing stuff)
  void SelectAllVertex();
  void CreatePolyFromVertices_Convex(); //create convex facet from selected vertices
  void CreatePolyFromVertices_Order(); //create facet from selected vertices following selection order
  void CreateDifference(); //creates the difference from 2 selected facets
  void SelectVertex(int x1,int y1,int x2,int y2,BOOL shiftDown,BOOL ctrlDown,BOOL circularSelection);
  void SelectVertex(int x,int y,BOOL shiftDown,BOOL ctrlDown);
  void SelectVertex(int facet);
  void UnselectAllVertex();
  int  GetNbSelectedVertex();
  void PaintSelectedVertices(BOOL hiddenVertex);
  //void RemoveSelectedVertex();
  void GetSelection(int **selection,int *nbSel);
  void SetSelection(int **selection,int *nbSel);

  // OpenGL Rendering/Initialisation
  void Render(GLfloat *matView,BOOL renderVolume,BOOL renderTexture,int showMode,BOOL filter,BOOL showHidden,BOOL showMesh,BOOL showDir);
  int  RestoreDeviceObjects();
  int  InvalidateDeviceObjects();

  // Geometry
  //BOOL     AskToReset_Geom(Worker *work);
  int      GetNbFacet();
  int      GetNbVertex();
  int      GetNbStructure();
  char     *GetStructureName(int idx);
  VERTEX3D GetCenter();
  Facet    *GetFacet(int facet);
  VERTEX3D *GetVertex(int idx);
  void     MoveVertexTo(int idx,double x,double y,double z);
  VERTEX3D GetFacetCenter(int facet);
  BOOL     IsInFacet(int facet,double u,double v);
  void     Collapse(double vT,double fT,double lT,BOOL doSelectedOnly,GLProgress *prg);
  void     SetFacetTexture(int facet,double ratio,BOOL corrMap);
  void     BuildPipe(double L,double R,double s,int step);
  void     BuildFacetList(Facet *f);
  AABB     GetBB();
  void     Rebuild();
  void	   MergecollinearSides(Facet *f,double fT);
  void     BuildTexture(BYTE *hits);
  void     ShiftVertex();
  int      HasIsolatedVertices();
  void     DeleteIsolatedVertices(BOOL selectedOnly);
  void	   SelectIsolatedVertices();
  void     SetNormeRatio(float r);
  float    GetNormeRatio();
  void     SetAutoNorme(BOOL enable);
  BOOL     GetAutoNorme();
  void     SetCenterNorme(BOOL enable);
  BOOL     GetCenterNorme();
  void	   RebuildLists();
  void	   CalcTotalOutGassing();
  void     InitializeGeometry(int facet_number=-1,BOOL noVertexShift=FALSE);           // Initialiase all geometry related variable
  void     LoadProfile(FileReader *file,Dataport *dpHit,int version);
 
  // Texture scaling
  int textureMode;                        // Pressure / Impingement rate / Density
  TEXTURE_SCALE_TYPE texture_limits[3];   // Min/max values for texture scaling: Pressure/Impingement rate/Density
  BOOL  texAutoScale;                     // Autoscale flag
  BOOL  texAutoScaleIncludeConstantFlow;  // Include constant flow when calculating autoscale values
  BOOL  texColormap;                      // Colormap flag
  BOOL  texLogScale;                      // Texture in LOG scale

  // Structure viewing (-1 => all)
  int viewStruct;

  // Temporary variable (used by LoadXXX)
  llong tNbHit;
  llong tNbDesorption;
  llong tNbDesorptionMax;
  llong tNbLeak;
  llong tNbAbsorption;
  double distTraveledTotal;

  // Memory usage (in bytes)
  DWORD GetGeometrySize(std::vector<double> *moments);
  DWORD GetHitsSize(std::vector<double> *moments);

  // Raw data buffer (geometry)
  void CopyGeometryBuffer(BYTE *buffer,std::vector<double> *moments);

  // AC matrix
  DWORD GetMaxElemNumber();
  void CopyElemBuffer(BYTE *buffer);

private:

  SHGEOM    sh;
  VERTEX3D  center;                     // Center (3D space)
  char      *strName[MAX_SUPERSTR];     // Structure name
  char      *strFileName[MAX_SUPERSTR]; // Structure file name
  char      strPath[512];               // Path were are stored files (super structure)

  // Geometry
  Facet    **facets;    // All facets of this geometry
  VERTEX3D  *vertices3; // Vertices (3D space)
  AABB bb;              // Global Axis Aligned Bounding Box (AABB)
  int nbSelected;       // Number of selected facets
  int nbSelectedVertex; // Number of selected vertex
  float normeRatio;     // Norme factor (direction field)
  BOOL  autoNorme;      // Auto normalize (direction field)
  BOOL  centerNorme;    // Center vector (direction field)

  void CalculateFacetParam(int facet); // Facet parameters
  void Merge(int nbV,int nbF,VERTEX3D *nV,Facet **nF); // Merge geometry
   void LoadTXTGeom(FileReader *file,int *nbV,int *nbF,VERTEX3D **V,Facet ***F,int strIdx=0);
  void InsertTXTGeom(FileReader *file,int *nbV,int *nbF,VERTEX3D **V,Facet ***F,int strIdx=0,BOOL newStruct=FALSE);
  void InsertGEOGeom(FileReader *file,int *nbV,int *nbF,VERTEX3D **V,Facet ***F,int strIdx=0,BOOL newStruct=FALSE);
  void InsertSYNGeom(FileReader *file,int *nbV,int *nbF,VERTEX3D **V,Facet ***F,int strIdx=0,BOOL newStruct=FALSE);
  void InsertSTLGeom(FileReader *file,int *nbV,int *nbF,VERTEX3D **V,Facet ***F,int strIdx=0,double scaleFactor=1.0,BOOL newStruct=FALSE);
  void RemoveLinkFacet();
  void UpdateName(FileReader *file);
  void SaveProfileGEO(FileWriter *file, Dataport *dpHit, int super = -1, BOOL saveSelected = FALSE, BOOL crashSave=FALSE);
  void SaveProfileTXT(FileWriter *file);
  void AdjustProfile();
  void SaveSuper(Dataport *dpHit,int s);

  BOOL isLoaded;  // Is loaded flag

  // Collapsing stuff
  int  AddRefVertex(VERTEX3D *p,VERTEX3D *refs,int *nbRef);
  void RemoveNullFacet();
  BOOL IsCoplanar(int i1,int i2);
  Facet *MergeFacet(Facet *f1,Facet *f2);
  BOOL GetCommonEdges(Facet *f1,Facet *f2,int *c1,int *c2,int *chainLength);
  void CollapseVertex(GLProgress *prg,double totalWork);
  double vThreshold;

  // Rendering/Selection stuff
  int selectHist[SEL_HISTORY];
  int nbSelectedHist;
  void AddToSelectionHist(int f);
  BOOL AlreadySelected(int f);

  int selectHistVertex[SEL_HISTORY];
  int nbSelectedHistVertex;
  void AddToSelectionHistVertex(int idx);
  BOOL AlreadySelectedVertex(int idx);

  std::vector<int> selectedVertexList;
  void EmptySelectedVertexList();
  void RemoveFromSelectedVertexList(int vertexId);
  void AddToSelectedVertexList(int vertexId);


  void DrawFacet(Facet *f,BOOL offset=FALSE,BOOL showHidden=FALSE,BOOL selOffset=FALSE);
  void FillFacet(Facet *f,BOOL addTextureCoord);
  void AddTextureCoord(Facet *f,VERTEX2D *p);
  void DrawPolys();
  void BuildGLList();
  void BuildShapeList();
  void RenderArrow(GLfloat *matView,float dx,float dy,float dz,float px,float py,float pz,float d);
  void DeleteGLLists(BOOL deletePoly=FALSE,BOOL deleteLine=FALSE);
  void BuildSelectList();
  
  void SetCullMode(int mode);
  GLMATERIAL fillMaterial;
  GLMATERIAL whiteMaterial;
  GLMATERIAL arrowMaterial;
  GLint lineList[MAX_SUPERSTR]; // Compiled geometry (wire frame)
  GLint polyList;               // Compiled geometry (polygon)
  GLint selectList;             // Compiled geometry (selection)
  GLint selectList2;            // Compiled geometry (selection with offset)
  GLint selectList3;            // Compiled geometry (no offset,hidden visible)
  GLint selectListVertex;       // Compiled geometry (selection)
  GLint selectList2Vertex;      // Compiled geometry (selection with offset)
  GLint selectList3Vertex;      // Compiled geometry (no offset,hidden visible)
  GLint arrowList;              // Compiled geometry of arrow used for direction field
  GLint sphereList;             // Compiled geometry of sphere used for direction field

  // Triangulation stuff
  int  FindEar(POLYGON *p);
  void Triangulate(Facet *f,BOOL addTextureCoord);
  void DrawEar(Facet *f,POLYGON *p,int ear,BOOL addTextureCoord);

};

#endif /* _GEOMETRYH_ */

