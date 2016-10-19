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
#include "SMP.h"
#include "Utils.h"
#include "GrahamScan.h"
#include "PugiXML/pugixml.hpp"
#include <vector>
#include <sstream>
#include <list>
#include "Clipper\clipper.hpp"

class Worker;

#define SEL_HISTORY  100
#define MAX_SUPERSTR 128
#define SYNVERSION 9
#define GEOVERSION   15
#define XMLVERSION 1

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
	void LoadTXT(FileReader *file, GLProgress *prg);
	void LoadSTR(FileReader *file, GLProgress *prg);
	void LoadSTL(FileReader *file, GLProgress *prg, double scaleFactor);
	void LoadASE(FileReader *file, GLProgress *prg);
	void LoadGEO(FileReader *file, GLProgress *prg, LEAK *pleak, int *nbleakLoad, HIT *pHits, int *nbHHitLoad, int *version, Worker *worker);
	void LoadSYN(FileReader *file, GLProgress *prg, LEAK *pleak, int *nbleakLoad, HIT *pHits, int *nbHHitLoad, int *version);
	bool LoadTextures(FileReader *file, GLProgress *prg, Dataport *dpHit, int version);
	void ImportDesorption_DES(FileReader *file);
	void ImportDesorption_SYN(FileReader *synFile, const size_t &source, const double &time,
		const size_t &mode, const double &eta0, const double &alpha, const double &cutoffdose,
		const std::vector<std::pair<double, double>> &convDistr,
		GLProgress *prg);
	void AnalyzeSYNfile(FileReader *f, GLProgress *progressDlg, int *nbFacet,
		int *nbTextured, int *nbDifferent, GLProgress *prg);
	BOOL IsLoaded();

	// Insert
	void InsertTXT(FileReader *file, GLProgress *prg, BOOL newStr);
	void InsertGEO(FileReader *file, GLProgress *prg, BOOL newStr);
	void InsertSYN(FileReader *file, GLProgress *prg, BOOL newStr);
	void InsertSTL(FileReader *file, GLProgress *prg, double scaleFactor, BOOL newStr);

	// Save
	void SaveTXT(FileWriter *file, Dataport *dhHit, BOOL saveSelected);
	void ExportTextures(FILE *file, int grouping, int mode, Dataport *dhHit, BOOL saveSelected);
	void ExportProfiles(FILE *file, int isTXT, Dataport *dhHit, Worker *worker);
	void SaveGEO(FileWriter *file, GLProgress *prg, Dataport *dpHit, std::vector<std::string> userMoments, Worker *worker,
		BOOL saveSelected, LEAK *pleak, int *nbleakSave, HIT *pHits, int *nbHHitSave, BOOL crashSave = FALSE);
	void SaveSTR(Dataport *dhHit, BOOL saveSelected);
	void SaveXML_geometry(pugi::xml_node saveDoc, Worker *work, GLProgress *prg, BOOL saveSelected);
	BOOL SaveXML_simustate(pugi::xml_node saveDoc, Worker *work, BYTE *buffer, SHGHITS *gHits, int nbLeakSave, int nbHHitSave,
		LEAK *pLeak, HIT *pHits, GLProgress *prg, BOOL saveSelected);
	void LoadXML_geom(pugi::xml_node loadXML, Worker *work, GLProgress *progressDlg);
	void InsertXML(pugi::xml_node loadXML, Worker *work, GLProgress *progressDlg, BOOL newStr);
	BOOL LoadXML_simustate(pugi::xml_node loadXML, Dataport *dpHit, Worker *work, GLProgress *progressDlg);

	void ClearFacetMeshLists();
	void BuildFacetMeshLists();
	void BuildFacetTextures(BYTE * hits);
	void ClearFacetTextures();

	// Selection (drawing stuff)
	void SelectAll();
	void UnSelectAll();
	void SelectArea(int x1, int y1, int x2, int y2, BOOL clear, BOOL unselect, BOOL vertexBound, BOOL circularSelection);
	void Select(int x, int y, BOOL clear, BOOL unselect, BOOL vertexBound, int width, int height);
	void Select(int facet);
	void Select(Facet *f);
	void UnselectAll();
	int  GetNbSelected();
	void UpdateSelection();

	// Vertex Selection (drawing stuff)
	void SelectAllVertex();
	void SelectVertex(int x1, int y1, int x2, int y2, BOOL shiftDown, BOOL ctrlDown, BOOL circularSelection);
	void SelectVertex(int x, int y, BOOL shiftDown, BOOL ctrlDown);
	void SelectVertex(int facet);
	void UnselectAllVertex();
	int  GetNbSelectedVertex();
	void PaintSelectedVertices(BOOL hiddenVertex);

	// OpenGL Rendering/Initialisation
	void Render(GLfloat *matView, BOOL renderVolume, BOOL renderTexture, int showMode, BOOL filter, BOOL showHidden, BOOL showMesh, BOOL showDir);
	int  RestoreDeviceObjects();
	int  InvalidateDeviceObjects();

	// Geometry
	//BOOL     AskToReset_Geom(Worker *work);
	BOOL     IsInFacet(int facet, const double &u, const double &v);
	void     SetFacetTexture(int facet, double ratio, BOOL corrMap);
	void     BuildPipe(double L, double R, double s, int step);
	void     BuildFacetList(Facet *f);
	void     Rebuild();
	void     SetNormeRatio(float r);
	float    GetNormeRatio();
	void     SetAutoNorme(BOOL enable);
	BOOL     GetAutoNorme();
	void     SetCenterNorme(BOOL enable);
	BOOL     GetCenterNorme();
	void	   RebuildLists();
	void	   CalcTotalOutGassing();
	void     InitializeGeometry(int facet_number = -1);           // Initialiase all geometry related variable
	void     LoadProfile(FileReader *file, Dataport *dpHit, int version);
	void UpdateName(FileReader *file);
	void UpdateName(const char *fileName);

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
	double distTraveledTotal_total;
	double distTraveledTotal_fullHitsOnly;

	// Memory usage (in bytes)
	size_t GetGeometrySize();
	size_t GetHitsSize(std::vector<double> *moments);

	// Raw data buffer (geometry)
	void CopyGeometryBuffer(BYTE *buffer);

	// AC matrix
	size_t GetMaxElemNumber();
	void CopyElemBuffer(BYTE *buffer);

	//Geometry Editor functions
	void CheckCollinear();
	void CheckNonSimple();
	void CheckIsolatedVertex();
	void CorrectNonSimple(int *nonSimpleList, int nbNonSimple);
	size_t AnalyzeNeighbors(Worker *work, GLProgress *prg);
	std::vector<size_t> GetConnectedFacets(size_t sourceFacetId, double maxAngleDiff);
	int GetNbVertex();
	VERTEX3D GetFacetCenter(int facet);
	int GetNbStructure();
	char* GetStructureName(int idx);
	void CreatePolyFromVertices_Convex();
	void CreatePolyFromVertices_Order();
	void CreateDifference();
	void ClipSelectedPolygons(ClipperLib::ClipType type, BOOL reverseOrder);
	void ClipPolygon(size_t id1, std::vector<std::vector<size_t>> clippingPaths, ClipperLib::ClipType type);
	void ClipPolygon(size_t id1, size_t id2, ClipperLib::ClipType type);
	void RegisterVertex(Facet *f, const VERTEX2D &vert, size_t id1, const std::vector<ProjectedPoint> &projectedPoints, std::vector<VERTEX3D> &newVertices, size_t registerLocation);
	void SelectCoplanar(int width, int height, double tolerance);
	VERTEX3D *GetVertex(int idx);
	Facet *GetFacet(int facet);
	int GetNbFacet();
	AABB GetBB();
	VERTEX3D GetCenter();
	int AddRefVertex(VERTEX3D *p, VERTEX3D *refs, int *nbRef, double vT);
	void CollapseVertex(Worker *work, GLProgress *prg, double totalWork, double vT);
	BOOL GetCommonEdges(Facet *f1, Facet *f2, int *c1, int *c2, int *chainLength);
	void MoveVertexTo(int idx, double x, double y, double z);
	void SwapNormal();
	void Extrude(int mode, VERTEX3D radiusBase, VERTEX3D offsetORradiusdir, BOOL againstNormal, double distanceORradius, double totalAngle, int steps);
	void ShiftVertex();
	void Merge(int nbV, int nbF, VERTEX3D *nV, Facet **nF);
	void RemoveLinkFacet();
	int HasIsolatedVertices();
	void  DeleteIsolatedVertices(BOOL selectedOnly);
	void  SelectIsolatedVertices();
	BOOL RemoveCollinear();
	void RemoveSelectedVertex();
	void RemoveSelected();
	void RemoveFacets(const std::vector<size_t> &facetIdList, BOOL doNotDestroy = FALSE);
	void RestoreFacets(std::vector<DeletedFacet> deletedFacetList, BOOL toEnd);
	int ExplodeSelected(BOOL toMap = FALSE, int desType = 1, double exponent = 0.0, double *values = NULL);
	BOOL RemoveNullFacet();
	void AlignFacets(int* selection, int nbSelected, int Facet_source, int Facet_dest, int Anchor_source, int Anchor_dest,
		int Aligner_source, int Aligner_dest, BOOL invertNormal, BOOL invertDir1, BOOL invertDir2, BOOL copy, Worker *worker);
	void MoveSelectedFacets(double dX, double dY, double dZ, BOOL copy, Worker *worker);
	void MirrorSelectedFacets(VERTEX3D P0, VERTEX3D N, BOOL copy, Worker *worker);
	void RotateSelectedFacets(const VERTEX3D &AXIS_P0, const VERTEX3D &AXIS_DIR, double theta, BOOL copy, Worker *worker);
	void CloneSelectedFacets();
	void MoveSelectedVertex(double dX, double dY, double dZ, BOOL copy, Worker *worker);
	void AddVertex(double X, double Y, double Z);
	void GetSelection(int **selection, int *nbSel);
	void SetSelection(int **selection, int *nbSel, BOOL isShiftDown, BOOL isCtrlDown);
	void AddStruct(const char *name);
	void DelStruct(int numToDel);
	void ScaleSelectedVertices(VERTEX3D invariant, double factor, BOOL copy, Worker *worker);
	void ScaleSelectedFacets(VERTEX3D invariant, double factorX, double factorY, double factorZ, BOOL copy, Worker *worker);
	BOOL IntersectingPlaneWithLine(const VERTEX3D &P0, const VERTEX3D &u, const VERTEX3D &V0, const VERTEX3D &n, VERTEX3D *intersectPoint, BOOL withinSection);
	std::vector<DeletedFacet> BuildIntersection(size_t *nbCreated);
	std::vector<DeletedFacet> SplitSelectedFacets(const VERTEX3D &base, const VERTEX3D &normal, size_t *nbCreated,/*Worker *worker,*/GLProgress *prg);
	Facet *MergeFacet(Facet *f1, Facet *f2);
	void Collapse(double vT, double fT, double lT, BOOL doSelectedOnly, Worker *work, GLProgress *prg);
	void RenumberNeighbors(const std::vector<int> &newRefs);
	void MergecollinearSides(Facet *f, double lT);
	void CalculateFacetParam_geometry(Facet *f);
	void RemoveFromStruct(int numToDel);
	void CreateLoft();

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

	void CalculateFacetParam(int facetId); // Facet parameters
	void LoadTXTGeom(FileReader *file, size_t *nbV, size_t *nbF, VERTEX3D **V, Facet ***F, size_t strIdx = 0);
	void InsertTXTGeom(FileReader *file, size_t *nbV, size_t *nbF, VERTEX3D **V, Facet ***F, size_t strIdx = 0, BOOL newStruct = FALSE);
	void InsertGEOGeom(FileReader *file, size_t *nbV, size_t *nbF, VERTEX3D **V, Facet ***F, size_t strIdx = 0, BOOL newStruct = FALSE);
	void InsertSYNGeom(FileReader *file, size_t *nbV, size_t *nbF, VERTEX3D **V, Facet ***F, size_t strIdx = 0, BOOL newStruct = FALSE);
	void InsertSTLGeom(FileReader *file, size_t *nbV, size_t *nbF, VERTEX3D **V, Facet ***F, size_t strIdx = 0, double scaleFactor = 1.0, BOOL newStruct = FALSE);
	void SaveProfileGEO(FileWriter *file, Dataport *dpHit, int super = -1, BOOL saveSelected = FALSE, BOOL crashSave = FALSE);
	void SaveProfileTXT(FileWriter *file);
	void AdjustProfile();
	void SaveSuper(Dataport *dpHit, int s);

	BOOL isLoaded;  // Is loaded flag

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


	void DrawFacet(Facet *f, BOOL offset = FALSE, BOOL showHidden = FALSE, BOOL selOffset = FALSE);
	void FillFacet(Facet *f, BOOL addTextureCoord);
	void AddTextureCoord(Facet *f, VERTEX2D *p);
	void DrawPolys();
	void BuildGLList();
	void BuildShapeList();
	void RenderArrow(GLfloat *matView, float dx, float dy, float dz, float px, float py, float pz, float d);
	void DeleteGLLists(BOOL deletePoly = FALSE, BOOL deleteLine = FALSE);
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
	void Triangulate(Facet *f, BOOL addTextureCoord);
	void DrawEar(Facet *f, POLYGON *p, int ear, BOOL addTextureCoord);

};

class PreviewVertex {
public:
	double x, y, z;  //3D coordinates
};

class PreviewFacet {
public:
	std::vector<PreviewVertex> indices;
};

class ClippingVertex {
public:

	ClippingVertex();
	VERTEX2D vertex; //Storing the actual vertex
	BOOL visited;
	BOOL inside;
	BOOL onClippingLine;
	BOOL isLink;
	double distance;
	std::list<ClippingVertex>::iterator link;
	size_t globalId;
};



BOOL operator<(const std::list<ClippingVertex>::iterator& a, const std::list<ClippingVertex>::iterator& b);

#endif /* _GEOMETRYH_ */

