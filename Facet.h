/*
  File:        Facet.h
  Description: Facet strucure
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

#ifndef FACETH
#define FACETH

#include "GLApp/GLApp.h"
#include "Shared.h"
#include "File.h"
#include "PugiXML/pugixml.hpp"
#include "Geometry.h"

class Facet {

public:

  typedef struct {
  
    int nbV;
    int nbF;
    Facet **facets;

  } FACETGROUP;

  typedef struct {
    int u;
    int v;
    int width;
    int height;
  } BOX;

  // Constructor/Desctructor/Initialisation
  Facet(int nbIndex);
  ~Facet();

  // Shared struct
  SHFACET sh;

  int      *indices;      // Indices (Reference to geometry vertex)
  Vector2d *vertices2;    // Vertices (2D plane space, UV coordinates)
  int     *cellPropertiesIds;      // -1 if full element, -2 if outside polygon, otherwise index in meshvector
  CellProperties* meshvector;
  size_t meshvectorsize;

  // Normalized plane equation (ax + by + cz + d = 0)
  double a;
  double b;
  double c;
  double d;
  double err;          // planeity error
  int texDimH;         // Texture dimension (a power of 2)
  int texDimW;         // Texture dimension (a power of 2)
  double tRatio;       // Texture sample per unit
  BOOL	textureVisible; //Draw the texture?
  BOOL  collinear;      //All vertices are on a line (non-simple)
  BOOL	volumeVisible;	//Draw volume?
  //SHELEM *mesh;        // Element mesh
  BOOL    hasMesh;     // Temporary flag (loading)
  
  double *outgassingMap; //outgassing map cell values (loaded from file)

  //Dynamic outgassing stuff
  BOOL textureError;   // Disable rendering if the texture has an error
  BOOL hasOutgassingFile; //true if a desorption file was loaded and had info about this facet
  double totalFlux;
  double totalDose;

  //Parametric stuff
  std::string userOutgassing;
  std::string userSticking;
  std::string userOpacity;

  //Smart selection stuff

  std::vector<NeighborFacet> neighbors;

  SHHITS counterCache; //Local copy of facet counter for the current moment, updated by Worker::Update() or moment change

  // GUI stuff
  BOOL  *visible;         // Edge visible flag
  BOOL   selected;        // Selected flag
  BOX    selectedElem;    // Selected mesh element
  GLint  glElem;          // Surface elements boundaries
  GLint  glSelElem;       // Selected surface elements boundaries
  GLint  glList;          // Geometry with texture
  GLuint glTex;           // Handle to OpenGL texture
  VHIT   *dirCache;       // Direction field cache (for rendering)

  //Facet methods

  void  ConvertOldDesorbType();
  BOOL  IsTXTLinkFacet();
  Vector3d GetRealCenter();
  void  LoadTXT(FileReader *file);
  void  SaveTXT(FileWriter *file);
  void  LoadGEO(FileReader *file,int version,int nbVertex);
  void  LoadSYN(FileReader *file,int version,int nbVertex);
  void  LoadXML(pugi::xml_node f,int nbVertex,BOOL isMolflowFile,int vertexOffset=0);
  void  SaveGEO(FileWriter *file,int idx);
  void  SaveXML_geom(pugi::xml_node f);
  BOOL  IsCoplanarAndEqual(Facet *f,double threshold);
  int   GetIndex(int idx);
  void  Copy(Facet *f,BOOL copyMesh=FALSE);
  void  SwapNormal();
  void  Explode(FACETGROUP *group);
  void  FillVertexArray(Vector3d *v);
  void  BuildMeshList();
  void  InitVisibleEdge();
  BOOL  SetTexture(double width,double height,BOOL useMesh);
  size_t GetGeometrySize();
  size_t GetHitsSize(size_t nbMoments);
  size_t GetTexSwapSize(BOOL useColormap);
  size_t GetTexRamSize(size_t nbMoments);
  size_t GetTexSwapSizeForRatio(double ratio, BOOL useColor);
  size_t GetTexRamSizeForRatio(double ratio, BOOL useMesh, BOOL countDir, size_t nbMoments);
  size_t GetNbCellForRatio(double ratio);
  size_t GetNbCell();
  void  UpdateFlags();
  void  BuildTexture(AHIT *texBuffer,int textureMode,double min,double max,BOOL useColorMap,double dCoeff1,double dCoeff2,double dCoeff3,BOOL doLog,size_t m);
  BOOL  BuildMesh();
  void  BuildSelElemList();
  int   RestoreDeviceObjects();
  int   InvalidateDeviceObjects();
  void  DetectOrientation();
  double GetSmooth(int i,int j,AHIT *texBuffer,int textureMode,double scaleF);
  void  glVertex2u(double u,double v);
  void  ShiftVertex();
  void  RenderSelectedElem();
  void  SelectElem(int u,int v,int width,int height);
  void  UnselectElem();
  float GetMeshArea(int index,BOOL correct2sides = FALSE);
  size_t GetMeshNbPoint(int index);
  Vector2d GetMeshPoint(int index, int pointId);
  Vector2d GetMeshCenter(int index);
  double GetArea();
};

class DeletedFacet {
public:
	Facet *f;
	size_t ori_pos;
	BOOL replaceOri;
};

#endif /* FACETH */
