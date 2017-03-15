#pragma once

#include "Geometry.h"


#define TEXTURE_MODE_PRESSURE 0
#define TEXTURE_MODE_IMPINGEMENT 1
#define TEXTURE_MODE_DENSITY 2
#define SYNVERSION 9

class Worker;

class MolflowGeometry: public Geometry {

public:

	// Constructor/Destructor
	MolflowGeometry();

	// Load
	void LoadGEO(FileReader *file, GLProgress *prg, LEAK *leakCache, size_t *leakCacheSize, HIT *hitCache, size_t *hitCacheSize, int *version, Worker *worker);
	void LoadSYN(FileReader *file, GLProgress *prg, int *version);
	BOOL LoadTextures(FileReader *file, GLProgress *prg, Dataport *dpHit, int version);
	void ImportDesorption_DES(FileReader *file);
	void ImportDesorption_SYN(FileReader *synFile, const size_t &source, const double &time,
		const size_t &mode, const double &eta0, const double &alpha, const double &cutoffdose,
		const std::vector<std::pair<double, double>> &convDistr,
		GLProgress *prg);
	void AnalyzeSYNfile(FileReader *f, GLProgress *progressDlg, int *nbFacet,
		int *nbTextured, int *nbDifferent, GLProgress *prg);

	// Insert
	void InsertSYN(FileReader *file, GLProgress *prg, BOOL newStr);

	// Save
	void SaveTXT(FileWriter *file, Dataport *dhHit, BOOL saveSelected);
	void ExportTextures(FILE *file, int grouping, int mode, Dataport *dhHit, BOOL saveSelected);
	void ExportProfiles(FILE *file, int isTXT, Dataport *dhHit, Worker *worker);
	void SaveGEO(FileWriter *file, GLProgress *prg, Dataport *dpHit, std::vector<std::string> userMoments, Worker *worker,
		BOOL saveSelected, LEAK *pleak, size_t *nbleakSave, HIT *hitCache, size_t *nbHHitSave, BOOL crashSave = FALSE);
	
	void SaveXML_geometry(pugi::xml_node saveDoc, Worker *work, GLProgress *prg, BOOL saveSelected);
	BOOL SaveXML_simustate(pugi::xml_node saveDoc, Worker *work, BYTE *buffer, SHGHITS *gHits, int nbLeakSave, int nbHHitSave,
		LEAK *leakCache, HIT *hitCache, GLProgress *prg, BOOL saveSelected);
	void LoadXML_geom(pugi::xml_node loadXML, Worker *work, GLProgress *progressDlg);
	void InsertXML(pugi::xml_node loadXML, Worker *work, GLProgress *progressDlg, BOOL newStr);
	BOOL LoadXML_simustate(pugi::xml_node loadXML, Dataport *dpHit, Worker *work, GLProgress *progressDlg);

	// Geometry
	void     BuildPipe(double L, double R, double s, int step);
	void     LoadProfile(FileReader *file, Dataport *dpHit, int version);

	// Memory usage (in bytes)
	size_t GetGeometrySize();
	size_t GetHitsSize(std::vector<double> *moments);

	// Raw data buffer (geometry)
	void CopyGeometryBuffer(BYTE *buffer);

	// AC matrix
	size_t GetMaxElemNumber();

	// Texture scaling
	TEXTURE_SCALE_TYPE texture_limits[3];   // Min/max values for texture scaling: Pressure/Impingement rate/Density
	BOOL  texAutoScaleIncludeConstantFlow;  // Include constant flow when calculating autoscale values

#pragma region GeometryRender.cpp
	void BuildFacetTextures(BYTE *hits,BOOL renderRegularTexture,BOOL renderDirectionTexture);
	void BuildFacetDirectionTextures(BYTE *hits);
#pragma endregion

	// Temporary variable (used by LoadXXX)
	double distTraveledTotal_total;
	double distTraveledTotal_fullHitsOnly;

private:

	void InsertSYNGeom(FileReader *file, size_t *nbV, size_t *nbF, InterfaceVertex **V, Facet ***F, size_t strIdx = 0, BOOL newStruct = FALSE);
	void SaveProfileGEO(FileWriter *file, Dataport *dpHit, int super = -1, BOOL saveSelected = FALSE, BOOL crashSave = FALSE);

};

