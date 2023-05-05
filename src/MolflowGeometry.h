/*
Program:     MolFlow+ / Synrad+
Description: Monte Carlo simulator for ultra-high vacuum and synchrotron radiation
Authors:     Jean-Luc PONS / Roberto KERSEVAN / Marton ADY / Pascal BAEHR
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
#pragma once

#include "Geometry_shared.h"
#include "PugiXML/pugixml.hpp"
#include "Simulation/MolflowSimGeom.h"
#include <cereal/archives/xml.hpp>

#define TEXTURE_MODE_PRESSURE 0
#define TEXTURE_MODE_IMPINGEMENT 1
#define TEXTURE_MODE_DENSITY 2



#define SYNVERSION 12

class Worker;

class MolflowGeometry: public Geometry {

public:

	// Constructor/Destructor
	MolflowGeometry();

	// Load
	void LoadGEO(FileReader *file, GLProgress_Abstract& prg, int *version, Worker *worker);
	void LoadSYN(FileReader *file, GLProgress_Abstract& prg, int *version, Worker *worker);
	bool LoadTexturesGEO(FileReader *file, GLProgress_Abstract& prg, GlobalSimuState &globState, int version);
	//void ImportDesorption_DES(FileReader *file); //Deprecated
	void ImportDesorption_SYN(FileReader *synFile, const size_t source, const double time,
		const size_t mode, const double eta0, const double alpha, const double cutoffdose,
		const std::vector<std::pair<double, double>> &convDistr,
		GLProgress_Abstract& prg);
	void AnalyzeSYNfile(FileReader *f, GLProgress_Abstract& prg, size_t *nbNewFacet,
		size_t *nbTextured, size_t *nbDifferent);

	// Insert
	void InsertSYN(FileReader *file, GLProgress_Abstract& prg, bool newStr);

	// Save
	void SaveTXT(FileWriter *file, GlobalSimuState &globState, bool saveSelected);
	void ExportTextures(FILE *file, int grouping, int mode, GlobalSimuState &globState, bool saveSelected);
	void ExportProfiles(FILE *file, int isTXT, Worker *worker);
	void SaveGEO(FileWriter *file, GLProgress_Abstract& prg, GlobalSimuState &globState, Worker *worker,
                 bool saveSelected, bool crashSave = false);
	
	void SaveXML_geometry(pugi::xml_node &saveDoc, Worker *work, GLProgress_Abstract& prg, bool saveSelected);
	bool SaveXML_simustate(pugi::xml_node saveDoc, Worker *work, GlobalSimuState &globState, GLProgress_Abstract& prg, bool saveSelected);
	void LoadXML_geom(pugi::xml_node loadXML, Worker *work, GLProgress_Abstract& prg);
	void InsertXML(pugi::xml_node loadXML, Worker *work, GLProgress_Abstract& prg, bool newStr);
	bool LoadXML_simustate(pugi::xml_node loadXML, GlobalSimuState &globState, Worker *work, GLProgress_Abstract& prg);
    bool CompareXML_simustate(const std::string &fileName_lhs, const std::string &fileName_rhs,
                              const std::string &fileName_out, double cmpThreshold) override;
	// Geometry
    void     BuildPipe(double L, double R, double s, int step);
    void     BuildPrisma(double L, double R, double angle, double s, int step);
	void     LoadProfileGEO(FileReader *file, GlobalSimuState &globState, int version);

	// Memory usage (in bytes)
	size_t GetGeometrySize();
	size_t GetHitsSize(size_t nbMoments);

	// Raw data buffer (geometry)
	void CopyGeometryBuffer(BYTE *buffer,const OntheflySimulationParams& ontheflyParams);

	// AC matrix
	size_t GetMaxElemNumber();

	// Texture scaling
	//TEXTURE_SCALE_TYPE texture_limits[3];   // Min/max values for texture scaling: Pressure/Impingement rate/Density
	short  texAutoScaleIncludeConstantFlow;  // Include constant flow when calculating autoscale values: 1 include, 0 moments only, 2 constant flow only

	

#pragma region GeometryRender.cpp
	void BuildFacetTextures(GlobalSimuState &globState, bool renderRegularTexture, bool renderDirectionTexture);
	void BuildFacetDirectionTextures(BYTE *texture);
#pragma endregion

	void SerializeForLoader(cereal::BinaryOutputArchive&);

    bool InitOldStruct(MolflowSimulationModel *model);
    void InitInterfaceFacets(std::vector<std::shared_ptr<SimulationFacet>> sFacets, Worker* work) override;

private:

	void InsertSYNGeom(FileReader *file, size_t strIdx = 0, bool newStruct = false);
	void SaveProfileGEO(FileWriter *file, GlobalSimuState &globState, int super = -1, bool saveSelected = false, bool crashSave = false);

};
