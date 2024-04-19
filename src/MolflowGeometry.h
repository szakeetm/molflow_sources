
#pragma once

#include "Geometry_shared.h"
#include <pugixml.hpp>
#include "Simulation/MolflowSimGeom.h"
#include <cereal/archives/xml.hpp>
#include "MolflowTypes.h"

#define TEXTURE_MODE_PRESSURE 0
#define TEXTURE_MODE_IMPINGEMENT 1
#define TEXTURE_MODE_DENSITY 2



#define SYNVERSION 12

class Worker;

class MolflowGeometry: public InterfaceGeometry {

public:

	// Constructor/Destructor
	MolflowGeometry();

	// Load
	void LoadGEO(FileReader& file, GLProgress_Abstract& prg, int *version, Worker *worker);
	void LoadSYN(FileReader& file, GLProgress_Abstract& prg, int *version, Worker *worker);
	bool LoadTexturesGEO(FileReader& file, GLProgress_Abstract& prg, const std::shared_ptr<GlobalSimuState> globalState, int version);
	//void ImportDesorption_DES(FileReader& file); //Deprecated
	void ImportDesorption_SYN(FileReader& synFile, const size_t source, const double time,
		const size_t mode, const double eta0, const double alpha, const double cutoffdose,
		const std::vector<std::pair<double, double>> &convDistr,
		GLProgress_Abstract& prg);
	void AnalyzeSYNfile(FileReader& f, GLProgress_Abstract& prg, size_t *nbNewFacet,
		size_t *nbTextured, size_t *nbDifferent);

	// Insert
	void InsertSYN(FileReader& file, GLProgress_Abstract& prg, bool newStr);

	// Save
	void SaveTXT(FileWriter& file, const std::shared_ptr<GlobalSimuState> globalState, bool saveSelected);
	void ExportTextures(FILE *file, int grouping, int mode, const std::shared_ptr<GlobalSimuState> globalState, bool saveSelected);
	void ExportProfiles(FILE *file, int isTXT, Worker *worker);
	void SaveGEO(FileWriter& file, GLProgress_Abstract& prg, const std::shared_ptr<GlobalSimuState> globalState, Worker *worker,
                 bool saveSelected, bool crashSave = false);
	
	void InsertModel(const std::shared_ptr<MolflowSimulationModel> loadedModel, const MolflowInterfaceSettings& interfaceSettings, Worker *work, GLProgress_Abstract& prg, bool newStr);
	bool CompareXML_simustate(const std::string &fileName_lhs, const std::string &fileName_rhs,
                              const std::string &fileName_out, double cmpThreshold) override;
	// Geometry
    void     BuildPipe(double L, double R, double s, int step);
    //void     BuildPrisma(double L, double R, double angle, double s, int step);
	void     LoadProfileGEO(FileReader& file, const std::shared_ptr<GlobalSimuState> globalState, int version);

	// Memory usage (in bytes)
	size_t GetGeometrySize();
	size_t GetHitsSize(size_t nbMoments);

	// AC matrix
	size_t GetMaxElemNumber();

	// Texture scaling
	//TEXTURE_SCALE_TYPE texture_limits[3];   // Min/max values for texture scaling: Pressure/Impingement rate/Density
	AutoScaleMode  texAutoScaleMode = AutoscaleMomentsAndConstFlow;  // Include constant flow when calculating autoscale values: 1 include, 0 moments only, 2 constant flow only
	std::tuple<double, double> GetTextureAutoscaleMinMax();

#pragma region GeometryRender.cpp
	void BuildFacetTextures(const std::shared_ptr<GlobalSimuState> globalState, bool renderRegularTexture, bool renderDirectionTexture);
#pragma endregion

    void SetInterfaceFacets(std::vector<std::shared_ptr<SimulationFacet>> sFacets, bool insert, size_t vertexOffset, int structOffset) override;

private:

	void InsertSYNGeom(FileReader& file, size_t strIdx = 0, bool newStruct = false);
	void SaveProfileGEO(FileWriter& file, const std::shared_ptr<GlobalSimuState> globalState, int super = -1, bool saveSelected = false, bool crashSave = false);

};
