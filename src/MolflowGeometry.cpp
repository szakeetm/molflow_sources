#include "Simulation/MolflowSimFacet.h"
#include "MolflowGeometry.h"
#include "MolFlow.h"
#include "Facet_shared.h"
#include "Helper/MathTools.h"
#include "ProfilePlotter.h"
#include "ProfileModes.h"
#include "ConvergencePlotter.h"
//#include "versionId.h"
#include <iomanip>
#include <cfloat> // DBL_EPSILON

#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <Helper/StringHelper.h>

using namespace pugi;

#if defined(MOLFLOW)
extern MolFlow* mApp;
#endif

#if defined(SYNRAD)
extern SynRad* mApp;
#endif

/**
* \brief Basic constructor that initializes a clean (none) geometry
*/
MolflowGeometry::MolflowGeometry() {
	Clear(); //Contains resettexturelimits
}

/**
* \brief Calculates the memory size for the whole geometry
* \return calculated memory usage for the whole geometry
*/
size_t MolflowGeometry::GetGeometrySize() {

	Worker* work = &mApp->worker;

	// Compute number of bytes allocated
	size_t memoryUsage = 0;
	memoryUsage += sizeof(GeomProperties);
	memoryUsage += sizeof(OntheflySimulationParams);
	memoryUsage += sh.nbVertex * sizeof(Vector3d);
	for (int i = 0; i < sh.nbFacet; i++)
		memoryUsage += facets[i]->GetGeometrySize();

	//Parameters
	memoryUsage += sizeof(size_t); //number of parameters
	for (auto& i : work->interfaceParameterCache) {

		memoryUsage += sizeof(size_t); //parameter size

		memoryUsage += i.GetSize() * 2 * sizeof(double);
	}
	memoryUsage += sizeof(size_t); //number of temperatures
	//memoryUsage += sizeof(double) * (int)(work->temperatures).size(); //temperatures

	//interfaceMomentCache size already passed
	memoryUsage += sizeof(double) * (int)(work->interfaceMomentCache).size(); //moments
	return memoryUsage;
}

/**
* \brief Compute number of bytes allocated from the hits size of all facets
* \param nbMoments vector containing all moments
* \return calculated size of memory usage from all facet hits in the geometry
*/
size_t MolflowGeometry::GetHitsSize(const size_t nbMoments) {

	// Compute number of bytes allocated
	size_t memoryUsage = 0;
	memoryUsage += sizeof(GlobalHitBuffer) + (1 + nbMoments) * mApp->worker.model->sp.globalHistogramParams.GetDataSize();
	for (int i = 0; i < sh.nbFacet; i++) {
		memoryUsage += facets[i]->GetHitsSize(nbMoments);
	}

	return memoryUsage;
}

/**
* \brief Compute the maximal (surface) element number (TODO: check if unused)
* \return max element number
*/
size_t MolflowGeometry::GetMaxElemNumber() {

	size_t nbElem = 0;
	for (size_t i = 0; i < sh.nbFacet; i++) {
		InterfaceFacet* f = facets[i];
		if (!f->cellPropertiesIds.empty()) nbElem += f->sh.texWidth * f->sh.texHeight;
		else          return 0;
	}
	return nbElem;

}

/**
* \brief Testing purpose function, construct a PIPE
* \param L length
* \param R radius
* \param s sticking value
* \param step number of facets used to construct the circular hull
*/
void  MolflowGeometry::BuildPipe(double L, double R, double s, int step) {
	Clear();

	sh.nbVertex = 2 * step;
	std::vector<InterfaceVertex>(sh.nbVertex).swap(vertices3);

	sh.nbFacet = step + 2;


	sh.nbSuper = 1;
	structNames = { "Pipe" };

	try {
		facets.resize(sh.nbFacet, nullptr);
	}
	catch (const std::exception&) {
		throw Error("Couldn't allocate memory for facets");
	}

	// Vertices
	for (int i = 0; i < step; i++) {
		double angle = (double)i / (double)step * 2 * PI;
		vertices3[2 * i].x = R * cos(angle);
		vertices3[2 * i].y = R * sin(angle);
		vertices3[2 * i].z = 0.0;
		vertices3[2 * i + 1].x = R * cos(angle);
		vertices3[2 * i + 1].y = R * sin(angle);
		vertices3[2 * i + 1].z = L;
	}

	try {
		// Cap facet
		facets[0] = new InterfaceFacet(step);
		facets[0]->sh.sticking = 1.0;
		facets[0]->sh.desorbType = DES_COSINE;
		facets[0]->sh.outgassing = 1.0;
		for (int i = 0; i < step; i++)
			facets[0]->indices[i] = 2 * i;

		facets[1] = new InterfaceFacet(step);
		facets[1]->sh.sticking = 1.0;
		facets[1]->sh.desorbType = DES_NONE;
		for (int i = 0; i < step; i++)
			facets[1]->indices[step - i - 1] = 2 * i + 1;

		// Wall facet
		for (int i = 0; i < step; i++) {
			facets[i + 2] = new InterfaceFacet(4);
			facets[i + 2]->sh.sticking = s;
			facets[i + 2]->indices[0] = 2 * i;
			facets[i + 2]->indices[1] = 2 * i + 1;
			if (i < step - 1) {
				facets[i + 2]->indices[2] = 2 * (i + 1) + 1;
				facets[i + 2]->indices[3] = 2 * (i + 1);
			}
			else {

				facets[i + 2]->indices[2] = 1;
				facets[i + 2]->indices[3] = 0;
			}
		}
	}
	catch (std::bad_alloc) {
		Clear();
		throw Error("Couldn't reserve memory for the facets");
	}
	catch (...) {
		Clear();
		throw Error("Unspecified Error while building pipe");
	}
	InitializeGeometry();
	
}

/**
* \brief Testing purpose function, construct an angled PRISMA
* \param L length
* \param R radius
* \param s sticking value
* \param step number of facets used to construct the circular hull
*/
/*
void  MolflowGeometry::BuildPrisma(double L, double R, double angle, double s, int step) {
	Clear();

	//mApp->ClearAllSelections();
	//mApp->ClearAllViews();

	int nbTF = 9 * nbDecade;
	int nbTV = 4 * nbTF;

	sh.nbVertex = 2 * step + nbTV;
	std::vector<InterfaceVertex>(sh.nbVertex).swap(vertices3);

	sh.nbFacet = step + 2 + nbTF;


	sh.nbSuper = 1;
	structNames[0]="Prisma";

	try {
		facets.resize(sh.nbFacet, nullptr);
	}
	catch (const std::exception&) {
		throw Error("Couldn't allocate memory for facets");
	}

	// Vertices
	for (int i = 0; i < step; i++) {
		double step_angle = (double)i / (double)step * 2 * PI;
		vertices3[2 * i + nbTV].x = R * cos(step_angle);
		vertices3[2 * i + nbTV].y = R * sin(step_angle);
		vertices3[2 * i + nbTV].z = 0;
		vertices3[2 * i + 1 + nbTV].x = R * cos(step_angle);
		vertices3[2 * i + 1 + nbTV].y = R * sin(step_angle) + L * cos(M_PI_2 - angle);
		vertices3[2 * i + 1 + nbTV].z = L * cos(angle);
	}

	try {
		// Cap facet
		facets[0 + nbTF] = new InterfaceFacet(step);
		facets[0 + nbTF]->sh.sticking = 1.0;
		facets[0 + nbTF]->sh.desorbType = DES_COSINE;
		facets[0 + nbTF]->sh.outgassing = 1.0;
		for (int i = 0; i < step; i++)
			facets[0 + nbTF]->indices[i] = 2 * i + nbTV;

		facets[1 + nbTF] = new InterfaceFacet(step);
		facets[1 + nbTF]->sh.sticking = 1.0;
		facets[1 + nbTF]->sh.desorbType = DES_NONE;
		for (int i = 0; i < step; i++)
			facets[1 + nbTF]->indices[step - i - 1] = 2 * i + 1 + nbTV;

		// Wall facet
		for (int i = 0; i < step; i++) {
			facets[i + 2 + nbTF] = new InterfaceFacet(4);
			//facets[i + 2 + nbTF]->sp.reflection.diffusePart = 1.0; //constructor does this already
			//facets[i + 2 + nbTF]->sp.reflection.specularPart = 0.0; //constructor does this already
			facets[i + 2 + nbTF]->sh.sticking = s;
			facets[i + 2 + nbTF]->indices[0] = 2 * i + nbTV;
			facets[i + 2 + nbTF]->indices[1] = 2 * i + 1 + nbTV;
			if (i < step - 1) {
				facets[i + 2 + nbTF]->indices[2] = 2 * (i + 1) + 1 + nbTV;
				facets[i + 2 + nbTF]->indices[3] = 2 * (i + 1) + nbTV;
			}
			else {

				facets[i + 2 + nbTF]->indices[2] = 1 + nbTV;
				facets[i + 2 + nbTF]->indices[3] = 0 + nbTV;
			}
		}
	}
	catch (std::bad_alloc) {
		Clear();
		throw Error("Couldn't reserve memory for the facets");
	}
	catch (...) {
		throw Error("Unspecified Error while building pipe");
	}
	InitializeGeometry();
	
}
*/

/**
* \brief File handling for inserting a SYN geometry + initialisation
* \param file name of the input file
* \param prg GLProgress_GUI (TODO: which is never used)
* \param newStr newStructure if a super structure will be used or not
*/
void MolflowGeometry::InsertSYN(FileReader& file, GLProgress_Abstract& prg, bool newStr) {

	int structId = viewStruct;
	if (structId == -1) structId = 0;
	InsertSYNGeom(file, structId, newStr);
	InitializeGeometry();
	
	//AdjustProfile();

}

/**
* \brief Inserting the SYN geometry
* \param file name of the input file
* \param strIdx struct ID
* \param newStruct if a super structure will be used or not
*/
void MolflowGeometry::InsertSYNGeom(FileReader& file, size_t strIdx, bool newStruct) {

	UnselectAll();

	file.ReadKeyword("version"); file.ReadKeyword(":");
	int version2;
	version2 = file.ReadInt();
	if (version2 > SYNVERSION) {
		char errMsg[512];
		sprintf(errMsg, "Unsupported SYN version V%d", version2);
		throw Error(errMsg);
	}

	if (version2 >= 12) {
		file.ReadKeyword("newReflectionModel");
		file.ReadKeyword(":");
		/*worker->sp.newReflectionModel =*/ file.ReadInt();
		file.ReadKeyword("lowFluxMode");
		file.ReadKeyword(":");
		/*worker->ontheflyParams.lowFluxMode =*/ file.ReadInt();
		file.ReadKeyword("lowFluxCutoff");
		file.ReadKeyword(":");
		/*worker->ontheflyParams.lowFluxCutoff =*/ file.ReadDouble();
	}

	file.ReadKeyword("totalHit"); file.ReadKeyword(":");
	file.ReadSizeT();
	if (version2 >= 10) {
		file.ReadKeyword("totalHitEquiv"); file.ReadKeyword(":");
		file.ReadDouble();
	}
	file.ReadKeyword("totalDes"); file.ReadKeyword(":");
	file.ReadSizeT();
	if (version2 >= 6) {
		file.ReadKeyword("no_scans"); file.ReadKeyword(":");
		/*loaded_no_scans =*/ file.ReadDouble();
	}

	file.ReadKeyword("totalLeak"); file.ReadKeyword(":");
	file.ReadSizeT();
	if (version2 > 2) {
		file.ReadKeyword("totalFlux"); file.ReadKeyword(":");
		file.ReadDouble();
		file.ReadKeyword("totalPower"); file.ReadKeyword(":");
		file.ReadDouble();
	}
	if (version2 >= 10) {
		file.ReadKeyword("totalAbsEquiv"); file.ReadKeyword(":");
		file.ReadDouble();
		file.ReadKeyword("totalDist"); file.ReadKeyword(":");
		file.ReadDouble();
	}
	file.ReadKeyword("maxDes"); file.ReadKeyword(":");
	file.ReadSizeT();
	file.ReadKeyword("nbVertex"); file.ReadKeyword(":");
	size_t nbNewVertex = file.ReadSizeT();
	file.ReadKeyword("nbFacet"); file.ReadKeyword(":");
	size_t nbNewFacets = file.ReadSizeT();
	file.ReadKeyword("nbSuper"); file.ReadKeyword(":");
	size_t nbNewSuper = file.ReadSizeT();
	file.ReadKeyword("nbFormula"); file.ReadKeyword(":");
	size_t nbF = file.ReadSizeT();
	file.ReadKeyword("nbView"); file.ReadKeyword(":");
	size_t nbV = file.ReadSizeT();
	file.ReadKeyword("nbSelection"); file.ReadKeyword(":");
	size_t nbS = file.ReadSizeT();
	if (version2 > 1) {
		file.ReadKeyword("nbRegions"); file.ReadKeyword(":");
		size_t nbR = file.ReadSizeT();
		file.ReadKeyword("PARfiles"); file.ReadKeyword("{");
		for (size_t i = 0; i < nbR; i++) {
			file.ReadString();
		}
		file.ReadKeyword("}");
	}

	file.ReadKeyword("formulas"); file.ReadKeyword("{");
	for (size_t i = 0; i < nbF; i++) {
		char tmpName[256];
		char tmpExpr[512];
		strcpy(tmpName, file.ReadString());
		strcpy(tmpExpr, file.ReadString());
		//mApp->AddFormula(tmpName,tmpExpr); //we don't add SynRad formulas to MolFlow
	}
	file.ReadKeyword("}");

	file.ReadKeyword("views"); file.ReadKeyword("{");
	for (size_t i = 0; i < nbV; i++) {
		CameraView v;
		v.name=file.ReadString();
		v.projMode = static_cast<ProjectionMode>(file.ReadInt());
		v.camAngleOx = file.ReadDouble();
		v.camAngleOy = file.ReadDouble();
		v.camAngleOz = 0.0; //No support for Z angle in current SYN version
		v.camDist = file.ReadDouble();
		v.camOffset.x = file.ReadDouble();
		v.camOffset.y = file.ReadDouble();
		v.camOffset.z = file.ReadDouble();
		v.performXY = static_cast<CameraPlaneMode>(file.ReadInt());
		v.lightAngleOx = v.lightAngleOy = 0.0;
		v.vLeft = file.ReadDouble();
		v.vRight = file.ReadDouble();
		v.vTop = file.ReadDouble();
		v.vBottom = file.ReadDouble();
		mApp->AddView(v);
	}
	file.ReadKeyword("}");

	file.ReadKeyword("selections"); file.ReadKeyword("{");
	for (size_t i = 0; i < nbS; i++) {
		SelectionGroup s;
		char tmpName[256];
		strcpy(tmpName, file.ReadString());
		s.name = strdup(tmpName);
		size_t nbSel = file.ReadSizeT();
		for (size_t j = 0; j < nbSel; j++) {
			s.facetIds.push_back(file.ReadInt() + sh.nbFacet);
		}
		mApp->AddSelection(s);
	}
	file.ReadKeyword("}");

	file.ReadKeyword("structures"); file.ReadKeyword("{");
	for (size_t i = 0; i < nbNewSuper; i++) {
		structNames.push_back(file.ReadString());
	}
	file.ReadKeyword("}");

	// Reallocate memory
	try {
		facets.resize(nbNewFacets + sh.nbFacet, nullptr);
	}
	catch (const std::exception&) {
		throw Error("Couldn't allocate memory for facets");
	}

	vertices3.resize(nbNewVertex + sh.nbVertex);

	// Read geometry vertices
	file.ReadKeyword("vertices"); file.ReadKeyword("{");
	for (size_t i = sh.nbVertex; i < (sh.nbVertex + nbNewVertex); i++) {
		// Check idx
		int idx = file.ReadInt();
		if (idx != i - sh.nbVertex + 1) throw Error(file.MakeError("Wrong vertex index !"));
		vertices3[i].x = file.ReadDouble();
		vertices3[i].y = file.ReadDouble();
		vertices3[i].z = file.ReadDouble();
		vertices3[i].selected = false;
	}
	file.ReadKeyword("}");

	// Read leaks
	file.ReadKeyword("leaks"); file.ReadKeyword("{");
	file.ReadKeyword("nbLeak"); file.ReadKeyword(":");
	size_t nbleak_local = file.ReadSizeT();
	for (size_t i = 0; i < nbleak_local; i++) {
		size_t idx = file.ReadSizeT();
		//if( idx != i ) throw Error(file.MakeError("Wrong leak index !"));
		file.ReadDouble();
		file.ReadDouble();
		file.ReadDouble();

		file.ReadDouble();
		file.ReadDouble();
		file.ReadDouble();
	}
	file.ReadKeyword("}");
	// Read hit cache

	file.ReadKeyword("hits"); file.ReadKeyword("{");
	file.ReadKeyword("nbHHit"); file.ReadKeyword(":");
	size_t nbHHit_local = file.ReadSizeT();
	for (size_t i = 0; i < nbHHit_local; i++) {
		size_t idx = file.ReadSizeT();
		//if( idx != i ) throw Error(file.MakeError("Wrong hit cache index !"));
		file.ReadDouble(); //x
		file.ReadDouble(); //y
		file.ReadDouble(); //z
		file.ReadDouble(); //dF
		file.ReadDouble(); //dP
		file.ReadInt();    //type
	}
	file.ReadKeyword("}");

	// Read geometry facets (indexed from 1)
	for (size_t i = sh.nbFacet; i < (sh.nbFacet + nbNewFacets); i++) {
		file.ReadKeyword("facet");
		// Check idx
		size_t idx = file.ReadSizeT();
		if (idx != i + 1 - sh.nbFacet) throw Error(file.MakeError("Wrong facet index !"));
		file.ReadKeyword("{");
		file.ReadKeyword("nbIndex");
		file.ReadKeyword(":");
		size_t nb = file.ReadSizeT();

		if (nb < 3) {
			char errMsg[512];
			sprintf(errMsg, "Facet %zd has only %zd vertices. ", i, nb);
			throw Error(errMsg);
		}

		facets[i] = new InterfaceFacet(nb);
		facets[i]->LoadSYN_facet(file, version2, nbNewVertex);
		facets[i]->selected = true;
		for (size_t j = 0; j < nb; j++)
			facets[i]->indices[j] += sh.nbVertex;
		file.ReadKeyword("}");
		if (newStruct) {
			if (facets[i]->sh.superIdx != -1) //-1 = facet member of all structures
				facets[i]->sh.superIdx += static_cast<int>(sh.nbSuper);
			if (facets[i]->sh.superDest > 0) facets[i]->sh.superDest += sh.nbSuper;
		}
		else {
			if (facets[i]->sh.superIdx != -1) //-1 = facet member of all structures
				facets[i]->sh.superIdx += static_cast<int>(strIdx);
			if (facets[i]->sh.superDest > 0) facets[i]->sh.superDest += strIdx;
		}
		if (facets[i]->sh.teleportDest > 0) facets[i]->sh.teleportDest += sh.nbFacet; //Offset teleport target
	}

	sh.nbVertex += nbNewVertex;
	sh.nbFacet += nbNewFacets;
	if (newStruct) sh.nbSuper += nbNewSuper;
	else if (sh.nbSuper < strIdx + nbNewSuper) sh.nbSuper = strIdx + nbNewSuper;
	//return result;
}

/**
* \brief For saving profile data (simulation) into GEO format
* \param file name of the output file
* \param results results from the simulation
* \param super TODO: check if truly unused
* \param saveSelected prevents profile from being saved?
* \param crashSave prevents profile from being saved?
* TODO: Function doesn't seem to cancel properly
*/
void MolflowGeometry::SaveProfileGEO(FileWriter& file, const std::shared_ptr<GlobalSimuState> globalState, int super, bool saveSelected, bool crashSave) {

	//if (!crashSave && !saveSelected) buffer = (BYTE *)buffer->buff;
	file.Write("profiles {\n");
	// Profiles
	int nbProfile = 0;
	int* profileFacet = (int*)malloc((sh.nbFacet) * sizeof(int));
	for (int i = 0; i < sh.nbFacet; i++)
		if ((!saveSelected && !crashSave) && facets[i]->sh.isProfile)
			profileFacet[nbProfile++] = i;

	file.Write(" number: "); file.Write(nbProfile, "\n");
	file.Write(" facets: ");
	for (int i = 0; i < nbProfile; i++) //doesn't execute when crashSave or saveSelected...
		file.Write(profileFacet[i], "\t");

	file.Write("\n");
	size_t facetHitsSize = (1 + mApp->worker.interfaceMomentCache.size()) * sizeof(FacetHitBuffer);
	for (size_t m = 0; (m <= mApp->worker.interfaceMomentCache.size()) || (m == 0); m++) {
		char tmp[128];
		sprintf(tmp, " moment %zd {\n", m);
		file.Write(tmp);

		for (int j = 0; j < PROFILE_SIZE; j++) {
			for (int i = 0; i < nbProfile; i++) { //doesn't execute when crashSave or saveSelected...
				InterfaceFacet* f = GetFacet(profileFacet[i]);
				const std::vector<ProfileSlice>& pr = globalState->facetStates[profileFacet[i]].momentResults[m].profile;
				//char tmp2[128];
				file.Write(static_cast<size_t>(pr[j].countEquiv), "\t"); //Backwards compatibility
				file.Write(pr[j].sum_1_per_ort_velocity, "\t");
				file.Write(pr[j].sum_v_ort);
				file.Write("\t");
			}

			if (nbProfile > 0) file.Write("\n");
		}
		file.Write(" }\n");
	}
	file.Write("}\n");
	SAFE_FREE(profileFacet);
}

/**
* \brief For loading profile data (simulation) from GEO format
* \param file name of the input file
* \param results results from the simulation
* \param version version of the GEO description
*/
void MolflowGeometry::LoadProfileGEO(FileReader& file, const std::shared_ptr<GlobalSimuState> globalState, int version) {

	file.ReadKeyword("profiles"); file.ReadKeyword("{");
	// Profiles
	int nbProfile;
	file.ReadKeyword("number"); file.ReadKeyword(":"); nbProfile = file.ReadInt();
	int* profileFacet = (int*)malloc((nbProfile) * sizeof(int));
	file.ReadKeyword("facets"); file.ReadKeyword(":");
	for (int i = 0; i < nbProfile; i++)
		profileFacet[i] = file.ReadInt();
	size_t facetHitsSize = (1 + mApp->worker.interfaceMomentCache.size()) * sizeof(FacetHitBuffer);
	for (size_t m = 0; m <= mApp->worker.interfaceMomentCache.size() || (version < 10 && m == 0); m++) {
		if (version >= 10) {
			file.ReadKeyword("moment");
			if (m != file.ReadInt()) {
				throw Error("Unexpected profile moment");
				break;
			}
			file.ReadKeyword("{");
		}

		for (int j = 0; j < PROFILE_SIZE; j++) {
			for (int i = 0; i < nbProfile; i++) {
				InterfaceFacet* f = GetFacet(profileFacet[i]);
				std::vector<ProfileSlice>& pr = globalState->facetStates[profileFacet[i]].momentResults[m].profile;
				pr[j].countEquiv = static_cast<double>(file.ReadSizeT());
				if (version >= 13) pr[j].sum_1_per_ort_velocity = file.ReadDouble();
				if (version >= 13) pr[j].sum_v_ort = file.ReadDouble();
			}
		}
		if (version >= 10) file.ReadKeyword("}");
	}
	SAFE_FREE(profileFacet);
}

/**
* \brief For loading geometry data from GEO format
* \param file name of the input file
* \param prg GLProgress_GUI window
* \param version version of the GEO description
* \param worker thread worker that executes the task
*/
void MolflowGeometry::LoadGEO(FileReader& file, GLProgress_Abstract& prg, int* version, Worker* worker) {

	//mApp->ClearAllSelections();
	//mApp->ClearAllViews();
	prg.SetMessage("Clearing current geometry...");
	Clear();
	//mApp->ClearFormulas();

	// Globals
	prg.SetMessage("Reading GEO file header...");
	file.ReadKeyword("version"); file.ReadKeyword(":");
	*version = file.ReadInt();
	if (*version > GEOVERSION) {
		char errMsg[512];
		sprintf(errMsg, "Unsupported GEO version V%d", *version);
		throw Error(errMsg);
	}

	file.ReadKeyword("totalHit"); file.ReadKeyword(":");
	worker->globalState->globalStats.globalHits.nbMCHit = file.ReadSizeT();
	worker->globalState->globalStats.globalHits.nbHitEquiv = static_cast<double>(worker->globalState->globalStats.globalHits.nbMCHit);

	file.ReadKeyword("totalDes"); file.ReadKeyword(":");
	worker->globalState->globalStats.globalHits.nbDesorbed = file.ReadSizeT();

	file.ReadKeyword("totalLeak"); file.ReadKeyword(":");
	worker->globalState->globalStats.nbLeakTotal = file.ReadSizeT();

	if (*version >= 12) {
		file.ReadKeyword("totalAbs"); file.ReadKeyword(":");
		worker->globalState->globalStats.globalHits.nbAbsEquiv = (double)file.ReadSizeT();
		if (*version >= 15) {
			file.ReadKeyword("totalDist_total");
		}
		else { //between versions 12 and 15
			file.ReadKeyword("totalDist");
		}
		file.ReadKeyword(":");
		worker->globalState->globalStats.distTraveled_total = file.ReadDouble();
		if (*version >= 15) {
			file.ReadKeyword("totalDist_fullHitsOnly"); file.ReadKeyword(":");
			worker->globalState->globalStats.distTraveledTotal_fullHitsOnly = file.ReadDouble();
		}
	}
	else {
		worker->globalState->globalStats.globalHits.nbAbsEquiv = 0.0;
		worker->globalState->globalStats.distTraveled_total = 0.0;
		worker->globalState->globalStats.distTraveledTotal_fullHitsOnly = 0.0;
	}
	file.ReadKeyword("maxDes"); file.ReadKeyword(":");
	worker->model->otfParams.desorptionLimit = file.ReadSizeT();
	file.ReadKeyword("nbVertex"); file.ReadKeyword(":");
	sh.nbVertex = file.ReadInt();
	file.ReadKeyword("nbFacet"); file.ReadKeyword(":");
	sh.nbFacet = file.ReadInt();
	file.ReadKeyword("nbSuper"); file.ReadKeyword(":");
	sh.nbSuper = file.ReadInt();
	int nbF = 0; std::vector<std::vector<std::string>> loadFormulas;
	int nbV = 0;
	if (*version >= 2) {
		file.ReadKeyword("nbFormula"); file.ReadKeyword(":");
		nbF = file.ReadInt(); loadFormulas.reserve(nbF);
		file.ReadKeyword("nbView"); file.ReadKeyword(":");
		nbV = file.ReadInt();
	}
	int nbS = 0;
	if (*version >= 8) {
		file.ReadKeyword("nbSelection"); file.ReadKeyword(":");
		nbS = file.ReadInt();
	}
	if (*version >= 7) {
		file.ReadKeyword("gasMass"); file.ReadKeyword(":");
		worker->model->sp.gasMass = file.ReadDouble();
	}
	if (*version >= 16) { //time-dependent version with variable time windows
		file.ReadKeyword("userMoments"); file.ReadKeyword("{");
		file.ReadKeyword("nb"); file.ReadKeyword(":");
		int nb = file.ReadInt();

		for (int i = 0; i < nb; i++) {
			UserMoment um;
			um.content=file.ReadString();
			um.timeWindow=file.ReadDouble();
			worker->userMoments.emplace_back(um);
		}
		file.ReadKeyword("}");
	}
	else if (*version >= 10) { //time-dependent version with fixed time window length
		file.ReadKeyword("userMoments"); file.ReadKeyword("{");
		file.ReadKeyword("nb"); file.ReadKeyword(":");
		int nb = file.ReadInt();

		for (int i = 0; i < nb; i++) {
			UserMoment um;
			um.content= file.ReadString();
			um.timeWindow=0.0; // Try to set a fixed time window later for GEO versions >= 11
			worker->userMoments.emplace_back(um);
		}
		file.ReadKeyword("}");
	}
	if (*version >= 11) { //pulse version
		file.ReadKeyword("desorptionStart"); file.ReadKeyword(":");
		/*worker->desorptionStartTime =*/ file.ReadDouble();

		file.ReadKeyword("desorptionStop"); file.ReadKeyword(":");
		/*worker->desorptionStopTime =*/ file.ReadDouble();

		file.ReadKeyword("timeWindow"); file.ReadKeyword(":");
		worker->model->sp.timeWindowSize = file.ReadDouble();

		if (*version < 16) { // use fixed time window for user moments
			for (auto& uMoment : worker->userMoments) {
				uMoment.timeWindow = worker->model->sp.timeWindowSize;
			}
		}
		file.ReadKeyword("useMaxwellian"); file.ReadKeyword(":");
		worker->model->sp.useMaxwellDistribution = file.ReadInt();
	}

	if (*version >= 12) { //2013.aug.22
		file.ReadKeyword("calcConstantFlow"); file.ReadKeyword(":");
		worker->model->sp.calcConstantFlow = file.ReadInt();

	}
	if (*version >= 2) {
		file.ReadKeyword("formulas"); file.ReadKeyword("{");
		for (int i = 0; i < nbF; i++) {
			char tmpName[256];
			char tmpExpr[512];
			strcpy(tmpName, file.ReadString());
			strcpy(tmpExpr, file.ReadString());
			//mApp->AddFormula(tmpName, tmpExpr); //parse after selection groups are loaded

			std::vector<std::string> newFormula;
			newFormula.emplace_back(tmpName);
			newFormula.emplace_back(tmpExpr);
			loadFormulas.push_back(newFormula);
		}
		file.ReadKeyword("}");

		file.ReadKeyword("views"); file.ReadKeyword("{");
		for (int i = 0; i < nbV; i++) {
			CameraView v;
			v.name=file.ReadString();
			v.projMode = static_cast<ProjectionMode>(file.ReadInt());
			v.camAngleOx = file.ReadDouble();
			v.camAngleOy = file.ReadDouble();
			v.camAngleOz = 0.0; //No support for Z angle in current GEO version
			v.camDist = file.ReadDouble();
			v.camOffset.x = file.ReadDouble();
			v.camOffset.y = file.ReadDouble();
			v.camOffset.z = file.ReadDouble();
			v.performXY = static_cast<CameraPlaneMode>(file.ReadInt());
			v.lightAngleOx = v.lightAngleOy = 0.0;
			v.vLeft = file.ReadDouble();
			v.vRight = file.ReadDouble();
			v.vTop = file.ReadDouble();
			v.vBottom = file.ReadDouble();
			mApp->AddView( v);
		}
		file.ReadKeyword("}");
	}

	if (*version >= 8) {
		file.ReadKeyword("selections"); file.ReadKeyword("{");
		for (int i = 0; i < nbS; i++) {
			SelectionGroup s;
			char tmpName[256];
			strcpy(tmpName, file.ReadString());
			s.name = strdup(tmpName);
			int nbSel = file.ReadInt();

			for (int j = 0; j < nbSel; j++) {
				s.facetIds.push_back(file.ReadInt());
			}
			mApp->AddSelection(s);
		}
		file.ReadKeyword("}");
	}

	for (int i = 0; i < nbF; i++) { //parse formulas now that selection groups are loaded
		mApp->appFormulas->AddFormula(loadFormulas[i][0], loadFormulas[i][1]);
	}

	file.ReadKeyword("structures"); file.ReadKeyword("{");
	for (int i = 0; i < sh.nbSuper; i++) {
		structNames.push_back(file.ReadString());
	}
	file.ReadKeyword("}");

	// Allocate memory
	try {
		facets.resize(sh.nbFacet, nullptr);
	}
	catch (const std::exception&) {
		throw Error("Couldn't allocate memory for facets");
	}
	std::vector<InterfaceVertex>(sh.nbVertex).swap(vertices3);

	// Read vertices
	prg.SetMessage("Reading vertices...");
	file.ReadKeyword("vertices"); file.ReadKeyword("{");
	for (int i = 0; i < sh.nbVertex; i++) {
		// Check idx
		int idx = file.ReadInt();
		if (idx != i + 1) throw Error(file.MakeError("Wrong vertex index !"));
		vertices3[i].x = file.ReadDouble();
		vertices3[i].y = file.ReadDouble();
		vertices3[i].z = file.ReadDouble();
		vertices3[i].selected = false;
	}
	file.ReadKeyword("}");

	if (*version >= 6) {
		prg.SetMessage("Reading leaks and hits...");
		// Read leaks
		file.ReadKeyword("leaks"); file.ReadKeyword("{");
		file.ReadKeyword("nbLeak"); file.ReadKeyword(":");
		worker->globalState->globalStats.leakCacheSize = file.ReadInt();
		for (int i = 0; i < worker->globalState->globalStats.leakCacheSize; i++) {
			int idx = file.ReadInt();
			if (idx != i) throw Error(file.MakeError("Wrong leak index !"));
			if (i < LEAKCACHESIZE) {
				worker->globalState->globalStats.leakCache[i].pos.x = file.ReadDouble();
				worker->globalState->globalStats.leakCache[i].pos.y = file.ReadDouble();
				worker->globalState->globalStats.leakCache[i].pos.z = file.ReadDouble();

				worker->globalState->globalStats.leakCache[i].dir.x = file.ReadDouble();
				worker->globalState->globalStats.leakCache[i].dir.y = file.ReadDouble();
				worker->globalState->globalStats.leakCache[i].dir.z = file.ReadDouble();
			}
			else { //Saved file has more leaks than we could load
				for (int skipIndex = 0; skipIndex < 6; skipIndex++)
					file.ReadDouble();
			}
		}
		file.ReadKeyword("}");

		// Read hit cache
		file.ReadKeyword("hits"); file.ReadKeyword("{");
		file.ReadKeyword("nbHHit"); file.ReadKeyword(":");
		worker->globalState->globalStats.hitCacheSize = file.ReadInt();
		for (int i = 0; i < worker->globalState->globalStats.hitCacheSize; i++) {
			int idx = file.ReadInt();
			if (idx != i) throw Error(file.MakeError("Wrong hit cache index !"));
			if (i < HITCACHESIZE) {
				worker->globalState->globalStats.hitCache[i].pos.x = file.ReadDouble();
				worker->globalState->globalStats.hitCache[i].pos.y = file.ReadDouble();
				worker->globalState->globalStats.hitCache[i].pos.z = file.ReadDouble();

				worker->globalState->globalStats.hitCache[i].type = file.ReadInt();
			}
			else { //Saved file has more hits than we could load
				for (int i = 0; i < 3; i++)
					file.ReadDouble();
				file.ReadInt();
			}
		}
		file.ReadKeyword("}");
	}

	// Read facets
	prg.SetMessage("Reading facets...");
	for (int i = 0; i < sh.nbFacet; i++) {
		file.ReadKeyword("facet");
		// Check idx
		int idx = file.ReadInt();
		if (idx != i + 1) throw Error(file.MakeError("Wrong facet index !"));
		file.ReadKeyword("{");
		file.ReadKeyword("nbIndex");
		file.ReadKeyword(":");
		int nbI = file.ReadInt();
		if (nbI < 3) {
			char errMsg[512];
			sprintf(errMsg, "Facet %d has only %d vertices. ", i + 1, nbI);
			throw Error(errMsg);
		}
		prg.SetProgress((double)i / (double)sh.nbFacet);
		facets[i] = new InterfaceFacet(nbI);
		facets[i]->LoadGEO(file, *version, sh.nbVertex);
		file.ReadKeyword("}");
	}

	InitializeGeometry();
	
	//AdjustProfile();
	//isLoaded = true; //InitializeGeometry() sets to true
	UpdateName(file);

	// Update mesh
	prg.SetMessage("Building mesh...");
	for (int i = 0; i < sh.nbFacet; i++) {
		double p = (double)i / (double)sh.nbFacet;
		prg.SetProgress(p);
		InterfaceFacet* f = facets[i];
		if (!f->SetTexture(f->sh.texWidth_precise, f->sh.texHeight_precise, f->hasMesh)) {
			char errMsg[512];
			sprintf(errMsg, "Not enough memory to build mesh on Facet %d. ", i + 1);
			throw Error(errMsg);
		}
		BuildFacetList(f);
	}

}

/**
* \brief For loading geometry data from SYN format
* \param file name of the input file
* \param prg GLProgress_GUI window
* \param version version of the SYN description
* \param worker thread worker that executes the task
*/
void MolflowGeometry::LoadSYN(FileReader& file, GLProgress_Abstract& prg, int* version, Worker* worker) {

	//mApp->ClearAllSelections();
	//mApp->ClearAllViews();
	prg.SetMessage("Clearing current geometry...");
	Clear();
	//mApp->ClearFormulas();

	// Globals
	prg.SetMessage("Reading SYN file header...");
	file.ReadKeyword("version"); file.ReadKeyword(":");
	*version = file.ReadInt();
	if (*version > SYNVERSION) {
		char errMsg[512];
		sprintf(errMsg, "Unsupported SYN version V%d", *version);
		throw Error(errMsg);
	}

	if (*version >= 12) {
		file.ReadKeyword("newReflectionModel");
		file.ReadKeyword(":");
		/*worker->sp.newReflectionModel =*/ file.ReadInt();
		file.ReadKeyword("lowFluxMode");
		file.ReadKeyword(":");
		/*worker->ontheflyParams.lowFluxMode =*/ file.ReadInt();
		file.ReadKeyword("lowFluxCutoff");
		file.ReadKeyword(":");
		/*worker->ontheflyParams.lowFluxCutoff =*/ file.ReadDouble();
	}

	file.ReadKeyword("totalHit"); file.ReadKeyword(":");
	worker->globalState->globalStats.globalHits.nbMCHit = 0;
	worker->globalState->globalStats.globalHits.nbHitEquiv = 0.0;
	file.ReadSizeT();
	if (*version >= 10) {
		file.ReadKeyword("totalHitEquiv"); file.ReadKeyword(":");
		file.ReadDouble();
	}
	file.ReadKeyword("totalDes"); file.ReadKeyword(":");
	worker->globalState->globalStats.globalHits.nbDesorbed = 0;
	file.ReadSizeT();
	if (*version >= 6) {
		file.ReadKeyword("no_scans"); file.ReadKeyword(":");
		/*loaded_no_scans = */file.ReadDouble();
	}
	file.ReadKeyword("totalLeak"); file.ReadKeyword(":");
	worker->globalState->globalStats.nbLeakTotal = 0; file.ReadSizeT();
	if (*version > 2) {
		file.ReadKeyword("totalFlux"); file.ReadKeyword(":");
		file.ReadDouble();
		file.ReadKeyword("totalPower"); file.ReadKeyword(":");
		file.ReadDouble();
	}
	if (*version >= 10) {
		file.ReadKeyword("totalAbsEquiv"); file.ReadKeyword(":");
		file.ReadDouble();
		file.ReadKeyword("totalDist"); file.ReadKeyword(":");
		file.ReadDouble();
	}
	file.ReadKeyword("maxDes"); file.ReadKeyword(":");
	worker->model->otfParams.desorptionLimit = 0; file.ReadSizeT();
	file.ReadKeyword("nbVertex"); file.ReadKeyword(":");
	sh.nbVertex = file.ReadInt();
	file.ReadKeyword("nbFacet"); file.ReadKeyword(":");
	sh.nbFacet = file.ReadInt();
	file.ReadKeyword("nbSuper"); file.ReadKeyword(":");
	sh.nbSuper = file.ReadInt();
	int nbF = 0;
	int nbV = 0;

	file.ReadKeyword("nbFormula"); file.ReadKeyword(":");
	nbF = file.ReadInt();
	file.ReadKeyword("nbView"); file.ReadKeyword(":");
	nbV = file.ReadInt();
	int nbS = 0;
	file.ReadKeyword("nbSelection"); file.ReadKeyword(":");
	nbS = file.ReadInt();

	if (*version > 1) {
		file.ReadKeyword("nbRegions"); file.ReadKeyword(":");
		int nbR = file.ReadInt();
		//result=PARfileList(nbR);

		file.ReadKeyword("PARfiles"); file.ReadKeyword("{");
		for (int i = 0; i < nbR; i++) {
			/*char tmp[512];
			strcpy(tmp,file.ReadString());
			result.fileNames[i]=strdup(tmp);*/
			file.ReadString();
		}
		file.ReadKeyword("}");
	}

	file.ReadKeyword("formulas"); file.ReadKeyword("{");
	for (int i = 0; i < nbF; i++) {
		char tmpName[256];
		char tmpExpr[512];
		strcpy(tmpName, file.ReadString());
		strcpy(tmpExpr, file.ReadString());
		//mApp->AddFormula(tmpName, tmpExpr);
	}
	file.ReadKeyword("}");

	file.ReadKeyword("views"); file.ReadKeyword("{");
	for (int i = 0; i < nbV; i++) {
		
		CameraView v;
		v.name = file.ReadString();
		v.projMode = static_cast<ProjectionMode>(file.ReadInt());
		v.camAngleOx = file.ReadDouble();
		v.camAngleOy = file.ReadDouble();
		v.camAngleOz = 0.0; //No support for Z angle in current SYN version
		v.camDist = file.ReadDouble();
		v.lightAngleOx = v.lightAngleOy = 0.0;
		v.camOffset.x = file.ReadDouble();
		v.camOffset.y = file.ReadDouble();
		v.camOffset.z = file.ReadDouble();
		v.performXY = static_cast<CameraPlaneMode>(file.ReadInt());
		v.vLeft = file.ReadDouble();
		v.vRight = file.ReadDouble();
		v.vTop = file.ReadDouble();
		v.vBottom = file.ReadDouble();
		mApp->AddView(v);
	}
	file.ReadKeyword("}");

	file.ReadKeyword("selections"); file.ReadKeyword("{");
	for (int i = 0; i < nbS; i++) {
		SelectionGroup s;
		char tmpName[256];
		strcpy(tmpName, file.ReadString());
		s.name = strdup(tmpName);
		int nbSel = file.ReadInt();

		for (int j = 0; j < nbSel; j++) {
			s.facetIds.push_back(file.ReadInt());
		}
		mApp->AddSelection(s);
	}
	file.ReadKeyword("}");

	file.ReadKeyword("structures"); file.ReadKeyword("{");
	for (int i = 0; i < sh.nbSuper; i++) {
		structNames.push_back(file.ReadString());
	}
	file.ReadKeyword("}");

	// Allocate memory
	try {
		facets.resize(sh.nbFacet, nullptr);
	}
	catch (const std::exception&) {
		throw Error("Couldn't allocate memory for facets");
	}

	vertices3.resize(sh.nbVertex); vertices3.shrink_to_fit();

	// Read vertices
	prg.SetMessage("Reading vertices...");
	file.ReadKeyword("vertices"); file.ReadKeyword("{");
	for (int i = 0; i < sh.nbVertex; i++) {
		// Check idx
		int idx = file.ReadInt();
		if (idx != i + 1) throw Error(file.MakeError("Wrong vertex index !"));
		vertices3[i].x = file.ReadDouble();
		vertices3[i].y = file.ReadDouble();
		vertices3[i].z = file.ReadDouble();
		vertices3[i].selected = false;
	}
	file.ReadKeyword("}");
	prg.SetMessage("Reading leaks and hits...");
	// Read leaks
	file.ReadKeyword("leaks"); file.ReadKeyword("{");
	file.ReadKeyword("nbLeak"); file.ReadKeyword(":");

	int nbleak_local = file.ReadInt();
	for (int i = 0; i < nbleak_local; i++) {
		int idx = file.ReadInt();
		if (idx != i) throw Error(file.MakeError("Wrong leak index !"));
		//(pleak+i)->pos.x = 
		file.ReadDouble();
		//(pleak+i)->pos.y = 
		file.ReadDouble();
		//(pleak+i)->pos.z = 
		file.ReadDouble();

		//(pleak+i)->dir.x = 
		file.ReadDouble();
		//(pleak+i)->dir.y = 
		file.ReadDouble();
		//(pleak+i)->dir.z = 
		file.ReadDouble();
	}
	file.ReadKeyword("}");

	// Read hit cache
	file.ReadKeyword("hits"); file.ReadKeyword("{");
	file.ReadKeyword("nbHHit"); file.ReadKeyword(":");

	int nbHHit_local = file.ReadInt();
	for (int i = 0; i < nbHHit_local; i++) {
		int idx = file.ReadInt();
		if (idx != i) throw Error(file.MakeError("Wrong hit cache index !"));
		//(hitCache+i)->pos.x = 
		file.ReadDouble();
		//(hitCache+i)->pos.y = 
		file.ReadDouble();
		//(hitCache+i)->pos.z = 
		file.ReadDouble();
		//(hitCache+i)->dF = 
		file.ReadDouble();
		//(hitCache+i)->dP = 
		file.ReadDouble();
		//(hitCache+i)->type = 
		file.ReadInt();
	}
	file.ReadKeyword("}");
	// Read facets
	prg.SetMessage("Reading facets...");
	for (int i = 0; i < sh.nbFacet; i++) {
		file.ReadKeyword("facet");
		// Check idx
		int idx = file.ReadInt();
		if (idx != i + 1) throw Error(file.MakeError("Wrong facet index !"));
		file.ReadKeyword("{");
		file.ReadKeyword("nbIndex");
		file.ReadKeyword(":");
		int nbI = file.ReadInt();
		if (nbI < 3) {
			char errMsg[512];
			sprintf(errMsg, "Facet %d has only %d vertices. ", i + 1, nbI);
			throw Error(errMsg);
		}
		prg.SetProgress((double)i / (double)sh.nbFacet);
		facets[i] = new InterfaceFacet(nbI);
		facets[i]->LoadSYN_facet(file, *version, sh.nbVertex);
		file.ReadKeyword("}");
	}

	prg.SetMessage("Initalizing geometry and building mesh...");
	InitializeGeometry();
	
	//AdjustProfile();
	//isLoaded = true; //InitializeGeometry() sets to true
	UpdateName(file);

	// Update mesh
	prg.SetMessage("Drawing textures...");
	for (int i = 0; i < sh.nbFacet; i++) {
		double p = (double)i / (double)sh.nbFacet;
		prg.SetProgress(p);
		InterfaceFacet* f = facets[i];
		BuildFacetList(f);
	}
}

/**
* \brief For loading texture data from GEO format
* \param file name of the input file
* \param prg GLProgress_GUI window
* \param results simulation results describing the texture
* \param version version of the GEO description
*/
bool MolflowGeometry::LoadTexturesGEO(FileReader& file, GLProgress_Abstract& prg, const std::shared_ptr<GlobalSimuState> globalState, int version) {

	if (file.SeekFor("{textures}")) {
		char tmp[256];
		//versions 3+
		// Block dpHit during the whole disc reading

		//TextureCell readVal;

		// Globals

		/*gHits->globalStats.hit.nbMCHit = loaded_nbMCHit;
		gHits->globalStats.hit.nbDesorbed = loaded_nbDesorption;
		gHits->globalStats.hit.nbAbsEquiv = loaded_nbAbsEquiv;
		gHits->nbLeakTotal = loaded_nbLeak;
		gHits->distTraveled_total = distTraveled_total;

		gHits->distTraveledTotal_fullHitsOnly = distTraveledTotal_fullHitsOnly;*/

		// Read facets
		if (version >= 13) {
			file.ReadKeyword("min_pressure_all"); file.ReadKeyword(":");
			texture_limits[0].autoscale.min.steady_state = file.ReadDouble();
			file.ReadKeyword("min_pressure_moments_only"); file.ReadKeyword(":");
			texture_limits[0].autoscale.min.moments_only = file.ReadDouble();
			file.ReadKeyword("max_pressure_all"); file.ReadKeyword(":");
			texture_limits[0].autoscale.max.steady_state = file.ReadDouble();
			file.ReadKeyword("max_pressure_moments_only"); file.ReadKeyword(":");
			texture_limits[0].autoscale.max.moments_only = file.ReadDouble();

			file.ReadKeyword("min_impingement_all"); file.ReadKeyword(":");
			texture_limits[1].autoscale.min.steady_state = file.ReadDouble();
			file.ReadKeyword("min_impingement_moments_only"); file.ReadKeyword(":");
			texture_limits[1].autoscale.min.moments_only = file.ReadDouble();
			file.ReadKeyword("max_impingement_all"); file.ReadKeyword(":");
			texture_limits[1].autoscale.max.steady_state = file.ReadDouble();
			file.ReadKeyword("max_impingement_moments_only"); file.ReadKeyword(":");
			texture_limits[1].autoscale.max.moments_only = file.ReadDouble();

			file.ReadKeyword("min_density_all"); file.ReadKeyword(":");

			texture_limits[2].autoscale.min.steady_state = file.ReadDouble();
			file.ReadKeyword("min_density_moments_only"); file.ReadKeyword(":");

			texture_limits[2].autoscale.min.moments_only = file.ReadDouble();
			file.ReadKeyword("max_density_all"); file.ReadKeyword(":");

			texture_limits[2].autoscale.max.steady_state = file.ReadDouble();
			file.ReadKeyword("max_density_moments_only"); file.ReadKeyword(":");

			texture_limits[2].autoscale.max.moments_only = file.ReadDouble();

			size_t facetHitsSize = (1 + mApp->worker.interfaceMomentCache.size()) * sizeof(FacetHitBuffer);
			for (size_t m = 0; m <= mApp->worker.interfaceMomentCache.size() || (m == 0 /*&& version<10*/); m++) {
				//if (version>=10) {
				file.ReadKeyword("moment");
				if (m != file.ReadInt()) {
					throw Error("Unexpected profile moment");
					break;
				}
				file.ReadKeyword("{");
				//}

				for (int i = 0; i < sh.nbFacet; i++) {
					InterfaceFacet* f = facets[i];
					if (f->hasMesh/* || version<8*/) {
						prg.SetProgress((double)(i + m * sh.nbFacet) / (double)(mApp->worker.interfaceMomentCache.size() * sh.nbFacet) * 0.33 + 0.66);
						file.ReadKeyword("texture_facet");
						// Check idx
						int idx = file.ReadInt();

						if (idx != i + 1) {
							sprintf(tmp, "Wrong facet index. Expected %d, read %d.", i + 1, idx);
							throw Error(file.MakeError(tmp));
						}

						file.ReadKeyword("{");

						int ix, iy;

						///Load textures, for GEO file version 3+

						size_t profSize = (f->sh.isProfile) ? ((1 + mApp->worker.interfaceMomentCache.size()) * (PROFILE_SIZE * sizeof(ProfileSlice))) : 0;
						size_t h = f->sh.texHeight;
						size_t w = f->sh.texWidth;

						std::vector<TextureCell>& texture = globalState->facetStates[i].momentResults[m].texture;

						size_t texWidth_file, texHeight_file;
						//In case of rounding errors, the file might contain different texture dimensions than expected.
						if (version >= 14) {
							file.ReadKeyword("width"); file.ReadKeyword(":"); texWidth_file = file.ReadInt();
							file.ReadKeyword("height"); file.ReadKeyword(":"); texHeight_file = file.ReadInt();
						}
						else {
							texWidth_file = f->sh.texWidth;
							texHeight_file = f->sh.texHeight;
						}

						for (iy = 0; iy < (std::min(f->sh.texHeight, texHeight_file)); iy++) { //MIN: If stored texture is larger, don't read extra cells
							for (ix = 0; ix < (std::min(f->sh.texWidth, texWidth_file)); ix++) { //MIN: If stored texture is larger, don't read extra cells
								size_t index = iy * f->sh.texWidth + ix;
								texture[index].countEquiv = static_cast<double>(file.ReadSizeT());
								texture[index].sum_1_per_ort_velocity = file.ReadDouble();
								texture[index].sum_v_ort_per_area = file.ReadDouble();

							}
							for (size_t ie = 0; ie < texWidth_file - f->sh.texWidth; ie++) {//Executed if file texture is bigger than expected texture
								//Read extra cells from file without doing anything
								file.ReadSizeT();
								file.ReadDouble();
								file.ReadDouble();
							}
						}
						for (size_t ie = 0; ie < texHeight_file - f->sh.texHeight; ie++) {//Executed if file texture is bigger than expected texture
							//Read extra cells ffrom file without doing anything
							for (int iw = 0; iw < texWidth_file; iw++) {
								file.ReadSizeT();

								file.ReadDouble();
								file.ReadDouble();
							}
						}
						file.ReadKeyword("}");
					}
				}
				/*if (version>=10)*/ file.ReadKeyword("}");
			}
		}

		return true;

	}
	else
	{
		//old versions
		return false;
	}

}

/**
* \brief For saving the geometry data into GEO format
* \param file name of the output file
* \param prg GLProgress_GUI window
* \param results simulation results describing the texture
* \param worker thread worker handling this task
* \param saveSelected if a selection is to be saved
* \param crashSave if crash save is enabled
*/
void MolflowGeometry::SaveGEO(FileWriter& file, GLProgress_Abstract& prg, const std::shared_ptr<GlobalSimuState> globalState, Worker* worker,
	bool saveSelected, bool crashSave) {

	prg.SetMessage("Counting hits...");
	if (!IsLoaded()) throw Error("Nothing to save !");

	// Block dpHit during the whole disc writing
	auto lock = GetHitLock(globalState.get(), 10000);
	if (!lock) return;

	double dCoef = 1.0;
	int ix, iy;

	prg.SetMessage("Writing geometry details...");
	file.Write("version:"); file.Write(GEOVERSION, "\n");
	file.Write("totalHit:"); file.Write((!crashSave && !saveSelected) ? globalState->globalStats.globalHits.nbMCHit : 0, "\n");
	file.Write("totalDes:"); file.Write((!crashSave && !saveSelected) ? globalState->globalStats.globalHits.nbDesorbed : 0, "\n");
	file.Write("totalLeak:"); file.Write((!crashSave && !saveSelected) ? globalState->globalStats.nbLeakTotal : 0, "\n");
	file.Write("totalAbs:"); file.Write((!crashSave && !saveSelected) ? (size_t)globalState->globalStats.globalHits.nbAbsEquiv : 0, "\n");
	file.Write("totalDist_total:"); file.Write((!crashSave && !saveSelected) ? globalState->globalStats.distTraveled_total : 0, "\n");
	file.Write("totalDist_fullHitsOnly:"); file.Write((!crashSave && !saveSelected) ? globalState->globalStats.distTraveledTotal_fullHitsOnly : 0, "\n");
	file.Write("maxDes:"); file.Write((!crashSave && !saveSelected) ? worker->model->otfParams.desorptionLimit : 0, "\n");

	auto selectedFacets = GetSelectedFacets();
	file.Write("nbVertex:"); file.Write(sh.nbVertex, "\n");
	file.Write("nbFacet:"); file.Write(saveSelected ? selectedFacets.size() : sh.nbFacet, "\n");
	file.Write("nbSuper:"); file.Write(sh.nbSuper, "\n");
	file.Write("nbFormula:"); file.Write((!saveSelected) ? mApp->appFormulas->formulas.size() : 0, "\n");

	file.Write("nbView:"); file.Write(mApp->views.size(), "\n");
	file.Write("nbSelection:"); file.Write((!saveSelected) ? mApp->selections.size() : 0, "\n");

	file.Write("gasMass:"); file.Write(worker->model->sp.gasMass, "\n");

	file.Write("userMoments {\n");
	file.Write(" nb:"); file.Write((int)worker->userMoments.size());
	for (size_t u = 0; u < worker->userMoments.size(); u++) {
		file.Write("\n \"");
		file.Write(worker->userMoments[u].content.c_str());
		file.Write("\" : ");
		file.Write(worker->userMoments[u].timeWindow);
	}
	file.Write("\n}\n");

	file.Write("desorptionStart:"); file.Write(/*worker->desorptionStartTime*/0.0, "\n");
	file.Write("desorptionStop:"); file.Write(/*worker->desorptionStopTime*/1.0, "\n");
	file.Write("timeWindow:"); file.Write(worker->model->sp.timeWindowSize, "\n");
	file.Write("useMaxwellian:"); file.Write(worker->model->sp.useMaxwellDistribution, "\n");
	file.Write("calcConstantFlow:"); file.Write(worker->model->sp.calcConstantFlow, "\n");

	file.Write("formulas {\n");
	if (!saveSelected) {
		for (auto& formula : mApp->appFormulas->formulas) {
			file.Write("  \"");
			file.Write(formula.GetName());
			file.Write("\" \"");
			file.Write(formula.GetExpression());
			file.Write("\"\n");
		}
	}
	file.Write("}\n");

	file.Write("views {\n");
	for (const auto& view : mApp->views) {
		file.Write("  \"");
		file.Write(view.name);
		file.Write("\"\n");
		file.Write(view.projMode, " ");
		file.Write(view.camAngleOx, " ");
		file.Write(view.camAngleOy, " ");
		file.Write(view.camDist, " ");
		file.Write(view.camOffset.x, " ");
		file.Write(view.camOffset.y, " ");
		file.Write(view.camOffset.z, " ");
		file.Write(view.performXY, " ");
		file.Write(view.vLeft, " ");
		file.Write(view.vRight, " ");
		file.Write(view.vTop, " ");
		file.Write(view.vBottom, "\n");
	}
	file.Write("}\n");

	file.Write("selections {\n");
	for (size_t i = 0; (i < mApp->selections.size()) && !saveSelected; i++) { //don't save selections when exporting part of the geometry (saveSelected)

		file.Write("  \"");
		file.Write(mApp->selections[i].name);
		file.Write("\"\n ");
		file.Write(mApp->selections[i].facetIds.size(), "\n");
		for (auto& sel : mApp->selections[i].facetIds) {
			file.Write("  ");
			file.Write(sel, "\n");
		}
		//file.Write("\n");
	}
	file.Write("}\n");

	file.Write("structures {\n");
	for (size_t i = 0; i < sh.nbSuper; i++) {
		file.Write("  \"");
		file.Write(structNames[i]);
		file.Write("\"\n");
	}
	file.Write("}\n");
	//vertices
	prg.SetMessage("Writing vertices...");
	file.Write("vertices {\n");
	for (size_t i = 0; i < sh.nbVertex; i++) {
		prg.SetProgress(0.33 * ((double)i / (double)sh.nbVertex));
		file.Write("  ");
		file.Write(i + 1, " ");
		file.Write(vertices3[i].x, " ");
		file.Write(vertices3[i].y, " ");
		file.Write(vertices3[i].z, "\n");
	}
	file.Write("}\n");

	//leaks
	prg.SetMessage("Writing leaks...");
	file.Write("leaks {\n");
	file.Write("  nbLeak:"); file.Write((!crashSave && !saveSelected) ? worker->globalState->globalStats.leakCacheSize : 0, "\n");
	for (int i = 0; (i < worker->globalState->globalStats.leakCacheSize) && (!crashSave && !saveSelected); i++) {

		file.Write("  ");
		file.Write(i, " ");
		file.Write(worker->globalState->globalStats.leakCache[i].pos.x, " ");
		file.Write(worker->globalState->globalStats.leakCache[i].pos.y, " ");
		file.Write(worker->globalState->globalStats.leakCache[i].pos.z, " ");

		file.Write(worker->globalState->globalStats.leakCache[i].dir.x, " ");
		file.Write(worker->globalState->globalStats.leakCache[i].dir.y, " ");
		file.Write(worker->globalState->globalStats.leakCache[i].dir.z, "\n");
	}

	file.Write("}\n");

	//hit cache (lines and dots)
	prg.SetMessage("Writing hit cache...");

	file.Write("hits {\n");
	file.Write("  nbHHit:"); file.Write((!crashSave && !saveSelected) ? worker->globalState->globalStats.hitCacheSize : 0, "\n");
	for (int i = 0; (i < worker->globalState->globalStats.hitCacheSize) && (!crashSave && !saveSelected); i++) {

		file.Write("  ");
		file.Write(i, " ");
		file.Write(worker->globalState->globalStats.hitCache[i].pos.x, " ");
		file.Write(worker->globalState->globalStats.hitCache[i].pos.y, " ");
		file.Write(worker->globalState->globalStats.hitCache[i].pos.z, " ");

		file.Write(worker->globalState->globalStats.hitCache[i].type, "\n");
	}

	file.Write("}\n");

	//facets

	prg.SetMessage("Writing facets...");

	for (int i = 0, k = 0; i < sh.nbFacet; i++) {
		prg.SetProgress(0.33 + ((double)i / (double)sh.nbFacet) * 0.33);

		if (!saveSelected || facets[i]->selected) facets[i]->SaveGEO(file, k++);

	}

	prg.SetMessage("Writing profiles...");
	SaveProfileGEO(file, globalState, -1, saveSelected, crashSave);

	///Save textures, for GEO file version 3+

	char tmp[256];
	file.Write("{textures}\n");

	file.Write("min_pressure_all:"); file.Write(
		(!crashSave && !saveSelected) ? texture_limits[0].autoscale.min.steady_state : 0, "\n");
	file.Write("min_pressure_moments_only:"); file.Write(
		(!crashSave && !saveSelected) ? texture_limits[0].autoscale.min.moments_only : 0, "\n");
	file.Write("max_pressure_all:"); file.Write(
		(!crashSave && !saveSelected) ? texture_limits[0].autoscale.max.steady_state : 1, "\n");
	file.Write("max_pressure_moments_only:"); file.Write(
		(!crashSave && !saveSelected) ? texture_limits[0].autoscale.max.moments_only : 1, "\n");

	file.Write("min_impingement_all:"); file.Write(
		(!crashSave && !saveSelected) ? texture_limits[1].autoscale.min.steady_state : 0, "\n");
	file.Write("min_impingement_moments_only:"); file.Write(
		(!crashSave && !saveSelected) ? texture_limits[1].autoscale.min.moments_only : 0, "\n");
	file.Write("max_impingement_all:"); file.Write(
		(!crashSave && !saveSelected) ? texture_limits[1].autoscale.max.steady_state : 1, "\n");
	file.Write("max_impingement_moments_only:"); file.Write(
		(!crashSave && !saveSelected) ? texture_limits[1].autoscale.max.moments_only : 1, "\n");

	file.Write("min_density_all:"); file.Write(
		(!crashSave && !saveSelected) ? texture_limits[2].autoscale.min.steady_state : 0, "\n");
	file.Write("min_density_moments_only:"); file.Write(
		(!crashSave && !saveSelected) ? texture_limits[2].autoscale.min.moments_only : 0, "\n");
	file.Write("max_density_all:"); file.Write(
		(!crashSave && !saveSelected) ? texture_limits[2].autoscale.max.steady_state : 1, "\n");
	file.Write("max_density_moments_only:"); file.Write(
		(!crashSave && !saveSelected) ? texture_limits[2].autoscale.max.moments_only : 1, "\n");

	//Selections
	//SaveSelections();

	prg.SetMessage("Writing textures...");
	size_t facetHitsSize = (1 + mApp->worker.interfaceMomentCache.size()) * sizeof(FacetHitBuffer);
	for (size_t m = 0; m <= mApp->worker.interfaceMomentCache.size(); m++) {
		sprintf(tmp, "moment %zd {\n", m);
		file.Write(tmp);
		for (size_t i = 0, k = 0; i < sh.nbFacet; i++) {
			if (!saveSelected || facets[i]->selected) {
				k++; //facet id in the group of selected facets
				prg.SetProgress((double)(i + m * sh.nbFacet) / (double)(mApp->worker.interfaceMomentCache.size() * sh.nbFacet) * 0.33 + 0.66);
				InterfaceFacet* f = facets[i];
				if (f->hasMesh) {
					size_t h = f->sh.texHeight;
					size_t w = f->sh.texWidth;
					size_t profSize = (f->sh.isProfile) ? (PROFILE_SIZE * sizeof(ProfileSlice) * (1 + (int)mApp->worker.interfaceMomentCache.size())) : 0;
					const std::vector<TextureCell>& texture = globalState->facetStates[i].momentResults[m].texture;

					//char tmp[256];
					sprintf(tmp, "texture_facet %zd {\n", k); //Starts from 1, not 0, first element is k=1

					file.Write(tmp);
					file.Write("width:"); file.Write(f->sh.texWidth); file.Write(" height:"); file.Write(f->sh.texHeight); file.Write("\n");
					for (iy = 0; iy < h; iy++) {
						for (ix = 0; ix < w; ix++) {
							file.Write((!crashSave && !saveSelected) ? static_cast<size_t>(texture[iy * f->sh.texWidth + ix].countEquiv) : 0, "\t");
							file.Write((!crashSave && !saveSelected) ? texture[iy * f->sh.texWidth + ix].sum_1_per_ort_velocity : 0, "\t");
							file.Write((!crashSave && !saveSelected) ? texture[iy * f->sh.texWidth + ix].sum_v_ort_per_area : 0, "\t");
						}
						file.Write("\n");
					}

					file.Write(" }\n"); //close facet
				}
			}
		} //end facet
		file.Write("}\n");
	}
	//if (!crashSave && !saveSelected) ReleaseDataport(buffer);

}

/**
* \brief For saving the geometry data into TXT format
* \param file name of the output file
* \param results simulation results describing the texture
* \param saveSelected if a selection is to be saved
*/
void MolflowGeometry::SaveTXT(FileWriter& file, const std::shared_ptr<GlobalSimuState> globalState, bool saveSelected) {

	if (!IsLoaded()) throw Error("Nothing to save !");

	// Unused
	file.Write(0, "\n");

	// Globals

	// Unused
	file.Write(globalState->globalStats.globalHits.nbMCHit, "\n");
	file.Write(globalState->globalStats.nbLeakTotal, "\n");

	file.Write(globalState->globalStats.globalHits.nbDesorbed, "\n");
	file.Write(0, "\n"); //Desorption limit

	file.Write(sh.nbVertex, "\n");
	file.Write(saveSelected ? GetNbSelectedFacets() : sh.nbFacet, "\n");

	// Read geometry vertices
	for (int i = 0; i < sh.nbVertex; i++) {
		file.Write(vertices3[i].x, " ");
		file.Write(vertices3[i].y, " ");
		file.Write(vertices3[i].z, "\n");
	}

	// Facets
	for (int i = 0; i < sh.nbFacet; i++) {
		InterfaceFacet* f = facets[i];
		int j;
		if (saveSelected) {
			if (f->selected) {
				file.Write(f->sh.nbIndex, " ");
				for (j = 0; j < f->sh.nbIndex - 1; j++)
					file.Write(f->indices[j] + 1, " ");
				file.Write(f->indices[j] + 1, "\n");
			}
		}
		else {
			file.Write(f->sh.nbIndex, " ");
			for (j = 0; j < f->sh.nbIndex - 1; j++)
				file.Write(f->indices[j] + 1, " ");
			file.Write(f->indices[j] + 1, "\n");
		}
	}

	// Params
	for (int i = 0; i < sh.nbFacet; i++) {

		// Update facet hits from shared mem
		InterfaceFacet* f = facets[i];
		/*FacetHitBuffer *shF = (FacetHitBuffer *)(buffer + f->sp.hitOffset);
		memcpy(&(f->sp.tmpCounter), shF, sizeof(FacetHitBuffer));*/
		if (saveSelected) {
			if (f->selected) f->SaveTXT(file);
		}
		else {

			f->SaveTXT(file);
		}

	}

	SaveProfileTXT(file);
}

/**
* \brief For exporting textures depending on the texture mode
* \param file name of the output file
* \param grouping if facets should be grouped for the output
* \param mode texture mode; which type of data describes it
* \param results simulation results describing the texture
* \param saveSelected if a selection is to be saved (TODO: chefk if actually used)
* \param sMode simulation mode
*/
void
MolflowGeometry::ExportTextures(FILE* file, int grouping, int mode, const std::shared_ptr<GlobalSimuState> globalState, bool saveSelected) {

	auto lock = GetHitLock(globalState.get(), 1000);
	if (!lock) return;
	if (grouping == 1) fprintf(file, "X_coord_cm\tY_coord_cm\tZ_coord_cm\tValue\t\n"); //mode 10: special ANSYS export

	size_t facetHitsSize = (1 + mApp->worker.interfaceMomentCache.size()) * sizeof(FacetHitBuffer);
	for (size_t m = 0; m <= mApp->worker.interfaceMomentCache.size(); m++) {
		if (m == 0) fprintf(file, " moment 0 (Constant Flow){\n");
		else fprintf(file, " moment %zd (%g s)[w=%g]{\n", m, mApp->worker.interfaceMomentCache[m - 1].time, mApp->worker.interfaceMomentCache[m - 1].window);
		// Facets
		for (size_t fInd = 0; fInd < sh.nbFacet; fInd++) {
			InterfaceFacet* f = facets[fInd];

			if (f->selected) {
				if (grouping == 0) fprintf(file, "FACET%zu\n", fInd + 1u); //mode 10: special ANSYS export

				if (!f->cellPropertiesIds.empty() || f->sh.countDirection) {

					char tmp[256];
					char out[512];
					if (!globalState->initialized) return;
					size_t nbMoments = mApp->worker.interfaceMomentCache.size();
					size_t profSize = (f->sh.isProfile) ? (PROFILE_SIZE * sizeof(ProfileSlice) * (1 + nbMoments)) : 0;
					size_t w = f->sh.texWidth;
					size_t h = f->sh.texHeight;
					size_t tSize = w * h * sizeof(TextureCell);
					size_t dSize = w * h * sizeof(DirectionCell);

					const auto& facetSnapshot = globalState->facetStates[fInd].momentResults[m];
					const std::vector<TextureCell>& texture = globalState->facetStates[fInd].momentResults[m].texture;
					const std::vector<DirectionCell>& dirs = globalState->facetStates[fInd].momentResults[m].direction;

					for (size_t r = 0; r < h; r++) {
						for (size_t c = 0; c < w; c++) {
							size_t index = c + r * w;
							tmp[0] = out[0] = 0;
							switch (mode) {

							case 0: // Element area
								sprintf(tmp, "%g", f->GetMeshArea(index));
								break;

							case 1: //MC Hits

								if (!grouping || texture[index].countEquiv > 0.0) {
									double val = GetPhysicalValue(f, PhysicalMode::MCHits, 1.0, 1.0, 1.0, (int)index, facetSnapshot).value;

									sprintf(tmp, "%g", val);
								}
								break;

							case 2: //Impingement rate

								if (!grouping || texture[index].countEquiv > 0.0) {
									double moleculesPerTP = mApp->worker.GetMoleculesPerTP(m);
									double val = GetPhysicalValue(f, PhysicalMode::ImpingementRate, moleculesPerTP, 1.0, mApp->worker.model->sp.gasMass, (int)index, facetSnapshot).value;

									sprintf(tmp, "%g", val);
								}
								break;

							case 3: //Particle density
							{

								if (!grouping || texture[index].countEquiv > 0.0) {
									double moleculesPerTP = mApp->worker.GetMoleculesPerTP(m);
									double densityCorrection = f->DensityCorrection();
									double rho = GetPhysicalValue(f, PhysicalMode::ParticleDensity, moleculesPerTP, densityCorrection, mApp->worker.model->sp.gasMass, (int)index, facetSnapshot).value;

									sprintf(tmp, "%g", rho);
								}
								break;
							}
							case 4: //Gas density
							{
								if (!grouping || texture[index].countEquiv > 0.0) {
									double moleculesPerTP = mApp->worker.GetMoleculesPerTP(m);
									double densityCorrection = f->DensityCorrection();
									double rho_mass = GetPhysicalValue(f, PhysicalMode::GasDensity, moleculesPerTP, densityCorrection, mApp->worker.model->sp.gasMass, (int)index, facetSnapshot).value;

									sprintf(tmp, "%g", rho_mass);
								}
								break;
							}
							case 5:  // Pressure [mbar]

								if (!grouping || texture[index].sum_v_ort_per_area != 0.0) {
									double moleculesPerTP = mApp->worker.GetMoleculesPerTP(m);
									double p = GetPhysicalValue(f, PhysicalMode::Pressure, moleculesPerTP, 1.0, mApp->worker.model->sp.gasMass, (int)index, facetSnapshot).value;
									sprintf(tmp, "%g", p);
								}
								break;

							case 6: // Average velocity

								if (!grouping || texture[index].countEquiv > 0.0) {
									double val = GetPhysicalValue(f, PhysicalMode::AvgGasVelocity, 1.0, 1.0, 1.0, (int)index, facetSnapshot).value;
									sprintf(tmp, "%g", val);
								};
								break;

							case 7: // Velocity vector
								if (f->sh.countDirection) {
									Vector3d v_vect = GetPhysicalValue(f, PhysicalMode::GasVelocityVector, 1.0, 1.0, 1.0, (int)index, facetSnapshot).vect;
									sprintf(tmp, "%g,%g,%g",
										v_vect.x, v_vect.y, v_vect.z);
								}
								else {
									sprintf(tmp, "Direction not recorded");
								}
								break;

							case 8: // Velocity vector Count
								if (f->sh.countDirection) {
									size_t count = GetPhysicalValue(f, PhysicalMode::NbVelocityVectors, 1.0, 1.0, 1.0, (int)index, facetSnapshot).count;
									sprintf(tmp, "%zd", count);
								}
								else {

									sprintf(tmp, "None");
								}
								break;
							default:
								sprintf(tmp, "Unknown mode");
								break;
							} //end switch

							if (grouping == 1 && tmp[0]) {
								Vector2d facetCenter = f->GetMeshCenter(index);
								sprintf(out, "%g\t%g\t%g\t%s\t\n",
									f->sh.O.x + facetCenter.u * f->sh.U.x + facetCenter.v * f->sh.V.x,
									f->sh.O.y + facetCenter.u * f->sh.U.y + facetCenter.v * f->sh.V.y,
									f->sh.O.z + facetCenter.u * f->sh.U.z + facetCenter.v * f->sh.V.z,
									tmp);
							}
							else sprintf(out, "%s", tmp);

							if (out) fprintf(file, "%s", out);
							if (c < w - 1 && grouping == 0)
								fprintf(file, "\t");
						} //h
						if (grouping == 0) fprintf(file, "\n");
					} //w
				} //if mesh
				else {
					fprintf(file, "No mesh.\n");
				}
				if (grouping == 0) fprintf(file, "\n"); //Current facet exported. 
			} //if selected

		} //end facet
		fprintf(file, " }\n");

	} //end moment

}

/**
* \brief For exporting profile data (simulation)
* \param file name of the output file
* \param isTXT if TXT output will be used
* \param results simulation results describing the texture
* \param worker thread worker handling the task
*/
void MolflowGeometry::ExportProfiles(FILE* file, int isTXT, Worker* worker) {

	char sep = isTXT ? '\t' : ',';
	//if(!IsLoaded()) throw Error("Nothing to save !");

	// Block dpHit during the whole disc writing
	/*BYTE *buffer = nullptr;
	if (buffer)
		if (AccessDataport(buffer))
			buffer = (BYTE *)buffer->buff;*/

			// Globals
			//BYTE *buffer = (BYTE *)dpHit->buff;
			//SHGHITS *gHits = (SHGHITS *)buffer;
	std::ostringstream header;
	header << "Facet number" << sep << "Profile_type" << sep << "O_x" << sep << "O_y" << sep << "O_z" << sep << "U_x" << sep << "U_y" << sep << "U_z" << sep;
	header << "V_x" << sep << "V_y" << sep << "V_z" << sep << "U_length" << sep << "V_length" << sep << "Center_x" << sep << "Center_y" << sep << "Center_z" << sep << "V_max" << sep << "MC Hits" << sep << "Hit equiv." << sep;
	for (int i = 0; i < PROFILE_SIZE; i++)
		header << i + 1 << sep;
	header << '\n';

	fputs(header.str().c_str(), file);

	size_t facetHitsSize = (1 + mApp->worker.interfaceMomentCache.size()) * sizeof(FacetHitBuffer);
	for (size_t m = 0; m <= mApp->worker.interfaceMomentCache.size(); m++) {
		if (m == 0) fputs(" moment 0 (Constant Flow){\n", file);
		else fprintf(file, " moment %zd (%g s)[w=%g]{\n", m, mApp->worker.interfaceMomentCache[m - 1].time, mApp->worker.interfaceMomentCache[m - 1].window);
		// Facets

		for (int i = 0; i < sh.nbFacet; i++) {
			InterfaceFacet* f = facets[i];

			if (f->selected) {
				std::ostringstream line;

				line << i + 1 << sep << profileRecordModeDescriptions[(ProfileRecordModes)f->sh.profileType].second << sep << f->sh.O.x << sep << f->sh.O.y << sep << f->sh.O.z << sep << f->sh.U.x << sep << f->sh.U.y << sep << f->sh.U.z << sep;
				line << f->sh.V.x << sep << f->sh.V.y << sep << f->sh.V.z << sep << f->sh.U.Norme() << sep << f->sh.V.Norme() << sep << f->sh.center.x << sep << f->sh.center.y << sep << f->sh.center.z << sep << f->sh.maxSpeed << sep << f->facetHitCache.nbMCHit << sep << f->facetHitCache.nbHitEquiv << sep;

				if (f->sh.isProfile) {

					double dCoef = 1.0;
					//double nbDes = (shGHit->globalStats.hit.nbDesorbed > 0) ? (double)shGHit->globalStats.hit.nbDesorbed : 1.0;
					size_t profOffset = PROFILE_SIZE * sizeof(ProfileSlice) * m;
					const std::vector<ProfileSlice>& prof = worker->globalState->facetStates[i].momentResults[m].profile;
					double scaleX, scaleY;
					switch (f->sh.profileType) {
					case PROFILE_U:
					case PROFILE_V:
						scaleY = 1.0 / (f->GetArea() / (double)PROFILE_SIZE * 1E-4) * worker->model->sp.gasMass / 1000 / 6E23 * 0.0100; //0.01: Pa->mbar
						scaleY *= worker->GetMoleculesPerTP(m);

						for (int j = 0; j < PROFILE_SIZE; j++)
							line << prof[j].sum_v_ort * scaleY << sep;
						break;
					case PROFILE_VELOCITY:
					case PROFILE_ORT_VELOCITY:
					case PROFILE_TAN_VELOCITY:
						scaleX = f->sh.maxSpeed / (double)PROFILE_SIZE;
						for (int j = 0; j < PROFILE_SIZE; j++)
							line << prof[j].countEquiv / (f->facetHitCache.nbHitEquiv + static_cast<double>(f->facetHitCache.nbDesorbed)) << sep;
						break;
					case PROFILE_ANGULAR:
						scaleX = 90.0 / (double)PROFILE_SIZE;
						for (int j = 0; j < PROFILE_SIZE; j++)
							line << prof[j].countEquiv / (f->facetHitCache.nbHitEquiv + static_cast<double>(f->facetHitCache.nbDesorbed)) << sep;
						break;
					}
				}

				else {
					line << "No profile.";
				}
				line << '\n';
				fputs(line.str().c_str(), file);

			}

		}
		fputs(" }\n", file);

	}
}

/* Commenting out as deprecated
void MolflowGeometry::ImportDesorption_DES(FileReader& file) {

	//if(!IsLoaded()) throw Error("Nothing to save !");

	// Block dpHit during the whole disc writing

	for (int i = 0; i < sp.nbFacet; i++) { //clear previous desorption maps
		facets[i]->hasOutgassingFile = false;
		facets[i]->sp.useOutgassingFile = false;
		facets[i]->sp.desorbType = DES_NONE; //clear previously set desorptions
		facets[i]->selected = false;
		facets[i]->UnselectElem();
	}
	// Facets
	while (!file.IsEof()) {
		file.ReadKeyword("facet");
		int facetId = file.ReadInt() - 1;
		file.ReadKeyword("{");
		if (!(facetId >= 0 && facetId < sp.nbFacet)) {
			file.MakeError("Invalid facet Id (loaded desorption file for a different geometry?)");
			return;
		}

		Facet *f = facets[facetId];
		f->hasOutgassingFile = true;
		f->sp.useOutgassingFile = true; //turn on file usage by default
		f->sp.desorbType = DES_COSINE; //auto-set to cosine
		SelectFacet(facetId);
		file.ReadKeyword("cell_size_cm"); file.ReadKeyword(":");
		double ratio = f->sp.outgassingFileRatio = file.ReadDouble();
		if (f->sp.outgassingFileRatio != 0.0) {
			f->sp.outgassingFileRatio = 1.0 / f->sp.outgassingFileRatio; //cell size -> samples per cm
			ratio = f->sp.outgassingFileRatio;
		}
		double nU = f->sp.U.Norme();
		double nV = f->sp.V.Norme();
		size_t w = f->sp.outgassingMapWidth = (size_t)ceil(nU*ratio); //double precision written to file
		size_t h = f->sp.outgassingMapHeight = (size_t)ceil(nV*ratio); //double precision written to file
		f->outgassingMapWindow = (double*)malloc(w*h*sizeof(double));
		if (!f->outgassingMapWindow) throw Error("Not enough memory to store outgassing map.");
		for (int i = 0; i < w; i++) {
			for (int j = 0; j < h; j++) {
				f->outgassingMapWindow[i + j*w] = file.ReadDouble();
			}
		}
		file.ReadKeyword("}");
	}
	UpdateSelection();
	//InitializeGeometry();

	//Debug memory check
	//_ASSERTE (!_CrtDumpMemoryLeaks());;
	_ASSERTE(_CrtCheckMemory());

}
*/

/**
* \brief For importing desorption data from a SYN file
* \param file name of the input file
* \param source what the source to calculate the dose is
* \param time time to calculate the dose
* \param mode mode used for outgassing calculation
* \param eta0 coefficient for outgassing calculation in mode==1
* \param alpha exponent for outgassing calculation in mode==1
* \param cutoffdose cutoff dose for outgassing calculation in mode==1
* \param convDistr distribution for outgassing calculation in mode==2
* \param prg GLProgress_GUI window where visualising of the import progress is shown
*/
void MolflowGeometry::ImportDesorption_SYN(
	FileReader& file, const size_t source, const double time,
	const size_t mode, const double eta0, const double alpha, const double cutoffdose,
	const std::vector<std::pair<double, double>>& convDistr,
	GLProgress_Abstract& prg) {

	//UnselectAll();
	char tmp[512];
	std::vector<double> xdims, ydims;
	double no_scans = 1.0;

	file.ReadKeyword("version"); file.ReadKeyword(":");
	int version = file.ReadInt();
	if (version > SYNVERSION) {
		char errMsg[512];
		sprintf(errMsg, "Unsupported SYN version V%d", version);
		throw Error(errMsg);

	}

	if (version >= 12) {
		file.ReadKeyword("newReflectionModel");
		file.ReadKeyword(":");
		/*worker->sp.newReflectionModel =*/ file.ReadInt();
		file.ReadKeyword("lowFluxMode");
		file.ReadKeyword(":");
		/*worker->ontheflyParams.lowFluxMode =*/ file.ReadInt();
		file.ReadKeyword("lowFluxCutoff");
		file.ReadKeyword(":");
		/*worker->ontheflyParams.lowFluxCutoff =*/ file.ReadDouble();
	}

	//now read number of facets
	file.ReadKeyword("totalHit"); file.ReadKeyword(":");
	file.ReadSizeT();
	if (version >= 10) {
		file.ReadKeyword("totalHitEquiv"); file.ReadKeyword(":");
		file.ReadDouble();
	}
	file.ReadKeyword("totalDes"); file.ReadKeyword(":");
	file.ReadSizeT();
	if (version >= 6) {
		file.ReadKeyword("no_scans"); file.ReadKeyword(":");
		no_scans = file.ReadDouble();
	}
	file.ReadKeyword("totalLeak"); file.ReadKeyword(":");
	file.ReadSizeT();
	if (version > 2) {
		file.ReadKeyword("totalFlux"); file.ReadKeyword(":");
		file.ReadDouble();
		file.ReadKeyword("totalPower"); file.ReadKeyword(":");
		file.ReadDouble();

	}
	if (version >= 10) {
		file.ReadKeyword("totalAbsEquiv"); file.ReadKeyword(":");
		file.ReadDouble();
		file.ReadKeyword("totalDist"); file.ReadKeyword(":");
		file.ReadDouble();
	}
	file.ReadKeyword("maxDes"); file.ReadKeyword(":");
	file.ReadSizeT();
	file.ReadKeyword("nbVertex"); file.ReadKeyword(":");
	file.ReadInt();
	file.ReadKeyword("nbFacet"); file.ReadKeyword(":");
	size_t nbNewFacet = file.ReadInt(); //gotcha! :)
	xdims.reserve(nbNewFacet);
	ydims.reserve(nbNewFacet);

	//now go for the facets to get their texture ratio
	for (size_t i = 0; i < std::min(nbNewFacet, GetNbFacet()); i++) {
		prg.SetProgress(0.5 * (double)i / (double)std::min(nbNewFacet, GetNbFacet()));
		file.JumpSection("facet");
		// Check idx
		int idx = file.ReadInt();
		if (idx != i + 1) throw Error(file.MakeError("Wrong facet index !"));

		file.JumpSection("texDimX"); file.ReadKeyword(":");
		xdims.push_back(file.ReadDouble());
		//if (!IsZero(xdims[i])) SelectFacet(GetFacet(i));
		file.ReadKeyword("texDimY"); file.ReadKeyword(":");
		ydims.push_back(file.ReadDouble());

	}

	//now read actual textures
	//read header
	file.SeekFor("{textures}");
	file.ReadKeyword("minHit_MC"); file.ReadKeyword(":");
	file.ReadSizeT();
	file.ReadKeyword("maxHit_MC"); file.ReadKeyword(":");
	file.ReadSizeT();
	file.ReadKeyword("minHit_flux"); file.ReadKeyword(":");
	file.ReadDouble();
	file.ReadKeyword("maxHit_flux"); file.ReadKeyword(":");
	file.ReadDouble();
	file.ReadKeyword("minHit_power"); file.ReadKeyword(":");
	file.ReadDouble();
	file.ReadKeyword("maxHit_power"); file.ReadKeyword(":");
	file.ReadDouble();

	//read texture values
	for (size_t i = 0; i < std::min(nbNewFacet, GetNbFacet()); i++) {
		prg.SetProgress(0.5 + 0.5 * (double)i / (double)std::min(nbNewFacet, GetNbFacet()));
		if (!IsZero(xdims[i])) { //has texture
			InterfaceFacet* f = GetFacet(i);

			file.ReadKeyword("texture_facet");
			// Check idx
			int idx = file.ReadInt();

			if (idx != i + 1) {
				sprintf(tmp, "Wrong facet index. Expected %zd, read %d.", i + 1, idx);
				throw Error(file.MakeError(tmp));
			}

			//Now load values
			file.ReadKeyword("{");

			size_t ix, iy;

			if (f->selected) {
				f->hasOutgassingFile = true;
				f->sh.useOutgassingFile = true; //turn on file usage by default
				f->sh.desorbType = DES_COSINE; //auto-set to cosine
				f->ogMap.outgassingMapWidth = (size_t)ceil(xdims[i] * 0.9999999);
				f->ogMap.outgassingMapHeight = (size_t)ceil(ydims[i] * 0.9999999);
				f->ogMap.outgassingFileRatioU = xdims[i] / f->sh.U.Norme();
				f->ogMap.outgassingFileRatioV = ydims[i] / f->sh.V.Norme();
				try {
					std::vector<double>(f->ogMap.outgassingMapWidth * f->ogMap.outgassingMapHeight).swap(f->ogMap.outgassingMap);
				}
				catch (...) {
					throw Error("Not enough memory to store outgassing map.");
				}
				f->ogMap.totalDose = f->sh.totalOutgassing = f->ogMap.totalFlux = 0.0;
			}

			size_t texWidth_file, texHeight_file;
			//In case of rounding errors, the file might contain different texture dimensions than expected.
			if (version >= 8) {
				file.ReadKeyword("width"); file.ReadKeyword(":"); texWidth_file = file.ReadInt();
				file.ReadKeyword("height"); file.ReadKeyword(":"); texHeight_file = file.ReadInt();
			}
			else {
				texWidth_file = f->ogMap.outgassingMapWidth;
				texHeight_file = f->ogMap.outgassingMapHeight;
			}

			for (iy = 0; iy < (std::min(f->ogMap.outgassingMapHeight, texHeight_file)); iy++) { //MIN: If stored texture is larger, don't read extra cells
				for (ix = 0; ix < (std::min(f->ogMap.outgassingMapWidth, texWidth_file)); ix++) { //MIN: If stored texture is larger, don't read extra cells
					size_t index = iy * f->ogMap.outgassingMapWidth + ix;
					//Read original values
					size_t MC = file.ReadSizeT();
					double cellArea = 1.0;
					if (version >= 7) cellArea = file.ReadDouble();
					if (cellArea < 1E-10) cellArea = 1.0; //to avoid division by zero
					double flux = file.ReadDouble() / no_scans; //not normalized by cell area
					double power = file.ReadDouble() / no_scans; //not normalized by cell area

					if (f->selected) {
						//Calculate dose
						double dose;
						if (source == 0) dose = (double)MC * time;
						else if (source == 1) dose = flux * time / cellArea;
						else if (source == 2) dose = power * time / cellArea;

						double outgassing;
						if (dose == 0) outgassing = 0; //to avoid division by zero later
						else {
							//Convert to outgassing
							if (mode == 0) {
								if (source == 0) outgassing = (double)MC * MBARLS_TO_PAM3S / 1.38E-23 / f->sh.temperature;
								else if (source == 1) outgassing = flux * MBARLS_TO_PAM3S / 1.38E-23 / f->sh.temperature; //Division by 10 because the user will want to see the same outgassing in mbar*l/s
								else if (source == 2) outgassing = power * MBARLS_TO_PAM3S / 1.38E-23 / f->sh.temperature; //(Outgassing is stored internally in Pa*m3/s, for consistent SI unit calculations)
							}
							else if (mode == 1) {
								double moleculePerPhoton = eta0 * pow(std::max(1.0, dose / cutoffdose), alpha);
								outgassing = flux * moleculePerPhoton;
							}
							else if (mode == 2) {
								double moleculePerPhoton = InterpolateY(dose, convDistr, true, true, true);
								outgassing = flux * moleculePerPhoton;
							}
						}
						//Apply outgassing
						f->ogMap.outgassingMap[index] = outgassing * 1.38E-23 * f->sh.temperature; //1[Pa*m3/s] = kT [particles/sec]

						//Facet diagnostic info
						f->ogMap.totalDose += flux * time;
						f->ogMap.totalFlux += flux;
						f->sh.totalOutgassing += f->ogMap.outgassingMap[index];

					} //if selected
				}
				for (size_t ie = 0; ie < texWidth_file - f->ogMap.outgassingMapWidth; ie++) {//Executed if file texture is bigger than expected texture
					//Read extra cells from file without doing anything
					//Read original values
					file.ReadSizeT(); //MC
					if (version >= 7) file.ReadDouble(); //area
					file.ReadDouble(); //flux
					file.ReadDouble(); //power
				}
			}
			for (size_t ie = 0; ie < texHeight_file - f->ogMap.outgassingMapHeight; ie++) {//Executed if file texture is bigger than expected texture
				//Read extra cells ffrom file without doing anything
				for (size_t iw = 0; iw < texWidth_file; iw++) {
					//Read original values
					file.ReadSizeT(); //MC
					if (version >= 7) file.ReadDouble(); //area
					file.ReadDouble(); //flux
					file.ReadDouble(); //power
				}
			}
			file.ReadKeyword("}");
		} //if has texture
	}
	//end
	//UpdateSelection();
}

/**
* \brief To analyze desorption data from a SYN file
* \param file name of the input file
* \param prg GLProgress_GUI dialog (TODO: but is it ever used?)
* \param nbNewFacet number of facets in the file
* \param nbTextured number of textured facets in the file
* \param nbDifferent number that is only set to 0 but never used (TODO: check usage)
* \param prg GLProgress_GUI window where visualising of the analysation progress is shown
*/
void MolflowGeometry::AnalyzeSYNfile(FileReader& file, GLProgress_Abstract& prg, size_t* nbNewFacet,
	size_t* nbTextured, size_t* nbDifferent) {
	//init
	*nbTextured = 0;
	*nbNewFacet = 0;
	*nbDifferent = 0;

	UnselectAll();
	//char tmp[512];

	file.ReadKeyword("version"); file.ReadKeyword(":");
	int version;
	version = file.ReadInt();
	if (version > SYNVERSION) {
		char errMsg[512];
		sprintf(errMsg, "Unsupported SYN version V%d", version);
		throw Error(errMsg);

	}

	if (version >= 12) {
		file.ReadKeyword("newReflectionModel");
		file.ReadKeyword(":");
		/*worker->sp.newReflectionModel =*/ file.ReadInt();
		file.ReadKeyword("lowFluxMode");
		file.ReadKeyword(":");
		/*worker->ontheflyParams.lowFluxMode =*/ file.ReadInt();
		file.ReadKeyword("lowFluxCutoff");
		file.ReadKeyword(":");
		/*worker->ontheflyParams.lowFluxCutoff =*/ file.ReadDouble();
	}

	//now read number of facets
	file.ReadKeyword("totalHit"); file.ReadKeyword(":");
	file.ReadSizeT();
	if (version >= 10) {
		file.ReadKeyword("totalHitEquiv"); file.ReadKeyword(":");
		file.ReadDouble();
	}
	file.ReadKeyword("totalDes"); file.ReadKeyword(":");
	file.ReadSizeT();
	if (version >= 6) {
		file.ReadKeyword("no_scans"); file.ReadKeyword(":");
		/*no_scans = */file.ReadDouble();

	}
	file.ReadKeyword("totalLeak"); file.ReadKeyword(":");
	file.ReadSizeT();
	if (version > 2) {
		file.ReadKeyword("totalFlux"); file.ReadKeyword(":");
		file.ReadDouble();
		file.ReadKeyword("totalPower"); file.ReadKeyword(":");
		file.ReadDouble();
	}
	if (version >= 10) {
		file.ReadKeyword("totalAbsEquiv"); file.ReadKeyword(":");
		file.ReadDouble();
		file.ReadKeyword("totalDist"); file.ReadKeyword(":");
		file.ReadDouble();
	}
	file.ReadKeyword("maxDes"); file.ReadKeyword(":");
	file.ReadSizeT();
	file.ReadKeyword("nbVertex"); file.ReadKeyword(":");
	file.ReadInt();
	file.ReadKeyword("nbFacet"); file.ReadKeyword(":");
	*nbNewFacet = file.ReadInt(); //gotcha! :)

	//now go for the facets to get their texture ratio, etc.
	for (size_t i = 0; i < *nbNewFacet && i < GetNbFacet(); i++) {
		prg.SetProgress((double)i / (double)std::min(*nbNewFacet, GetNbFacet()));
		file.JumpSection("facet");
		// Check idx
		int idx = file.ReadInt();
		if (idx != i + 1) throw Error(file.MakeError("Wrong facet index !"));

		file.JumpSection("mesh"); file.ReadKeyword(":");
		if (file.ReadInt()) { //has mesh
			(*nbTextured)++;
			SelectFacet(i);
			/*file.ReadKeyword("texDimX");file.ReadKeyword(":");
			if ((this->GetFacet(i)->sp.texWidth_precise-file.ReadDouble())>1E-8) {
			(*nbDifferent)++;
			continue;

			}
			file.ReadKeyword("texDimY");file.ReadKeyword(":");
			if ((this->GetFacet(i)->sp.texHeight_precise-file.ReadDouble())>1E-8) {
			(*nbDifferent)++;
			}*/

		}
	}
	UpdateSelection();

}

/**
* \brief To save geometry data into a XML file
* \param saveDoc xml output file
* \param work thread worker handling the task
* \param prg GLProgress_GUI window where visualising of the export progress is shown
* \param saveSelected saveSelected if a selection is to be saved
*/
/*
void MolflowGeometry::SaveXML_geometry(xml_node& saveDoc, Worker* work, GLProgress_Abstract& prg, bool saveSelected) { //scheduled to be removed
	//TiXmlDeclaration* decl = new TiXmlDeclaration("1.0")="")="");
	//saveDoc->LinkEndChild(decl);

	xml_node rootNode;
	if (mApp->useOldXMLFormat) {
		rootNode = saveDoc;
	}
	else {
		rootNode = saveDoc.append_child("SimulationEnvironment");
		rootNode.attribute("type") = "molflow";
		rootNode.append_attribute("version") = appVersionId;
	}
	xml_node geomNode = rootNode.append_child("Geometry");

	prg.SetMessage("Writing vertices...");
	geomNode.append_child("Vertices").append_attribute("nb") = sh.nbVertex; //creates Vertices node, adds nb attribute and sets its value to sp.nbVertex
	for (int i = 0; i < sh.nbVertex; i++) {
		prg.SetProgress(0.166 * ((double)i / (double)sh.nbVertex));
		xml_node v = geomNode.child("Vertices").append_child("Vertex");
		v.append_attribute("id") = i;
		v.append_attribute("x") = vertices3[i].x;
		v.append_attribute("y") = vertices3[i].y;
		v.append_attribute("z") = vertices3[i].z;
	}

	prg.SetMessage("Writing facets...");
	geomNode.append_child("Facets");
	geomNode.child("Facets").append_attribute("nb") = sh.nbFacet;
	for (int i = 0, k = 0; i < sh.nbFacet; i++) {
		prg.SetProgress(0.166 + ((double)i / (double)sh.nbFacet) * 0.166);
		if (!saveSelected || facets[i]->selected) {
			xml_node f = geomNode.child("Facets").append_child("Facet");
			f.append_attribute("id") = k++;
			facets[i]->SaveXML_geom(f);
		}
	}

	prg.SetMessage("Writing model details...");
	geomNode.append_child("Structures").append_attribute("nb") = sh.nbSuper;

	for (int i = 0, k = 0; i < sh.nbSuper; i++) {
		xml_node s = geomNode.child("Structures").append_child("Structure");
		s.append_attribute("id") = i;
		s.append_attribute("name") = (strName) ? strName[i] : "";

	}
	xml_node interfNode = rootNode.append_child("Interface");

	xml_node selNode = interfNode.append_child("Selections");
	selNode.append_attribute("nb") = (!saveSelected) * (mApp->selections.size());
	for (size_t i = 0; (i < mApp->selections.size()) && !saveSelected; i++) { //don't save selections when exporting part of the geometry (saveSelected)
		xml_node newSel = selNode.append_child("Selection");
		newSel.append_attribute("id") = i;
		newSel.append_attribute("name") = mApp->selections[i].name.c_str();
		newSel.append_attribute("nb") = mApp->selections[i].selection.size();
		for (int j = 0; j < mApp->selections[i].selection.size(); j++) {
			xml_node newItem = newSel.append_child("selItem");
			newItem.append_attribute("id") = j;
			newItem.append_attribute("facet") = mApp->selections[i].selection[j];
		}
	}

	xml_node viewNode = interfNode.append_child("Views");
	viewNode.append_attribute("nb") = (!saveSelected) * (mApp->nbView);
	for (int i = 0; (i < mApp->nbView) && !saveSelected; i++) { //don't save views when exporting part of the geometry (saveSelected)
		xml_node newView = viewNode.append_child("View");
		newView.append_attribute("id") = i;
		newView.append_attribute("name") = mApp->views[i].name.c_str();
		newView.append_attribute("projMode") = mApp->views[i].projMode;
		newView.append_attribute("camAngleOx") = mApp->views[i].camAngleOx;
		newView.append_attribute("camAngleOy") = mApp->views[i].camAngleOy;
		newView.append_attribute("camAngleOz") = mApp->views[i].camAngleOz;
		newView.append_attribute("camDist") = mApp->views[i].camDist;
		newView.append_attribute("lightAngleOx") = mApp->views[i].lightAngleOx;
		newView.append_attribute("lightAngleOy") = mApp->views[i].lightAngleOy;
		newView.append_attribute("camOffset.x") = mApp->views[i].camOffset.x;
		newView.append_attribute("camOffset.y") = mApp->views[i].camOffset.y;
		newView.append_attribute("camOffset.z") = mApp->views[i].camOffset.z;
		newView.append_attribute("performXY") = mApp->views[i].performXY;
		newView.append_attribute("vLeft") = mApp->views[i].vLeft;
		newView.append_attribute("vRight") = mApp->views[i].vRight;
		newView.append_attribute("vTop") = mApp->views[i].vTop;
		newView.append_attribute("vBottom") = mApp->views[i].vBottom;
	}

	xml_node formulaNode = interfNode.append_child("Formulas");
	formulaNode.append_attribute("nb") = (!saveSelected) * (mApp->appFormulas->formulas.size());
	if (!saveSelected) { //don't save formulas when exporting part of the geometry (saveSelected)
		for (size_t i = 0; i < mApp->appFormulas->formulas.size(); i++) {
			xml_node newFormula = formulaNode.append_child("Formula");
			newFormula.append_attribute("id") = i;
			newFormula.append_attribute("name") = mApp->appFormulas->formulas[i].GetName();
			newFormula.append_attribute("expression") = mApp->appFormulas->formulas[i].GetExpression();
		}
	}

	if (mApp->profilePlotter) {
		std::vector<int> ppViews = mApp->profilePlotter->GetViews();
		xml_node profilePlotterNode = interfNode.append_child("ProfilePlotter");
		profilePlotterNode.append_child("Parameters").append_attribute("logScale") = (int)mApp->profilePlotter->IsLogScaled(); //backward compatibility: 0 or 1
		xml_node viewsNode = profilePlotterNode.append_child("Views");
		for (int v : ppViews) {
			xml_node view = viewsNode.append_child("View");
			view.append_attribute("facetId") = v;
		}
	}

	if (mApp->convergencePlotter) {
		std::vector<int> cpViews = mApp->convergencePlotter->GetViews();
		xml_node convergencePlotterNode = interfNode.append_child("ConvergencePlotter");
		convergencePlotterNode.append_child("Parameters").append_attribute("logScale") = (int)mApp->convergencePlotter->IsLogScaled(); //backward compatibility: 0 or 1
		xml_node viewsNode = convergencePlotterNode.append_child("Views");
		for (int v : cpViews) {
			xml_node view = viewsNode.append_child("View");
			view.append_attribute("formulaHash") = v;
		}
	}

	xml_node simuParamNode = rootNode.append_child("MolflowSimuSettings");

	simuParamNode.append_child("Gas").append_attribute("mass") = work->model->sp.gasMass;
	simuParamNode.child("Gas").append_attribute("enableDecay") = (int)work->model->sp.enableDecay; //backward compatibility: 0 or 1
	simuParamNode.child("Gas").append_attribute("halfLife") = work->model->sp.halfLife;

	xml_node timeSettingsNode = simuParamNode.append_child("TimeSettings");

	xml_node userMomentsNode = timeSettingsNode.append_child("UserMoments");
	userMomentsNode.append_attribute("nb") = work->userMoments.size();
	for (size_t i = 0; i < work->userMoments.size(); i++) {
		xml_node newUserEntry = userMomentsNode.append_child("UserEntry");
		newUserEntry.append_attribute("id") = i;
		newUserEntry.append_attribute("content") = work->userMoments[i].first.c_str();
		newUserEntry.append_attribute("window") = work->userMoments[i].second;
	}

	timeSettingsNode.append_attribute("timeWindow") = work->model->sp.timeWindowSize;
	timeSettingsNode.append_attribute("useMaxwellDistr") = (int)work->model->sp.useMaxwellDistribution; //backward compatibility: 0 or 1
	timeSettingsNode.append_attribute("calcConstFlow") = (int)work->model->sp.calcConstantFlow; //backward compatibility: 0 or 1

	xml_node motionNode = simuParamNode.append_child("Motion");
	motionNode.append_attribute("type") = work->model->sp.motionType;
	if (work->model->sp.motionType == 1) { //fixed motion
		xml_node v = motionNode.append_child("VelocityVector");
		v.append_attribute("vx") = work->model->sp.motionVector2.x;
		v.append_attribute("vy") = work->model->sp.motionVector2.y;
		v.append_attribute("vz") = work->model->sp.motionVector2.z;
	}
	else if (work->model->sp.motionType == 2) { //rotation
		xml_node v = motionNode.append_child("AxisBasePoint");
		v.append_attribute("x") = work->model->sp.motionVector1.x;
		v.append_attribute("y") = work->model->sp.motionVector1.y;
		v.append_attribute("z") = work->model->sp.motionVector1.z;
		xml_node v2 = motionNode.append_child("RotationVector");
		v2.append_attribute("x") = work->model->sp.motionVector2.x;
		v2.append_attribute("y") = work->model->sp.motionVector2.y;
		v2.append_attribute("z") = work->model->sp.motionVector2.z;
	}

	auto forcesNode = simuParamNode.append_child("MeasureForces");
	forcesNode.append_attribute("enabled") = work->model->sp.enableForceMeasurement;
	auto torqueNode = forcesNode.append_child("Torque");
	auto v = torqueNode.append_child("refPoint");
	v.append_attribute("x") = work->model->sp.torqueRefPoint.x;
	v.append_attribute("y") = work->model->sp.torqueRefPoint.y;
	v.append_attribute("z") = work->model->sp.torqueRefPoint.z;

	xml_node paramNode = simuParamNode.append_child("Parameters");
	size_t nonCatalogParameters = 0;

	for (size_t i = 0; i < work->parameters.size(); i++) {
		if (!work->parameters[i].fromCatalog) { //Don't save catalog parameters
			xml_node newParameter = paramNode.append_child("Parameter");
			newParameter.append_attribute("id") = nonCatalogParameters;
			newParameter.append_attribute("name") = work->parameters[i].name.c_str();
			newParameter.append_attribute("nbMoments") = (int)work->parameters[i].GetSize();
			newParameter.append_attribute("logXinterp") = work->parameters[i].logXinterp;
			newParameter.append_attribute("logYinterp") = work->parameters[i].logYinterp;
			for (size_t m = 0; m < work->parameters[i].GetSize(); m++) {
				xml_node newMoment = newParameter.append_child("Moment");
				newMoment.append_attribute("id") = m;
				newMoment.append_attribute("t") = work->parameters[i].GetX(m);
				newMoment.append_attribute("value") = work->parameters[i].GetY(m);
			}
			nonCatalogParameters++;
		}
	}
	paramNode.append_attribute("nb") = nonCatalogParameters;
	xml_node globalHistNode = simuParamNode.append_child("Global_histograms");
	if (work->model->sp.globalHistogramParams.recordBounce) {
		xml_node nbBounceNode = globalHistNode.append_child("Bounces");
		nbBounceNode.append_attribute("binSize") = work->model->sp.globalHistogramParams.nbBounceBinsize;
		nbBounceNode.append_attribute("max") = work->model->sp.globalHistogramParams.nbBounceMax;
	}
	if (work->model->sp.globalHistogramParams.recordDistance) {
		xml_node distanceNode = globalHistNode.append_child("Distance");
		distanceNode.append_attribute("binSize") = work->model->sp.globalHistogramParams.distanceBinsize;
		distanceNode.append_attribute("max") = work->model->sp.globalHistogramParams.distanceMax;
	}
#ifdef MOLFLOW
	if (work->model->sp.globalHistogramParams.recordTime) {
		xml_node timeNode = globalHistNode.append_child("Time");
		timeNode.append_attribute("binSize") = work->model->sp.globalHistogramParams.timeBinsize;
		timeNode.append_attribute("max") = work->model->sp.globalHistogramParams.timeMax;
	}
#endif
}
*/

/**
* \brief To save simulation data into a XML file
* \param saveDoc xml output file
* \param work thread worker handling the task
* \param results simulation results
* \param prg GLProgress_GUI window where visualising of the export progress is shown
* \param saveSelected saveSelected if a selection is to be saved (TODO: check if necessary)
* \return bool if saving is successfull (always is here)
*/
/*
bool MolflowGeometry::SaveXML_simustate(xml_node saveDoc, Worker* work, const std::shared_ptr<GlobalSimuState> globalState, GLProgress_Abstract& prg, bool saveSelected) { //scheduled to be removed
	xml_node rootNode;
	//Lock during writing
	std::unique_lock<std::timed_mutex> lock(globalState->particleLogMutex, std::chrono::seconds(10));
	if (!lock.owns_lock()) {
		return false;
	}
	if (mApp->useOldXMLFormat) {
		rootNode = saveDoc;
	}
	else {
		rootNode = saveDoc.child("SimulationEnvironment");
	}
	xml_node resultNode = rootNode.append_child("MolflowResults");
	prg.SetMessage("Saving simulation results...");
	xml_node momentsNode = resultNode.append_child("Moments");
	momentsNode.append_attribute("nb") = work->interfaceMomentCache.size() + 1;
	size_t facetHitsSize = (1 + mApp->worker.interfaceMomentCache.size()) * sizeof(FacetHitBuffer);
	for (size_t m = 0; m <= mApp->worker.interfaceMomentCache.size(); m++) {

		std::ostringstream msg;
		msg << "Saving moment " << (m + 1) << "/" << (mApp->worker.interfaceMomentCache.size() + 1) << "...";
		prg.SetMessage(msg.str(), false);
		prg.SetProgress(0.5 + 0.5 * (double)m / (1.0 + (double)mApp->worker.interfaceMomentCache.size()));

		xml_node newMoment = momentsNode.append_child("Moment");
		newMoment.append_attribute("id") = m;
		if (m == 0) {
			newMoment.append_attribute("time") = "Constant flow";
			newMoment.append_attribute("timeWindow") = 0;
		}
		else {
			newMoment.append_attribute("time") = work->interfaceMomentCache[m - 1].first;
			newMoment.append_attribute("timeWindow") = work->interfaceMomentCache[m - 1].second;
		}

		if (m == 0) { //Write global results. Later these results will probably be time-dependent as well.
			xml_node globalNode = newMoment.append_child("Global");

			xml_node hitsNode = globalNode.append_child("Hits");
			hitsNode.append_attribute("totalHit") = globalState->globalStats.globalHits.nbMCHit;
			hitsNode.append_attribute("totalHitEquiv") = globalState->globalStats.globalHits.nbHitEquiv;
			hitsNode.append_attribute("totalDes") = globalState->globalStats.globalHits.nbDesorbed;
			hitsNode.append_attribute("totalAbsEquiv") = globalState->globalStats.globalHits.nbAbsEquiv;
			hitsNode.append_attribute("totalDist_total") = globalState->globalStats.distTraveled_total;
			hitsNode.append_attribute("totalDist_fullHitsOnly") = globalState->globalStats.distTraveledTotal_fullHitsOnly;
			hitsNode.append_attribute("totalLeak") = globalState->globalStats.nbLeakTotal;
			hitsNode.append_attribute("maxDesorption") = work->model->otfParams.desorptionLimit;

			xml_node hitCacheNode = globalNode.append_child("Hit_Cache");
			hitCacheNode.append_attribute("nb") = work->globalState->globalStats.hitCacheSize;

			for (int i = 0; i < work->globalState->globalStats.hitCacheSize; i++) {
				xml_node newHit = hitCacheNode.append_child("Hit");
				newHit.append_attribute("id") = i;
				newHit.append_attribute("posX") = work->globalState->globalStats.hitCache[i].pos.x;
				newHit.append_attribute("posY") = work->globalState->globalStats.hitCache[i].pos.y;
				newHit.append_attribute("posZ") = work->globalState->globalStats.hitCache[i].pos.z;
				newHit.append_attribute("type") = work->globalState->globalStats.hitCache[i].type;
			}

			xml_node leakCacheNode = globalNode.append_child("Leak_Cache");
			leakCacheNode.append_attribute("nb") = work->globalState->globalStats.leakCacheSize;
			for (int i = 0; i < work->globalState->globalStats.leakCacheSize; i++) {
				xml_node newLeak = leakCacheNode.append_child("Leak");
				newLeak.append_attribute("id") = i;
				newLeak.append_attribute("posX") = work->globalState->globalStats.leakCache[i].pos.x;
				newLeak.append_attribute("posY") = work->globalState->globalStats.leakCache[i].pos.y;
				newLeak.append_attribute("posZ") = work->globalState->globalStats.leakCache[i].pos.z;
				newLeak.append_attribute("dirX") = work->globalState->globalStats.leakCache[i].dir.x;
				newLeak.append_attribute("dirY") = work->globalState->globalStats.leakCache[i].dir.y;
				newLeak.append_attribute("dirZ") = work->globalState->globalStats.leakCache[i].dir.z;
			}
		} //end global node


		bool hasHistogram = work->model->sp.globalHistogramParams.recordBounce || work->model->sp.globalHistogramParams.recordDistance;
#ifdef MOLFLOW
		hasHistogram = hasHistogram || work->model->sp.globalHistogramParams.recordTime;
#endif										
		if (hasHistogram) {
			xml_node histNode = newMoment.append_child("Histograms");
			//Retrieve histogram map from hits dp
			auto& globalHist = work->globalState->globalHistograms[m];
			if (work->model->sp.globalHistogramParams.recordBounce) {
				auto& nbHitsHistogram = globalHist.nbHitsHistogram;
				xml_node hist = histNode.append_child("Bounces");
				size_t histSize = work->model->sp.globalHistogramParams.GetBounceHistogramSize();
				hist.append_attribute("size") = histSize;
				hist.append_attribute("binSize") = work->model->sp.globalHistogramParams.nbBounceBinsize; //redundancy for human-reading or export
				hist.append_attribute("max") = work->model->sp.globalHistogramParams.nbBounceMax; //redundancy for human-reading or export
				for (size_t h = 0; h < histSize; h++) {
					xml_node bin = hist.append_child("Bin");
					auto value = bin.append_attribute("start");
					if (h == histSize - 1) value = "overRange";
					else value = h * work->model->sp.globalHistogramParams.nbBounceBinsize;
					bin.append_attribute("count") = nbHitsHistogram[h];
				}
			}
			if (work->model->sp.globalHistogramParams.recordDistance) {
				auto& distanceHistogram = globalHist.distanceHistogram;
				xml_node hist = histNode.append_child("Distance");
				size_t histSize = work->model->sp.globalHistogramParams.GetDistanceHistogramSize();
				hist.append_attribute("size") = histSize;
				hist.append_attribute("binSize") = work->model->sp.globalHistogramParams.distanceBinsize; //redundancy for human-reading or export
				hist.append_attribute("max") = work->model->sp.globalHistogramParams.distanceMax; //redundancy for human-reading or export
				for (size_t h = 0; h < histSize; h++) {
					xml_node bin = hist.append_child("Bin");
					auto value = bin.append_attribute("start");
					if (h == histSize - 1) value = "overRange";
					else value = h * work->model->sp.globalHistogramParams.distanceBinsize;
					bin.append_attribute("count") = distanceHistogram[h];
				}
			}
			if (work->model->sp.globalHistogramParams.recordTime) {
				auto& timeHistogram = globalHist.timeHistogram;
				xml_node hist = histNode.append_child("Time");
				size_t histSize = work->model->sp.globalHistogramParams.GetTimeHistogramSize();
				hist.append_attribute("size") = histSize;
				hist.append_attribute("binSize") = work->model->sp.globalHistogramParams.timeBinsize; //redundancy for human-reading or export
				hist.append_attribute("max") = work->model->sp.globalHistogramParams.timeMax; //redundancy for human-reading or export
				for (size_t h = 0; h < histSize; h++) {
					xml_node bin = hist.append_child("Bin");
					auto value = bin.append_attribute("start");
					if (h == histSize - 1) value = "overRange";
					else value = h * work->model->sp.globalHistogramParams.timeBinsize;
					bin.append_attribute("count") = timeHistogram[h];
				}
			}
		}

		xml_node facetResultsNode = newMoment.append_child("FacetResults");

		for (size_t i = 0; i < sh.nbFacet; i++) {
			InterfaceFacet* f = GetFacet(i);
			xml_node newFacetResult = facetResultsNode.append_child("Facet");
			newFacetResult.append_attribute("id") = i;

			xml_node facetHitNode = newFacetResult.append_child("Hits");
			const FacetHitBuffer& facetCounter = globalState->facetStates[i].momentResults[m].hits;
			facetHitNode.append_attribute("nbHit") = facetCounter.nbMCHit;
			facetHitNode.append_attribute("nbHitEquiv") = facetCounter.nbHitEquiv;
			facetHitNode.append_attribute("nbDes") = facetCounter.nbDesorbed;
			facetHitNode.append_attribute("nbAbsEquiv") = facetCounter.nbAbsEquiv;
			facetHitNode.append_attribute("sum_v_ort") = facetCounter.sum_v_ort;
			facetHitNode.append_attribute("sum_1_per_v") = facetCounter.sum_1_per_ort_velocity;
			facetHitNode.append_attribute("sum_v") = facetCounter.sum_1_per_velocity;

			if (work->model->sp.enableForceMeasurement) { //don't save all-zero quantities if not measured

				auto forcesNode = newFacetResult.append_child("Forces");

				auto impulseNode = forcesNode.append_child("Impulse");
				impulseNode.append_attribute("x") = facetCounter.impulse.x;
				impulseNode.append_attribute("y") = facetCounter.impulse.y;
				impulseNode.append_attribute("z") = facetCounter.impulse.z;

				auto impulse_square_Node = forcesNode.append_child("Impulse_square");
				impulse_square_Node.append_attribute("x") = facetCounter.impulse_square.x;
				impulse_square_Node.append_attribute("y") = facetCounter.impulse_square.y;
				impulse_square_Node.append_attribute("z") = facetCounter.impulse_square.z;

				auto impulse_momentum_Node = forcesNode.append_child("Impulse_momentum");
				impulse_momentum_Node.append_attribute("x") = facetCounter.impulse_momentum.x;
				impulse_momentum_Node.append_attribute("y") = facetCounter.impulse_momentum.y;
				impulse_momentum_Node.append_attribute("z") = facetCounter.impulse_momentum.z;

			}

			if (f->sh.isProfile) {
				xml_node profileNode = newFacetResult.append_child("Profile");
				profileNode.append_attribute("size") = PROFILE_SIZE;
				const std::vector<ProfileSlice>& pr = globalState->facetStates[i].momentResults[m].profile;
				for (size_t p = 0; p < std::min(PROFILE_SIZE, globalState->facetStates[i].momentResults[m].profile.size()); p++) {
					xml_node slice = profileNode.append_child("Slice");
					slice.append_attribute("id") = p;
					slice.append_attribute("countEquiv") = pr[p].countEquiv;
					slice.append_attribute("sum_1_per_v") = pr[p].sum_1_per_ort_velocity;
					slice.append_attribute("sum_v_ort") = pr[p].sum_v_ort;
				}
			}

			size_t profSize = (f->sh.isProfile) ? (PROFILE_SIZE * sizeof(ProfileSlice) * (1 + mApp->worker.interfaceMomentCache.size())) : 0;

			//Textures
			size_t h = f->sh.texHeight;
			size_t w = f->sh.texWidth;

			if (f->hasMesh) {
				xml_node textureNode = newFacetResult.append_child("Texture");
				textureNode.append_attribute("width") = f->sh.texWidth;
				textureNode.append_attribute("height") = f->sh.texHeight;

				const std::vector<TextureCell>& texture = globalState->facetStates[i].momentResults[m].texture;
				std::stringstream countText, sum1perText, sumvortText;
				countText << '\n'; //better readability in file
				sum1perText << std::setprecision(8) << '\n';
				sumvortText << std::setprecision(8) << '\n';

				for (size_t iy = 0; iy < h; iy++) {
					for (size_t ix = 0; ix < w; ix++) {
						countText << texture[iy * f->sh.texWidth + ix].countEquiv << '\t';
						sum1perText << texture[iy * f->sh.texWidth + ix].sum_1_per_ort_velocity << '\t';
						sumvortText << texture[iy * f->sh.texWidth + ix].sum_v_ort_per_area << '\t';

					}
					countText << '\n';
					sum1perText << '\n';
					sumvortText << '\n';
				}
				textureNode.append_child("count").append_child(node_cdata).set_value(countText.str().c_str());
				textureNode.append_child("sum_1_per_v").append_child(node_cdata).set_value(sum1perText.str().c_str());
				textureNode.append_child("sum_v_ort").append_child(node_cdata).set_value(sumvortText.str().c_str());

			} //end texture
			size_t textureSize = (1 + (int)work->interfaceMomentCache.size()) * w * h * sizeof(TextureCell);

			if (f->sh.countDirection && f->dirCache) {
				xml_node dirNode = newFacetResult.append_child("Directions");
				dirNode.append_attribute("width") = f->sh.texWidth;
				dirNode.append_attribute("height") = f->sh.texHeight;

				const std::vector<DirectionCell>& dirs = globalState->facetStates[i].momentResults[m].direction;

				std::stringstream dirText, dirCountText;
				dirText << std::setprecision(8) << '\n'; //better readability in file
				dirCountText << '\n';

				for (size_t iy = 0; iy < h; iy++) {
					for (size_t ix = 0; ix < w; ix++) {
						dirText << dirs[iy * f->sh.texWidth + ix].dir.x << ",";
						dirText << dirs[iy * f->sh.texWidth + ix].dir.y << ",";
						dirText << dirs[iy * f->sh.texWidth + ix].dir.z << "\t";
						dirCountText << dirs[iy * f->sh.texWidth + ix].count << "\t";

					}
					dirText << "\n";
					dirCountText << "\n";
				}
				dirNode.append_child("vel.vectors").append_child(node_cdata).set_value(dirText.str().c_str());
				dirNode.append_child("count").append_child(node_cdata).set_value(dirCountText.str().c_str());
			} //end directions
			size_t dirSize = f->sh.countDirection ? (1 + (int)work->interfaceMomentCache.size()) * w * h * sizeof(DirectionCell) : 0;

			size_t angleMapRecordedDataSize = sizeof(size_t) * (f->sh.anglemapParams.phiWidth *
				(f->sh.anglemapParams.thetaLowerRes +
					f->sh.anglemapParams.thetaHigherRes));

			//Facet histograms (1 per moment) comes here
			bool hasHistogram = f->sh.facetHistogramParams.recordBounce || f->sh.facetHistogramParams.recordDistance;
#ifdef MOLFLOW
			hasHistogram = hasHistogram || f->sh.facetHistogramParams.recordTime;
#endif
			if (hasHistogram) {
				xml_node histNode = newFacetResult.append_child("Histograms");
				//Retrieve histogram map from hits dp
				auto& histogram = globalState->facetStates[i].momentResults[m].histogram;
				if (f->sh.facetHistogramParams.recordBounce) {
					auto& nbHitsHistogram = histogram.nbHitsHistogram;
					xml_node hist = histNode.append_child("Bounces");
					size_t histSize = f->sh.facetHistogramParams.GetBounceHistogramSize();
					hist.append_attribute("size") = histSize;
					hist.append_attribute("binSize") = f->sh.facetHistogramParams.nbBounceBinsize; //redundancy for human-reading or export
					hist.append_attribute("max") = f->sh.facetHistogramParams.nbBounceMax; //redundancy for human-reading or export
					for (size_t h = 0; h < histSize; h++) {
						xml_node bin = hist.append_child("Bin");
						auto value = bin.append_attribute("start");
						if (h == histSize - 1) value = "overRange";
						else value = h * f->sh.facetHistogramParams.nbBounceBinsize;
						bin.append_attribute("count") = nbHitsHistogram[h];
					}
				}
				if (f->sh.facetHistogramParams.recordDistance) {
					auto& distanceHistogram = histogram.distanceHistogram;
					xml_node hist = histNode.append_child("Distance");
					size_t histSize = f->sh.facetHistogramParams.GetDistanceHistogramSize();
					hist.append_attribute("size") = histSize;
					hist.append_attribute("binSize") = f->sh.facetHistogramParams.distanceBinsize; //redundancy for human-reading or export
					hist.append_attribute("max") = f->sh.facetHistogramParams.distanceMax; //redundancy for human-reading or export
					for (size_t h = 0; h < histSize; h++) {
						xml_node bin = hist.append_child("Bin");
						auto value = bin.append_attribute("start");
						if (h == histSize - 1) value = "overRange";
						else value = h * f->sh.facetHistogramParams.distanceBinsize;
						bin.append_attribute("count") = distanceHistogram[h];
					}
				}
				if (f->sh.facetHistogramParams.recordTime) {
					auto& timeHistogram = histogram.timeHistogram;
					xml_node hist = histNode.append_child("Time");
					size_t histSize = f->sh.facetHistogramParams.GetTimeHistogramSize();
					hist.append_attribute("size") = histSize;
					hist.append_attribute("binSize") = f->sh.facetHistogramParams.timeBinsize; //redundancy for human-reading or export
					hist.append_attribute("max") = f->sh.facetHistogramParams.timeMax; //redundancy for human-reading or export
					for (size_t h = 0; h < histSize; h++) {
						xml_node bin = hist.append_child("Bin");
						auto value = bin.append_attribute("start");
						if (h == histSize - 1) value = "overRange";
						else value = h * f->sh.facetHistogramParams.timeBinsize;
						bin.append_attribute("count") = timeHistogram[h];
					}
				}
			}
		}
	}

	//Texture Min/Max
	xml_node minMaxNode = resultNode.append_child("TextureMinMax");
	minMaxNode.append_child("With_constant_flow").append_child("Pressure").append_attribute("min") = texture_limits[0].autoscale.min.steady_state;
	minMaxNode.child("With_constant_flow").child("Pressure").append_attribute("max") = texture_limits[0].autoscale.max.steady_state;
	minMaxNode.child("With_constant_flow").append_child("Density").append_attribute("min") = texture_limits[1].autoscale.min.steady_state;
	minMaxNode.child("With_constant_flow").child("Density").append_attribute("max") = texture_limits[1].autoscale.max.steady_state;
	minMaxNode.child("With_constant_flow").append_child("Imp.rate").append_attribute("min") = texture_limits[2].autoscale.min.steady_state;
	minMaxNode.child("With_constant_flow").child("Imp.rate").append_attribute("max") = texture_limits[2].autoscale.max.steady_state;

	minMaxNode.append_child("Moments_only").append_child("Pressure").append_attribute("min") = texture_limits[0].autoscale.min.moments_only;
	minMaxNode.child("Moments_only").child("Pressure").append_attribute("max") = texture_limits[0].autoscale.max.moments_only;
	minMaxNode.child("Moments_only").append_child("Density").append_attribute("min") = texture_limits[1].autoscale.min.moments_only;
	minMaxNode.child("Moments_only").child("Density").append_attribute("max") = texture_limits[1].autoscale.max.moments_only;
	minMaxNode.child("Moments_only").append_child("Imp.rate").append_attribute("min") = texture_limits[2].autoscale.min.moments_only;
	minMaxNode.child("Moments_only").child("Imp.rate").append_attribute("max") = texture_limits[2].autoscale.max.moments_only;

	//Convergence results
	xml_node convNode = resultNode.append_child("Convergence");

	int formulaId = 0;
	for (const auto& formulaVec : mApp->appFormulas->convergenceValues) {
		std::stringstream convText;
		convText << std::setprecision(10) << '\n';
		convText << std::scientific;
		for (const auto& convVal : formulaVec.conv_vec) {
			convText << convVal.first << "\t" << convVal.second << "\n";
		}
		xml_node newFormulaNode = convNode.append_child("ConvData");
		newFormulaNode.append_attribute("Formula") = mApp->appFormulas->formulas[formulaId]->GetExpression();
		xml_node newConv = newFormulaNode.append_child(node_cdata);
		newConv.set_value(convText.str().c_str());
		formulaId++;
	}

	return true;
}
*/

/**
* \brief To load geometry data from a XML file
* \param loadXML xml input file
* \param work thread worker handling the task
* \param prg GLProgress_GUI window where visualising of the import progress is shown
*/

/*

/**
* \brief To load geometry data from a XML file and insert into an existing structure
* \param loadXML xml input file
* \param work thread worker handling the task
* \param prg GLProgress_GUI window where visualising of the insert progress is shown
* \param newStr if a new super structure is to be used
*/
void MolflowGeometry::InsertModel(const std::shared_ptr<MolflowSimulationModel> loadedModel, const MolflowInterfaceSettings& interfaceSettings, Worker* work, GLProgress_Abstract& prg, bool newStr) {

	int targetStructId = viewStruct;
	if (targetStructId == -1) targetStructId = 0;
	UnselectAll();

	size_t oldNbVertex = sh.nbVertex;
	SetInterfaceVertices(loadedModel->vertices3, true);
	int oldNbStruct = (int)sh.nbSuper;
	SetInterfaceStructures(loadedModel->structures, true,newStr, targetStructId);
	//Parameters (needs to precede facets)
	TimeDependentParameters::InsertParametersBeforeCatalog(work->interfaceParameterCache, loadedModel->tdParams.parameters);
	
	for (const auto& selection : interfaceSettings.selections) {
		SelectionGroup s;
		s.name = selection.name;
		for (const auto& id_original : selection.facetIds)
			s.facetIds.push_back(id_original + sh.nbFacet); //offset selection numbers
		mApp->AddSelection(s);
	}

	for (const auto& newView : interfaceSettings.userViews) {
		mApp->AddView(newView);
	}

	for (const auto & newFormula : interfaceSettings.userFormulas) {
		auto expr_cpy = newFormula.expression;
		mApp->OffsetFormula(expr_cpy, (int)sh.nbFacet);
		mApp->appFormulas->AddFormula(newFormula.name, expr_cpy);
	}

	SetInterfaceFacets(loadedModel->facets, true, oldNbVertex, newStr ? oldNbStruct : targetStructId); //increases sh.nbFacet, must be after selection/formula offset

	InitializeGeometry();
}


/**
* \brief Compare two XML input files by comparing simulation results against each other with a threshold
* \param fileName_lhs xml input file 1
* \param fileName_rhs xml input file 2
* \param fileName_out xml output file
* \param cmpThreshold threshold for equality checks
* \return bool true if compare quasi equal
*/
bool MolflowGeometry::CompareXML_simustate(const std::string& fileName_lhs, const std::string& fileName_rhs,
	const std::string& fileName_out, double cmpThreshold) {
	return false;
}

void MolflowGeometry::SetInterfaceFacets(std::vector<std::shared_ptr<SimulationFacet>> sFacets, bool insert, size_t vertexOffset, int structOffset) {
	//General Facets
	size_t index = insert ? facets.size() : 0;
	InterfaceGeometry::SetInterfaceFacets(sFacets, insert,vertexOffset,structOffset);

	// Init Molflow properties
	for (auto& sFac : sFacets) {
		auto mfFac = std::static_pointer_cast<MolflowSimFacet>(sFac);
		auto& intFacet = facets[index];

		// Molflow
		intFacet->ogMap = mfFac->ogMap;
		intFacet->angleMapCache = mfFac->angleMap.pdf;

		if (intFacet->ogMap.outgassingMapWidth > 0 || intFacet->ogMap.outgassingMapHeight > 0
			|| intFacet->ogMap.outgassingFileRatioU > 0.0 || intFacet->ogMap.outgassingFileRatioV > 0.0) {
			intFacet->hasOutgassingFile = true;
		}
		index++;
	}
}

PhysicalValue InterfaceGeometry::GetPhysicalValue(InterfaceFacet* f, const PhysicalMode& mode, const double moleculesPerTP, const double densityCorrection, const double gasMass, const int index, const FacetMomentSnapshot& facetSnap) {

	//if x==y==-1 and buffer=NULL then returns facet value, otherwise texture cell [x,y] value
	//buff is either NULL or a (BYTE*) pointer to texture or direction buffer, must be locked by AccessDataport before call
	//texture position must be calculated before call (index=i+j*w)


	PhysicalValue result;

	if (index == -1) { //Facet mode
		//Not implemented yet, calculation in formula editor and Facet Details windows
	}
	else { //Texture cell mode
		switch (mode) {
		case PhysicalMode::CellArea:
			result.value = f->GetMeshArea(index);
			break;
		case PhysicalMode::MCHits:
			result.value = facetSnap.texture[index].countEquiv;
			break;
		case PhysicalMode::ImpingementRate:
		{
			double area = (f->GetMeshArea(index, true));
			if (area == 0.0) area = 1.0;
			result.value = facetSnap.texture[index].countEquiv / (area * 1E-4) * moleculesPerTP;
			break;
		}
		case PhysicalMode::ParticleDensity:
		{
			double density = facetSnap.texture[index].sum_1_per_ort_velocity / (f->GetMeshArea(index, true) * 1E-4) * moleculesPerTP * densityCorrection;
			result.value = density;
			break;
		}
		case PhysicalMode::GasDensity:
		{
			double density = facetSnap.texture[index].sum_1_per_ort_velocity / (f->GetMeshArea(index, true) * 1E-4) * moleculesPerTP * densityCorrection;
			result.value = density * gasMass / 1000.0 / 6E23;
			break;
		}
		case PhysicalMode::Pressure:
			result.value = facetSnap.texture[index].sum_v_ort_per_area * 1E4 * (gasMass / 1000 / 6E23) * 0.0100 * moleculesPerTP;  //1E4 is conversion from m2 to cm2; 0.01 is Pa->mbar
			break;
		case PhysicalMode::AvgGasVelocity:
			result.value = 4.0 * facetSnap.texture[index].countEquiv / facetSnap.texture[index].sum_1_per_ort_velocity; //Different from FacetDetails, since textures don't record sum_1_per_v to save memory
			break;
		case PhysicalMode::GasVelocityVector:
		{
			double denominator = (facetSnap.direction[index].count > 0) ? (1.0 / facetSnap.direction[index].count) : 1.0;
			result.vect = facetSnap.direction[index].dir * denominator;
			break;
		}
		case PhysicalMode::NbVelocityVectors:
			result.count = facetSnap.direction[index].count;
			break;
		}
	}

	return result;
}

std::tuple<double, double> MolflowGeometry::GetTextureAutoscaleMinMax() {
	double autoscaleMin, autoscaleMax;
	if (texAutoScaleMode == AutoscaleMomentsOnly) {
		autoscaleMin = texture_limits[textureMode].autoscale.min.moments_only;
		autoscaleMax = texture_limits[textureMode].autoscale.max.moments_only;
	}
	else if (texAutoScaleMode == AutoscaleMomentsAndConstFlow) {
		autoscaleMin = std::min(
			texture_limits[textureMode].autoscale.min.steady_state,
			texture_limits[textureMode].autoscale.min.moments_only);
		autoscaleMax = std::max(
			texture_limits[textureMode].autoscale.max.steady_state,
			texture_limits[textureMode].autoscale.max.moments_only);
	}
	else if (texAutoScaleMode == AutoscaleConstFlow) {
		autoscaleMin = texture_limits[textureMode].autoscale.min.steady_state;
		autoscaleMax = texture_limits[textureMode].autoscale.max.steady_state;
	}
	else {
		throw Error("Unknown texture autoscale mode (enum index {})", static_cast<int>(texAutoScaleMode));
	}
	return { autoscaleMin, autoscaleMax };
}