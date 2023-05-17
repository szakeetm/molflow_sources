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

#ifdef _WIN32
#define NOMINMAX
//#include <Windows.h>
#include <direct.h>
#include <process.h>
#include "SMP.h"

#else

#endif

#include <future>
#include <cmath>
//#include <cstdlib>
#include <fstream>
//#include <istream>
#include <filesystem>
#include <cereal/archives/binary.hpp>
#include <cereal/types/utility.hpp>
#include <stdexcept>
/*//#include <cereal/archives/xml.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>*/
#include <IO/InterfaceXML.h>
#include <Buffer_shared.h>
#include <Simulation/IDGeneration.h>
#include <Simulation/CDFGeneration.h>

#include "MolflowGeometry.h"
#include "Worker.h"
#include "GLApp/GLApp.h"
#include "GLApp/GLMessageBox.h"
#include "GLApp/GLProgress_GUI.hpp"
#include "GLApp/GLUnitDialog.h"
#include "Helper/MathTools.h"
#include "Helper/StringHelper.h"
#include "Facet_shared.h"
//#include "Simulation.h" //SHELEM
#include "Interface/GlobalSettings.h"
#include "Interface/FacetAdvParams.h"
#include "Interface/ProfilePlotter.h"


#include "ConvergencePlotter.h"

#if defined(MOLFLOW)

#include "MolFlow.h"

#endif

#if defined(SYNRAD)
#include "SynRad.h"
#endif

#include <ZipLib/ZipArchive.h>
//#include "ziplib/ZipArchiveEntry.h"
#include <ZipLib/ZipFile.h>
//#include "File.h" //File utils (Get extension, etc)
//#include "ProcessControl.h"
#include "versionId.h"
#include "TimeMoments.h"

/*
//Leak detection
#if defined(_DEBUG)
#define _CRTDBG_MAP_ALLOC
#define DEBUG_NEW new(_NORMAL_BLOCK, __FILE__, __LINE__)
#define new DEBUG_NEW
#endif
*/
using namespace pugi;

#if defined(MOLFLOW)
extern MolFlow* mApp;
#endif

#if defined(SYNRAD)
extern SynRad* mApp;
#endif

/**
* \brief Default constructor for a worker
*/
Worker::Worker() : simManager(0) {

	model = std::make_shared<MolflowSimulationModel>();
	//Molflow specific
	temperatures = std::vector<double>();
	desorptionParameterIDs = std::set<size_t>();
	moments = std::vector<Moment>();
	userMoments = std::vector<UserMoment>(); //strings describing moments, to be parsed
	CDFs = std::vector<std::vector<CDF_p>>();
	IDs = std::vector<IntegratedDesorption>();
	parameters = std::vector<Parameter>();
	needsReload = true;  //When main and subprocess have different geometries, needs to reload (synchronize)
	abortRequested = false;
	displayedMoment = 0; //By default, steady-state is displayed
	fullFileName = "";

	ResetWorkerStats();
	geom = new MolflowGeometry();

	simuTimer.ReInit();
	//startTime = 0.0f;
	//stopTime = 0.0f;
	//simuTime = 0.0f;

	fullFileName = "";
}

/**
* \brief Getter to retrieve the geometry on this worker
* \return pointer of the geometry object
*/
MolflowGeometry* Worker::GetMolflowGeometry() {
	return geom;
}

/**
* \brief Function for saving geometry to a set file
* \param fileName output file name with extension
* \param prg GLProgress_GUI window where loading is visualised
* \param askConfirm if a window should be opened to ask if the file should be overwritten
* \param saveSelected if a selection is to be saved
* \param autoSave if automatic saving is enabled
* \param crashSave if save on crash is enabled
*/
void Worker::SaveGeometry(std::string fileName, GLProgress_Abstract& prg, bool askConfirm, bool saveSelected, bool autoSave,
	bool crashSave) {

	try {
		if (needsReload && (!crashSave && !saveSelected)) RealReload();
		// Update persistent anglemap
		SendAngleMaps();
	}
	catch (const std::exception& e) {
		char errMsg[512];
		sprintf(errMsg, "Error reloading worker. Trying safe save (geometry only):\n%s", e.what());
		GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
		crashSave = true;
	}
	std::string compressCommandLine;
	std::string fileNameWithGeo; //file name with .geo extension (instead of .geo7z)
	std::string fileNameWithGeo7z;
	std::string fileNameWithXML;
	std::string fileNameWithZIP;
	std::string fileNameWithoutExtension; //file name without extension

	std::string ext = FileUtils::GetExtension(fileName);

	bool ok = true;
	if (ext.empty()) {
		fileName = fileName + (mApp->compressSavedFiles ? ".zip" : ".xml");
		ext = FileUtils::GetExtension(fileName);
	}

	// Read a file
	bool isTXT = Contains({ "txt", "TXT" }, ext);
	bool isSTR = Contains({ "str", "STR" }, ext);
	bool isGEO = ext == "geo";
	bool isGEO7Z = ext == "geo7z";
	bool isXML = ext == "xml";
	bool isXMLzip = ext == "zip";
	bool isSTL = ext == "stl";

	if (isTXT || isGEO || isGEO7Z || isSTR || isXML || isXMLzip || isSTL) {
#ifdef _WIN32
		//Check (using native handle) if background compressor is still alive
		if ((isGEO7Z) && WAIT_TIMEOUT == WaitForSingleObject(mApp->compressProcessHandle, 0)) {
			GLMessageBox::Display("Compressing a previous save file is in progress. Wait until that finishes "
				"or close process \"compress.exe\"\nIf this was an autosave attempt,"
				"you have to lower the autosave frequency.", "Can't save right now.", GLDLG_OK,
				GLDLG_ICONERROR);
			return;
		}
#endif
		if (isGEO) {
			fileNameWithoutExtension = fileName.substr(0, fileName.length() - 4);
			fileNameWithGeo7z = fileName + "7z";
			fileNameWithGeo = fileName;
		}
		else if (isGEO7Z) {
			fileNameWithoutExtension = fileName.substr(0, fileName.length() - 6);
			fileNameWithGeo = fileName.substr(0, fileName.length() - 2);
			fileNameWithGeo7z = fileName;
			std::ostringstream tmp;
			tmp << "A .geo file of the same name exists. Overwrite that file ?\n" << fileNameWithGeo;
			if (!autoSave && FileUtils::Exist(fileNameWithGeo)) {
				ok = (GLMessageBox::Display(tmp.str().c_str(), "Question", GLDLG_OK | GLDLG_CANCEL,
					GLDLG_ICONWARNING) == GLDLG_OK);
			}
		}

		if (isXML || isXMLzip) {
			fileNameWithoutExtension = FileUtils::StripExtension(fileName);
			fileNameWithXML = fileNameWithoutExtension + ".xml";
			fileNameWithZIP = fileNameWithoutExtension + ".zip";
		}

		if (isSTL) {
			//Nothing to prepare
			ok = true;
		}

		if (ok) {

			if (isSTR) {
				geom->SaveSTR(saveSelected);
			}
			else {
				try {
					std::string toOpen;

					if (isGEO7Z) {
						toOpen = fileNameWithGeo;
					}
					else if (!(isXML || isXMLzip)) {
						toOpen = fileName; //Txt, stl, geo, etc...
						auto file = FileWriter(toOpen);//We first write a GEO file, then compress it to GEO7Z later
						if (isTXT) {
							geom->SaveTXT(file, interfaceGlobalState, saveSelected);
						}
						else if (isGEO || isGEO7Z) {
							geom->SaveGEO(file, prg, interfaceGlobalState, this, saveSelected, crashSave);
						}
						else if (isSTL) {
							geom->SaveSTL(file, prg);
						}
					}
					else if (isXML || isXMLzip) {
						auto mf_model = std::dynamic_pointer_cast<MolflowSimulationModel>(model);
						std::stringstream xmlStream; //Will store XML file content
						{ //Scope to store XML tree
							xml_document saveDoc;
							FlowIO::WriterInterfaceXML writer(mApp->useOldXMLFormat, false);
							this->uInput.facetViewSettings.clear();
							for (size_t facetId = 0; facetId < geom->GetNbFacet(); facetId++) {
								auto facet = geom->GetFacet(facetId);
								bool textureVisible = facet->textureVisible;
								bool volumeVisible = facet->volumeVisible;
								this->uInput.facetViewSettings.emplace_back(textureVisible, volumeVisible);
							}
							this->uInput.userMoments = userMoments;
							this->uInput.parameters = parameters;
							this->uInput.selections = mApp->selections;

							writer.uInput = this->uInput;

							if (saveSelected)
								writer.SaveGeometry(saveDoc, mf_model, GetGeometry()->GetSelectedFacets());
							else
								writer.SaveGeometry(saveDoc, mf_model);
							FlowIO::WriterInterfaceXML::WriteInterface(saveDoc, mApp, saveSelected);

							xml_document geom_only;
							geom_only.reset(saveDoc);
							bool success = false; //success: simulation state could be saved
							if (!crashSave && !saveSelected) {
								try {
									//success = geom->SaveXML_simustate(saveDoc, this, interfaceGlobalState, prg, saveSelected);
									success = writer.SaveSimulationState(saveDoc, mf_model, interfaceGlobalState);
								}
								catch (const std::exception& e) {
									GLMessageBox::Display(e.what(), "Error saving simulation state.", GLDLG_OK,
										GLDLG_ICONERROR);
									return;
								}
							}

							prg.SetMessage("Assembling xml file in memory...");

							if (success) {
								saveDoc.save(xmlStream);
							}
							else {
								geom_only.save(xmlStream);
							}
						}//saveDoc goes out of scope for the compress duration

						if (isXMLzip) {
							prg.SetProgress(0.75);
							prg.SetMessage("Compressing xml to zip...");

							if (FileUtils::Exist(fileNameWithZIP)) {
								try {
									std::filesystem::remove(fileNameWithZIP);
								}
								catch (const std::exception& e) {
									std::string msg =
										"Can't overwrite \n" + fileNameWithZIP + "\nMaybe file is in use.\n" + e.what();
									GLMessageBox::Display(msg.c_str(), "File save error", GLDLG_OK, GLDLG_ICONERROR);
									return;
								}
							}
							try {
								auto zipFile = ZipFile::Open(fileNameWithZIP);
								auto entry = zipFile->CreateEntry(FileUtils::GetFilename(fileNameWithXML));
								auto method = DeflateMethod::Create();
								method->SetCompressionLevel(DeflateMethod::CompressionLevel::Fastest); //+15% size is acceptable
								entry->SetCompressionStream(   // data from contentStream are pumped here (into the memory)
									xmlStream,
									method,
									ZipArchiveEntry::CompressionMode::Deferred
								);
								prg.SetProgress(0.9);
								prg.SetMessage("Writing zip file...");
								ZipFile::SaveAndClose(zipFile, fileNameWithZIP);
							}
							catch (const std::exception& e) {
								std::string msg =
									"Can't write \n" + fileNameWithZIP + "\n" + e.what();
								GLMessageBox::Display(msg.c_str(), "File save error", GLDLG_OK, GLDLG_ICONERROR);
								return;
							}
						}
						else { //simple XML
							prg.SetProgress(0.9);
							prg.SetMessage("Writing xml file...");

							if (FileUtils::Exist(fileNameWithXML)) {
								try {
									std::filesystem::remove(fileNameWithXML);
								}
								catch (const std::exception& e) {
									std::string msg =
										"Can't overwrite \n" + fileNameWithXML + "\nMaybe file is in use.\n" + e.what();
									GLMessageBox::Display(msg.c_str(), "File save error", GLDLG_OK, GLDLG_ICONERROR);
									return;
								}
							}
							try {
								std::ofstream outFile;
								outFile.open(fileNameWithXML);
								outFile << xmlStream.rdbuf();
							} //outFile goes out of scope: file flushed
							catch (const std::exception& e) {
								std::string msg =
									"Can't write \n" + fileNameWithXML + "\n" + e.what();
								GLMessageBox::Display(msg.c_str(), "File save error", GLDLG_OK, GLDLG_ICONERROR);
								return;
							}
						} //end xml
					} //end xml or zip

				}
				catch (const std::exception& e) {
					GLMessageBox::Display(e.what(), "Error writing file.", GLDLG_OK, GLDLG_ICONERROR);
					return;
				}
			} //end not str
		} //end "ok"
	} //end valid extension
	else
	{
		throw std::runtime_error("SaveGeometry(): Invalid file extension [only xml,zip,geo,geo7z,txt,stl or str]");
	}

	//File written, compress it if the user wanted to
	if (ok && isGEO7Z) {

#ifdef _WIN32
		std::string compressorName = "compress.exe";
#else
		std::string compressorName = "./compress";
#endif

		if (FileUtils::Exist(compressorName)) { //compress GEO file to GEO7Z using 7-zip launcher "compress.exe"
			std::ostringstream tmp;
			tmp << compressorName << " \"" << fileNameWithGeo << "\" Geometry.geo";
#ifdef _WIN32
			char* command[1];
			command[0] = new char[512];
			sprintf(command[0], "%s", tmp.str().c_str());
			size_t procId = StartProc(command, STARTPROC_BACKGROUND);
			mApp->compressProcessHandle = OpenProcess(PROCESS_ALL_ACCESS, true, (unsigned long)procId);

			delete[] command[0];
#else
			//In Linux, compressing to old format will be blocking
			system(tmp.str().c_str());
#endif

			fileName = fileNameWithGeo7z;
		}
		else {
			GLMessageBox::Display("compress.exe (part of Molfow) not found.\n Will save as uncompressed GEO file.",
				"Compressor not found", GLDLG_OK, GLDLG_ICONERROR);
			fileName = fileNameWithGeo;
		}
	}
	else if (ok && isGEO) {
		fileName = fileNameWithGeo;
	}
	if (!autoSave && !saveSelected && !isSTL) { //STL file is just a copy
		SetCurrentFileName(fileName.c_str());
		mApp->UpdateTitle();
	}

	// Save angle maps together with Molflow_Autosave.zip; otherwise data is lost after crash
	/*if(autoSave){
		std::vector<std::string> listOfFiles = ExportAngleMaps(fileNameWithoutExtension, true);
		for(auto angleFile : listOfFiles){
			ZipFile::AddFile(fileNameWithZIP, angleFile,FileUtils::GetFilename(angleFile));
			//At this point, if no error was thrown, the compression is successful
			remove(angleFile.c_str());
		}
	}*/
	}

/**
* \brief Function for saving of the profile data (simulation)
* \param fileName output file name
*/
void Worker::ExportProfiles(const char* fn) {

	std::string fileName = fn;
	char tmp[512];
	// Read a file
	FILE* f = nullptr;

	if (FileUtils::GetExtension(fileName).empty()) {
		fileName += ".csv"; //set to default CSV format
		if (FileUtils::Exist(fileName)) {
			sprintf(tmp, "Overwrite existing file ?\n%s", fileName.c_str());
			if (GLMessageBox::Display(tmp, "Question", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONWARNING) != GLDLG_OK) return;
		}
	}
	bool isTXT = FileUtils::GetExtension(fileName) == "txt";


	f = fopen(fileName.c_str(), "w");
	if (!f) {
		sprintf(tmp, "Cannot open file for writing %s", fileName.c_str());
		throw std::runtime_error(tmp);
	}
	geom->ExportProfiles(f, isTXT, this);
	fclose(f);

}

/**
* \brief Function for exportation of the angle maps
* \param fileName output file name
 * \param saveAll true if all files -- otherwise only selected -- should be saved
 * \return Vector with strings containing the file names of all angle map files
*/
std::optional<std::vector<std::string>> Worker::ExportAngleMaps(const std::string& fileName, bool saveAll) {
	//returns false if cancelled or error, vector of file name to export otherwise, empty vector if no sleected facets have angle map
	bool overwriteAll = false;

	//Geometry *geom = GetGeometry();
	std::vector<size_t> angleMapFacetIndices;
	for (size_t i = 0; i < geom->GetNbFacet(); i++) {
		InterfaceFacet* f = geom->GetFacet(i);
		// saveAll facets e.g. when auto saving or just when selected
		if ((saveAll || f->selected) && !f->angleMapCache.empty()) {
			angleMapFacetIndices.push_back(i);
		}
	}

	std::string extension = FileUtils::GetExtension(fileName);
	bool isTXT = Contains({ "txt","TXT" }, extension);
	std::vector<std::string> listOfFiles;
	for (size_t facetIndex : angleMapFacetIndices) {
		std::string saveFileName;
		if (angleMapFacetIndices.size() == 1) {
			saveFileName = fileName; //as user specified
		}
		else {
			std::stringstream tmp;
			tmp << FileUtils::StripExtension(fileName) << "_facet" << facetIndex + 1 << "." << extension; //for example anglemap_facet22.csv
			saveFileName = tmp.str();
		}

		if (FileUtils::Exist(saveFileName)) {
			if (!overwriteAll) {
				std::vector<std::string> buttons = { "Cancel", "Overwrite" };
				if (angleMapFacetIndices.size() > 1) buttons.emplace_back("Overwrite All");
				int answer = GLMessageBox::Display("Overwrite existing file ?\n" + saveFileName, "Question", buttons,
					GLDLG_ICONWARNING);
				if (answer == 0) return std::nullopt; //User cancel, return empty vector, resulting in parent function cancel
				overwriteAll = (answer == 2);
			}
		}

		std::ofstream file;
		file.open(saveFileName);
		if (!file.is_open()) {
			std::string tmp = "Cannot open file for writing " + saveFileName;
			throw std::runtime_error(tmp.c_str());
		}
		file << geom->GetFacet(facetIndex)->GetAngleMap(isTXT ? 2 : 1);
		file.close();
		listOfFiles.push_back(saveFileName);
	}

	return listOfFiles; // false if angleMapFacetIndices.size() == 0
}

/**
* \brief Function for loading of the geometry
* \param fileName input file name
* \param insert if geometry will be inserted into an existing geometry view
* \param newStr if a new structure needs to be created
*/
void Worker::LoadGeometry(const std::string& fileName, bool insert, bool newStr) {
	if (!insert) {
		needsReload = true;
	}
	else {
		try {
			RealReload();
		}
		catch (const std::exception& e) {
			GLMessageBox::Display(e.what(), "Error (Reloading Geometry)", GLDLG_OK, GLDLG_ICONERROR);
		}
	}
	//char CWD[MAX_PATH];
	//_getcwd(CWD, MAX_PATH);

	std::string ext = FileUtils::GetExtension(fileName);

	if (ext.empty())

		throw std::runtime_error("LoadGeometry(): No file extension, can't determine type");


	auto prg = GLProgress_GUI("Reading file...", "Please wait");
	prg.SetVisible(true);

	ResetWorkerStats();

	if (!insert) {
		//Clear hits and leaks cache
		ResetMoments();
		model->wp = WorkerParams(); //reset to default
	}

	if (ext == "txt" || ext == "TXT") {

		try {

			if (insert) mApp->changedSinceSave = true;

			auto file = FileReader(fileName);

			if (!insert) {
				geom->LoadTXT(file, prg, this);
				//RealReload();
				fullFileName = fileName;
				RealReload();
				SendToHitBuffer(); //Global hit counters and hit/leak cache
				SendFacetHitCounts(); // From facetHitCache to dpHit's const.flow counter
			}
			else { //insert

				geom->InsertTXT(file, prg, newStr);
				MarkToReload();
			}
		}

		catch (const std::exception& e) {
			if (!insert) geom->Clear();
			throw;
		}

	}
	else if (ext == "stl" || ext == "STL") {
		try {
			int ret = GLUnitDialog::Display("", "Choose STL file units:",
				GLDLG_MM | GLDLG_CM | GLDLG_M | GLDLG_INCH | GLDLG_FOOT | GLDLG_CANCEL_U,
				GLDLG_ICONNONE);
			double scaleFactor = 1.0;
			switch (ret) {
			case GLDLG_MM:
				scaleFactor = 0.1;
				break;
			case GLDLG_CM:
				scaleFactor = 1.0;
				break;
			case GLDLG_M:
				scaleFactor = 100;
				break;
			case GLDLG_INCH:
				scaleFactor = 2.54;
				break;
			case GLDLG_FOOT:
				scaleFactor = 30.48;
				break;
			default:
				break;
			}
			if (ret != GLDLG_CANCEL_U) {
				prg.SetMessage("Resetting worker...");
				prg.SetVisible(true);
				prg.SetMessage("Reading geometry...");
				auto file = FileReader(fileName);
				if (!insert) {
					geom->LoadSTL(file, prg, scaleFactor);
					fullFileName = fileName;
					mApp->DisplayCollapseDialog();
				}
				else { //insert
					mApp->changedSinceSave = true;
					geom->InsertSTL(file, prg, scaleFactor, newStr);
					MarkToReload();
				}
			}
		}
		catch (const std::exception& e) {
			if (!insert) geom->Clear();
			throw;

		}

	}
	else if (ext == "str" || ext == "STR") {
		if (insert) throw std::runtime_error("STR file inserting is not supported.");
		try {
			auto file = FileReader(fileName);
			prg.SetVisible(true);
			geom->LoadSTR(file, prg);
			//RealReload();

			fullFileName = fileName;
		}

		catch (const std::exception& e) {
			geom->Clear();
			throw;
		}

	}
	else if (ext == "syn" || ext == "syn7z") { //Synrad file
		int version;
		prg.SetVisible(true);
		try {
			std::unique_ptr<FileReader> file;
			if (ext == "syn7z") {
				//decompress file
				prg.SetMessage("Decompressing file...");
				file.reset(FileReader::ExtractFrom7zAndOpen(fileName, "Geometry.syn"));
			}
			else {
				file.reset(new FileReader(fileName));  //original file opened
			}

			if (!insert) {
				geom->LoadSYN(*file, prg, &version, this);
				model->otfParams.desorptionLimit = 0;
			}
			else { //insert
				geom->InsertSYN(*file, prg, newStr);
			}

			prg.SetMessage("Reloading worker with new geometry...");
			MarkToReload();
			if (!insert) fullFileName = fileName;
		}

		catch (const std::exception& e) {
			if (!insert) geom->Clear();
			throw;
		}

	}
	else if (ext == "geo" || ext == "geo7z") {
		int version;
		prg.SetVisible(true);
		try {
			std::unique_ptr<FileReader> file;
			if (ext == "geo7z") {
				//decompress file
				prg.SetMessage("Decompressing file...");
				file.reset(FileReader::ExtractFrom7zAndOpen(fileName, "Geometry.geo"));
			}
			else { //not geo7z
				file.reset(new FileReader(fileName)); //geo file, open it directly
			}

			if (!insert) {

				geom->LoadGEO(*file, prg, &version, this);
				// Add moments only after user Moments are completely initialized
				if (TimeMoments::ParseAndCheckUserMoments(&moments, &userMoments, prg)) {
					GLMessageBox::Display("Overlap in time moments detected! Check in Moments Editor!", "Warning",
						GLDLG_OK, GLDLG_ICONWARNING);
					return;
				}

				prg.SetMessage("Reloading worker with new geometry...");
				RealReload(); //for the loading of textures

				if (version >= 8)
					geom->LoadProfileGEO(*file, interfaceGlobalState, version);

				SendToHitBuffer(); //Global hit counters and hit/leak cache
				SendFacetHitCounts(); // From facetHitCache to dpHit's const.flow counter
				SendAngleMaps();

				prg.SetMessage("Loading textures...");
				LoadTexturesGEO(*file, version);
				fullFileName = fileName;
			}
			else { //insert
				mApp->changedSinceSave = true;
				geom->InsertGEO(*file, prg, newStr);
				MarkToReload();
			}
		}

		catch (const std::exception& e) {
			if (!insert) geom->Clear();
			throw;
		}

	}
	else if (ext == "xml" || ext == "zip") { //XML file, optionally in ZIP container
		xml_document loadXML;
		xml_parse_result parseResult;
		std::string parseFileName = fileName;
		prg.SetVisible(true);
		try {
			if (ext == "zip") { //compressed in ZIP container
				//decompress file
				prg.SetMessage("Decompressing zip file...");

				ZipArchive::Ptr zip = ZipFile::Open(fileName);
				if (zip == nullptr) {
					throw std::runtime_error("Can't open ZIP file");
				}
				size_t numitems = zip->GetEntriesCount();
				bool notFoundYet = true;
				for (int i = 0; i < numitems && notFoundYet; i++) { //extract first xml file found in ZIP archive
					auto zipItem = zip->GetEntry(i);
					std::string zipFileName = zipItem->GetName();

					if (FileUtils::GetExtension(zipFileName) == "xml") { //if it's an .xml file
						notFoundYet = false;

						FileUtils::CreateDir("tmp");// If doesn't exist yet

						parseFileName = "tmp/" + zipFileName;
						ZipFile::ExtractFile(fileName, zipFileName, parseFileName);
					}
					/*else if(FileUtils::GetExtension(zipFileName) == "csv"){ // otherwise extract angle maps
						ZipFile::ExtractFile(fileName, zipFileName, zipFileName);
					}*/
				}
				/*if (!notFoundYet) {
					for (int i = 0; i < numitems; i++) { //extract first xml file found in ZIP archive
						auto zipItem = zip->GetEntry(i);
						std::string zipFileName = zipItem->GetName();
						if(FileUtils::GetExtension(zipFileName) == "csv"){ // extract angle maps
							ZipFile::ExtractFile(fileName, zipFileName, zipFileName);
							std::cout << "Import fileName: "<< fileName<<" - "<< zipFileName<<std::endl;
							ImportAngleMaps(fileName);
						}
					}
				}*/
				if (notFoundYet) {
					throw std::runtime_error("Didn't find any XML file in the ZIP file.");
				}

			}

			prg.SetMessage("Reading and parsing XML file...");
			// Parse zip file name or original
			parseResult = loadXML.load_file(parseFileName.c_str()); //parse xml file directly

			ResetWorkerStats();
			if (!parseResult) {
				//Parse error
				std::stringstream err;
				err << "XML parsed with errors.\n";
				err << "Error description: " << parseResult.description() << "\n";
				err << "Error offset: " << parseResult.offset << "\n";
				throw std::runtime_error(err.str().c_str());
			}

			prg.SetMessage("Building geometry...");
			xml_node rootNode = loadXML.root();
			if (appVersionId >= 2680) {
				xml_node envNode = loadXML.child("SimulationEnvironment");
				if (!envNode.empty())
					rootNode = envNode;
			}

			if (!insert) {
				//geom->LoadXML_geom(rootNode, this, prg);

				geom->Clear();
				FlowIO::LoaderInterfaceXML loader;
				auto mf_model = std::dynamic_pointer_cast<MolflowSimulationModel>(model);
				{
					try {
						loader.LoadGeometry(parseFileName, mf_model, prg);
					}
					catch (const std::exception& ex) {
						std::string msg = "There was an error loading this file:\n" + std::string(ex.what());
						throw std::runtime_error(msg);
					}
					/*
					auto future = std::async(std::launch::async, &FlowIO::LoaderInterfaceXML::LoadGeometry, &loader,
											 parseFileName, mf_model, load_progress);
					do {
						prg.SetProgress(load_progress);
						ProcessSleep(100);
					} while (future.wait_for(std::chrono::seconds(0)) != std::future_status::ready);

					if (future.get()) {
						prg.SetVisible(false);
						throw std::runtime_error("There was an error loading this file, check console for details.");
					}
					*/
				}
				geom->InitOldStruct(mf_model.get());

				//std::future<int> resultFromDB = std::async(std::launch::async, &FlowIO::LoaderInterfaceXML::LoadGeometryPtr, ldr, "Data", model.get(), load_progress);
				//loader.LoadGeometry(parseFileName, model, load_progress);
				//if (allowUpdateCheck) updateThread = std::thread(&AppUpdater::PerformUpdateCheck, (AppUpdater*)this, false); //Launch parallel update-checking thread
				prg.SetMessage("Loading interface settings...");

				xml_node interfNode = rootNode.child("Interface");
				userMoments = loader.uInput.userMoments;
				uInput = loader.uInput;
				InsertParametersBeforeCatalog(loader.uInput.parameters);

				*geom->GetGeomProperties() = model->sh;

				// Move actual geom to interface geom
				geom->InitInterfaceVertices(model->vertices3);
				geom->InitInterfaceFacets(model->facets, this);

				// Needs to be called after Interface Facets are loaded, as these are used e.g. when updating the ProfilePlotter state
				FlowIO::LoaderInterfaceXML::LoadInterface(interfNode, mApp);

				if (loader.uInput.facetViewSettings.size() == geom->GetNbFacet()) {
					for (size_t facetId = 0; facetId < geom->GetNbFacet(); facetId++) {
						auto facet = geom->GetFacet(facetId);
						facet->textureVisible = std::get<0>(loader.uInput.facetViewSettings[facetId]);
						facet->volumeVisible = std::get<1>(loader.uInput.facetViewSettings[facetId]);
					}
				}
				else {
					std::cerr << "Amount of view settings doesn't equal number of facets: "
						<< loader.uInput.facetViewSettings.size() << " <> " << GetGeometry()->GetNbFacet()
						<< std::endl;
				}

				prg.SetMessage("Parsing user moments...");
				// Add moments only after user Moments are completely initialized

				{
					try {
						TimeMoments::ParseAndCheckUserMoments(&moments, &userMoments, prg);
					}
					catch (...) {
						GLMessageBox::Display("Overlap in time moments detected! Check in Moments Editor!", "Warning",
							GLDLG_OK, GLDLG_ICONWARNING);
						return;
					}
					/*
					auto future = std::async(std::launch::async, TimeMoments::ParseAndCheckUserMoments, &moments, &userMoments, load_progress);
					do {
						prg.SetProgress(load_progress);
						ProcessSleep(100);
					} while (future.wait_for(std::chrono::seconds(0)) != std::future_status::ready);

					if (future.get()) {
						GLMessageBox::Display("Overlap in time moments detected! Check in Moments Editor!", "Warning",
											  GLDLG_OK, GLDLG_ICONWARNING);
						return;
					}
					*/
				}



				this->uInput = loader.uInput;
				// Init after load stage
				//geom->InitializeMesh();
				prg.SetMessage("Initializing geometry...");

				geom->InitializeGeometry();
				prg.SetMessage("Building mesh...");
				auto nbFacet = geom->GetNbFacet();


				for (size_t i = 0; i < nbFacet; i++) {
					double p = (double)i / (double)nbFacet;

					prg.SetProgress(p);
					auto f = geom->GetFacet(i);
					f->InitVisibleEdge();
					if (!f->SetTexture(f->sh.texWidth_precise, f->sh.texHeight_precise, f->hasMesh)) {
						char errMsg[512];
						sprintf(errMsg, "Not enough memory to build mesh on Facet %zd. ", i + 1);
						throw Error(errMsg);
					}
					geom->BuildFacetList(f);
					const double nU = f->sh.U.Norme();
					const double nV = f->sh.V.Norme();

					f->tRatioU = f->sh.texWidth_precise / nU;
					f->tRatioV = f->sh.texHeight_precise / nV;

					if (std::abs(f->tRatioU - f->tRatioV) <= DBL_EPSILON) {
						f->tRatioV = f->tRatioU;
					}

				}
				prg.SetMessage("Calculating OpenGL render data...");
				geom->InitializeInterfaceGeometry();

				geom->UpdateName(fileName.c_str());

				prg.SetMessage("Reloading worker with new geometry...");
				try {
					fullFileName = fileName;

					if (ext == "xml" || ext == "zip")
						prg.SetMessage("Restoring simulation state...");

					simManager.ForwardGlobalCounter(&interfaceGlobalState, &particleLog);
					RealReload(); //To create the dpHit dataport for the loading of textures, profiles, etc...
					{
						/*
						auto future = std::async(std::launch::async, FlowIO::LoadSimulationState,
												 parseFileName, mf_model, &interfaceGlobalState, load_progress);
						do {
							prg.SetProgress(load_progress);
							ProcessSleep(100);
						} while (future.wait_for(std::chrono::seconds(0)) != std::future_status::ready);

						try{
							future.get(); //exception thrown if it was stored
						}
						catch (const std::exception &e){
							throw;
						}
						*/

						FlowIO::LoaderXML::LoadSimulationState(parseFileName, mf_model, &interfaceGlobalState, prg);

						/*
						future = std::async(std::launch::async, FlowIO::LoaderInterfaceXML::LoadConvergenceValues,
											parseFileName, &mApp->formula_ptr->convergenceValues, load_progress);
						do {
							prg.SetProgress(load_progress);
							ProcessSleep(100);
						} while (future.wait_for(std::chrono::seconds(0)) != std::future_status::ready);
						*/
						FlowIO::LoaderInterfaceXML::LoadConvergenceValues(parseFileName, &mApp->formula_ptr->convergenceValues, prg);

					}
					simManager.simulationChanged = true; //mark for loading


					CalculateTextureLimits(); // Load texture limits on init

					// actually loads all caches
					RetrieveHistogramCache(); //So interface gets histogram data for disp.moment right after loadin
					SendToHitBuffer(); //Send global hits without sending facet counters, as they are directly written during the load process (mutiple moments)
					SendFacetHitCounts(); //Send hits without sending facet counters, as they are directly written during the load process (mutiple moments)
					SendAngleMaps();

					RebuildTextures();
				}
				catch (const std::exception& e) {
					if (!mApp->profilePlotter) {
						mApp->profilePlotter = new ProfilePlotter();
						mApp->profilePlotter->SetWorker(this);
					}
					mApp->profilePlotter->Reset(); //To avoid trying to display non-loaded simulation results
					if (!mApp->convergencePlotter) {
						mApp->convergencePlotter = new ConvergencePlotter(this, mApp->formula_ptr);
						mApp->convergencePlotter->SetWorker(this);
					}
					mApp->convergencePlotter->Reset(); //To avoid trying to display non-loaded simulation results
					GLMessageBox::Display(e.what(), "Error while loading simulation state", GLDLG_CANCEL,
						GLDLG_ICONWARNING);
					throw;
				}
			}
			else { //insert
				geom->InsertXML(rootNode, this, prg, newStr);
				model->sh = *geom->GetGeomProperties();
				mApp->changedSinceSave = true;
				ResetWorkerStats();

				MarkToReload();
			}
		}
		catch (const std::exception& e) {
			if (!insert) geom->Clear();
			throw;
		}

	}
	else {
		throw std::runtime_error(
			"LoadGeometry(): Invalid file extension [Only xml,zip,geo,geo7z,syn.syn7z,txt,stl or str]");
	}

	// Readers that load the geometry directly into the sim model
	// need to update the interface geometry afterwards
	if (ext == "xml" || ext == "zip") {
		if (!insert) {
			SimModelToInterfaceGeom();
		}
		else {
			InterfaceGeomToSimModel();
		}
	}
	else if (insert)
		InterfaceGeomToSimModel();

	if (!insert) {
		CalcTotalOutgassing();
		/*
		//Refresh is already done by the caller Molflow::LoadFile()

		if (mApp->timeSettings) mApp->timeSettings->RefreshMoments(); //Sets displayed moment to 0
		if (mApp->momentsEditor) mApp->momentsEditor->Refresh();
		if (mApp->parameterEditor) mApp->parameterEditor->UpdateCombo();
		if (mApp->timewisePlotter) mApp->timewisePlotter->Refresh();
		*/
	}

	if (insert) {
		mApp->UpdateFacetlistSelected();
		mApp->UpdateViewers();
	}
}

bool Worker::SimModelToInterfaceGeom() {

	*geom->GetGeomProperties() = model->sh;

	bool hasVolatile = false;

	for (size_t facIdx = 0; facIdx < model->sh.nbFacet; facIdx++) {
		assert(model->sh.nbFacet == model->facets.size());
		SimulationFacet* sFac = model->facets[facIdx].get();
		{
			InterfaceFacet* facet = geom->GetFacet(facIdx);

			facet->sh = sFac->sh;
			/*sFac.indices = facet->indices;
			sFac.vertices2 = facet->vertices2;
			sFac.ogMap = facet->ogMap;
			sFac.angleMap.pdf = angleMapVector;
			sFac.textureCellIncrements = textIncVector;*/
		}
	}

	return true;
}


/**
* \brief Function for loading textures from a GEO file
* \param f input file handle
* \param version version of the GEO data description
*/
void Worker::LoadTexturesGEO(FileReader& f, int version) {
	auto prg = GLProgress_GUI("Loading textures", "Please wait");
	try {
		prg.SetVisible(true);
		geom->LoadTexturesGEO(f, prg, interfaceGlobalState, version);
		RebuildTextures();
	}
	catch (const std::exception& e) {
		char tmp[256];
		sprintf(tmp,
			"Couldn't load some textures. To avoid continuing a partially loaded state, it is recommended to reset the simulation.\n%s",
			e.what());
		GLMessageBox::Display(tmp, "Error while loading textures.", GLDLG_OK, GLDLG_ICONWARNING);
	}
}

/**
* \brief Function that inserts a list of new paramters at the beginning of the catalog parameters
* \param newParams vector containing new parameters to be inserted
* \return index to insert position
*/
size_t Worker::InsertParametersBeforeCatalog(const std::vector<Parameter>& newParams) {
	size_t index = 0;
	for (; index != parameters.size() && parameters[index].fromCatalog == false; index++);
	parameters.insert(parameters.begin() + index, newParams.begin(),
		newParams.end()); //Insert to front (before catalog parameters)
	return index; //returns insert position
}

/* //Moved to worker_shared.cpp

void Worker::Update(float appTime) {
	if (needsReload) RealReload();
	//if (!needsReload) {
		// Check calculation ending
		bool done = true;
		bool error = true;
		if (dpControl) {
			if (AccessDataport(dpControl)) {
				int i = 0;
				SHCONTROL *master = (SHCONTROL *)dpControl->buff;
				for (int proc = 0; proc < nbProcess && done; proc++) {
					done = done && (master->states[proc] == PROCESS_DONE);
					error = error && (master->states[proc] == PROCESS_ERROR);
#if defined(MOLFLOW)
					if (master->states[proc] == PROCESS_RUNAC) calcACprg = master->cmdParam[proc];
#endif
				}
				ReleaseDataport(dpControl);
			}
		}

		// End of simulation reached (Stop GUI)
		if ((error || done) && running && appTime != 0.0f) {
			InnerStop(appTime);
			if (error) ThrowSubProcError();
		}

		// Retrieve hit count recording from the shared memory
		if (dpHit) {

			if (AccessDataport(dpHit)) {
				BYTE *buffer = (BYTE *)dpHit->buff;

				mApp->changedSinceSave = true;
				// Globals
				SHGHITS *gHits = (SHGHITS *)buffer;

				// Copy Global hits and leaks
				nbMCHit = gHits->globalHits.hit.nbMCHit;
				nbAbsEquiv = gHits->globalHits.hit.nbAbsEquiv;
				nbDesorption = gHits->globalHits.hit.nbDesorbed;
				distTraveled_total = gHits->distTraveled_total;
				distTraveledTotal_fullHitsOnly = gHits->distTraveledTotal_fullHitsOnly;


				nbLeakTotal = gHits->nbLeakTotal;
				hitCacheSize = gHits->hitCacheSize;
				memcpy(hitCache, gHits->hitCache, sizeof(HIT)*hitCacheSize);
				leakCacheSize = gHits->leakCacheSize;
				memcpy(leakCache, gHits->leakCache, sizeof(LEAK)*leakCacheSize); //will display only first leakCacheSize leaks

				// Refresh local facet hit cache for the displayed moment
				int nbFacet = geom->GetNbFacet();
				for (int i = 0; i < nbFacet; i++) {
					Facet *f = geom->GetFacet(i);
					f->facetHitCache=(*((FacetHitBuffer*)(buffer + f->model->wp.hitOffset+displayedMoment*sizeof(FacetHitBuffer))));
				}
				try {
					if (mApp->needsTexture || mApp->needsDirection) geom->BuildFacetTextures(buffer,mApp->needsTexture,mApp->needsDirection);
				}
				catch (const std::exception &e) {
					GLMessageBox::Display(e.what(), "Error building texture", GLDLG_OK, GLDLG_ICONERROR);
					ReleaseDataport(dpHit);
					return;
				}
				ReleaseDataport(dpHit);
			}

		}
	//}

}
*/

/**
* \brief Saves current AngleMap from cache to results
*/
int Worker::SendAngleMaps() {
	size_t nbFacet = geom->GetNbFacet();
	std::vector<std::vector<size_t>> angleMapCaches;
	for (size_t i = 0; i < nbFacet; i++) {
		InterfaceFacet* f = geom->GetFacet(i);
		angleMapCaches.push_back(f->angleMapCache);
	}

	if (interfaceGlobalState.facetStates.size() != angleMapCaches.size())
		return 1;
	if (!interfaceGlobalState.tMutex.try_lock_for(std::chrono::seconds(10)))
		return 1;
	for (size_t i = 0; i < angleMapCaches.size(); i++) {
		auto* sFac = (MolflowSimFacet*)model->facets[i].get();
		if (sFac->sh.anglemapParams.record)
			interfaceGlobalState.facetStates[i].recordedAngleMapPdf = angleMapCaches[i];
		//else if(sFac->sh.desorbType == DES_ANGLEMAP)
		sFac->angleMap.pdf = angleMapCaches[i];
	}
	interfaceGlobalState.tMutex.unlock();
	return 0;
}

bool Worker::InterfaceGeomToSimModel() {
	//Converts iterface-extended geometry to subprocess-formatted geometry, but doesn't forward or copy anything
	//result is stored in worker::model
	//InterfaceVertex -> Vector3d
	//InterfaceFacet-> MolflowSimFacet
	//cellPropertiesIds->textIncVector
	//etc.   

	//auto geom = GetMolflowGeometry();
	// TODO: Proper clear call before for Real reload?
	auto mf_model = std::dynamic_pointer_cast<MolflowSimulationModel>(model);
	mf_model->structures.clear();
	mf_model->facets.clear();
	mf_model->vertices3.clear();

	for (size_t nbV = 0; nbV < geom->GetNbVertex(); ++nbV) {
		mf_model->vertices3.emplace_back(*geom->GetVertex(nbV)); //InterfaceVertex->Vertex3d conversion
	}
	// Parse usermoments to regular moment intervals
	//mf_model->tdParams.moments = this->moments;


	mf_model->tdParams.CDFs.clear();
	mf_model->tdParams.IDs.clear();

	// we create it directly on the Sim side
	//mf_model->tdParams.CDFs = this->CDFs;
	//        mf_model->tdParams.IDs = this->IDs;
	{
		//mf_model->tdParams.IDs.clear();
		// we create it directly on the Sim side
		/*for (auto &id : this->IDs) {
			mf_model->tdParams.IDs.push_back(id.GetValues());
		}*/
	}

	mf_model->tdParams.parameters.clear();
	for (auto& param : this->parameters)
		mf_model->tdParams.parameters.emplace_back(param); //parameter->Distribution2D conversion

	std::vector<Moment> momentIntervals;
	momentIntervals.reserve(this->moments.size());
	for (auto& moment : this->moments) {
		momentIntervals.emplace_back(
			std::make_pair(moment.first - (0.5 * moment.second), moment.first + (0.5 * moment.second)));
	}
	mf_model->tdParams.moments = momentIntervals;

	mf_model->sh = *geom->GetGeomProperties();

	mf_model->structures.resize(mf_model->sh.nbSuper); //Create structures

	bool hasVolatile = false;

	for (size_t facIdx = 0; facIdx < mf_model->sh.nbFacet; facIdx++) {
		MolflowSimFacet sFac;
		{
			InterfaceFacet* facet = geom->GetFacet(facIdx);

			//std::vector<double> outgMapVector(sh.useOutgassingFile ? sh.outgassingMapWidth*sh.outgassingMapHeight : 0);
			//memcpy(outgMapVector.data(), outgassingMapWindow, sizeof(double)*(sh.useOutgassingFile ? sh.outgassingMapWidth*sh.outgassingMapHeight : 0));
			size_t mapSize = facet->sh.anglemapParams.GetMapSize();
			if (facet->angleMapCache.size() != facet->sh.anglemapParams.GetRecordedMapSize()) {
				// on mismatch between cached values, check if interface just got out of sync (record) or interface and simulation side are out of sync (no record)
				if (facet->sh.anglemapParams.record) {
					facet->angleMapCache.clear();
					facet->angleMapCache.resize(mapSize);
					needsAngleMapStatusRefresh = true; //Mark that facet adv. parameters needs to be updated
				}
				else {
					/*auto errString = fmt::format(
							"[Facet #{}] Recorded Data Size is different from actual size: {} / {}\n",
							facIdx + 1,
							facet->angleMapCache.size(),
							facet->sh.anglemapParams.GetRecordedMapSize());
					fmt::print(stderr, errString);
					throw std::runtime_error(errString);*/
				}
			}

			std::vector<double> textIncVector;

			// Add surface elements area (reciprocal)
			if (facet->sh.isTextured) {
				textIncVector.resize(facet->sh.texHeight * facet->sh.texWidth);
				if (!facet->cellPropertiesIds.empty()) {
					size_t add = 0;
					for (size_t j = 0; j < facet->sh.texHeight; j++) {
						for (size_t i = 0; i < facet->sh.texWidth; i++) {
							double area = facet->GetMeshArea(add, true);

							if (area > 0.0) {
								// Use the sign bit to store isFull flag
								textIncVector[add] = 1.0 / area;
							}
							else {
								textIncVector[add] = 0.0;
							}
							add++;
						}
					}
				}
				else {

					double rw = facet->sh.U.Norme() / facet->sh.texWidth_precise;
					double rh = facet->sh.V.Norme() / facet->sh.texHeight_precise;
					double area = rw * rh;
					size_t add = 0;
					for (int j = 0; j < facet->sh.texHeight; j++) {
						for (int i = 0; i < facet->sh.texWidth; i++) {
							if (area > 0.0) {
								textIncVector[add] = 1.0 / area;
							}
							else {
								textIncVector[add] = 0.0;
							}
							add++;
						}
					}
				}
			}


			sFac.sh = facet->sh;
			sFac.indices = facet->indices;
			sFac.vertices2 = facet->vertices2;
			sFac.ogMap = facet->ogMap;
			sFac.angleMap.pdf = facet->angleMapCache;
			sFac.textureCellIncrements = textIncVector;
		}

		//Some initialization
		try {
			if (!sFac.InitializeOnLoad(facIdx, mf_model->tdParams.moments.size()))
				return false;
		}
		catch (const std::exception& err) {
			//Log::console_error("Failed to initialize facet (F#%d)\n", facIdx + 1);
			throw;
		}

		hasVolatile |= sFac.sh.isVolatile;

		if ((sFac.sh.superDest || sFac.sh.isVolatile) &&
			((sFac.sh.superDest - 1) >= mf_model->sh.nbSuper || sFac.sh.superDest < 0)) {
			// Geometry error
			//ClearSimulation();
			//ReleaseDataport(loader);
			std::ostringstream err;
			err << "Invalid structure (wrong link on F#" << facIdx + 1 << ")";
			//SetErrorSub(err.str().c_str());
			std::cerr << err.str() << std::endl;
			return false;
		}

		mf_model->facets.push_back(std::make_shared<MolflowSimFacet>(sFac));
	}

	if (!mf_model->facets.empty() && !mf_model->vertices3.empty())
		mf_model->initialized = true;
	return true;
}

/**
* \brief Function that reloads the whole simulation (resets simulation, rebuilds ray tracing etc) and synchronises subprocesses to main process
* \param sendOnly if only certain parts should be reloaded (geometry reloading / ray tracing tree)
*/
void Worker::RealReload(bool sendOnly) { //Sharing geometry with workers

	GLProgress_GUI prg("Performing preliminary calculations on geometry...",
		"Passing Geometry to workers");
	prg.SetVisible(true);

	if (!sendOnly) {
		if (model->otfParams.nbProcess == 0 && !geom->IsLoaded()) {
			return;
		}

		try {
			PrepareToRun();
		}
		catch (const std::exception& e) {
			GLMessageBox::Display(e.what(), "Error (Full reload)", GLDLG_OK, GLDLG_ICONWARNING);
			std::stringstream err;
			err << "Error (Full reload) " << e.what();
			throw std::runtime_error(err.str());
		}
	}

	prg.SetMessage("Reloading structures for simulation unit...");
	ReloadSim(sendOnly, prg); //Convert interf. geom to worker::mode and construct global counters, then copy to simManagar.simulations[0]
	
	auto mf_model = std::dynamic_pointer_cast<MolflowSimulationModel>(model);
	mf_model->CalcTotalOutgassing(); // needs IDs

	if (!sendOnly) {
		try {
			prg.SetMessage("Asking subprocesses to clear geometry...");
			simManager.ResetSimulations();
			prg.SetMessage("Clearing Logger...");
			particleLog.clear();
			prg.SetMessage("Creating hit buffer...");
			simManager.ResetHits();
		}
		catch (const std::exception& e) {
			GLMessageBox::Display(e.what(), "Error (Full reload)", GLDLG_OK, GLDLG_ICONWARNING);
			std::stringstream err;
			err << "Error (Full reload) " << e.what();
			throw std::runtime_error(err.str());
		}

		if (model->otfParams.enableLogging) {
			try {
				particleLog.resize(model->otfParams.logLimit);
			}
			catch (...) {
				throw Error(
					"Failed to create 'Particle Log' vector.\nMost probably out of memory.\nReduce number of logged particles in Particle Logger.");
			}
		}
	}

	// Send and Load geometry on simulation side
	simManager.simulationChanged = true; //mark for loading by threads

	prg.SetMessage("Finishing reload...");
	needsReload = false;
}

/**
* \brief Serialization function for a binary cereal archive for the worker attributes
* \return output string stream containing the result of the archiving
*/
std::ostringstream Worker::SerializeParamsForLoader() {
	std::ostringstream result;
	cereal::BinaryOutputArchive outputArchive(result);

	outputArchive(
		CEREAL_NVP(model->otfParams)
	);
	return result;
}

/**
* \brief Resets workers global hit cache
* Reset function mainly used for initialisation / reload procedures
*/
void Worker::ResetWorkerStats() {

	interfaceGlobalState.Reset();
	particleLog.clear();
	//memset(&interfaceGlobalState.globalHits, 0, sizeof(GlobalHitBuffer));


}

/**
* \brief Starts the simulation process
*/
void Worker::Start() {
	// Sanity checks
	// Is there some desorption in the system? (depends on pre calculation)
	if (model->wp.finalOutgassingRate_Pa_m3_sec <= 0.0) {
		// Do another check for existing desorp facets, needed in case a desorp parameter's final value is 0
		bool found = false;
		size_t nbF = geom->GetNbFacet();
		size_t i = 0;
		while (i < nbF && !found) {
			found = (geom->GetFacet(i)->sh.desorbType != DES_NONE);
			if (!found) i++;
		}

		if (!found)
			throw Error("No desorption facet found");
	}
	if (model->wp.totalDesorbedMolecules <= 0.0)
		throw std::runtime_error("Total outgassing is zero.");

	if (model->otfParams.desorptionLimit > 0 && model->otfParams.desorptionLimit <= interfaceGlobalState.globalHits.globalHits.nbDesorbed)
		throw std::runtime_error("Desorption limit has already been reached.");

	try {
		if (simManager.StartSimulation()) {
			throw std::logic_error("Processes are already done!");
		}
	}
	catch (const std::exception& e) {
		throw Error(e.what()); //convert to runtime error
	}
}

/**
* \brief Resets/clears all moment variables
*/
void Worker::ResetMoments() {
	displayedMoment = 0;
	moments.clear();
	userMoments.clear();
}

/**
* \brief Returns how many physical molecules one test particle represents
* \param moment if there is constant flow or time-dependent mode
* \return amount of physical molecules represented by one test particle
*/
double Worker::GetMoleculesPerTP(size_t moment) const {
	if (interfaceGlobalState.globalHits.globalHits.nbDesorbed == 0) return 0; //avoid division by 0
	if (moment == 0) {
		//Constant flow
		//Each test particle represents a certain real molecule influx per second
		return model->wp.finalOutgassingRate / interfaceGlobalState.globalHits.globalHits.nbDesorbed;
	}
	else {
		//Time-dependent mode
		//Each test particle represents a certain absolute number of real molecules. Since Molflow displays per-second values (imp.rate, etc.), the sampled time window length is only a fraction of a second.
		//For example, if dt=0.1s, we have collected only 1/10th of what would happen during a second. Hence we DIVIDE by the time window length, even if it's uninuitional.
		return (model->wp.totalDesorbedMolecules / mApp->worker.moments[moment - 1].second) /
			interfaceGlobalState.globalHits.globalHits.nbDesorbed;
	}
}

/**
* \brief Importing desorption data from a SYN file
* \param fileName name of the input file
* \param source what the source to calculate the dose is
* \param time time to calculate the dose
* \param mode mode used for outgassing calculation
* \param eta0 coefficient for outgassing calculation in mode==1
* \param alpha exponent for outgassing calculation in mode==1
* \param cutoffdose cutoff dose for outgassing calculation in mode==1
* \param convDistr distribution for outgassing calculation in mode==2
* \param prg GLProgress_GUI window where visualising of the import progress is shown
*/
void Worker::ImportDesorption_SYN(const char* fileName, const size_t source, const double time,
	const size_t mode, const double eta0, const double alpha, const double cutoffdose,
	const std::vector<std::pair<double, double>>& convDistr,
	GLProgress_Abstract& prg) {
	std::string ext = FileUtils::GetExtension(fileName);
	if (!Contains({ "syn7z", "syn" }, ext))
		throw std::runtime_error("ImportDesorption_SYN(): Invalid file extension [Only syn, syn7z]");

	// Read a file

	std::unique_ptr<FileReader> file;

	bool isSYN7Z = (iequals(ext, "syn7z"));
	bool isSYN = (iequals(ext, "syn"));

	if (isSYN || isSYN7Z) {
		if (isSYN7Z) {
			file.reset(FileReader::ExtractFrom7zAndOpen(fileName, "Geometry.syn"));
		}
		else {
			file.reset(new FileReader(fileName));  //original file opened
		}

		geom->ImportDesorption_SYN(*file, source, time, mode, eta0, alpha, cutoffdose, convDistr, prg);
		CalcTotalOutgassing();
	}
}

/**
* \brief To analyze desorption data from a SYN file
* \param fileName name of the input file
* \param nbFacet number of facets in the file
* \param nbTextured number of textured facets in the file
* \param nbDifferent (TODO: check usage)
*/
void Worker::AnalyzeSYNfile(const char* fileName, size_t* nbFacet, size_t* nbTextured, size_t* nbDifferent) {
	std::string ext = FileUtils::GetExtension(fileName);
	bool isSYN = (ext == "syn") || (ext == "SYN");
	bool isSYN7Z = (ext == "syn7z") || (ext == "SYN7Z");

	if (!(isSYN || isSYN7Z))
		throw std::runtime_error("AnalyzeSYNfile(): Invalid file extension [Only syn, syn7z]");

	// Read a file
	std::unique_ptr<FileReader> file;

	auto prg = GLProgress_GUI("Analyzing SYN file...", "Please wait");
	prg.SetVisible(true);

	if (isSYN || isSYN7Z) {
		prg.SetVisible(true);
		if (isSYN7Z) {
			//decompress file
			prg.SetMessage("Decompressing file...");
			file.reset(FileReader::ExtractFrom7zAndOpen(fileName, "Geometry.syn"));
		}
		else { //syn
			file.reset(new FileReader(fileName));  //original file opened
		}
		geom->AnalyzeSYNfile(*file, prg, nbFacet, nbTextured, nbDifferent);
	}
}

/**
* \brief Do calculations necessary before launching simulation
* determine latest moment
* Generate integrated desorption functions
* match parameters
* Generate speed distribution functions
* Angle map
*/
void Worker::PrepareToRun() {

	//determine latest moment
	if (!moments.empty())
		model->wp.latestMoment = (moments.end() - 1)->first + (moments.end() - 1)->second / 2.0;
	else {
		model->wp.latestMoment = model->wp.timeWindowSize * .5;
	}

	Geometry* g = GetGeometry();

	//Reset Maxwell-Boltzmann precalculated distributions
	temperatures = std::vector<double>();
	desorptionParameterIDs = std::set<size_t>();
	CDFs = std::vector<std::vector<CDF_p>>();
	IDs = std::vector<IntegratedDesorption>();

	bool needsAngleMapStatusRefresh = false;

	for (size_t i = 0; i < g->GetNbFacet(); i++) {
		InterfaceFacet* f = g->GetFacet(i);

		//match time-dependent parameters
		if (f->userOutgassing.length() > 0) {
			int id = GetParamId(f->userOutgassing);
			if (id == -1) { //parameter not found
				char tmp[256];
				sprintf(tmp, "Facet #%zd: Outgassing parameter \"%s\" isn't defined.", i + 1,
					f->userOutgassing.c_str());
				throw std::runtime_error(tmp);
			}
			else f->sh.outgassing_paramId = id;
		}
		else f->sh.outgassing_paramId = -1;

		if (f->userOpacity.length() > 0) {
			int id = GetParamId(f->userOpacity);
			if (id == -1) { //parameter not found
				char tmp[256];
				sprintf(tmp, "Facet #%zd: Opacity parameter \"%s\" isn't defined.", i + 1, f->userOpacity.c_str());
				throw std::runtime_error(tmp);
			}
			else f->sh.opacity_paramId = id;
		}
		else f->sh.opacity_paramId = -1;

		if (f->userSticking.length() > 0) {
			int id = GetParamId(f->userSticking);
			if (id == -1) { //parameter not found
				char tmp[256];
				sprintf(tmp, "Facet #%zd: Sticking parameter \"%s\" isn't defined.", i + 1, f->userSticking.c_str());
				throw std::runtime_error(tmp);
			}
			else f->sh.sticking_paramId = id;
		}
		else f->sh.sticking_paramId = -1;

	}
}

/**
* \brief Get ID (if it exists) of the Commulative Distribution Function (CFD) for a particular temperature (bin)
* \param temperature temperature for the CFD
* \return ID of the CFD
*/
int Worker::GetCDFId(double temperature) const {
	return CDFGeneration::GetCDFId(temperatures, temperature);
}

/**
* \brief Generate a new Commulative Distribution Function (CFD) for a particular temperature (bin)
* \param temperature for the CFD
* \return Previous size of temperatures vector, which determines new ID
*/
int Worker::GenerateNewCDF(double temperature) {
	/*size_t i = temperatures.size();
	temperatures.push_back(temperature);
	CDFs.push_back(Generate_CDF(temperature, model->wp.gasMass, CDF_SIZE));
	return (int) i;*/

	auto [id_new, cdf_vec] = CDFGeneration::GenerateNewCDF(temperatures, temperature, model->wp.gasMass);
	CDFs.emplace_back(std::move(cdf_vec));

	return (int)id_new;
}

/**
* \brief Generate a new ID (integrated desorption) for desorption parameter for time-dependent simulations
* \param paramId parameter ID
* \return Previous size of IDs vector, which determines new id in the vector
*/
int Worker::GenerateNewID(size_t paramId) {
	//This function is called if parameter with index paramId doesn't yet have a cumulative des. function
	auto mf_model = std::dynamic_pointer_cast<MolflowSimulationModel>(model);
	auto [id_new, id_vec] = IDGeneration::GenerateNewID(desorptionParameterIDs, paramId, mf_model.get());
	Parameter& par = parameters[paramId]; //we'll reference it a lot
	IntegratedDesorption result;
	result.logXinterp = par.logXinterp;
	result.logYinterp = par.logYinterp;
	result.SetValues(std::move(id_vec), false);
	IDs.emplace_back(std::move(result));

	return (int)id_new;
}

/**
* \brief Get ID (if it exists) of the integrated desorption (ID) function for a particular paramId
* \param paramId parameter ID
* \return Id of the integrated desorption function
*/
int Worker::GetIDId(size_t paramId) const {
	return IDGeneration::GetIDId(desorptionParameterIDs, paramId);
}

/**
* \brief Compute the outgassing of all source facet depending on the mode (file, regular, time-dependent) and set it to the global settings
*/
void Worker::CalcTotalOutgassing() {
	// Compute the outgassing of all source facet
	auto mf_model = std::dynamic_pointer_cast<MolflowSimulationModel>(model);
	mf_model->wp.totalDesorbedMolecules = mf_model->wp.finalOutgassingRate_Pa_m3_sec = mf_model->wp.finalOutgassingRate = 0.0;
	Geometry* g = GetGeometry();

	for (size_t i = 0; i < g->GetNbFacet(); i++) {
		InterfaceFacet* f = g->GetFacet(i);
		if (f->sh.desorbType != DES_NONE) { //there is a kind of desorption
			if (f->sh.useOutgassingFile) { //outgassing file
				for (size_t l = 0; l < (f->ogMap.outgassingMapWidth * f->ogMap.outgassingMapHeight); l++) {
					mf_model->wp.totalDesorbedMolecules +=
						mf_model->wp.latestMoment * f->ogMap.outgassingMap[l] / (1.38E-23 * f->sh.temperature);
					mf_model->wp.finalOutgassingRate += f->ogMap.outgassingMap[l] / (1.38E-23 * f->sh.temperature);
					mf_model->wp.finalOutgassingRate_Pa_m3_sec += f->ogMap.outgassingMap[l];
				}
			}
			else { //regular outgassing
				if (f->sh.outgassing_paramId == -1) { //constant outgassing
					mf_model->wp.totalDesorbedMolecules +=
						mf_model->wp.latestMoment * f->sh.outgassing / (1.38E-23 * f->sh.temperature);
					mf_model->wp.finalOutgassingRate +=
						f->sh.outgassing / (1.38E-23 * f->sh.temperature);  //Outgassing molecules/sec
					mf_model->wp.finalOutgassingRate_Pa_m3_sec += f->sh.outgassing;
				}
				else { //time-dependent outgassing
					if (f->sh.IDid >= mf_model->tdParams.IDs.size())
						throw std::runtime_error(fmt::format("Trying to access Integrated Desorption {} of {} for facet #{}", f->sh.IDid, mf_model->tdParams.IDs.size(), i));

					double lastValue = mf_model->tdParams.IDs[f->sh.IDid].back().second;
					mf_model->wp.totalDesorbedMolecules += lastValue / (1.38E-23 * f->sh.temperature);
					size_t lastIndex = parameters[f->sh.outgassing_paramId].GetSize() - 1;
					double finalRate_mbar_l_s = parameters[f->sh.outgassing_paramId].GetY(lastIndex);
					mf_model->wp.finalOutgassingRate +=
						finalRate_mbar_l_s * MBARLS_TO_PAM3S /
						(1.38E-23 * f->sh.temperature); //0.1: mbar*l/s->Pa*m3/s
					mf_model->wp.finalOutgassingRate_Pa_m3_sec += finalRate_mbar_l_s * MBARLS_TO_PAM3S;
				}
			}
		}
	}
	if (mApp->globalSettings) mApp->globalSettings->UpdateOutgassing();

}

/**
* \brief Get ID of a parameter (if it exists) for a corresponding name
* \param name name of the parameter that shall be looked up
* \return ID corresponding to the found parameter
*/
int Worker::GetParamId(const std::string& name) {
	int foundId = -1;
	for (int i = 0; foundId == -1 && i < (int)parameters.size(); i++)
		if (name == parameters[i].name) foundId = i;
	return foundId;
}