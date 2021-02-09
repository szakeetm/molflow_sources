/*
Program:     MolFlow+ / Synrad+
Description: Monte Carlo simulator for ultra-high vacuum and synchrotron radiation
Authors:     Jean-Luc PONS / Roberto KERSEVAN / Marton ADY
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

#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
#define NOMINMAX
//#include <Windows.h>
#include <direct.h>
#include <process.h>
#include "SMP.h"
#else

#endif

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <istream>
#include <filesystem>
#include <cereal/archives/binary.hpp>
#include <cereal/types/utility.hpp>
//#include <cereal/archives/xml.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>

#include "MolflowGeometry.h"
#include "Worker.h"
#include "GLApp/GLApp.h"
#include "GLApp/GLMessageBox.h"

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

#include "ziplib/ZipArchive.h"
#include "ziplib/ZipArchiveEntry.h"
#include "ziplib/ZipFile.h"
#include "File.h" //File utils (Get extension, etc)
#include "ProcessControl.h"
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
extern MolFlow *mApp;
#endif

#if defined(SYNRAD)
extern SynRad*mApp;
#endif

/**
* \brief Default constructor for a worker
*/
Worker::Worker() : simManager("molflow", "MFLW"), model{} {

    //Molflow specific
    temperatures = std::vector<double>();
    desorptionParameterIDs = std::vector<size_t>();
    moments = std::vector<Moment>();
    userMoments = std::vector<UserMoment>(); //strings describing moments, to be parsed
    CDFs = std::vector<std::vector<CDF_p>>();
    IDs = std::vector<IntegratedDesorption>();
    parameters = std::vector<Parameter>();
    needsReload = true;  //When main and subprocess have different geometries, needs to reload (synchronize)
    displayedMoment = 0; //By default, steady-state is displayed

    ResetWorkerStats();
    geom = new MolflowGeometry();

    simuTimer.ReInit();
    //startTime = 0.0f;
    //stopTime = 0.0f;
    //simuTime = 0.0f;

    strcpy(fullFileName, "");
}

/**
* \brief Getter to retrieve the geometry on this worker
* \return pointer of the geometry object
*/
MolflowGeometry *Worker::GetMolflowGeometry() {
    return geom;
}

/**
* \brief Function for saving geometry to a set file
* \param fileName output file name with extension
* \param prg GLProgress window where loading is visualised
* \param askConfirm if a window should be opened to ask if the file should be overwritten
* \param saveSelected if a selection is to be saved
* \param autoSave if automatic saving is enabled
* \param crashSave if save on crash is enabled
*/
void Worker::SaveGeometry(std::string fileName, GLProgress *prg, bool askConfirm, bool saveSelected, bool autoSave,
                          bool crashSave) {

    try {
        if (needsReload && (!crashSave && !saveSelected)) RealReload();
    }
    catch (std::exception &e) {
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
    std::string path = FileUtils::GetPath(fileName);

    bool ok = true;
    if (ext.empty()) {
        fileName = fileName + (mApp->compressSavedFiles ? ".zip" : ".xml");
        ext = FileUtils::GetExtension(fileName);
        if (!autoSave && FileUtils::Exist(fileName)) {
            char tmp[1024];
            sprintf(tmp, "Overwrite existing file ?\n%s", fileName.c_str());
            if (askConfirm)
                ok = (GLMessageBox::Display(tmp, "Question", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONWARNING) == GLDLG_OK);
        }
    }

    // Read a file

    FileWriter *f = NULL;
    bool isTXT = Contains({"txt", "TXT"}, ext);
    bool isSTR = Contains({"str", "STR"}, ext);
    bool isGEO = ext == "geo";
    bool isGEO7Z = ext == "geo7z";
    bool isXML = ext == "xml";
    bool isXMLzip = ext == "zip";
    bool isSTL = ext == "stl";

    if (isTXT || isGEO || isGEO7Z || isSTR || isXML || isXMLzip || isSTL) {
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
        //Check (using native handle) if background compressor is still alive
        if ((isGEO7Z) && WAIT_TIMEOUT == WaitForSingleObject(mApp->compressProcessHandle, 0)) {
            GLMessageBox::Display("Compressing a previous save file is in progress. Wait until that finishes "
                "or close process \"compress.exe\"\nIf this was an autosave attempt,"
                "you have to lower the autosave frequency.", "Can't save right now.", GLDLG_OK, GLDLG_ICONERROR);
            return;
        }
#endif
        if (isGEO) {
            fileNameWithoutExtension = fileName.substr(0, fileName.length() - 4);
            fileNameWithGeo7z = fileName + "7z";
            fileNameWithGeo = fileName;
        } else if (isGEO7Z) {
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
            fileNameWithoutExtension = fileName.substr(0, fileName.length() - 4);
            fileNameWithXML = fileNameWithoutExtension + ".xml";
            fileNameWithZIP = fileNameWithoutExtension + ".zip";
        }
        if (isXMLzip) {
            char tmp[1024];
            sprintf(tmp, "An .xml file of the same name exists. Overwrite that file ?\n%s", fileNameWithZIP.c_str());
            if (!autoSave && FileUtils::Exist(fileNameWithXML)) {

                ok = (GLMessageBox::Display(tmp, "Question", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONWARNING) == GLDLG_OK);
            }
        }
        if (isSTL) {
            //Nothing to prepare
            ok = true;
        }

        if (!autoSave && ok && FileUtils::Exist(fileName)) {
            char tmp[1024];
            sprintf(tmp, "Overwrite existing file ?\n%s", fileName.c_str());
            if (askConfirm)
                ok = (GLMessageBox::Display(tmp, "Question", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONWARNING) == GLDLG_OK);
        }

        if (ok) {
            // Get copy of hit buffer, once a load can be initiated
            bool buffer_old;
            try {
                buffer_old = simManager.GetLockedHitBuffer();
                if(!buffer_old){
                    throw std::runtime_error("Error getting access to hit buffer.");
                }
            }
            catch (std::exception &e) {
                GLMessageBox::Display(e.what(), "Error getting access to hit buffer.", GLDLG_OK, GLDLG_ICONERROR);
                return;
            }

            if (isSTR) {
                geom->SaveSTR(saveSelected);
            } else {
                try {
                    if (isGEO7Z) {
                        f = new FileWriter(fileNameWithGeo); //We first write a GEO file, then compress it to GEO7Z later
                    } else if (!(isXML || isXMLzip))
                        f = new FileWriter(fileName); //Txt, stl, geo, etc...
                }
                catch (std::exception &e) {
                    SAFE_DELETE(f);
                    GLMessageBox::Display(e.what(), "Error writing file.", GLDLG_OK, GLDLG_ICONERROR);
                    return;
                }

                if (isTXT) {
                    geom->SaveTXT(f, globState, saveSelected);
                }
                else if (isGEO || isGEO7Z) {
                    /*
                    // Retrieve leak cache
                    int nbLeakSave, nbHHitSave;
                    LEAK leakCache[LEAKCACHESIZE];
                    if (!crashSave && !saveSelected) GetLeak(leakCache, &nbLeakSave);

                    // Retrieve hit cache (lines and dots)
                    HIT hitCache[HITCACHESIZE];
                    if (!crashSave && !saveSelected) GetHHit(hitCache, &nbHHitSave);
                    */
                    geom->SaveGEO(f, prg, globState, this, saveSelected, crashSave);
                } else if (isSTL) {
                    geom->SaveSTL(f, prg);
                } else if (isXML || isXMLzip) {
                    xml_document saveDoc;
                    geom->SaveXML_geometry(saveDoc, this, prg, saveSelected);
                    xml_document geom_only;
                    geom_only.reset(saveDoc);
                    bool success = false; //success: simulation state could be saved
                    if (!crashSave && !saveSelected) {
                        try {
                            //BYTE *buffer = simManager.GetLockedHitBuffer();
                            //GlobalHitBuffer *gHits = (GlobalHitBuffer *) buffer;
                            /*
                            int nbLeakSave, nbHHitSave;
                            LEAK leakCache[LEAKCACHESIZE];
                            GetLeak(leakCache, &nbLeakSave);
                            HIT hitCache[HITCACHESIZE];
                            GetHHit(hitCache, &nbHHitSave);
                            */

                            success = geom->SaveXML_simustate(saveDoc, this, globState, prg, saveSelected);
                            //SAFE_DELETE(buffer);
                        }
                        catch (std::exception &e) {
                            SAFE_DELETE(f);
                            simManager.UnlockHitBuffer();
                            GLMessageBox::Display(e.what(), "Error saving simulation state.", GLDLG_OK,
                                                  GLDLG_ICONERROR);
                            return;
                        }
                    }

                    prg->SetMessage("Writing xml file to disk...");
                    if (success) {
                        if (!saveDoc.save_file(fileNameWithXML.c_str()))
                            throw std::runtime_error("Error writing XML file."); //successful save
                    } else {
                        if (!geom_only.save_file(fileNameWithXML.c_str()))
                            throw std::runtime_error("Error writing XML file."); //simu state error
                    } //saveDoc goes out of scope for the compress duration
                    if (isXMLzip) {
                        prg->SetProgress(0.75);
                        prg->SetMessage("Compressing xml to zip...");
                        //mApp->compressProcessHandle=CreateThread(0, 0, ZipThreadProc, 0, 0, 0);

                        //Zipper library
                        if (FileUtils::Exist(fileNameWithZIP)) {
                            try {
                                remove(fileNameWithZIP.c_str());
                            }
                            catch (std::exception& e) {
                                SAFE_DELETE(f);
                                simManager.UnlockHitBuffer();
                                std::string msg = "Error compressing to \n" + fileNameWithZIP + "\nMaybe file is in use.";
                                GLMessageBox::Display(e.what(), msg.c_str(), GLDLG_OK, GLDLG_ICONERROR);
                                return;
                            }
                        }
                        ZipFile::AddFile(fileNameWithZIP, fileNameWithXML, FileUtils::GetFilename(fileNameWithXML));
                        //At this point, if no error was thrown, the compression is successful
                        try {
                            remove(fileNameWithXML.c_str());
                        }
                        catch (std::exception& e) {
                            SAFE_DELETE(f);
                            simManager.UnlockHitBuffer();
                            std::string msg = "Error removing\n" + fileNameWithXML + "\nMaybe file is in use.";
                            GLMessageBox::Display(e.what(), msg.c_str(), GLDLG_OK, GLDLG_ICONERROR);
                            return;
                        }
                    }
                }
                simManager.UnlockHitBuffer();
            }
            /*if (!autoSave && !saveSelected) {
                strcpy(fullFileName, fileName);
                remove("Molflow_AutoSave.zip");
                }*/
        }
    } else {
        SAFE_DELETE(f);
        throw std::runtime_error("SaveGeometry(): Invalid file extension [only xml,zip,geo,geo7z,txt,stl or str]");
    }

    SAFE_DELETE(f);

    //File written, compress it if the user wanted to
    if (ok && isGEO7Z) {

#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
        std::string compressorName = "compress.exe";
#else
        std::string compressorName = "./compress";
#endif

        if (FileUtils::Exist(compressorName)) { //compress GEO file to GEO7Z using 7-zip launcher "compress.exe"
            std::ostringstream tmp;
            tmp << compressorName << " \"" << fileNameWithGeo << "\" Geometry.geo";
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
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
        } else {
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
void Worker::ExportProfiles(const char *fn) {

    std::string fileName = fn;
    char tmp[512];

    // Read a file
    FILE *f = NULL;

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
        char tmp[256];
        sprintf(tmp, "Cannot open file for writing %s", fileName.c_str());
        throw std::runtime_error(tmp);
    }
    bool buffer_old = simManager.GetLockedHitBuffer();
    if(!buffer_old)
        throw std::runtime_error("Cannot access shared hit buffer");
    geom->ExportProfiles(f, isTXT, this);
    simManager.UnlockHitBuffer();
    fclose(f);

}

/**
* \brief Function for exportation of the angle maps
* \param fileName output file name
 * \param saveAll true if all files -- otherwise only selected -- should be saved
 * \return Vector with strings containing the file names of all angle map files
*/
std::vector<std::string> Worker::ExportAngleMaps(std::string fileName, bool saveAll) {
    bool overwriteAll = false;

    Geometry *geom = GetGeometry();
    std::vector<size_t> angleMapFacetIndices;
    for (size_t i = 0; i < geom->GetNbFacet(); i++) {
        InterfaceFacet *f = geom->GetFacet(i);
        // saveAll facets e.g. when auto saving or just when selected
        if ((saveAll || f->selected) && f->sh.anglemapParams.hasRecorded){
            angleMapFacetIndices.push_back(i);
        }
    }

    std::vector<std::string> listOfFiles;
    for (size_t facetIndex : angleMapFacetIndices) {
        std::string saveFileName;
        if (angleMapFacetIndices.size() == 1) {
            saveFileName = FileUtils::StripExtension(fileName) + ".csv";
        } else {
            std::stringstream tmp;
            tmp << FileUtils::StripExtension(fileName) << "_facet" << facetIndex + 1 << ".csv";
            saveFileName = tmp.str();
        }

        if (FileUtils::Exist(saveFileName)) {
            if (!overwriteAll) {
                std::vector<std::string> buttons = {"Cancel", "Overwrite"};
                if (angleMapFacetIndices.size() > 1) buttons.push_back("Overwrite All");
                int answer = GLMessageBox::Display("Overwrite existing file ?\n" + saveFileName, "Question", buttons,
                                                   GLDLG_ICONWARNING);
                if (answer == 0) break; //User cancel
                overwriteAll = (answer == 2);
            }
        }

        std::ofstream file;
        file.open(saveFileName);
        if (!file.is_open()) {
            std::string tmp = "Cannot open file for writing " + saveFileName;
            throw std::runtime_error(tmp.c_str());
        }
        file << geom->GetFacet(facetIndex)->GetAngleMap(1);
        file.close();
        listOfFiles.push_back(saveFileName);
    }

    return listOfFiles; // false if angleMapFacetIndices.size() == 0
}

[[maybe_unused]] bool Worker::ImportAngleMaps(std::string fileName) {

    for (auto &p : std::filesystem::directory_iterator("")) {
        std::stringstream fileName;
        fileName << p.path().string();
        if (FileUtils::GetExtension(fileName.str()) == "csv") std::cout << p.path() << '\n';
    }

    return true; // false if angleMapFacetIndices.size() == 0
}

/*void Worker::ImportDesorption(const char *fileName) {
	//if (needsReload) RealReload();

	// Read a file
	FileReader *f=new FileReader(fileName);
	geom->ImportDesorption(f,dpHit);
	SAFE_DELETE(f);
	changedSinceSave=true;
	Reload();
	}*/

/**
* \brief Function for loading of the geometry
* \param fileName input file name
* \param insert if geometry will be inserted into an existing geometry view
* \param newStr if a new structure needs to be created
*/
void Worker::LoadGeometry(const std::string &fileName, bool insert, bool newStr) {
    if (!insert) {
        needsReload = true;
    } else {
        try {
            RealReload();
        } catch(std::exception& e){
            GLMessageBox::Display(e.what(), "Error (Reloading Geometry)", GLDLG_OK, GLDLG_ICONERROR);
        }
    }
    //char CWD[MAX_PATH];
    //_getcwd(CWD, MAX_PATH);

    std::string ext = FileUtils::GetExtension(fileName);

    if (ext == "")

        throw std::runtime_error("LoadGeometry(): No file extension, can't determine type");

    // Read a file
    FileReader *f = NULL;
    GLProgress *progressDlg = new GLProgress("Reading file...", "Please wait");
    progressDlg->SetVisible(true);
    progressDlg->SetProgress(0.0);

    ResetWorkerStats();

    if (!insert) {
        //Clear hits and leaks cache
        ResetMoments();
        model.wp.globalHistogramParams = HistogramParams();

        //default values
        model.wp.enableDecay = false;
        model.wp.gasMass = 28;
    }

    /*
    bool isASE = (stricmp(ext, "ase") == 0);
    bool isSTR = (stricmp(ext, "str") == 0);
    bool isSTL = (stricmp(ext, "stl") == 0);
    bool isTXT = (stricmp(ext, "txt") == 0);
    bool isGEO = (stricmp(ext, "geo") == 0);
    bool isGEO7Z = (stricmp(ext, "geo7z") == 0);
    bool isSYN = (stricmp(ext, "syn") == 0);
    bool isSYN7Z = (stricmp(ext, "syn7z") == 0);
    bool isXML = (stricmp(ext, "xml") == 0);
    bool isXMLzip = (stricmp(ext, "zip") == 0);
    */

    if (ext == "txt" || ext == "TXT") {

        try {

            if (insert) mApp->changedSinceSave = true;

            f = new FileReader(fileName);

            if (!insert) {
                geom->LoadTXT(f, progressDlg, this);
                SAFE_DELETE(f);
                //RealReload();
                strcpy(fullFileName, fileName.c_str());
            } else { //insert

                geom->InsertTXT(f, progressDlg, newStr);
                SAFE_DELETE(f);
                Reload();
            }
        }

        catch (std::exception &e) {
            if (!insert) geom->Clear();
            SAFE_DELETE(f);
            progressDlg->SetVisible(false);
            SAFE_DELETE(progressDlg);
            throw e;
        }

    } else if (ext == "stl" || ext == "STL") {
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
            }
            if (ret != GLDLG_CANCEL_U) {
                progressDlg->SetMessage("Resetting worker...");
                progressDlg->SetVisible(true);
                progressDlg->SetMessage("Reading geometry...");
                f = new FileReader(fileName);
                if (!insert) {
                    geom->LoadSTL(f, progressDlg, scaleFactor);
                    SAFE_DELETE(f);
                    strcpy(fullFileName, fileName.c_str());
                    mApp->DisplayCollapseDialog();
                } else { //insert
                    mApp->changedSinceSave = true;
                    geom->InsertSTL(f, progressDlg, scaleFactor, newStr);
                    SAFE_DELETE(f);
                    Reload();
                }
            }
        }
        catch (std::exception &e) {
            if (!insert) geom->Clear();
            SAFE_DELETE(f);
            progressDlg->SetVisible(false);
            SAFE_DELETE(progressDlg);
            throw e;

        }

    } else if (ext == "str" || ext == "STR") {
        if (insert) throw std::runtime_error("STR file inserting is not supported.");
        try {
            f = new FileReader(fileName);
            progressDlg->SetVisible(true);
            geom->LoadSTR(f, progressDlg);
            SAFE_DELETE(f);
            //RealReload();

            strcpy(fullFileName, fileName.c_str());
        }

        catch (std::exception &e) {
            geom->Clear();
            SAFE_DELETE(f);
            progressDlg->SetVisible(false);
            SAFE_DELETE(progressDlg);
            throw e;
        }

    } else if (ext == "syn" || ext == "syn7z") { //Synrad file
        int version;
        progressDlg->SetVisible(true);
        try {
            if (ext == "syn7z") {
                //decompress file
                progressDlg->SetMessage("Decompressing file...");
                f = ExtractFrom7zAndOpen(fileName, "Geometry.syn");
            } else {
                f = new FileReader(fileName);  //original file opened
            }

            if (!insert) {
                progressDlg->SetMessage("Resetting worker...");

                geom->LoadSYN(f, progressDlg, &version, this);
                SAFE_DELETE(f);
                model.otfParams.desorptionLimit = 0;
            } else { //insert
                geom->InsertSYN(f, progressDlg, newStr);
                SAFE_DELETE(f);
            }

            progressDlg->SetMessage("Reloading worker with new geometry...");
            Reload();
            if (!insert) strcpy(fullFileName, fileName.c_str());
        }

        catch (std::exception &e) {
            if (!insert) geom->Clear();

            SAFE_DELETE(f);
            //if (isSYN7Z) remove(tmp2);

            progressDlg->SetVisible(false);
            SAFE_DELETE(progressDlg);
            throw e;
        }

    } else if (ext == "geo" || ext == "geo7z") {
        int version;
        progressDlg->SetVisible(true);
        try {
            if (ext == "geo7z") {
                //decompress file
                progressDlg->SetMessage("Decompressing file...");
                f = ExtractFrom7zAndOpen(fileName, "Geometry.geo");
            } else { //not geo7z
                f = new FileReader(fileName); //geo file, open it directly
            }
            if (insert) mApp->changedSinceSave = true;

            if (!insert) {

                geom->LoadGEO(f, progressDlg, &version, this);

                // Add moments only after user Moments are completely initialized
                if(TimeMoments::ParseAndCheckUserMoments(&moments, userMoments)){
                    GLMessageBox::Display("Overlap in time moments detected! Check in Moments Editor!", "Warning", GLDLG_OK, GLDLG_ICONWARNING);
                    return;
                }

                progressDlg->SetMessage("Reloading worker with new geometry...");
                RealReload(); //for the loading of textures

                bool buffer_old = simManager.GetLockedHitBuffer();
                if(!buffer_old)
                    throw std::runtime_error("Cannot access shared hit buffer");
                if (version >= 8)
                    geom->LoadProfileGEO(f, globState, version);
                simManager.UnlockHitBuffer();

                SendToHitBuffer(); //Global hit counters and hit/leak cache
                SendFacetHitCounts(); // From facetHitCache to dpHit's const.flow counter

                progressDlg->SetMessage("Loading textures...");
                LoadTexturesGEO(f, version);
                strcpy(fullFileName, fileName.c_str());
            } else { //insert
                mApp->changedSinceSave = true;
                geom->InsertGEO(f, progressDlg, newStr);
                Reload();
            }
            SAFE_DELETE(f);
        }

        catch (std::exception &e) {
            if (!insert) geom->Clear();
            SAFE_DELETE(f);
            //if (isGEO7Z) remove(tmp2);
            progressDlg->SetVisible(false);
            SAFE_DELETE(progressDlg);
            throw e;
        }

    } else if (ext == "xml" || ext == "zip") { //XML file, optionally in ZIP container
        xml_document loadXML;
        xml_parse_result parseResult;
        progressDlg->SetVisible(true);
        try {
            if (ext == "zip") { //compressed in ZIP container
                //decompress file
                progressDlg->SetMessage("Decompressing zip file...");

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

                        std::string tmpFileName = "tmp/" + zipFileName;
                        ZipFile::ExtractFile(fileName, zipFileName, tmpFileName);
                        progressDlg->SetMessage("Reading and parsing XML file...");

                        parseResult = loadXML.load_file(tmpFileName.c_str()); //load and parse it
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
            else {
                parseResult = loadXML.load_file(fileName.c_str()); //parse xml file directly
            }

            ResetWorkerStats();
            if (!parseResult) {
                //Parse error
                std::stringstream err;
                err << "XML parsed with errors.\n";
                err << "Error description: " << parseResult.description() << "\n";
                err << "Error offset: " << parseResult.offset << "\n";
                throw std::runtime_error(err.str().c_str());
            }

            progressDlg->SetMessage("Building geometry...");
            xml_node rootNode = loadXML;
            if(appVersionId >= 2680){
                xml_node envNode = loadXML.child("SimulationEnvironment");
                if(!envNode.empty())
                    rootNode = envNode;
            }

            if (!insert) {
                geom->LoadXML_geom(rootNode, this, progressDlg);
                // Add moments only after user Moments are completely initialized
                if(TimeMoments::ParseAndCheckUserMoments(&moments, userMoments)){
                    GLMessageBox::Display("Overlap in time moments detected! Check in Moments Editor!", "Warning", GLDLG_OK, GLDLG_ICONWARNING);
                    return;
                }

                geom->UpdateName(fileName.c_str());

                progressDlg->SetMessage("Reloading worker with new geometry...");
                try {
                    RealReload(); //To create the dpHit dataport for the loading of textures, profiles, etc...
                    strcpy(fullFileName, fileName.c_str());

                    if (ext == "xml" || ext == "zip")
                        progressDlg->SetMessage("Restoring simulation state...");
                    bool buffer_old = simManager.GetLockedHitBuffer();
                    if(!buffer_old)
                        throw std::runtime_error("Cannot access shared hit buffer");
                    geom->LoadXML_simustate(rootNode, globState, this, progressDlg);
                    RetrieveHistogramCache(); //So interface gets histogram data for disp.moment right after loading
                    simManager.UnlockHitBuffer();
                    SendToHitBuffer(); //Send global hits without sending facet counters, as they are directly written during the load process (mutiple moments)
                    SendFacetHitCounts(); //Send hits without sending facet counters, as they are directly written during the load process (mutiple moments)
                    RebuildTextures();
                }
                catch (std::exception &e) {
                    if (!mApp->profilePlotter) {
                        mApp->profilePlotter = new ProfilePlotter();
                        mApp->profilePlotter->SetWorker(this);
                    }
                    mApp->profilePlotter->Reset(); //To avoid trying to display non-loaded simulation results
                    if (!mApp->convergencePlotter) {
                        mApp->convergencePlotter = new ConvergencePlotter(this,mApp->formula_ptr);
                        mApp->convergencePlotter->SetWorker(this);
                    }
                    mApp->convergencePlotter->Reset(); //To avoid trying to display non-loaded simulation results
                    GLMessageBox::Display(e.what(), "Error while loading simulation state", GLDLG_CANCEL,
                                          GLDLG_ICONWARNING);
                }
            } else { //insert
                geom->InsertXML(rootNode, this, progressDlg, newStr);
                mApp->changedSinceSave = true;
                ResetWorkerStats();

                Reload();
            }
        }
        catch (std::exception &e) {
            if (!insert) geom->Clear();
            progressDlg->SetVisible(false);
            SAFE_DELETE(progressDlg);
            throw e;
        }

    } else if (ext == "ase" || ext == "ASE") {
        if (insert) throw std::runtime_error("ASE file inserting is not supported.");
        try {
            ResetWorkerStats();
            f = new FileReader(fileName);
            progressDlg->SetVisible(true);
            geom->LoadASE(f, progressDlg);
            SAFE_DELETE(f);
            //RealReload();
            strcpy(fullFileName, fileName.c_str());

        }
        catch (std::exception &e) {
            geom->Clear();
            SAFE_DELETE(f);
            progressDlg->SetVisible(false);
            SAFE_DELETE(progressDlg);
            throw e;
        }
    } else {
        progressDlg->SetVisible(false);
        SAFE_DELETE(progressDlg);
        throw std::runtime_error("LoadGeometry(): Invalid file extension [Only xml,zip,geo,geo7z,syn.syn7z,txt,ase,stl or str]");
    }
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

    globalHitCache = globState.globalHits;
    progressDlg->SetVisible(false);
    SAFE_DELETE(progressDlg);
    if (insert) {
        mApp->UpdateFacetlistSelected();
        mApp->UpdateViewers();
    }
}

/**
* \brief Function for loading textures from a GEO file
* \param f input file handle
* \param version version of the GEO data description
*/
void Worker::LoadTexturesGEO(FileReader *f, int version) {
    GLProgress *progressDlg = new GLProgress("Loading textures", "Please wait");
    progressDlg->SetProgress(0.0);
    try {
        bool buffer_old = simManager.GetLockedHitBuffer();
        if(!buffer_old)
            throw std::runtime_error("Cannot access shared hit buffer");
        progressDlg->SetVisible(true);
        geom->LoadTexturesGEO(f, progressDlg, globState, version);
        simManager.UnlockHitBuffer();
        RebuildTextures();
    }
    catch (std::exception &e) {
        char tmp[256];
        sprintf(tmp,
                "Couldn't load some textures. To avoid continuing a partially loaded state, it is recommended to reset the simulation.\n%s",
                e.what());
        GLMessageBox::Display(tmp, "Error while loading textures.", GLDLG_OK, GLDLG_ICONWARNING);
    }
    progressDlg->SetVisible(false);
    SAFE_DELETE(progressDlg);
}

/**
* \brief Function that updates various variables when stopping a simulation
* \param appTime current time of the application
*/
void Worker::InnerStop(float appTime) {
    simuTimer.Stop();
}

/**
* \brief Function that handles starting and stopping of the simulation
* \param appTime current time of the application
* \param sMode simulation mode (MC/AC)
*/
void Worker::StartStop(float appTime) {

    if (IsRunning()) {

        // Stop
        InnerStop(appTime);
        try {
            Stop();
            Update(appTime);
        }

        catch (std::exception &e) {
            GLMessageBox::Display(e.what(), "Error (Stop)", GLDLG_OK, GLDLG_ICONERROR);
            return;
        }
    } else {

        // Start
        try {
            if (needsReload) RealReload(); //Synchronize subprocesses to main process
            simuTimer.Start();
            Start();
        }
        catch (std::exception &e) {
            //isRunning = false;
            GLMessageBox::Display(e.what(), "Error (Start)", GLDLG_OK, GLDLG_ICONERROR);
            return;
        }

        // Particular case when simulation ends before getting RUN state
        if (simManager.allProcsDone) {
            Update(appTime);
            GLMessageBox::Display("Max desorption reached", "Information (Start)", GLDLG_OK, GLDLG_ICONINFO);
        }

    }
}

/**
* \brief Function that inserts a list of new paramters at the beginning of the catalog parameters
* \param newParams vector containing new parameters to be inserted
* \return index to insert position
*/
size_t Worker::InsertParametersBeforeCatalog(const std::vector<Parameter> &newParams) {
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
					f->facetHitCache=(*((FacetHitBuffer*)(buffer + f->model.wp.hitOffset+displayedMoment*sizeof(FacetHitBuffer))));
				}
				try {
					if (mApp->needsTexture || mApp->needsDirection) geom->BuildFacetTextures(buffer,mApp->needsTexture,mApp->needsDirection);
				}
				catch (Error &e) {
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

bool Worker::MolflowGeomToSimModel(){
    //auto geom = GetMolflowGeometry();
    // TODO: Proper clear call before for Real reload?
    model.structures.clear();
    model.vertices3.clear();

    for(int nbV = 0; nbV < geom->GetNbVertex(); ++nbV) {
        model.vertices3.emplace_back(*geom->GetVertex(nbV));
    }
    // Parse usermoments to regular moment intervals
    //model.tdParams.moments = this->moments;
    model.tdParams.CDFs = this->CDFs;
    //        model.tdParams.IDs = this->IDs;
    {
        model.tdParams.IDs.clear();
        for(auto& id : this->IDs){
            model.tdParams.IDs.push_back(id.GetValues());
        }
    }
    for(auto& param : this->parameters)
        model.tdParams.parameters.emplace_back(param);

    std::vector<Moment> momentIntervals;
    momentIntervals.reserve(this->moments.size());
    for(auto& moment : this->moments){
        momentIntervals.emplace_back(std::make_pair(moment.first - (0.5 * moment.second), moment.first + (0.5 * moment.second)));
    }
    model.tdParams.moments = momentIntervals;

    model.sh = *geom->GetGeomProperties();

    model.structures.resize(model.sh.nbSuper); //Create structures
    //TODO: Globalize Size values
    size_t angleMapTotalSize = 0;
    size_t dirTotalSize = 0;
    size_t profTotalSize = 0;
    size_t textTotalSize = 0;
    size_t histogramTotalSize = 0;

    bool hasVolatile = false;

    for (size_t facIdx = 0; facIdx < model.sh.nbFacet; facIdx++) {
        SubprocessFacet sFac;
        {
            InterfaceFacet* facet = geom->GetFacet(facIdx);

            //std::vector<double> outgMapVector(sh.useOutgassingFile ? sh.outgassingMapWidth*sh.outgassingMapHeight : 0);
            //memcpy(outgMapVector.data(), outgassingMapWindow, sizeof(double)*(sh.useOutgassingFile ? sh.outgassingMapWidth*sh.outgassingMapHeight : 0));
            size_t mapSize = facet->sh.anglemapParams.GetMapSize();
            std::vector<size_t> angleMapVector(mapSize);
            if(facet->angleMapCache.size() != facet->sh.anglemapParams.GetRecordedDataSize()){
                std::cerr << "Recorded Data Size is different from actual size: "<<facet->angleMapCache.size() << " / " << facet->sh.anglemapParams.GetRecordedDataSize() << std::endl;
                return false;
            }
            memcpy(angleMapVector.data(), facet->angleMapCache.data(), facet->sh.anglemapParams.GetRecordedDataSize());

            std::vector<double> textIncVector;

            // Add surface elements area (reciprocal)
            if (facet->sh.isTextured) {
                textIncVector.resize(facet->sh.texHeight * facet->sh.texWidth);
                if (facet->cellPropertiesIds) {
                    size_t add = 0;
                    for (size_t j = 0; j < facet->sh.texHeight; j++) {
                        for (size_t i = 0; i < facet->sh.texWidth; i++) {
                            double area = facet->GetMeshArea(add, true);

                            if (area > 0.0) {
                                // Use the sign bit to store isFull flag
                                textIncVector[add] = 1.0 / area;
                            } else {
                                textIncVector[add] = 0.0;
                            }
                            add++;
                        }
                    }
                } else {

                    double rw = facet->sh.U.Norme() / (double) (facet->sh.texWidthD);
                    double rh = facet->sh.V.Norme() / (double) (facet->sh.texHeightD);
                    double area = rw * rh;
                    size_t add = 0;
                    for (int j = 0; j < facet->sh.texHeight; j++) {
                        for (int i = 0; i < facet->sh.texWidth; i++) {
                            if (area > 0.0) {
                                textIncVector[add] = 1.0 / area;
                            } else {
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
            sFac.angleMap.pdf = angleMapVector;
            sFac.textureCellIncrements = textIncVector;
        }

        //Some initialization
        if (!sFac.InitializeOnLoad(facIdx, model.tdParams.moments.size())) return false;
        // Increase size counters
        //histogramTotalSize += 0;
        /*angleMapTotalSize += sFac.angleMapSize;
        dirTotalSize += sFac.directionSize* (1 + model.tdParams.moments.size());
        profTotalSize += sFac.profileSize* (1 + model.tdParams.moments.size());
        textTotalSize += sFac.textureSize* (1 + model.tdParams.moments.size());
*/
        hasVolatile |= sFac.sh.isVolatile;

        if ((sFac.sh.superDest || sFac.sh.isVolatile) && ((sFac.sh.superDest - 1) >= model.sh.nbSuper || sFac.sh.superDest < 0)) {
            // Geometry error
            //ClearSimulation();
            //ReleaseDataport(loader);
            std::ostringstream err;
            err << "Invalid structure (wrong link on F#" << facIdx + 1 << ")";
            //SetErrorSub(err.str().c_str());
            std::cerr << err.str() << std::endl;
            return false;
        }

        if (sFac.sh.superIdx == -1) { //Facet in all structures
            for (auto& s : model.structures) {
                s.facets.push_back(sFac);
            }
        }
        else {
            model.structures[sFac.sh.superIdx].facets.push_back(sFac); //Assign to structure
        }
    }
    return true;
}

/**
* \brief Function that reloads the whole simulation (resets simulation, rebuilds ray tracing etc) and synchronises subprocesses to main process
* \param sendOnly if only certain parts should be reloaded (geometry reloading / ray tracing tree)
*/
void Worker::RealReload(bool sendOnly) { //Sharing geometry with workers
    GLProgress *progressDlg = new GLProgress("Performing preliminary calculations on geometry...",
                                             "Passing Geometry to workers");
    progressDlg->SetVisible(true);
    progressDlg->SetProgress(0.0);

    if (!sendOnly) {
        if (model.otfParams.nbProcess == 0 && !geom->IsLoaded()) {
            progressDlg->SetVisible(false);
            SAFE_DELETE(progressDlg);
            return;
        }

        try {
            progressDlg->SetMessage("Do preliminary calculations...");
            PrepareToRun();

            size_t logDpSize = 0;
            if (model.otfParams.enableLogging) {
                logDpSize = sizeof(size_t) + model.otfParams.logLimit * sizeof(ParticleLoggerItem);
            }
            size_t hitSize = geom->GetHitsSize(moments.size());

            progressDlg->SetMessage("Asking subprocesses to clear geometry...");
            simManager.ResetSimulations();
            progressDlg->SetMessage("Creating Logger...");
            simManager.ReloadLogBuffer(logDpSize, true);
            progressDlg->SetMessage("Creating hit buffer...");
            simManager.ReloadHitBuffer(hitSize);
        }
        catch (std::exception &e) {
            GLMessageBox::Display(e.what(), "Error (Full reload)", GLDLG_OK, GLDLG_ICONWARNING);
            progressDlg->SetVisible(false);
            SAFE_DELETE(progressDlg);
            std::stringstream err;
            err << "Error (Full reload) " << e.what();
            throw std::runtime_error(err.str());
        }
    }

    // Send and Load geometry
    //std::string loaderString = SerializeForLoader().str();
    progressDlg->SetMessage("Waiting for subprocesses to load geometry...");
    try {
        if(!MolflowGeomToSimModel()){
        //if (simManager.ShareWithSimUnits((BYTE *) loaderString.c_str(), loaderString.size(), LoadType::LOADGEOM)) {
            std::string errString = "Failed to send geometry to sub process!\n";
            GLMessageBox::Display(errString.c_str(), "Warning (LoadGeom)", GLDLG_OK, GLDLG_ICONWARNING);

            progressDlg->SetVisible(false);
            SAFE_DELETE(progressDlg);
            return;
        }
        simManager.ForwardSimModel(&model);

        if (simManager.ExecuteAndWait(COMMAND_LOAD, PROCESS_READY, 0, 0)) {
            //CloseLoaderDP();
            std::string errString = "Failed to send geometry to sub process:\n";
            errString.append(GetErrorDetails());
            throw std::runtime_error(errString);
        }
    }
    catch (std::exception& e) {
        GLMessageBox::Display(e.what(), "Error (LoadGeom)", GLDLG_OK, GLDLG_ICONERROR);
    }

    progressDlg->SetMessage("Constructing memory structure to store results...");
    if (!sendOnly) {
        globState.Resize(model);
    }

    //Old send hits location
    progressDlg->SetMessage("Closing dataport...");
    needsReload = false;
    progressDlg->SetVisible(false);
    SAFE_DELETE(progressDlg);
}

/**
* \brief Serialization function for a binary cereal archive for the worker attributes
* \return output string stream containing the result of the archiving
*/
std::ostringstream Worker::SerializeParamsForLoader() {
    std::ostringstream result;
    cereal::BinaryOutputArchive outputArchive(result);

    outputArchive(
            CEREAL_NVP(model.otfParams)
    );
    return result;
}

/**
* \brief Resets workers global hit cache
* Reset function mainly used for initialisation / reload procedures
*/
void Worker::ResetWorkerStats() {

    globState.Reset();
    //memset(&globState.globalHits, 0, sizeof(GlobalHitBuffer));


}

/**
* \brief Starts the simulation process
*/
void Worker::Start() {
    // Sanity checks
    // Is there some desorption in the system? (depends on pre calculation)
    if(model.wp.finalOutgassingRate_Pa_m3_sec <= 0.0){
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
    if (model.wp.totalDesorbedMolecules <= 0.0)
        throw std::runtime_error("Total outgassing is zero.");

    try {
        simManager.ForwardGlobalCounter(&globState);

        if (simManager.StartSimulation()) {
            throw std::logic_error("Processes are already done!");
        }
    }
    catch (std::exception& e) {
        throw e;
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
double Worker::GetMoleculesPerTP(size_t moment) {
    if (globalHitCache.globalHits.hit.nbDesorbed == 0) return 0; //avoid division by 0
    if (moment == 0) {
        //Constant flow
        //Each test particle represents a certain real molecule influx per second
        return model.wp.finalOutgassingRate / globalHitCache.globalHits.hit.nbDesorbed;
    } else {
        //Time-dependent mode
        //Each test particle represents a certain absolute number of real molecules
        return (model.wp.totalDesorbedMolecules / mApp->worker.moments[moment - 1].second) / globalHitCache.globalHits.hit.nbDesorbed;
    }
}

/* //Commenting out as deprecated
void Worker::ImportDesorption_DES(const char *fileName) {
	//if (needsReload) RealReload();
	// Read a file
	FileReader *f = new FileReader(fileName);
	geom->ImportDesorption_DES(f);
	SAFE_DELETE(f);
	mApp->changedSinceSave = true;
	Reload();
}
*/

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
* \param prg GLProgress window where visualising of the import progress is shown
*/
void Worker::ImportDesorption_SYN(const char *fileName, const size_t &source, const double &time,
                                  const size_t &mode, const double &eta0, const double &alpha, const double &cutoffdose,
                                  const std::vector<std::pair<double, double>> &convDistr,
                                  GLProgress *prg) {
    std::string ext = FileUtils::GetExtension(fileName);
    if (!Contains({"syn7z", "syn"}, ext))
        throw std::runtime_error("ImportDesorption_SYN(): Invalid file extension [Only syn, syn7z]");

    // Read a file

    FileReader *f = nullptr;

    GLProgress *progressDlg = new GLProgress("Analyzing SYN file...", "Please wait");
    progressDlg->SetProgress(0.0);
    progressDlg->SetVisible(true);
    bool isSYN7Z = (iequals(ext, "syn7z"));
    bool isSYN = (iequals(ext, "syn"));

    if (isSYN || isSYN7Z) {
        progressDlg->SetVisible(true);
        try {
            if (isSYN7Z) {
                f = ExtractFrom7zAndOpen(fileName, "Geometry.syn");
            } else {
                f = new FileReader(fileName);  //original file opened
            }

            geom->ImportDesorption_SYN(f, source, time, mode, eta0, alpha, cutoffdose, convDistr, prg);
            CalcTotalOutgassing();
            SAFE_DELETE(f);

        }
        catch (std::exception &e) {

            SAFE_DELETE(f);
            progressDlg->SetVisible(false);
            SAFE_DELETE(progressDlg);
            throw e;
        }
        progressDlg->SetVisible(false);
        SAFE_DELETE(progressDlg);
    }
}

/**
* \brief To analyse desorption data from a SYN file
* \param fileName name of the input file
* \param nbFacet number of facets in the file
* \param nbTextured number of textured facets in the file
* \param nbDifferent (TODO: check usage)
*/
void Worker::AnalyzeSYNfile(const char *fileName, size_t *nbFacet, size_t *nbTextured, size_t *nbDifferent) {
    std::string ext = FileUtils::GetExtension(fileName);
    bool isSYN = (ext == "syn") || (ext == "SYN");
    bool isSYN7Z = (ext == "syn7z") || (ext == "SYN7Z");

    if (!(isSYN || isSYN7Z))
        throw std::runtime_error("AnalyzeSYNfile(): Invalid file extension [Only syn, syn7z]");

    // Read a file
    FileReader *f = nullptr;

    GLProgress *progressDlg = new GLProgress("Analyzing SYN file...", "Please wait");
    progressDlg->SetProgress(0.0);
    progressDlg->SetVisible(true);

    if (isSYN || isSYN7Z) {
        progressDlg->SetVisible(true);
        try {
            if (isSYN7Z) {
                //decompress file
                progressDlg->SetMessage("Decompressing file...");
                f = ExtractFrom7zAndOpen(fileName, "Geometry.syn");
            } else {
                f = new FileReader(fileName);  //original file opened
            }

            geom->AnalyzeSYNfile(f, progressDlg, nbFacet, nbTextured, nbDifferent, progressDlg);

            SAFE_DELETE(f);

        }
        catch (std::exception &e) {
            SAFE_DELETE(f);
            progressDlg->SetVisible(false);
            SAFE_DELETE(progressDlg);
            throw e;
        }
        progressDlg->SetVisible(false);
        SAFE_DELETE(progressDlg);

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
    model.wp.latestMoment = model.wp.timeWindowSize * .5;

    if(!moments.empty())
        model.wp.latestMoment = (moments.end()-1)->first + (moments.end()-1)->second / 2.0;
    //model.wp.latestMoment += model.wp.timeWindowSize / 2.0;

    Geometry *g = GetGeometry();
    //Generate integrated desorption functions

    temperatures = std::vector<double>();
    desorptionParameterIDs = std::vector<size_t>();
    CDFs = std::vector<std::vector<CDF_p>>();
    IDs = std::vector<IntegratedDesorption>();

    bool needsAngleMapStatusRefresh = false;

    for (size_t i = 0; i < g->GetNbFacet(); i++) {
        InterfaceFacet *f = g->GetFacet(i);

        //match parameters
        if (f->userOutgassing.length() > 0) {
            int id = GetParamId(f->userOutgassing);
            if (id == -1) { //parameter not found
                char tmp[256];
                sprintf(tmp, "Facet #%zd: Outgassing parameter \"%s\" isn't defined.", i + 1,
                        f->userOutgassing.c_str());
                throw std::runtime_error(tmp);
            } else f->sh.outgassing_paramId = id;
        } else f->sh.outgassing_paramId = -1;

        if (f->userOpacity.length() > 0) {
            int id = GetParamId(f->userOpacity);
            if (id == -1) { //parameter not found
                char tmp[256];
                sprintf(tmp, "Facet #%zd: Opacity parameter \"%s\" isn't defined.", i + 1, f->userOpacity.c_str());
                throw std::runtime_error(tmp);
            } else f->sh.opacity_paramId = id;
        } else f->sh.opacity_paramId = -1;

        if (f->userSticking.length() > 0) {
            int id = GetParamId(f->userSticking);
            if (id == -1) { //parameter not found
                char tmp[256];
                sprintf(tmp, "Facet #%zd: Sticking parameter \"%s\" isn't defined.", i + 1, f->userSticking.c_str());
                throw std::runtime_error(tmp);
            } else f->sh.sticking_paramId = id;
        } else f->sh.sticking_paramId = -1;

        if (f->sh.outgassing_paramId >= 0) { //if time-dependent desorption
            int id = GetIDId(f->sh.outgassing_paramId);
            if (id >= 0)
                f->sh.IDid = id; //we've already generated an integrated des. for this time-dep. outgassing
            else
                f->sh.IDid = GenerateNewID(f->sh.outgassing_paramId); //Convert timedep outg. (PDF) to CDF
        }

        //Generate speed distribution functions
        int id = GetCDFId(f->sh.temperature);
        if (id >= 0)
            f->sh.CDFid = id; //we've already generated a CDF for this temperature
        else
            f->sh.CDFid = GenerateNewCDF(f->sh.temperature);

        //Angle map
        if (f->sh.desorbType == DES_ANGLEMAP) {
			if (!f->sh.anglemapParams.hasRecorded) {
                char tmp[256];
                sprintf(tmp, "Facet #%zd: Uses angle map desorption but doesn't have a recorded angle map.", i + 1);
                throw std::runtime_error(tmp);
            }
            if (f->sh.anglemapParams.record) {
                char tmp[256];
                sprintf(tmp, "Facet #%zd: Can't RECORD and USE angle map desorption at the same time.", i + 1);
                throw std::runtime_error(tmp);
            }
        }

        //First worker::update will do it
		if (f->sh.anglemapParams.record) {
			if (!f->sh.anglemapParams.hasRecorded) {
                //Initialize angle map and Set values to zero
                try {
                    f->angleMapCache.resize(f->sh.anglemapParams.GetMapSize(), 0);
                }
                catch(...) {
                    std::stringstream tmp;
                    tmp << "Not enough memory for incident angle map on facet " << i + 1;
                    throw std::runtime_error(tmp.str().c_str());
                }
                f->sh.anglemapParams.hasRecorded = true;
                if (f->selected) needsAngleMapStatusRefresh = true;
            }
		}

    }

    if (mApp->facetAdvParams && mApp->facetAdvParams->IsVisible() && needsAngleMapStatusRefresh)
        mApp->facetAdvParams->Refresh(geom->GetSelectedFacets());

    CalcTotalOutgassing();

}

/**
* \brief Get ID (if it exists) of the Commulative Distribution Function (CFD) for a particular temperature (bin)
* \param temperature temperature for the CFD
* \return ID of the CFD
*/
int Worker::GetCDFId(double temperature) {

    int i;
    for (i = 0; i < (int) temperatures.size() &&
                (std::abs(temperature - temperatures[i]) > 1E-5); i++); //check if we already had this temperature
    if (i >= (int) temperatures.size()) i = -1; //not found
    return i;
}

/**
* \brief Generate a new Commulative Distribution Function (CFD) for a particular temperature (bin)
* \param temperature for the CFD
* \return Previous size of temperatures vector, which determines new ID
*/
int Worker::GenerateNewCDF(double temperature) {
    size_t i = temperatures.size();
    temperatures.push_back(temperature);
    CDFs.push_back(Generate_CDF(temperature, model.wp.gasMass, CDF_SIZE));
    return (int) i;
}

/**
* \brief Generate a new ID (integrated desorption) for desorption parameter for time-dependent simulations
* \param paramId parameter ID
* \return Previous size of IDs vector, which determines new id in the vector
*/
int Worker::GenerateNewID(int paramId) {
    //This function is called if parameter with index paramId doesn't yet have a cumulative des. function
    size_t i = desorptionParameterIDs.size();
    desorptionParameterIDs.push_back(paramId); //mark that i.th integrated des. belongs to paramId
    IDs.push_back(Generate_ID(paramId)); //actually convert PDF to CDF
    return (int) i; //return own index
}

/**
* \brief Get ID (if it exists) of the integrated desorption (ID) function for a particular paramId
* \param paramId parameter ID
* \return Id of the integrated desorption function
*/
int Worker::GetIDId(int paramId) {

    int i;
    for (i = 0; i < (int) desorptionParameterIDs.size() &&
                (paramId != desorptionParameterIDs[i]); i++); //check if we already had this parameter Id
    if (i >= (int) desorptionParameterIDs.size()) i = -1; //not found
    return i;

}

/**
* \brief Compute the outgassing of all source facet depending on the mode (file, regular, time-dependent) and set it to the global settings
*/
void Worker::CalcTotalOutgassing() {
    // Compute the outgassing of all source facet
    model.wp.totalDesorbedMolecules = model.wp.finalOutgassingRate_Pa_m3_sec = model.wp.finalOutgassingRate = 0.0;
    Geometry *g = GetGeometry();

    for (int i = 0; i < g->GetNbFacet(); i++) {
        InterfaceFacet *f = g->GetFacet(i);
        if (f->sh.desorbType != DES_NONE) { //there is a kind of desorption
            if (f->sh.useOutgassingFile) { //outgassing file
                for (int l = 0; l < (f->ogMap.outgassingMapWidth * f->ogMap.outgassingMapHeight); l++) {
                    model.wp.totalDesorbedMolecules += model.wp.latestMoment * f->ogMap.outgassingMap[l] / (1.38E-23 * f->sh.temperature);
                    model.wp.finalOutgassingRate += f->ogMap.outgassingMap[l] / (1.38E-23 * f->sh.temperature);
                    model.wp.finalOutgassingRate_Pa_m3_sec += f->ogMap.outgassingMap[l];
                }
            } else { //regular outgassing
                if (f->sh.outgassing_paramId == -1) { //constant outgassing
                    model.wp.totalDesorbedMolecules += model.wp.latestMoment * f->sh.outgassing / (1.38E-23 * f->sh.temperature);
                    model.wp.finalOutgassingRate +=
                            f->sh.outgassing / (1.38E-23 * f->sh.temperature);  //Outgassing molecules/sec
                    model.wp.finalOutgassingRate_Pa_m3_sec += f->sh.outgassing;
                } else { //time-dependent outgassing
                    double lastValue = IDs[f->sh.IDid].GetValues().back().second;
                    model.wp.totalDesorbedMolecules += lastValue / (1.38E-23 * f->sh.temperature);
                    size_t lastIndex = parameters[f->sh.outgassing_paramId].GetSize() - 1;
                    double finalRate_mbar_l_s = parameters[f->sh.outgassing_paramId].GetY(lastIndex);
                    model.wp.finalOutgassingRate +=
                            finalRate_mbar_l_s * MBARLS_TO_PAM3S / (1.38E-23 * f->sh.temperature); //0.1: mbar*l/s->Pa*m3/s
                    model.wp.finalOutgassingRate_Pa_m3_sec += finalRate_mbar_l_s * MBARLS_TO_PAM3S;
                }
            }
        }
    }
    if (mApp->globalSettings) mApp->globalSettings->UpdateOutgassing();

}

/**
* \brief Generate cumulative distribution function (CFD) for the velocity
* \param gasTempKelvins gas temperature in Kelvin
* \param gasMassGramsPerMol molar gas mass in grams per mol
* \param size amount of points/bins of the CFD
* \return CFD as a Vector containing a pair of double values (x value = speed_bin, y value = cumulated value)
*/
std::vector<std::pair<double, double>>
Worker::Generate_CDF(double gasTempKelvins, double gasMassGramsPerMol, size_t size) {
    std::vector<std::pair<double, double>> cdf;
    cdf.reserve(size);
    constexpr double Kb = 1.38E-23;
    constexpr double R = 8.3144621;
    const double a = sqrt(Kb * gasTempKelvins /
                    (gasMassGramsPerMol * 1.67E-27)); //distribution a parameter. Converting molar mass to atomic mass

    //Generate cumulative distribution function
    double mostProbableSpeed = sqrt(2 * R * gasTempKelvins / (gasMassGramsPerMol / 1000.0));
    double binSize = 4.0 * mostProbableSpeed / (double) size; //distribution generated between 0 and 4*V_prob
    /*double coeff1=1.0/sqrt(2.0)/a;
    double coeff2=sqrt(2.0/PI)/a;
    double coeff3=1.0/(2.0*pow(a,2));

    for (size_t i=0;i<size;i++) {
    double x=(double)i*binSize;
    cdf.push_back(std::make_pair(x,erf(x*coeff1)-coeff2*x*exp(-pow(x,2)*coeff3)));
    }*/
    for (size_t i = 0; i < size; i++) {
        double x = (double) i * binSize;
        double x_square_per_2_a_square = pow(x, 2) / (2 * pow(a, 2));
        cdf.emplace_back(std::make_pair(x, 1 - exp(-x_square_per_2_a_square) * (x_square_per_2_a_square + 1)));

    }

    /* //UPDATE: not generating inverse since it was introducing sampling problems at the large tail for high speeds
    //CDF created, let's generate its inverse
    std::vector<std::pair<double,double>> inverseCDF;inverseCDF.reserve(size);
    binSize=1.0/(double)size; //Divide probability to bins
    for (size_t i=0;i<size;i++) {
    double p=(double)i*binSize;
    //inverseCDF.push_back(std::make_pair(p,InterpolateX(p,cdf,true)));
    inverseCDF.push_back(std::make_pair(p, InterpolateX(p, cdf, false)));

    }
    return inverseCDF;
    */
    return cdf;
}

/**
* \brief Generate integrated (cumulative) desorption (ID) function from time-dependent outgassing, for inverse lookup at RNG
* \param paramId outgassing parameter index (parameters can be lin/lin, log-lin, lin-log or log-log, each requiring different integration method)
* \return Integrated desorption: class containing a Vector containing a pair of double values (x value = time moment, y value = cumulative desorption value)
*/
IntegratedDesorption Worker::Generate_ID(int paramId) {

    std::vector<std::pair<double, double>> ID; //time-cumulative desorption pairs. Can have more points than the corresponding outgassing (some sections are divided to subsections)

    //We need to slightly modify the original outgassing:
    //Beginning: Add a point at t=0, with the first outgassing value (assuming constant outgassing from 0 to t1 - const. extrap.)
    //End: if latestMoment is after the last user-defined moment, copy last user value to latestMoment (constant extrapolation)
    //     if latestMoment is before, create a new point at latestMoment with interpolated outgassing value and throw away the rest (don't generate particles after latestMoment)
    std::vector<std::pair<double,double>> myOutgassing; //modified parameter

    //First, let's check at which index is the latest moment
    size_t indexAfterLatestMoment;
    Parameter& par = parameters[paramId]; //we'll reference it a lot
    for (indexAfterLatestMoment = 0; indexAfterLatestMoment < par.GetSize() &&
                                    (par.GetX(indexAfterLatestMoment) <
                                     model.wp.latestMoment); indexAfterLatestMoment++); //loop exits after first index after latestMoment
    bool lastUserMomentBeforeLatestMoment;
    if (indexAfterLatestMoment >= par.GetSize()) {
        indexAfterLatestMoment = par.GetSize() - 1; //not found, set as last moment
        lastUserMomentBeforeLatestMoment = true;
    } else {
        lastUserMomentBeforeLatestMoment = false;
    }

    //Construct integral from 0 to the simulation's latest moment
    //First point: t=0, Q(0)=Q(t0)
    ID.push_back(std::make_pair(0.0, 0.0)); //At t=0 no particles have desorbed yet
    if (par.GetX(0)>0.0) { //If the user starts later than t=0, copy first user value to t=0(const extrapolation)
        myOutgassing.push_back(std::make_pair(0.0,par.GetY(0)));
    } else { //user has set an outgassing at t=0, so just copy it to myOutgassing
        myOutgassing.push_back(par.GetValues()[0]);
    }
    //Consecutive points: user-defined points that are before latestMoment
    //We throw away user-defined moments after latestMoment:
    //Example: we sample the system until t=10s but outgassing is defined until t=1000s -> ignore values after 10s
    {
        const auto& valuesCopy = par.GetValues();
        if (lastUserMomentBeforeLatestMoment) {
            myOutgassing.insert(myOutgassing.end(),valuesCopy.begin(),valuesCopy.end()); //copy all values...
            myOutgassing.push_back(std::make_pair(model.wp.latestMoment,myOutgassing.back().second)); //...and create last point equal to last outgassing (const extrapolation)
        } else {
             if (indexAfterLatestMoment>0) {
                myOutgassing.insert(myOutgassing.end(),valuesCopy.begin(),valuesCopy.begin()+indexAfterLatestMoment-1); //copy values that are before latestMoment
             }
             if (!IsEqual(myOutgassing.back().first,model.wp.latestMoment)) { //if interpolation is needed
                myOutgassing.push_back(std::make_pair(model.wp.latestMoment,
                InterpolateY(model.wp.latestMoment,valuesCopy,par.logXinterp,par.logYinterp))); //interpolate outgassing value to t=latestMoment
            }
        }

    } //valuesCopy goes out of scope

    //Intermediate moments, from first to t=latestMoment
    for (size_t i = 1; i < myOutgassing.size(); i++) { //myOutgassing[0] is always at t=0, skipping
        if (IsEqual(myOutgassing[i].second, myOutgassing[i - 1].second)) {
            //easy case of two equal y0=y1 values, simple integration by multiplying, reducing number of points
            ID.push_back(std::make_pair(myOutgassing[i].first,
                ID.back().second +
                (myOutgassing[i].first - myOutgassing[i - 1].first) *
                myOutgassing[i].second * MBARLS_TO_PAM3S)); //integral = y1*(x1-x0)
        }
        else { //we need to split the user-defined section to 10 subsections, and integrate each
               //(terminology: section is between two consecutive user-defined time-value pairs)
            const int nbSteps = 10; //change here
            double sectionStartTime = myOutgassing[i - 1].first;
            double sectionEndTime = myOutgassing[i].first;
            double sectionTimeInterval = sectionEndTime - sectionStartTime;
            double sectionDelta = 1.0 / nbSteps;
            double subsectionTimeInterval,subsectionLogTimeInterval;
            double logSectionStartTime, logSectionEndTime, logSectionTimeInterval;
            if (par.logXinterp) { //time is logarithmic
                if (sectionStartTime == 0) logSectionStartTime = -99;
                else logSectionStartTime = log10(sectionStartTime);
                logSectionEndTime = log10(sectionEndTime);
                logSectionTimeInterval = logSectionEndTime - logSectionStartTime;
                subsectionLogTimeInterval = logSectionTimeInterval * sectionDelta;
            }
            else {
                subsectionTimeInterval = sectionTimeInterval * sectionDelta;
            }
            double previousSubsectionValue = myOutgassing[i - 1].second; //Points to start of section, will be updated to point to start of subsections
            double previousSubsectionTime = sectionStartTime;
            for (double sectionFraction = sectionDelta; sectionFraction < 1.0001; sectionFraction += sectionDelta) { //sectionFraction: where we are within the section [0..1]
                double subSectionEndTime;
                if (!par.logXinterp) {
                    //linX: Distribute sampled points evenly
                    subSectionEndTime = sectionStartTime + sectionFraction * sectionTimeInterval;
                }
                else {
                    //logX: Distribute sampled points logarithmically
                    subSectionEndTime = Pow10(logSectionStartTime + sectionFraction * logSectionTimeInterval);
                }
                double subSectionEndValue = InterpolateY(subSectionEndTime,myOutgassing,par.logXinterp,par.logYinterp);
                double subsectionDesorbedGas; //desorbed gas in this subsection, in mbar*l/s, will be converted to Pa*m3/s when adding to ID
                if (!par.logXinterp && !par.logYinterp) { //lin-lin interpolation
                    //Area under a straight section from (x0,y0) to (x1,y1) on a lin-lin plot: I = (x1-x0) * (y0+y1)/2
                    subsectionDesorbedGas = subsectionTimeInterval * 0.5 * (previousSubsectionValue + subSectionEndValue);
                }
                else if (par.logXinterp && !par.logYinterp) { //log-lin: time (X) is logarithmic, outgassing (Y) is linear
                    //From Mathematica/WolframAlpha: integral of a straight section from (x0,y0) to (x1,y1) on a log10-lin plot: I = (y1-y0)(x0-x1)/ln(x1/x0) - x0y0 + x1y1
                    double x0=previousSubsectionTime;
                    double x1=subSectionEndTime;
                    double y0=previousSubsectionValue;
                    double y1=subSectionEndValue;
                    subsectionDesorbedGas = (y1-y0)*(x0-x1)/log(x1/x0) - x0*y0 + x1*y1;
                }
                else if (!par.logXinterp && par.logYinterp) { //lin-log: time (X) is linear, outgassing (Y) is logarithmic
                    //Area under a straight section from (x0,y0) to (x1,y1) on a lin-log plot: I = 1/m * (y1-y0) where m=log10(y1/y0)/(x1-x0)
                    double logSubSectionEndValue = (subSectionEndValue>0.0) ? log10(subSectionEndValue) : -99;
                    double logPreviousSubSectionValue = (previousSubsectionTime>0.0) ? log10(previousSubsectionValue) : -99;
                    //I = (x2-x1)*(y2-y1)*log10(exp(1))/(log10(y2)-log10(y1)) from https://fr.mathworks.com/matlabcentral/answers/404930-finding-the-area-under-a-semilogy-plot-with-straight-line-segments
                    subsectionDesorbedGas = subsectionTimeInterval*(subSectionEndValue-previousSubsectionValue)
                        * log10(exp(1))/(logSubSectionEndValue-logPreviousSubSectionValue);
                }
                else { //log-log
                    //Area under a straight section from (x0,y0) to (x1,y1) on a log-log plot:
                    //I = (y0/(x0^m))/(m+1) * (x1^(m+1)-x0^(m+1)) if m!=-1
                    //I = x0*y0*log10(x1/x0) if m==-1
                    //where m= log10(y1/y0) / log10(x1/x0)
                    double logSubSectionEndValue = (subSectionEndValue>0.0) ? log10(subSectionEndValue) : -99;
                    double logPreviousSubSectionValue = (previousSubsectionTime>0.0) ? log10(previousSubsectionValue) : -99;
                    double m = (logSubSectionEndValue-logPreviousSubSectionValue) / subsectionLogTimeInterval; //slope
                    if (m != -1.0) {
                        subsectionDesorbedGas = previousSubsectionValue / pow(previousSubsectionTime, m) / (m + 1.0) * (pow(subSectionEndTime, m + 1.0) - pow(previousSubsectionTime, m + 1.0));
                    }
                    else { //m==-1
                        subsectionDesorbedGas = previousSubsectionTime * previousSubsectionValue * subsectionLogTimeInterval;
                    }
                }

                ID.push_back(std::make_pair(subSectionEndTime,ID.back().second + subsectionDesorbedGas * MBARLS_TO_PAM3S ));

                //Cache end values for next iteration
                previousSubsectionValue = subSectionEndValue;
                previousSubsectionTime = subSectionEndTime;
            } //end subsection loop
        }
    } //end section loop

    IntegratedDesorption result;
    result.logXinterp=par.logXinterp;
    result.logYinterp=par.logYinterp;
    result.SetValues(ID, false);
    return result;

}

/**
* \brief Get ID of a parameter (if it exists) for a corresponding name
* \param name name of the parameter that shall be looked up
* \return ID corresponding to the found parameter
*/
int Worker::GetParamId(const std::string name) {
    int foundId = -1;
    for (int i = 0; foundId == -1 && i < (int) parameters.size(); i++)
        if (name.compare(parameters[i].name) == 0) foundId = i;
    return foundId;
}