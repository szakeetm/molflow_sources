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

#include "Interface/Interface.h"

class Worker;
class ImportDesorption;
class TimeSettings;
class Movement;
class MeasureForce;
class FacetAdvParams;
class FacetDetails;
class Viewer3DSettings;
class TextureScaling;
class GlobalSettings;
class ProfilePlotter;
class PressureEvolution;
class TimewisePlotter;
class TexturePlotter;
class OutgassingMapWindow;
class MomentsEditor;
class ParameterEditor;

struct Error;

class MolFlow : public Interface {
public:
    MolFlow();

	//Public textfields so we can disable them from "Advanced facet parameters":
	GLTextField   *facetOutgassingText;
	GLTextField   *facetOutgPerAreaText;

    void LoadFile(const std::string &fileName) override;
	void InsertGeometry(bool newStr, const std::string &fileName) override;
	void SaveFile() override;
    
	//void ImportDesorption_DES(); //Deprecated
	void ExportProfiles();
	void ExportAngleMaps();
	void ImportAngleMaps();
	void CopyAngleMapToClipboard();
	void ClearAngleMapsOnSelection();
    void ClearFacetParams() override;
	
    void ApplyFacetParams();
	void UpdateFacetParams(bool updateSelection) override;
    void StartStopSimulation();
    void SaveConfig() override;
    void LoadConfig() override;
    
	void PlaceComponents() override;
    void UpdateFacetHits(bool allRows) override;
	void UpdatePlotters() override;
	void RefreshPlotterCombos();

	//Flow/sticking coeff. conversion
	void CalcPumpingSpeed();
	void CalcSticking();

    GLTextField   *facetStickingText;
	
    GLCombo       *facetDesType;
	GLTextField   *facetDesTypeN;
    GLCombo       *facetProfileCombo;

	GLLabel       *facetPumpingLabel;
	GLTextField   *facetPumpingSpeedText;	
    GLLabel       *facetStickingLabel;
    
	GLLabel       *facetTempLabel;
	GLTextField   *facetTemperatureText;
    GLLabel       *facetDesTypeLabel;
    GLLabel       *facetProfileLabel;
    GLToggle      *facetOutgToggleLabel;
	GLToggle      *facetOutgPerAreaToggleLabel;

	GLButton      *profilePlotterBtn;
	GLButton      *texturePlotterBtn;
	GLButton      *textureScalingBtn;

	GLTitledPanel *inputPanel;
	GLTitledPanel *outputPanel;

    //Dialog
	ImportDesorption*		importDesorption = nullptr;
	TimeSettings*			timeSettings = nullptr;
	Movement*				movement = nullptr;
	MeasureForce*			measureForces = nullptr;
    FacetAdvParams*			facetAdvParams = nullptr;
    PressureEvolution*		pressureEvolution = nullptr;
	TimewisePlotter*		timewisePlotter = nullptr;
	OutgassingMapWindow*	outgassingMapWindow = nullptr;
	MomentsEditor			*momentsEditor = nullptr;
	ParameterEditor			*parameterEditor = nullptr;

    // Testing
    //int     nbSt;
    //void LogProfile();
    void BuildPipe(double ratio,int steps) override;
	void EmptyGeometry() override;
	void CrashHandler(const std::exception &e);
    int  FrameMove() override;
	
	void SelectViewer(int s) override; //to call Molflow-specialized Viewer3dSettings

protected:
    int  OneTimeSceneInit() override;
    int  RestoreDeviceObjects() override;
	int  InvalidateDeviceObjects() override;
    void ProcessMessage(GLComponent *src,int message) override;
};