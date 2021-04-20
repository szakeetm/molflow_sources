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
#include <filesystem>
#include <sstream>
#include "Parameter.h"
#include "File.h"


/*
StringClass::StringClass()
{
	fromCatalog = false;
	name = "";
}*/
int Parameter::LoadParameterCatalog(std::vector<Parameter> &vec_param) {
    int res = 0;

    std::filesystem::path catalogPath = "parameter_catalog"; //string (POSIX) or wstring (Windows)
    if (!std::filesystem::exists(catalogPath)) return 1; //No param_catalog directory
    for (const auto & p : std::filesystem::directory_iterator(catalogPath)) {
        if (p.path().extension() == ".csv") {

            std::string csvPath = p.path().u8string();
            std::string csvName = p.path().filename().u8string();

            Parameter newParam{};
            newParam.fromCatalog = true;
            newParam.name = "[catalog] " + csvName;

            std::vector<std::vector<std::string>> table;
            try {

                FileReader f(csvPath);
                table = f.ImportCSV_string();
            }
            catch(std::exception &e) {
                ++res;
                char errMsg[512];
                sprintf(errMsg, "Failed to load CSV file.\n%s", e.what());
                //GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR); //Can't display dialog window: interface not yet initialized
                std::cerr << errMsg << "\n";
                continue;
            }
            //Parse
            for (size_t i = 0; i < table.size(); i++) {
                std::vector<std::string> row = table[i];
                if (row.size() != 2) {
                    ++res;
                    std::stringstream errMsg;
                    errMsg << p.path().filename() << " Row " << i + 1 << "has " << row.size() << " values instead of 2.";
                    //GLMessageBox::Display(errMsg.str().c_str(), "Error", GLDLG_OK, GLDLG_ICONERROR);
                    std::cerr << errMsg.str().c_str() << "\n";
                    break;
                }
                else {
                    double valueX, valueY;
                    try {
                        valueX = ::atof(row[0].c_str());
                    }
                    catch (std::exception& err) {
                        ++res;
                        char tmp[256];
                        sprintf(tmp, "Can't parse value \"%s\" in row %zd, first column:\n%s", row[0].c_str(), i + 1, err.what());
                        //GLMessageBox::Display(tmp, "Invalid parameter definition", GLDLG_OK, GLDLG_ICONWARNING);
                        std::cerr << tmp << "\n";
                        break;
                    }
                    try {
                        valueY = ::atof(row[1].c_str());
                    }
                    catch (std::exception& err) {
                        ++res;
                        char tmp[256];
                        sprintf(tmp, "Can't parse value \"%s\" in row %zd, second column:\n%s", row[1].c_str(), i + 1, err.what());
                        //GLMessageBox::Display(tmp, "Invalid parameter definition", GLDLG_OK, GLDLG_ICONWARNING);
                        std::cerr << tmp << "\n";
                        break;
                    }
                    newParam.AddPair(valueX, valueY, true); //insert in correct position
                }
            }
            vec_param.push_back(newParam);
        }
    }
}
