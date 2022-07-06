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

#ifndef MOLFLOW_PROJ_TIMEMOMENTS_H
#define MOLFLOW_PROJ_TIMEMOMENTS_H

#include <vector>
#include "MolflowTypes.h"

class TimeMoments {
public:
    static int CheckIntervalOverlap(const std::vector<Moment>& vecA, const std::vector<Moment>& vecB);
    static std::pair<int, int> CheckIntervalOverlap(const std::vector<std::vector<Moment>>& vecParsedMoments);
    static std::vector<Moment> ParseMoment(const std::string& userInput, double timeWindow);
    static int
    ParseAndCheckUserMoments(std::vector<Moment> *moments, std::vector<UserMoment> *userMoments,
                             double *progress);

    static int AddMoment(std::vector<Moment> *moments, std::vector<Moment> newMoments); //Adds a time serie to moments and returns the number of elements
};


#endif //MOLFLOW_PROJ_TIMEMOMENTS_H
