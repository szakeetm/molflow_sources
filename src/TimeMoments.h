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
#include <optional>
class GLProgress_Abstract;

class TimeMoments {
public:
    static std::optional<std::pair<int, int>> HasIntervalOverlap(const std::vector<std::vector<Moment>>& vecParsedMoments);
    static std::vector<Moment> ParseUserMoment(const UserMoment& input);
    static void ParseAndCheckUserMoments(std::vector<Moment>& moments, const std::vector<UserMoment>& userMoments,
        GLProgress_Abstract& prg);
};


#endif //MOLFLOW_PROJ_TIMEMOMENTS_H
