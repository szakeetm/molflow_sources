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

#include <vector>
#include "MolflowSimGeom.h"

//! Class that implements advanced methods that are related to the SimulationModel or an individual Particle
class Physics {
public:
    static double
    GenerateDesorptionTime(const std::vector<std::vector<IntegratedDesorptionEntry>> &IDs,
                           const SimulationFacet *src, const double rndVal, double latestMoment);

    static double GenerateRandomVelocity(const std::vector<IntegratedVelocityEntry>& maxwell_CDF_1K, const double sqrt_temperature, const double rndVal);

    static void TreatMovingFacet(std::shared_ptr<MolflowSimulationModel> model, const Vector3d &position, Vector3d &direction, double &velocity);

    static double GenerateImpactParameter(const double rho, const double rndVal);

    static double GetScatteringAngle(const double b, const double rho, const double massRatio);
    
    static double GetPostScatteringVelocity(const double pre_collision_velocity, const double rho, const double b, const double theta, const double massRatio);
};
