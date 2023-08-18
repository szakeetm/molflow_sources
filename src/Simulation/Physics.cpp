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

#include "Physics.h"
#include "Helper/MathTools.h"

double Physics::GenerateRandomVelocity(const std::vector<IntegratedVelocityEntry>& maxwell_CDF_1K, const double sqrt_temperature, const double rndVal) {
    double v_1K = InterpolateX(rndVal, maxwell_CDF_1K, false, false, true); //Allow extrapolate
    return sqrt_temperature*v_1K;
}

double Physics::GenerateDesorptionTime(const std::vector<std::vector<IntegratedDesorptionEntry>> &IDs,
                                       const SimulationFacet *src, double rndVal, double latestMoment) {
    if (!src->sh.outgassingParam.empty()) { //time-dependent desorption
        return InterpolateX(rndVal * IDs[src->sh.IDid].back().cumulativeDesValue, IDs[src->sh.IDid],
                            false, false, true); //allow extrapolate
    } else {
        return rndVal * latestMoment; //continous desorption between 0 and latestMoment
    }
}

/**
* \brief Updates particle direction and velocity if we are dealing with a moving facet (translated or rotated)
*/
void
Physics::TreatMovingFacet(std::shared_ptr<MolflowSimulationModel> model, const Vector3d &position, Vector3d &direction, double &velocity) {
    Vector3d localVelocityToAdd;
    if (model->sp.motionType == 1) { //Translation
        localVelocityToAdd = model->sp.motionVector2; //Fixed translational vector
    } else if (model->sp.motionType == 2) { //Rotation
        Vector3d distanceVector = 0.01 * (position -
                                          model->sp.motionVector1); //distance from base, with cm->m conversion, motionVector1 is rotation base point
        localVelocityToAdd = CrossProduct(model->sp.motionVector2, distanceVector); //motionVector2 is rotation axis
    }
    Vector3d oldVelocity, newVelocity;
    oldVelocity = direction * velocity;
    newVelocity = oldVelocity + localVelocityToAdd;
    direction = newVelocity.Normalized();
    velocity = newVelocity.Norme();
}