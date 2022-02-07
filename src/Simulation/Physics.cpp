//
// Created by pascal on 2/5/21.
//

#include "Physics.h"
#include "Helper/MathTools.h"

double Physics::GenerateRandomVelocity(const std::vector<std::vector<std::pair<double, double>>>& CDFs, int CDFId, const double rndVal) {
    //return FastLookupY(randomGenerator.rnd(),CDFs[CDFId],false);
    //double r = randomGenerator.rnd();
    double v = InterpolateX(rndVal, CDFs[CDFId], false, false, true); //Allow extrapolate
    return v;
}

double Physics::GenerateDesorptionTime(const std::vector<std::vector<std::pair<double, double>>> &IDs,
                                       const SubprocessFacet *src, double rndVal, double latestMoment) {
    if (src->sh.outgassing_paramId >= 0) { //time-dependent desorption
        return InterpolateX(rndVal * IDs[src->sh.IDid].back().second, IDs[src->sh.IDid],
                            false, false, true); //allow extrapolate
    } else {
        return rndVal * latestMoment; //continous desorption between 0 and latestMoment
    }
}

/**
* \brief Updates particle direction and velocity if we are dealing with a moving facet (translated or rotated)
*/
void
Physics::TreatMovingFacet(SimulationModel *model, const Vector3d &position, Vector3d &direction, double &velocity) {
    Vector3d localVelocityToAdd;
    if (model->wp.motionType == 1) { //Translation
        localVelocityToAdd = model->wp.motionVector2; //Fixed translational vector
    } else if (model->wp.motionType == 2) { //Rotation
        Vector3d distanceVector = 0.01 * (position -
                                          model->wp.motionVector1); //distance from base, with cm->m conversion, motionVector1 is rotation base point
        localVelocityToAdd = CrossProduct(model->wp.motionVector2, distanceVector); //motionVector2 is rotation axis
    }
    Vector3d oldVelocity, newVelocity;
    oldVelocity = direction * velocity;
    newVelocity = oldVelocity + localVelocityToAdd;
    direction = newVelocity.Normalized();
    velocity = newVelocity.Length();
}