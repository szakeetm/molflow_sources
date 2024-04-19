

#include "Physics.h"
#include "Helper/MathTools.h"

double Physics::GenerateRandomVelocity(const std::vector<IntegratedVelocityEntry>& maxwell_CDF_1K, const double sqrt_temperature, const double rndVal) {
    double v_1K = InterpolateX(rndVal, maxwell_CDF_1K, false, false, true); //Allow extrapolate
    return sqrt_temperature*v_1K;
}

double Physics::GenerateDesorptionTime(const std::vector<std::vector<IntegratedDesorptionEntry>> &IDs,
                                       const SimulationFacet *src, const double rndVal, double latestMoment) {
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

double Physics::GenerateImpactParameter(const double rho, const double rndVal) {
    //PDF: P(b) = 2b / rho^2    (if b<rho, 0 otherwise)
    //CDF: 0 below b=0, b^2/rho^2 between 0 and rho, 1 over rho
    //inversion method: b=rho*sqrt(rndVal)
    return rho * sqrt(rndVal);
}

double Physics::GetScatteringAngle(const double b, const double rho, const double massRatio) {
    //ctg(theta) = 1/2ab * (rho^2 * m_tp/m_bg + b^2 - a^2)
    //massratio: test particle mass / background gas mass
    double rho_sqr = pow(rho, 2.0);
    double b_sqr = pow(b, 2.0);
    double a_sqr = rho_sqr - b_sqr;
    double a = sqrt(a_sqr);

    double ctg_theta = 1.0 / (2.0 * a * b) * (rho_sqr * massRatio+ b_sqr - a_sqr);
    double tan_theta = 1.0 / ctg_theta;
    return atan(tan_theta);
}

double Physics::GetPostScatteringVelocity(const double pre_collision_velocity, const double rho, const double b, const double theta, const double massRatio) {
    //v_new^2 = v_old^2 / (1 + (m_tp / m_bg * rho / b * sin(theta))^2
    return sqrt(pow(pre_collision_velocity, 2.0) / (1.0 + pow(massRatio*rho/b*sin(theta), 2.0)));
}