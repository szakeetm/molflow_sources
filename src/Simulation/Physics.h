

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
