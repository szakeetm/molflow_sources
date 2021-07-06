//
// Created by pascal on 2/5/21.
//

#ifndef MOLFLOW_PROJ_PHYSICS_H
#define MOLFLOW_PROJ_PHYSICS_H

#include <vector>
#include <GeometrySimu.h>

class Physics {
public:
    static double
    GenerateDesorptionTime(const std::vector<std::vector<std::pair<double, double>>> &IDs,
                           const SubprocessFacet *src, double rndVal, double latestMoment);

    static double GenerateRandomVelocity(const std::vector<std::vector<std::pair<double, double>>>& CDFs, int CDFId, double rndVal);

    static void TreatMovingFacet(SimulationModel *model, const Vector3_t<FLOAT> &position, Vector3_t<FLOAT> &direction, double &velocity);
};


#endif //MOLFLOW_PROJ_PHYSICS_H
