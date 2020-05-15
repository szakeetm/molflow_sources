//
// Created by Pascal Baehr on 22.04.20.
//

#ifndef MOLFLOW_PROJ_MOLFLOWBUFFER_H
#define MOLFLOW_PROJ_MOLFLOWBUFFER_H

#include "Buffer_shared.h"

struct MolflowWorkerParams : public WorkerParams {
    double latestMoment;
    double totalDesorbedMolecules; //Number of molecules desorbed between t=0 and latest_moment
    double finalOutgassingRate; //Number of outgassing molecules / second at latest_moment (constant flow)
    double finalOutgassingRate_Pa_m3_sec;
    double gasMass;
    bool enableDecay;
    double halfLife;
    double timeWindowSize;
    bool useMaxwellDistribution; //true: Maxwell-Boltzmann distribution, false: All molecules have the same (V_avg) speed
    bool calcConstantFlow;

    int motionType;
    Vector3d motionVector1; //base point for rotation
    Vector3d motionVector2; //rotation vector or velocity vector
    size_t    sMode;                // Simu mode (MC_MODE or AC_MODE)
};

#endif //MOLFLOW_PROJ_MOLFLOWBUFFER_H
