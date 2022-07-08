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
};

#endif //MOLFLOW_PROJ_MOLFLOWBUFFER_H
