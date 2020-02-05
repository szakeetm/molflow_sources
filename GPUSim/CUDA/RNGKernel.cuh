//
// Created by pbahr on 27/11/2019.
//

#ifndef MOLFLOW_PROJ_RNGKERNEL_CUH
#define MOLFLOW_PROJ_RNGKERNEL_CUH

#pragma once
#include "cuda_runtime.h"
#include "device_launch_parameters.h"


namespace RandWrapper {
    void getRand(float* randNumbers, short nbRand);
}

#endif //MOLFLOW_PROJ_RNGKERNEL_CUH
