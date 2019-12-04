//
// Created by pbahr on 05/11/2019.
//

#ifndef MOLFLOW_PROJ_RANDOMWRAPPER_H
#define MOLFLOW_PROJ_RANDOMWRAPPER_H

#pragma once

#include <cstdint>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"


namespace RandWrapper {
    void getRand(float* randNumbers, short nbRand);
    //float* getCuRand(float* randNubers, uint32_t nbRand);
}

#endif //MOLFLOW_PROJ_RANDOMWRAPPER_H
