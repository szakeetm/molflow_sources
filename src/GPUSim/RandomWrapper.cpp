//
// Created by pbahr on 05/11/2019.
//

#ifndef MOLFLOW_PROJ_RANDOMWRAPPER_CPP
#define MOLFLOW_PROJ_RANDOMWRAPPER_CPP
#pragma once

#include <ctime>
#include <cstdlib>
#include "RandomWrapper.h"
// PRNG
#include <curand.h>
#include <curand_kernel.h>
#include <thrust/device_vector.h>
#include <thrust/random.h>

namespace RandWrapper {
    /*float *getRand(float* randNubers, uint32_t nbRand){
        //float *randomNumbers = new float[nbRand];
        srand(time(0));
        for(int i = 0; i < nbRand; i++)
            randNubers[i] = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

        return randNubers;
    }*/

    float *getCuRand(float* randNubers, uint32_t nbRand){

        //thrust::default_random_engine rng(ix+iy*2048);
        //rng.discard(10*(ix+iy*2048));
        //thrust::uniform_real_distribution<float> dist(0,1); // [0,1]

        //float *randomNumbers = new float[nbRand];
        srand(time(0));
        for(int i = 0; i < nbRand; i++)
            randNubers[i] = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

        return randNubers;
    }
    //thrust::default_random_engine rng(ix+iy*2048);
    //rng.discard(10*(ix+iy*2048));
    //thrust::uniform_real_distribution<float> dist(0,1); // [0,1]
}

#endif //MOLFLOW_PROJ_RANDOMWRAPPER_CPP