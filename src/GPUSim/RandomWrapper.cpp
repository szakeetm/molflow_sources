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