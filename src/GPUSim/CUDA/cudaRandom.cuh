//
// Created by pbahr on 04/11/2019.
//

#ifndef MOLFLOW_PROJ_CUDARANDOM_CUH
#define MOLFLOW_PROJ_CUDARANDOM_CUH

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

// PRNG
#pragma once
#include <curand_kernel.h>
#include "cuda_runtime.h"
#include "GPUDefines.h"

namespace crng {
    //__host__ int  initializeRandDevice(unsigned int kernelSize, void* states, unsigned int seed = 1234ULL);
    __host__ int  initializeRandDevice(unsigned int kernelSize, curandState_t **states, unsigned int seed = 1234ULL);
    __host__ int  initializeRandDevice_ref(unsigned int kernelSize, void *&states, unsigned int seed = 1234ULL);

    __host__ void generateRandDevice(unsigned int kernelSize, curandState_t *states);
    __host__ void generateRand(unsigned int kernelSize, curandState_t *states, float *randomNumbers);
    __host__ void destroyRand(curandState_t *states, float *randomNumbers);
    __host__ int  testRand(void** devData, size_t n);
    __host__ int  printDevDataAtHost(void* devData, size_t n);

    // float
    __host__ int initializeRandHost(unsigned int kernelSize, RN_T **randomNumbersPtr, unsigned int seed = 1234ULL);
    __host__ int  generateRandHost(unsigned int kernelSize, RN_T *randomNumbers);
    __host__ int  destroyRandHost(RN_T **randomNumbersPtr);
    __host__ int  destroyRandDevice(curandState_t **states);

    // offset for random number buffer
    __host__ int  offsetBufferZeroInit(unsigned int kernelSize, void *randomOffsets);

/*    class cudaRandom {
        void initializeRand();

    public:
        CUDA_CALLABLE_MEMBER cudaRandom(unsigned int N);

        CUDA_CALLABLE_MEMBER ~cudaRandom();

        CUDA_CALLABLE_MEMBER void generateRand();

        curandState_t *states;
        //unsigned int* randomNumbers;
        float *randomNumbers;
        unsigned int kernelSize;

    };*/
}

#endif //MOLFLOW_PROJ_CUDARANDOM_CUH
