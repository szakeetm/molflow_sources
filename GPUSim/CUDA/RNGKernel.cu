//based on http://ianfinlayson.net/class/cpsc425/notes/cuda-random
//

#include <math_constants.h>
// PRNG
#include <curand.h>
#include <curand_kernel.h>
//#include <thrust/device_vector.h>
//#include <thrust/random.h>

#include "RandomWrapper.h"
#include "stdio.h"
#define NB_RAND 10

/*thrust::host_vector<int> random_vector(size_t N)
{
    thrust::host_vector<int> vec(N);
    static thrust::default_random_engine rng;
    static thrust::uniform_real_distribution<float> dist(0,1);

    for (size_t i = 0; i < N; i++)
        vec[i] = dist(rng);

    return vec;
}*/

/* this GPU kernel function is used to initialize the random states */
/*__global__ void helloKernel(float* randNubers, short nbRand) {
    printf("Hello Cuworld");
    *//* we have to initialize the state *//*
    *//*curand_init(seed, *//**//* the seed can be the same for each core, here we pass the time in from the CPU *//**//*
                blockIdx.x, *//**//* the sequence number should be different for each core (unless you want all
                             cores to get the same sequence of numbers for some reason - use thread id! *//**//*
                0, *//**//* the offset is how much extra we advance in the sequence for each call, can be 0 *//**//*
                &states[blockIdx.x]);*//*
}

*//* this GPU kernel function is used to initialize the random states *//*
__global__ void init(unsigned int seed, curandState_t* states) {

    *//* we have to initialize the state *//*
    curand_init(seed, *//* the seed can be the same for each core, here we pass the time in from the CPU *//*
                blockIdx.x, *//* the sequence number should be different for each core (unless you want all
                             cores to get the same sequence of numbers for some reason - use thread id! *//*
                0, *//* the offset is how much extra we advance in the sequence for each call, can be 0 *//*
                &states[blockIdx.x]);
}

*//* this GPU kernel takes an array of states, and an array of ints, and puts a random int into each *//*
__global__ void randoms_bits(curandState_t* states, unsigned int* numbers) {
    *//* curand works like rand - except that it takes a state as a parameter *//*
    numbers[blockIdx.x] = curand(&states[blockIdx.x]) % 100;
}

*//* this GPU kernel takes an array of states, and an array of ints, and puts a random int into each *//*
__global__ void random_floats(curandState_t* states, float* numbers) {
    *//* curand works like rand - except that it takes a state as a parameter *//*
    for(int offset = 0; offset < NB_RAND; offset++)
        numbers[blockIdx.x + offset] = curand_uniform(&states[blockIdx.x]);
}*/

__global__ void test_kernel(void) {
}

namespace RandWrapper {
    void getRand(float* randNumbers, short nbRand)
    {
        test_kernel <<<1, 1>>> ();
        printf("Hello, world!");
    }
}