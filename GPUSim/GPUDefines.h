//
// Created by pbahr on 29/11/2019.
//

#ifndef MOLFLOW_PROJ_GPUDEFINES_H
#define MOLFLOW_PROJ_GPUDEFINES_H

// some GPU info we are making use of
#define CORESPERMP 1664 // 13MP * 128 CUDA Cores
#define WARPSCHEDULERS 4

#define N_LOOPS 200
#ifdef DEBUG
#define N_LOOPS 20
#endif

#define MAX_DEPTH 0 // 0=no recursion
#define SCENE_EPSILON 1.e-6f // 0=no recursion
#define NB_RAND 32


#ifndef DEBUG
#define LAUNCH_SIZE_X 1664
#define LAUNCH_SIZE_Y 64
#define LAUNCH_SIZE_Z 16
#else
#define LAUNCH_SIZE_X 1664
#define LAUNCH_SIZE_Y 64
#define LAUNCH_SIZE_Z 16
#endif

#define NB_INPOLYCHECKS 0 //0=always use middle point
#ifdef DEBUG
#define BOUND_CHECK
#endif

#ifdef DEBUGCOUNT
//#define NCOUNTBINS 100
#endif

#ifdef DEBUGMISS
#define NMISSES 10
#endif

#endif //MOLFLOW_PROJ_GPUDEFINES_H
