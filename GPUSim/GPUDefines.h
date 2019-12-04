//
// Created by pbahr on 29/11/2019.
//

#ifndef MOLFLOW_PROJ_GPUDEFINES_H
#define MOLFLOW_PROJ_GPUDEFINES_H

#define N_LOOPS 200
#ifdef DEBUG
#define N_LOOPS 10
#endif

//#define BOUND_CHECK
#define MAX_DEPTH 2 // 0=no recursion
#define NB_RAND 100
#define LAUNCH_SIZE_X 1024
#define LAUNCH_SIZE_Y 128
#define LAUNCH_SIZE_Z 4
#ifdef DEBUG
    #define LAUNCH_SIZE_X 4
    #define LAUNCH_SIZE_Y 1
    #define LAUNCH_SIZE_Z 1
#endif

#define NB_INPOLYCHECKS 3 //0=always use middle point

#endif //MOLFLOW_PROJ_GPUDEFINES_H
