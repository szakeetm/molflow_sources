//
// Created by pbahr on 29/11/2019.
//

#ifndef MOLFLOW_PROJ_GPUDEFINES_H
#define MOLFLOW_PROJ_GPUDEFINES_H

#define RNG64

// some GPU info we are making use of
/*#define CORESPERMP 1664 // 13MP * 128 CUDA Cores
#define WARPSCHEDULERS 4*/
#define CORESPERSM 1 // 30SM * 64 CUDA Cores
#define WARPSCHEDULERS 1 // 4WS ???* 3 GPC???

#define MAX_DEPTH 0 // 0=no recursion
#define SCENE_EPSILON 1.e-6f // 0=no recursion
#define NB_RAND 8 // for triangle, make it a factor of 8

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

// Random Number typedef
#ifdef RNG64
using RN_T = double;
#else
using RN_T = float;
#endif


#endif //MOLFLOW_PROJ_GPUDEFINES_H
