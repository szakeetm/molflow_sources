//
// Created by pbahr on 29/11/2019.
//

#ifndef MOLFLOW_PROJ_GPUDEFINES_H
#define MOLFLOW_PROJ_GPUDEFINES_H

#define PROFILE_SIZE 100
// some GPU info we are making use of
/*#define CORESPERMP 1664 // 13MP * 128 CUDA Cores
#define WARPSCHEDULERS 4*/
#define CORESPERSM 1 // 30SM * 64 CUDA Cores
#define WARPSCHEDULERS 1 // 4WS ???* 3 GPC???
#define EXTRAFACETCOUNTERS CORESPERSM*WARPSCHEDULERS

//#define MAX_DEPTH 0 // 0=no recursion
#define SCENE_EPSILON 1.e-6f // 0=no recursion

#define NB_INPOLYCHECKS 5u //0=always use middle point
// for triangle, make it a factor of 8
#ifdef WITHTRIANGLES
//#define NB_RAND_PER_STEP (8u + MAX_DEPTH*2u)
//#define NB_RAND_PER_STEP(d) (8u + (d)*2u) // without velocity
#define NB_RAND_PER_STEP(d) (9u + (d)*3u) // +1 for velocity
#else
//#define NB_RAND_PER_STEP (8u + NB_INPOLYCHECKS*2u + MAX_DEPTH*2u)
#define NB_RAND_PER_STEP(d) (8u + NB_INPOLYCHECKS*2u + (d)*2u)
#endif
//#define RAND_GEN_STEP (1u)
#define NB_RAND(g,r) ((g)*NB_RAND_PER_STEP(r))

#ifdef DEBUGCOUNT
//#define NCOUNTBINS 100
#endif

#ifdef DEBUGMISS
#define NMISSES 10
#endif

// Whether to use Ad hoc RNG or Bulk generation for GPU random numbers
//#define RNG_BULKED 1

// defined = use direct payloads, else use ptr
//#define PAYLOAD_DIRECT 1

// Random Number typedef
#ifdef RNG64
using RN_T = double;
#else
using RN_T = float;
#endif

#ifdef HIT64
using FLOAT_T = double;
#else
using FLOAT_T = float;
#endif

#endif //MOLFLOW_PROJ_GPUDEFINES_H
