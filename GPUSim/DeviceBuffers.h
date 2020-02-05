//
// Created by pbahr on 05/02/2020.
//

#ifndef MOLFLOW_PROJ_DEVICEBUFFERS_H
#define MOLFLOW_PROJ_DEVICEBUFFERS_H

#include "CUDABuffer.h"

namespace flowgpu {
    struct SBTMemory{
        CUDABuffer raygenRecordsBuffer;
        CUDABuffer missRecordsBuffer;
        CUDABuffer hitgroupRecordsBuffer;
    };


    struct SimulationMemory_perThread{
        CUDABuffer moleculeBuffer;
        CUDABuffer randBuffer;
        CUDABuffer randOffsetBuffer;
    };

    struct SimulationMemory_perFacet{
        CUDABuffer hitCounterBuffer;
        CUDABuffer missCounterBuffer;
    };

    class DeviceBuffers {

    };
}

#endif //MOLFLOW_PROJ_DEVICEBUFFERS_H
