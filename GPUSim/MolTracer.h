//
// Created by pbahr on 01/11/2019.
//

#ifndef MOLFLOW_PROJ_MOLTRACER_H
#define MOLFLOW_PROJ_MOLTRACER_H

#include "OptixControl.h"

namespace flowgpu {

    struct MolTracer {
        MolTracer(const Model* model) : optixHandle(model){}
        ~MolTracer();
        void setup();

        void launchRT();

        unsigned long long int fetchData();

        void run();

        void resize(const vec2i &newSize);

        vec2i                 kernelDimensions; // blocks and threads per block
        OptixControl        optixHandle;
        std::vector<uint32_t> pixels;
    };

}

#endif //MOLFLOW_PROJ_MOLTRACER_H
