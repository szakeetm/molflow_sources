//
// Created by pbahr on 01/11/2019.
//

#include "MolTracer.h"
#include "GPUDefines.h"

flowgpu::MolTracer::~MolTracer() {
    optixHandle.cleanup();
}

// Init and start RT launch
void flowgpu::MolTracer::launchRT(){
    optixHandle.render();
}

// Fetch data
unsigned long long int flowgpu::MolTracer::fetchData(){
    return optixHandle.downloadDataFromDevice(/*pixels.data()*/);
}

static uint32_t runCount = 0;

void flowgpu::MolTracer::run()
{

    // for testing only generate and upload random numbers once
    // generate new numbers whenever necessary, recursion = TraceProcessing only, poly checks only for ray generation with polygons
    if(runCount%(NB_RAND/(8+MAX_DEPTH*2+NB_INPOLYCHECKS*2))==0){
        std::cout << "#flowgpu: generating random numbers at run #" << runCount << std::endl;
        optixHandle.generateRand();
    }
    launchRT();
    ++runCount;

    // only fetch once at the end for now
    /*if(runCount==N_LOOPS)
        fetchData();*/

    //std::cout << "#flowgpu: camera at " << cameraFrame.position << std::endl;
    //std::cout << "#flowgpu: camera lookat " << cameraFrame.get_at() << std::endl;

}

void flowgpu::MolTracer::setup(){

    int launchX = LAUNCH_SIZE_X*LAUNCH_SIZE_Y*LAUNCH_SIZE_Z;//1024*1024*1;
    int launchY = 1;

    //int nBlocks = 1024;
    //int blockSize = 1024;

    vec2i newSize = vec2i(launchX,launchY);
    if(newSize != kernelDimensions)
        resize(vec2i(launchX,launchY));

    optixHandle.initSimulation();
}

void flowgpu::MolTracer::resize(const vec2i &newSize){
    kernelDimensions = newSize;
    optixHandle.resize(newSize);
    //pixels.resize(newSize.x*newSize.y);
}