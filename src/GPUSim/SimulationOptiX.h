// ======================================================================== //
// Copyright 2018-2019 Ingo Wald                                            //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

#pragma once

// our own classes, partly shared between host and device
#include "DeviceBuffers.h"
#include "LaunchParams.h"
#include "Model.h"
#include "HostData.h"

/*! \namespace flowgpu - Molflow GPU code */
namespace flowgpu {



    struct OptixState
    {

        /*! @{ CUDA device context and stream that optix pipeline will run
            on, as well as device properties for this device */
        CUcontext          cudaContext;
        cudaDeviceProp     deviceProps;

        CUstream                    stream                    = 0;
#ifdef MULTI_STREAMS
        CUstream                    stream2                    = 0;
        std::vector<CUstream> cuStreams;
#endif

        /*! @} */

        OptixDeviceContext          context                   = 0;        //! the optix context that our pipeline will run in.

        OptixTraversableHandle      asHandle                = {};
        CUDABuffer asBuffer;                                        //! buffer that keeps the (final, compacted) accel structure

        /*! @{ the module that contains our device programs */
        struct Modules {
            OptixModule                 geometryModule        = 0;
            OptixModule                 rayModule             = 0;
            OptixModule                 traceModule           = 0;
            OptixModule                 exceptionModule       = 0;
        } modules;

        OptixModuleCompileOptions   moduleCompileOptions;
        /* @} */

        /*!  all our program(group)s, for now only for one molecule and surface type */
        OptixProgramGroup           raygenPG            = 0;
        OptixProgramGroup           missPG              = 0;
        std::vector<OptixProgramGroup>           hitgroupPG;
        OptixProgramGroup           exceptionPG         = 0;


        /*!  and the SBT built around them */
        OptixShaderBindingTable     sbt                       = {};

        /*! @{ the pipeline we're building */
        OptixPipeline               pipeline                  = 0;
        OptixPipelineCompileOptions pipelineCompileOptions  = {};
        OptixPipelineLinkOptions    pipelineLinkOptions;
        /*! @} */

        /*! @{ our launch parameters, on the host, and the buffer to store
            them on the device */
        LaunchParams                      launchParams;
        CUDABuffer   launchParamsBuffer;
/*! @} */
    };

    /*! a sample OptiX-7 renderer that demonstrates how to set up
        context, module, programs, pipeline, SBT, etc, and perform a
        valid launch that renders some pixel (using a simple test
        pattern, in this case */
    class SimulationOptiX
    {
        // ------------------------------------------------------------------
        // publicly accessible interface
        // ------------------------------------------------------------------
    public:
        /*! constructor - performs all setup, including initializing
          optix, creates module, pipeline, programs, SBT, etc. */
        SimulationOptiX(const Model *model, const uint2 &launchSize);
        ~SimulationOptiX();

        void resetDeviceData(const uint2 &newSize);
            /*! upload some parts only on start */
        void initSimulation();

        /*! launch a bunch of molecules to trace */
        void launchMolecules();

        /*! resize thread buffer to given resolution */
        void resize(const uint2 &newSize);

        /*! resize t buffer to given resolution */
        void generateRand();

        /*! download the rendered color buffer */
        void downloadDataFromDevice(HostData* hostData);
        void resetDeviceBuffers();
        void updateHostData(HostData* tempData);

        void cleanup();

    protected:
        // ------------------------------------------------------------------
        // internal helper functions
        // ------------------------------------------------------------------

        /*! helper function that initializes optix and checks for errors */
        void initOptix();

        void initLaunchParams(const uint2 &newSize);

        /*! creates and configures a optix device context (in this simple
          example, only for the primary GPU device) */
        void createContext();

        /*! creates the module that contains all the programs we are going
          to use. in this simple example, we use a single module from a
          single .cu file, using a single embedded ptx string */
        void createModule();

        /*! does all setup for the raygen program(s) we are going to use */
        void createRaygenPrograms(std::vector<OptixProgramGroup> &programGroups);

        /*! does all setup for the miss program(s) we are going to use */
        void createMissPrograms(std::vector<OptixProgramGroup> &programGroups);

        /*! does all setup for the hitgroup program(s) we are going to use */
        void createHitgroupPrograms(std::vector<OptixProgramGroup> &programGroups);

        void createExceptionPrograms(std::vector<OptixProgramGroup> &programGroups);

        /*! assembles the full pipeline of all programs */
        void createPipeline(std::vector<OptixProgramGroup> &programGroups);

        /*! constructs the shader binding table */
        void buildSBTPolygon();
        void buildSBTTriangle();

        /*! build an acceleration structure for the given triangle mesh */
        OptixTraversableHandle buildAccelPolygon();
        OptixTraversableHandle buildAccelTriangle();

    protected:

        OptixState state;

        SBTMemory sbt_memory;
        SimulationMemory_perThread sim_memory;
        SimulationMemory_perFacet facet_memory;
        DeviceTriangleMemory tri_memory;
        DevicePolygonMemory poly_memory;
#ifdef DEBUG
        DeviceMemoryDebug memory_debug;
#endif

        /*! the model we are going to trace rays against */
        const Model *model;
    };

} // ::flowgpu
