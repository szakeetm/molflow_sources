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

#include <curand_kernel.h>
#include "OptixControl.h"
#include "LaunchParams.h"
#include "cudaRandom.cuh"
#include "GPUDefines.h"

extern void initializeRand(unsigned int kernelSize, curandState_t *states, float *randomNumbers);
extern void generateRand(unsigned int kernelSize, curandState_t *states, float *randomNumbers);
extern void destroyRand(curandState_t *states, float *randomNumbers);

// this include may only appear in a single source file:
#include <optix_function_table_definition.h>
#include <fstream>
#include <ctime>

/*! \namespace flowgpu - Molflow GPU code */
namespace flowgpu {

    extern "C" char geometry_ptx_code[];
    extern "C" char trace_ptx_code[];
    extern "C" char ray_ptx_code[];

    /*! SBT record for a raygen program */
    struct __align__( OPTIX_SBT_RECORD_ALIGNMENT ) RaygenRecord
    {
        __align__( OPTIX_SBT_RECORD_ALIGNMENT ) char header[OPTIX_SBT_RECORD_HEADER_SIZE];
        // just a dummy value - later examples will use more interesting
        // data here
        //void *data;
        PolygonRayGenData data_poly;
    };

    /*! SBT record for a miss program */
    struct __align__( OPTIX_SBT_RECORD_ALIGNMENT ) MissRecord
    {
        __align__( OPTIX_SBT_RECORD_ALIGNMENT ) char header[OPTIX_SBT_RECORD_HEADER_SIZE];
        // just a dummy value - later examples will use more interesting
        // data here
        //void *data;
        //void *data_poly;
        PolygonMeshSBTData data_poly;
    };

    /*! SBT record for a hitgroup program */
    struct __align__( OPTIX_SBT_RECORD_ALIGNMENT ) HitgroupRecord
    {
        __align__( OPTIX_SBT_RECORD_ALIGNMENT ) char header[OPTIX_SBT_RECORD_HEADER_SIZE];
        PolygonMeshSBTData data_poly;
    };

    /*! SBT record for a raygen program */
    struct __align__( OPTIX_SBT_RECORD_ALIGNMENT ) RaygenRecordTri
    {
        __align__( OPTIX_SBT_RECORD_ALIGNMENT ) char header[OPTIX_SBT_RECORD_HEADER_SIZE];
        // just a dummy value - later examples will use more interesting
        // data here
        //void *data;
        TriangleRayGenData data_tri;
    };

    /*! SBT record for a miss program */
    struct __align__( OPTIX_SBT_RECORD_ALIGNMENT ) MissRecordTri
    {
        __align__( OPTIX_SBT_RECORD_ALIGNMENT ) char header[OPTIX_SBT_RECORD_HEADER_SIZE];
        // just a dummy value - later examples will use more interesting
        // data here
        //void *data;
        //void *data_poly;
        TriangleMeshSBTData data_tri;
    };

    /*! SBT record for a hitgroup program */
    struct __align__( OPTIX_SBT_RECORD_ALIGNMENT ) HitgroupRecordTri
    {
        __align__( OPTIX_SBT_RECORD_ALIGNMENT ) char header[OPTIX_SBT_RECORD_HEADER_SIZE];
        TriangleMeshSBTData data_tri;
    };

    static void
    polygon_bound(std::vector<uint32_t> polyIndicies, uint32_t indexOffset, std::vector<vec3f> polyVertices,
                  uint32_t nbIndices, float *result)
    {

        float3 m_max{static_cast<float>(-1e100),static_cast<float>(-1e100),static_cast<float>(-1e100)};
        float3 m_min{static_cast<float>(1e100),static_cast<float>(1e100),static_cast<float>(1e100)};

        for(uint32_t ind = indexOffset; ind < indexOffset + nbIndices; ind++){
            auto polyIndex = polyIndicies[ind];
            const auto& vert = polyVertices[polyIndex];
            //auto newVert = vert / dot(vert,vert); // scale back
            m_min.x = std::min(vert.x,m_min.x);
            m_min.y = std::min(vert.y, m_min.y);
            m_min.z = std::min(vert.z, m_min.z);
            m_max.x = std::max(vert.x, m_max.x);
            m_max.y = std::max(vert.y, m_max.y);
            m_max.z = std::max(vert.z, m_max.z);
        }

        OptixAabb* aabb = reinterpret_cast<OptixAabb*>(result);

        *aabb = {
                m_min.x, m_min.y, m_min.z,
                m_max.x, m_max.y, m_max.z
        };
    }

    /*! constructor - performs all setup, including initializing
      optix, creates module, pipeline, programs, SBT, etc. */
    OptixControl::OptixControl(const Model *model)
            : model(model)
    {
        initOptix();

        std::cout << "#flowgpu: creating optix context ..." << std::endl;
        createContext();

        std::cout << "#flowgpu: setting up module ..." << std::endl;
        createModule();

        std::cout << "#flowgpu: creating raygen programs ..." << std::endl;
        createRaygenPrograms();
        std::cout << "#flowgpu: creating miss programs ..." << std::endl;
        createMissPrograms();
        std::cout << "#flowgpu: creating hitgroup programs ..." << std::endl;
        createHitgroupPrograms();

#ifdef WITHTRIANGLES
        launchParams.traversable = buildAccelTriangle();
#else
        launchParams.traversable = buildAccelPolygon();
#endif
        std::cout << "#flowgpu: setting up optix pipeline ..." << std::endl;
        createPipeline();

        std::cout << "#flowgpu: building SBT ..." << std::endl;
#ifdef WITHTRIANGLES
        buildSBTTriangle();
#else
        buildSBTPolygon();
#endif
        launchParamsBuffer.alloc(sizeof(launchParams));
        std::cout << "#flowgpu: context, module, pipeline, etc, all set up ..." << std::endl;

        std::cout << GDT_TERMINAL_GREEN;
        std::cout << "#flowgpu: Optix 7 Sample fully set up" << std::endl;
        std::cout << GDT_TERMINAL_DEFAULT;
    }

    OptixTraversableHandle OptixControl::buildAccelPolygon()
    {
        PING;
        PRINT(model->poly_meshes.size());

        // All buffers that should be uploaded to device memory
        poly_memory.aabbBuffer.resize(model->poly_meshes.size());
        poly_memory.vertexBuffer.resize(model->poly_meshes.size());
        poly_memory.vertex2Buffer.resize(model->poly_meshes.size());
        poly_memory.indexBuffer.resize(model->poly_meshes.size());
        poly_memory.polyBuffer.resize(model->poly_meshes.size());
        poly_memory.facprobBuffer.resize(model->poly_meshes.size());
        poly_memory.cdfBuffer.resize(model->poly_meshes.size());

        OptixTraversableHandle asHandle { 0 };

        // ==================================================================
        // triangle inputs
        // ==================================================================

        //std::vector<OptixBuildInput> triangleInput(model->poly_meshes.size());

        std::vector<OptixBuildInput> polygonInput(model->poly_meshes.size());
        std::vector<CUdeviceptr> d_aabb(model->poly_meshes.size());
        /*std::vector<CUdeviceptr> d_vertices(model->poly_meshes.size());
        std::vector<CUdeviceptr> d_vertices2(model->poly_meshes.size());
        std::vector<CUdeviceptr> d_indices(model->poly_meshes.size());
        std::vector<CUdeviceptr> d_polygons(model->poly_meshes.size());*/

        for (int meshID=0;meshID<model->poly_meshes.size();meshID++) {

            // upload the model to the device: the builder
            PolygonMesh &mesh = *model->poly_meshes[meshID];

            std::vector<OptixAabb> aabb(mesh.poly.size());
            int bbCount = 0;
            for(Polygon& poly : mesh.poly){
                polygon_bound(mesh.indices, poly.vertOffset, mesh.vertices3d, poly.nbVertices,
                              reinterpret_cast<float *>(&aabb[bbCount]));
                //std::cout << bbCount<<"# poly box: " << "("<<aabb[bbCount].minX <<","<<aabb[bbCount].minY <<","<<aabb[bbCount].minZ <<")-("<<aabb[bbCount].maxX <<","<<aabb[bbCount].maxY <<","<<aabb[bbCount].maxZ <<")"<<std::endl;
                //std::cout << bbCount<<"# poly uv: " << "("<<poly.U <<","<<poly.V <<")"<<std::endl;

                bbCount++;
            }

            poly_memory.aabbBuffer[meshID].alloc_and_upload(aabb);

            polygonInput[meshID] = {};
            polygonInput[meshID].type
                    = OPTIX_BUILD_INPUT_TYPE_CUSTOM_PRIMITIVES;
            uint32_t aabb_input_flags[1]       = {OPTIX_GEOMETRY_FLAG_NONE
                                                    | OPTIX_GEOMETRY_FLAG_DISABLE_ANYHIT};
            // create local variables, because we need a *pointer* to the
            // device pointers
            d_aabb[meshID] =  poly_memory.aabbBuffer[meshID].d_pointer();

            // in this example we have one SBT entry, and no per-primitive
            // materials:
            polygonInput[meshID].aabbArray.aabbBuffers   = &(d_aabb[meshID]);
            polygonInput[meshID].aabbArray.flags         = aabb_input_flags;
            polygonInput[meshID].aabbArray.numSbtRecords = 1;
            polygonInput[meshID].aabbArray.numPrimitives = mesh.poly.size();
            polygonInput[meshID].aabbArray.sbtIndexOffsetBuffer         = 0;
            polygonInput[meshID].aabbArray.sbtIndexOffsetSizeInBytes    = 0;
            polygonInput[meshID].aabbArray.primitiveIndexOffset         = 0;
        }
        // ==================================================================
        // BLAS setup
        // ==================================================================

        OptixAccelBuildOptions accelOptions = {};
        accelOptions.buildFlags             = OPTIX_BUILD_FLAG_NONE
                                            | OPTIX_BUILD_FLAG_ALLOW_COMPACTION
                                            | OPTIX_BUILD_FLAG_PREFER_FAST_TRACE
                ;
        accelOptions.operation              = OPTIX_BUILD_OPERATION_BUILD;

        OptixAccelBufferSizes blasBufferSizes;
        OPTIX_CHECK(optixAccelComputeMemoryUsage
                            (optixContext,
                             &accelOptions,
                             polygonInput.data(),
                             (int)model->poly_meshes.size(),  // num_build_inputs
                             &blasBufferSizes
                            ));

        // ==================================================================
        // prepare compaction
        // ==================================================================

        CUDABuffer compactedSizeBuffer;
        compactedSizeBuffer.alloc(sizeof(uint64_t));

        OptixAccelEmitDesc emitDesc;
        emitDesc.type   = OPTIX_PROPERTY_TYPE_COMPACTED_SIZE;
        emitDesc.result = compactedSizeBuffer.d_pointer();

        // ==================================================================
        // execute build (main stage)
        // ==================================================================

        CUDABuffer tempBuffer;
        tempBuffer.alloc(blasBufferSizes.tempSizeInBytes);

        CUDABuffer outputBuffer;
        outputBuffer.alloc(blasBufferSizes.outputSizeInBytes);

        OPTIX_CHECK(optixAccelBuild(optixContext,
                0,
                                    &accelOptions,
                                    polygonInput.data(),
                                    (int)model->poly_meshes.size(),
                                    tempBuffer.d_pointer(),
                                    tempBuffer.sizeInBytes,

                                    outputBuffer.d_pointer(),
                                    outputBuffer.sizeInBytes,

                                    &asHandle,

                                    &emitDesc,1
        ));
        CUDA_SYNC_CHECK();

        // ==================================================================
        // perform compaction
        // ==================================================================
        uint64_t compactedSize;
        compactedSizeBuffer.download(&compactedSize,1);

        asBuffer.alloc(compactedSize);
        OPTIX_CHECK(optixAccelCompact(optixContext,
                0,
                                      asHandle,
                                      asBuffer.d_pointer(),
                                      asBuffer.sizeInBytes,
                                      &asHandle));
        CUDA_SYNC_CHECK();

        // ==================================================================
        // aaaaaand .... clean up
        // ==================================================================
        outputBuffer.free(); // << the UNcompacted, temporary output buffer
        tempBuffer.free();
        compactedSizeBuffer.free();

        return asHandle;
    }

    OptixTraversableHandle OptixControl::buildAccelTriangle()
    {
        PING;
        PRINT(model->triangle_meshes.size());

        // All buffers that should be uploaded to device memory
        //memory.aabbBuffer.resize(model->triangle_meshes.size());
        tri_memory.vertexBuffer.resize(model->triangle_meshes.size());
        //memory.vertex2Buffer.resize(model->triangle_meshes.size());
        tri_memory.indexBuffer.resize(model->triangle_meshes.size());
        tri_memory.polyBuffer.resize(model->triangle_meshes.size());
        tri_memory.facprobBuffer.resize(model->triangle_meshes.size());
        tri_memory.cdfBuffer.resize(model->triangle_meshes.size());

        OptixTraversableHandle asHandle { 0 };

        // ==================================================================
        // triangle inputs
        // ==================================================================
        std::vector<OptixBuildInput> triangleInput(model->triangle_meshes.size());
        std::vector<CUdeviceptr> d_vertices(model->triangle_meshes.size());
        std::vector<CUdeviceptr> d_indices(model->triangle_meshes.size());
        std::vector<uint32_t> triangleInputFlags(model->triangle_meshes.size());

        //std::vector<OptixBuildInput> polygonInput(model->triangle_meshes.size());
        //std::vector<CUdeviceptr> d_aabb(model->triangle_meshes.size());
        /*std::vector<CUdeviceptr> d_vertices(model->triangle_meshes.size());
        std::vector<CUdeviceptr> d_vertices2(model->triangle_meshes.size());
        std::vector<CUdeviceptr> d_indices(model->triangle_meshes.size());
        std::vector<CUdeviceptr> d_polygons(model->triangle_meshes.size());*/

        for (int meshID=0;meshID<model->triangle_meshes.size();meshID++) {

            // upload the model to the device: the builder
            TriangleMesh &mesh = *model->triangle_meshes[meshID];

            /*for(vec3i& vert : mesh.indices){

                std::cout << vert << std::endl;
                //std::cout << bbCount<<"# poly uv: " << "("<<poly.U <<","<<poly.V <<")"<<std::endl;
            }
            //exit(0);*/

            tri_memory.vertexBuffer[meshID].alloc_and_upload(mesh.vertices3d);
            tri_memory.indexBuffer[meshID].alloc_and_upload(mesh.indices);

            triangleInput[meshID] = {};
            triangleInput[meshID].type
                    = OPTIX_BUILD_INPUT_TYPE_TRIANGLES;

            // create local variables, because we need a *pointer* to the
            // device pointers
            d_vertices[meshID] = tri_memory.vertexBuffer[meshID].d_pointer();
            d_indices[meshID]  = tri_memory.indexBuffer[meshID].d_pointer();

            triangleInput[meshID].triangleArray.vertexFormat        = OPTIX_VERTEX_FORMAT_FLOAT3;
            triangleInput[meshID].triangleArray.vertexStrideInBytes = sizeof(vec3f);
            triangleInput[meshID].triangleArray.numVertices         = (int)mesh.vertices3d.size();
            triangleInput[meshID].triangleArray.vertexBuffers       = &d_vertices[meshID];

            triangleInput[meshID].triangleArray.indexFormat         = OPTIX_INDICES_FORMAT_UNSIGNED_INT3;
            triangleInput[meshID].triangleArray.indexStrideInBytes  = sizeof(vec3i);
            triangleInput[meshID].triangleArray.numIndexTriplets    = (int)mesh.indices.size();
            triangleInput[meshID].triangleArray.indexBuffer         = d_indices[meshID];

            triangleInputFlags[meshID] = OPTIX_GEOMETRY_FLAG_NONE | OPTIX_GEOMETRY_FLAG_DISABLE_ANYHIT  ;

            // in this example we have one SBT entry, and no per-primitive
            // materials:
            triangleInput[meshID].triangleArray.flags               = &triangleInputFlags[meshID];
            triangleInput[meshID].triangleArray.numSbtRecords               = 1;
            triangleInput[meshID].triangleArray.sbtIndexOffsetBuffer        = 0;
            triangleInput[meshID].triangleArray.sbtIndexOffsetSizeInBytes   = 0;
            triangleInput[meshID].triangleArray.sbtIndexOffsetStrideInBytes = 0;
        }
        // ==================================================================
        // BLAS setup
        // ==================================================================

        OptixAccelBuildOptions accelOptions = {};
        accelOptions.buildFlags             = OPTIX_BUILD_FLAG_NONE
                                              | OPTIX_BUILD_FLAG_ALLOW_COMPACTION
                                              | OPTIX_BUILD_FLAG_PREFER_FAST_TRACE
                ;
        accelOptions.motionOptions.numKeys  = 1;
        accelOptions.operation              = OPTIX_BUILD_OPERATION_BUILD;

        OptixAccelBufferSizes blasBufferSizes;
        OPTIX_CHECK(optixAccelComputeMemoryUsage
                            (optixContext,
                             &accelOptions,
                             triangleInput.data(),
                             (int)model->triangle_meshes.size(),  // num_build_inputs
                             &blasBufferSizes
                            ));

        // ==================================================================
        // prepare compaction
        // ==================================================================

        CUDABuffer compactedSizeBuffer;
        compactedSizeBuffer.alloc(sizeof(uint64_t));

        OptixAccelEmitDesc emitDesc;
        emitDesc.type   = OPTIX_PROPERTY_TYPE_COMPACTED_SIZE;
        emitDesc.result = compactedSizeBuffer.d_pointer();

        // ==================================================================
        // execute build (main stage)
        // ==================================================================

        CUDABuffer tempBuffer;
        tempBuffer.alloc(blasBufferSizes.tempSizeInBytes);

        CUDABuffer outputBuffer;
        outputBuffer.alloc(blasBufferSizes.outputSizeInBytes);

        OPTIX_CHECK(optixAccelBuild(optixContext,
                                    0,
                                    &accelOptions,
                                    triangleInput.data(),
                                    (int)model->triangle_meshes.size(),
                                    tempBuffer.d_pointer(),
                                    tempBuffer.sizeInBytes,

                                    outputBuffer.d_pointer(),
                                    outputBuffer.sizeInBytes,

                                    &asHandle,

                                    &emitDesc,1
        ));
        CUDA_SYNC_CHECK();

        // ==================================================================
        // perform compaction
        // ==================================================================
        uint64_t compactedSize;
        compactedSizeBuffer.download(&compactedSize,1);

        asBuffer.alloc(compactedSize);
        OPTIX_CHECK(optixAccelCompact(optixContext,
                                      0,
                                      asHandle,
                                      asBuffer.d_pointer(),
                                      asBuffer.sizeInBytes,
                                      &asHandle));
        CUDA_SYNC_CHECK();

        // ==================================================================
        // aaaaaand .... clean up
        // ==================================================================
        outputBuffer.free(); // << the UNcompacted, temporary output buffer
        tempBuffer.free();
        compactedSizeBuffer.free();

        return asHandle;
    }

    /*! helper function that initializes optix and checks for errors */
    void OptixControl::initOptix()
    {
        std::cout << "#flowgpu: initializing optix..." << std::endl;

        // -------------------------------------------------------
        // check for available optix7 capable devices
        // -------------------------------------------------------
        cudaFree(0);
        int numDevices;
        cudaGetDeviceCount(&numDevices);
        if (numDevices == 0)
            throw std::runtime_error("#flowgpu: no CUDA capable devices found!");
        std::cout << "#flowgpu: found " << numDevices << " CUDA devices" << std::endl;

        // -------------------------------------------------------
        // initialize optix
        // -------------------------------------------------------
        OPTIX_CHECK( optixInit() );
        std::cout << GDT_TERMINAL_GREEN
                  << "#flowgpu: successfully initialized optix... yay!"
                  << GDT_TERMINAL_DEFAULT << std::endl;
    }

    static void context_log_cb(unsigned int level,
                               const char *tag,
                               const char *message,
                               void *)
    {
        fprintf( stderr, "[%2d][%12s]: %s\n", (int)level, tag, message );
    }

    /*! creates and configures a optix device context (in this simple
      example, only for the primary GPU device) */
    void OptixControl::createContext()
    {
        // for this sample, do everything on one device
        const int deviceID = 0;
        CUDA_CHECK(SetDevice(deviceID));
        CUDA_CHECK(StreamCreate(&stream));

        cudaGetDeviceProperties(&deviceProps, deviceID);
        std::cout << "#flowgpu: running on device: " << deviceProps.name << std::endl;

        CUresult  cuRes = cuCtxGetCurrent(&cudaContext);
        if( cuRes != CUDA_SUCCESS )
            fprintf( stderr, "Error querying current context: error code %d\n", cuRes );

        OPTIX_CHECK(optixDeviceContextCreate(cudaContext, 0, &optixContext));
        OPTIX_CHECK(optixDeviceContextSetLogCallback
                            (optixContext,context_log_cb,nullptr,4));

    }



    /*! creates the module that contains all the programs we are going
      to use. in this simple example, we use a single module from a
      single .cu file, using a single embedded ptx string */
    void OptixControl::createModule()
    {
        moduleCompileOptions.maxRegisterCount  = 50;
#ifdef DEBUG
        moduleCompileOptions.optLevel          = OPTIX_COMPILE_OPTIMIZATION_LEVEL_0;
        moduleCompileOptions.debugLevel        = OPTIX_COMPILE_DEBUG_LEVEL_LINEINFO;
#else
        moduleCompileOptions.optLevel          = OPTIX_COMPILE_OPTIMIZATION_DEFAULT;
        moduleCompileOptions.debugLevel        = OPTIX_COMPILE_DEBUG_LEVEL_NONE;
#endif
        pipelineCompileOptions = {};
        pipelineCompileOptions.traversableGraphFlags = OPTIX_TRAVERSABLE_GRAPH_FLAG_ALLOW_SINGLE_GAS;
        pipelineCompileOptions.usesMotionBlur     = false;
        pipelineCompileOptions.numPayloadValues   = 5; // values that get send as PerRayData
        pipelineCompileOptions.numAttributeValues = 5; // ret values e.g. by optixReportIntersection
        pipelineCompileOptions.exceptionFlags     = OPTIX_EXCEPTION_FLAG_NONE;
        pipelineCompileOptions.pipelineLaunchParamsVariableName = "optixLaunchParams";

        pipelineLinkOptions.overrideUsesMotionBlur = false;
        pipelineLinkOptions.maxTraceDepth          = 2;

        char log[2048];
        size_t sizeof_log = sizeof( log );

        {
            const std::string ptxCode = geometry_ptx_code;
            OPTIX_CHECK(optixModuleCreateFromPTX(optixContext,
                                                 &moduleCompileOptions,
                                                 &pipelineCompileOptions,
                                                 ptxCode.c_str(),
                                                 ptxCode.size(),
                                                 log,&sizeof_log,
                                                 &modules.geometryModule
            ));
            if (sizeof_log > 1) PRINT(log);
        }

        {
            const std::string ptxCode = ray_ptx_code;
            OPTIX_CHECK(optixModuleCreateFromPTX(optixContext,
                                                 &moduleCompileOptions,
                                                 &pipelineCompileOptions,
                                                 ptxCode.c_str(),
                                                 ptxCode.size(),
                                                 log,&sizeof_log,
                                                 &modules.rayModule
            ));
            if (sizeof_log > 1) PRINT(log);
        }

        {
            const std::string ptxCode = trace_ptx_code;
            OPTIX_CHECK(optixModuleCreateFromPTX(optixContext,
                                                 &moduleCompileOptions,
                                                 &pipelineCompileOptions,
                                                 ptxCode.c_str(),
                                                 ptxCode.size(),
                                                 log,&sizeof_log,
                                                 &modules.traceModule
            ));
            if (sizeof_log > 1) PRINT(log);
        }
    }



    /*! does all setup for the raygen program(s) we are going to use */
    void OptixControl::createRaygenPrograms()
    {
        // we do a single ray gen program in this example:
        raygenPGs.resize(1);

        OptixProgramGroupOptions pgOptions = {};
        OptixProgramGroupDesc pgDesc    = {};
        pgDesc.kind                     = OPTIX_PROGRAM_GROUP_KIND_RAYGEN;
        pgDesc.raygen.module            = modules.rayModule;
        pgDesc.raygen.entryFunctionName = "__raygen__startFromSource";

        // OptixProgramGroup raypg;
        char log[2048];
        size_t sizeof_log = sizeof( log );
        OPTIX_CHECK(optixProgramGroupCreate(optixContext,
                                            &pgDesc,
                                            1,
                                            &pgOptions,
                                            log,&sizeof_log,
                                            &raygenPGs[0]
        ));
        if (sizeof_log > 1) PRINT(log);
    }

    /*! does all setup for the miss program(s) we are going to use */
    void OptixControl::createMissPrograms()
    {
        // we do a single ray gen program in this example:
        missPGs.resize(1);

        OptixProgramGroupOptions pgOptions = {};
        OptixProgramGroupDesc pgDesc    = {};
        pgDesc.kind                     = OPTIX_PROGRAM_GROUP_KIND_MISS;
        pgDesc.miss.module            = modules.traceModule;
        pgDesc.miss.entryFunctionName = "__miss__molecule";

        // OptixProgramGroup raypg;
        char log[2048];
        size_t sizeof_log = sizeof( log );
        OPTIX_CHECK(optixProgramGroupCreate(optixContext,
                                            &pgDesc,
                                            1,
                                            &pgOptions,
                                            log,&sizeof_log,
                                            &missPGs[0]
        ));
        if (sizeof_log > 1) PRINT(log);
    }

    /*! does all setup for the hitgroup program(s) we are going to use */
    void OptixControl::createHitgroupPrograms()
    {
        // for this simple example, we set up a single hit group
        hitgroupPGs.resize(1);

        OptixProgramGroupOptions pgOptions = {};
        OptixProgramGroupDesc pgDesc    = {};
        pgDesc.kind                     = OPTIX_PROGRAM_GROUP_KIND_HITGROUP;
#ifdef WITHTRIANGLES
        pgDesc.hitgroup.moduleCH            = modules.traceModule;
        pgDesc.hitgroup.entryFunctionNameCH = "__closesthit__molecule_triangle";
#else
        pgDesc.hitgroup.moduleIS            = modules.geometryModule;
        pgDesc.hitgroup.entryFunctionNameIS = "__intersection__polygon";
        pgDesc.hitgroup.moduleCH            = modules.traceModule;
        pgDesc.hitgroup.entryFunctionNameCH = "__closesthit__molecule";
#endif


        pgDesc.hitgroup.moduleAH            = modules.traceModule;
        pgDesc.hitgroup.entryFunctionNameAH = "__anyhit__molecule";

        char log[2048];
        size_t sizeof_log = sizeof( log );
        OPTIX_CHECK(optixProgramGroupCreate(optixContext,
                                            &pgDesc,
                                            1,
                                            &pgOptions,
                                            log,&sizeof_log,
                                            &hitgroupPGs[0]
        ));
        if (sizeof_log > 1) PRINT(log);
    }


    /*! assembles the full pipeline of all programs */
    void OptixControl::createPipeline()
    {
        std::vector<OptixProgramGroup> programGroups;
        for (auto pg : raygenPGs)
            programGroups.push_back(pg);
        for (auto pg : missPGs)
            programGroups.push_back(pg);
        for (auto pg : hitgroupPGs)
            programGroups.push_back(pg);

        char log[2048];
        size_t sizeof_log = sizeof( log );
        OPTIX_CHECK(optixPipelineCreate(optixContext,
                                        &pipelineCompileOptions,
                                        &pipelineLinkOptions,
                                        programGroups.data(),
                                        (int)programGroups.size(),
                                        log,&sizeof_log,
                                        &pipeline
        ));
        if (sizeof_log > 1) PRINT(log);

        OPTIX_CHECK(optixPipelineSetStackSize
                            (/* [in] The pipeline to configure the stack size for */
                                    pipeline,
                                    /* [in] The direct stack size requirement for direct
                                       callables invoked from IS or AH. */
                                    2*1024,
                                    /* [in] The direct stack size requirement for direct
                                       callables invoked from RG, MS, or CH.  */
                                    2*1024,
                                    /* [in] The continuation stack requirement. */
                                    2*1024,
                                    /* [in] The maximum depth of a traversable graph
                                       passed to trace. */
                                    1));
        if (sizeof_log > 1) PRINT(log);
    }


    /*! constructs the shader binding table */
    void OptixControl::buildSBTPolygon()
    {
        // first allocate device memory and upload data
        for (int meshID=0;meshID<(int)model->poly_meshes.size();meshID++) {
            PolygonMesh &mesh = *model->poly_meshes[meshID];
            poly_memory.vertexBuffer[meshID].alloc_and_upload(mesh.vertices3d);
            poly_memory.vertex2Buffer[meshID].alloc_and_upload(mesh.vertices2d);
            poly_memory.indexBuffer[meshID].alloc_and_upload(mesh.indices);
            poly_memory.polyBuffer[meshID].alloc_and_upload(mesh.poly);
            poly_memory.cdfBuffer[meshID].alloc_and_upload(mesh.cdfs);
            poly_memory.facprobBuffer[meshID].alloc_and_upload(mesh.facetProbabilities);
        }

        // ------------------------------------------------------------------
        // build raygen records
        // ------------------------------------------------------------------
        std::vector<RaygenRecord> raygenRecords;
        for (int i=0;i<raygenPGs.size();i++) {
            RaygenRecord rec;
            OPTIX_CHECK(optixSbtRecordPackHeader(raygenPGs[i],&rec));
            rec.data_poly.vertex = (vec3f*)poly_memory.vertexBuffer[0].d_pointer();
            rec.data_poly.vertex2 = (vec2f*)poly_memory.vertex2Buffer[0].d_pointer();
            rec.data_poly.index = (uint32_t*)poly_memory.indexBuffer[0].d_pointer();
            rec.data_poly.poly = (Polygon*)poly_memory.polyBuffer[0].d_pointer();
            rec.data_poly.cdfs = (float*)poly_memory.cdfBuffer[0].d_pointer();
            rec.data_poly.facetProbabilities = (vec2f*)poly_memory.facprobBuffer[0].d_pointer();
            raygenRecords.push_back(rec);
        }
        raygenRecordsBuffer.alloc_and_upload(raygenRecords);
        sbt.raygenRecord = raygenRecordsBuffer.d_pointer();

        // ------------------------------------------------------------------
        // build miss records
        // ------------------------------------------------------------------
        int numObjects = (int)model->poly_meshes.size();
        std::vector<MissRecord> missRecords;
        for (int meshID=0;meshID<numObjects;meshID++) {
            MissRecord rec;
            OPTIX_CHECK(optixSbtRecordPackHeader(missPGs[0],&rec));
            //rec.data = nullptr; /* for now ... */
            //rec.data_poly = nullptr; /* for now ... */
            rec.data_poly.vertex = (vec3f*)poly_memory.vertexBuffer[meshID].d_pointer();
            rec.data_poly.vertex2 = (vec2f*)poly_memory.vertex2Buffer[meshID].d_pointer();
            rec.data_poly.index = (uint32_t*)poly_memory.indexBuffer[meshID].d_pointer();
            rec.data_poly.poly = (Polygon*)poly_memory.polyBuffer[meshID].d_pointer();
            missRecords.push_back(rec);
        }
        missRecordsBuffer.alloc_and_upload(missRecords);
        sbt.missRecordBase          = missRecordsBuffer.d_pointer();
        sbt.missRecordStrideInBytes = sizeof(MissRecord);
        sbt.missRecordCount         = (int)missRecords.size();

        // ------------------------------------------------------------------
        // build hitgroup records
        // ------------------------------------------------------------------
        std::vector<HitgroupRecord> hitgroupRecords;
        for (int meshID=0;meshID<numObjects;meshID++) {
            HitgroupRecord rec;
            // all meshes use the same code, so all same hit group
            OPTIX_CHECK(optixSbtRecordPackHeader(hitgroupPGs[0],&rec));
            //rec.data.color  = gdt::vec3f(0.5,1.0,0.5);
            //rec.data.vertex = (vec3f*)vertexBuffer[meshID].d_pointer();
            rec.data_poly.vertex = (vec3f*)poly_memory.vertexBuffer[meshID].d_pointer();
            rec.data_poly.vertex2 = (vec2f*)poly_memory.vertex2Buffer[meshID].d_pointer();
            rec.data_poly.index = (uint32_t*)poly_memory.indexBuffer[meshID].d_pointer();
            rec.data_poly.poly = (Polygon*)poly_memory.polyBuffer[meshID].d_pointer();
            hitgroupRecords.push_back(rec);
        }
        hitgroupRecordsBuffer.alloc_and_upload(hitgroupRecords);
        sbt.hitgroupRecordBase          = hitgroupRecordsBuffer.d_pointer();
        sbt.hitgroupRecordStrideInBytes = sizeof(HitgroupRecord);
        sbt.hitgroupRecordCount         = (int)hitgroupRecords.size();
    }

    /*! constructs the shader binding table */
    void OptixControl::buildSBTTriangle()
    {
        // first allocate device memory and upload data
        for (int meshID=0;meshID<(int)model->triangle_meshes.size();meshID++) {
            TriangleMesh &mesh = *model->triangle_meshes[meshID];
            tri_memory.polyBuffer[meshID].alloc_and_upload(mesh.poly);
            tri_memory.cdfBuffer[meshID].alloc_and_upload(mesh.cdfs);
            tri_memory.facprobBuffer[meshID].alloc_and_upload(mesh.facetProbabilities);
        }

        // ------------------------------------------------------------------
        // build raygen records
        // ------------------------------------------------------------------
        std::vector<RaygenRecordTri> raygenRecords;
        for (int i=0;i<raygenPGs.size();i++) {
            RaygenRecordTri rec;
            OPTIX_CHECK(optixSbtRecordPackHeader(raygenPGs[i],&rec));
            rec.data_tri.vertex = (vec3f*)tri_memory.vertexBuffer[0].d_pointer();
            rec.data_tri.index = (vec3i*)tri_memory.indexBuffer[0].d_pointer();
            rec.data_tri.poly = (Polygon*)tri_memory.polyBuffer[0].d_pointer();
            rec.data_tri.cdfs = (float*)tri_memory.cdfBuffer[0].d_pointer();
            rec.data_tri.facetProbabilities = (vec2f*)tri_memory.facprobBuffer[0].d_pointer();
            raygenRecords.push_back(rec);
        }
        raygenRecordsBuffer.alloc_and_upload(raygenRecords);
        sbt.raygenRecord = raygenRecordsBuffer.d_pointer();

        // ------------------------------------------------------------------
        // build miss records
        // ------------------------------------------------------------------
        int numObjects = (int)model->triangle_meshes.size();
        std::vector<MissRecordTri> missRecords;
        for (int meshID=0;meshID<numObjects;meshID++) {
            MissRecordTri rec;
            OPTIX_CHECK(optixSbtRecordPackHeader(missPGs[0],&rec));
            //rec.data = nullptr; /* for now ... */
            //rec.data_poly = nullptr; /* for now ... */
            rec.data_tri.vertex = (vec3f*)tri_memory.vertexBuffer[meshID].d_pointer();
            rec.data_tri.index = (vec3i*)tri_memory.indexBuffer[meshID].d_pointer();
            rec.data_tri.poly = (Polygon*)tri_memory.polyBuffer[meshID].d_pointer();
            missRecords.push_back(rec);
        }
        missRecordsBuffer.alloc_and_upload(missRecords);
        sbt.missRecordBase          = missRecordsBuffer.d_pointer();
        sbt.missRecordStrideInBytes = sizeof(MissRecordTri);
        sbt.missRecordCount         = (int)missRecords.size();

        // ------------------------------------------------------------------
        // build hitgroup records
        // ------------------------------------------------------------------
        std::vector<HitgroupRecordTri> hitgroupRecords;
        for (int meshID=0;meshID<numObjects;meshID++) {
            HitgroupRecordTri rec;
            // all meshes use the same code, so all same hit group
            OPTIX_CHECK(optixSbtRecordPackHeader(hitgroupPGs[0],&rec));
            //rec.data.color  = gdt::vec3f(0.5,1.0,0.5);
            //rec.data.vertex = (vec3f*)vertexBuffer[meshID].d_pointer();
            rec.data_tri.vertex = (vec3f*)tri_memory.vertexBuffer[meshID].d_pointer();
            rec.data_tri.index = (vec3i*)tri_memory.indexBuffer[meshID].d_pointer();
            rec.data_tri.poly = (Polygon*)tri_memory.polyBuffer[meshID].d_pointer();
            hitgroupRecords.push_back(rec);
        }
        hitgroupRecordsBuffer.alloc_and_upload(hitgroupRecords);
        sbt.hitgroupRecordBase          = hitgroupRecordsBuffer.d_pointer();
        sbt.hitgroupRecordStrideInBytes = sizeof(HitgroupRecordTri);
        sbt.hitgroupRecordCount         = (int)hitgroupRecords.size();
    }


    /*! upload some parts only on start */
    void OptixControl::initSimulation()
    {
        const uint32_t launchSize = launchParams.simConstants.size.x * launchParams.simConstants.size.y;
        const uint32_t nbHCBins = CORESPERMP*WARPSCHEDULERS;

        CuFacetHitCounter *hitCounter = new CuFacetHitCounter[nbHCBins*launchParams.simConstants.nbFacets]();
        uint32_t *missCounter = new uint32_t(0);

        hitCounterBuffer.upload(hitCounter,nbHCBins*launchParams.simConstants.nbFacets);
        missCounterBuffer.upload(missCounter,1);

        delete[] hitCounter;
        delete missCounter;

#ifdef DEBUGCOUNT
        uint32_t *detVal = new uint32_t[NCOUNTBINS]();
        uint32_t *uVal = new uint32_t[NCOUNTBINS]();
        uint32_t *vVal = new uint32_t[NCOUNTBINS]();

        memory_debug.detBuffer.upload(detVal,NCOUNTBINS*1);
        memory_debug.uBuffer.upload(uVal,NCOUNTBINS*1);
        memory_debug.vBuffer.upload(vVal,NCOUNTBINS*1);

        delete[] detVal;
        delete[] uVal;
        delete[] vVal;
#endif
#ifdef DEBUGPOS
        vec3f *pos = new vec3f[NBCOUNTS*LAUNCH_SIZE_X*LAUNCH_SIZE_Y]();
        uint32_t *offset = new uint32_t[LAUNCH_SIZE_X*LAUNCH_SIZE_Y]();
        memory_debug.posBuffer.upload(pos,NBCOUNTS*LAUNCH_SIZE_X*LAUNCH_SIZE_Y*1);
        memory_debug.posOffsetBuffer.upload(offset,LAUNCH_SIZE_X*LAUNCH_SIZE_Y*1);

        delete[] pos;
        delete[] offset;
#endif
#ifdef DEBUGMISS
        uint32_t *miss = new uint32_t[NMISSES*LAUNCH_SIZE_X*LAUNCH_SIZE_Y]();
        memory_debug.missBuffer.upload(miss,NMISSES*LAUNCH_SIZE_X*LAUNCH_SIZE_Y*1);

        delete[] miss;
#endif
        launchParamsBuffer.upload(&launchParams,1);
    }

    /*! render one frame */
    void OptixControl::render()
    {
        // sanity check: make sure we launch only after first resize is
        // already done:
        if (launchParams.simConstants.size.x == 0) return;

        // Initialize RNG
        //const uint32_t launchSize = launchParams.simConstants.size.x * launchParams.simConstants.size.y;
        //const uint32_t nbRand  = 10*launchSize;
        //float *randomNumbers = new float[nbRand];
        //CuFacetHitCounter *hitCounter = new CuFacetHitCounter[launchSize*launchParams.simConstants.nbFacets]();

        //RandWrapper::getRand(randomNumbers, nbRand);
        //curandState_t *states;
        //crng::initializeRand(launchParams.frame.size.x*launchParams.frame.size.y, states, randomNumbers);
        //crng::generateRand(launchParams.frame.size.x*launchParams.frame.size.y, states, randomNumbers);
        //launchParams.randomNumbers = randomNumbers;
        //randBuffer.upload(randomNumbers,nbRand);
        //hitCounterBuffer.upload(hitCounter,launchSize*launchParams.simConstants.nbFacets);
        //launchParamsBuffer.upload(&launchParams,1);

        OPTIX_CHECK(optixLaunch(/*! pipeline we're launching launch: */
                pipeline,stream,
                /*! parameters and SBT */
                launchParamsBuffer.d_pointer(),
                launchParamsBuffer.sizeInBytes,
                &sbt,
                /*! dimensions of the launch: */
                launchParams.simConstants.size.x,
                launchParams.simConstants.size.y,
                1
        ));

        // sync - make sure the frame is rendered before we download and
        // display (obviously, for a high-performance application you
        // want to use streams and double-buffering, but for this simple
        // example, this will have to do)
        CUDA_SYNC_CHECK();

        //delete[] randomNumbers;
    }

    void OptixControl::generateRand()
    {
        /*const uint32_t launchSize = launchParams.simConstants.size.x * launchParams.simConstants.size.y;
        const short nbRand  = NB_RAND*launchSize;
        float *randomNumbers = new float[nbRand];
        RandWrapper::getRand(randomNumbers, nbRand);
        uint32_t *randomOffset = new uint32_t[launchSize]();

        //curandState_t *states;
        //crng::initializeRand(launchParams.frame.size.x*launchParams.frame.size.y, states, randomNumbers);
        //crng::generateRand(launchParams.frame.size.x*launchParams.frame.size.y, states, randomNumbers);
        //launchParams.randomNumbers = randomNumbers;
        randBuffer.upload(randomNumbers,nbRand);
        delete[] randomNumbers;
        delete[] randomOffset;*/

        const unsigned int launchSize = launchParams.simConstants.size.x * launchParams.simConstants.size.y;
        const unsigned int nbRand  = NB_RAND*launchSize;

        //CUDABuffer stateBuff, randomBuff;
        //curandState_t *states;
        //TODO: Try upload or host API
        //printf( "Pre : %p --> %p \n", (void*)randBuffer.d_ptr, (void*)&randBuffer.d_ptr);
        //crng::testRand(&randomBuff.d_ptr, launchSize);
        //crng::initializeRandHost(launchSize, (float**) &randomBuff.d_ptr);
        crng::generateRandHost(launchSize, (float*) randBuffer.d_ptr);
        //printf( "Post: %p --> %p --> %p\n", (void*)randBuffer.d_ptr, (void*)&randBuffer.d_ptr, randBuffer.d_pointer());

        //std::cout<< " --- print Rand --- " << std::endl;
        //crng::printDevDataAtHost(randBuffer.d_ptr, nbRand);
        //std::cout<< " --- ---------- --- " << std::endl;
        //crng::initializeRand(launchSize, stateBuff.d_ptr, randomBuff.d_ptr);
        //crng::generateRand(launchSize, (curandState_t *)stateBuff.d_ptr, (float*)randomBuff.d_ptr);
        launchParams.randomNumbers = (float *)randBuffer.d_ptr;

        //randBuffer.upload(randomNumbers,nbRand);

        /*uint32_t *randomOffset = new uint32_t[launchSize]();
        randOffsetBuffer.upload(randomOffset,launchSize);
        delete[] randomOffset;*/
        crng::offsetBufferZeroInit(launchSize, (void*) randOffsetBuffer.d_ptr);

        //launchParamsBuffer.upload(&launchParams,1);

    }


    /*! resize buffers to given amount of threads */
    // initlaunchparams
    void OptixControl::resize(const vec2i &newSize)
    {

        const uint32_t nbRand = NB_RAND;

        // resize our cuda frame buffer
        // TODO: one counter per thread is a problem for memory
        moleculeBuffer.resize(newSize.x * newSize.y * sizeof(MolPRD));
        randBuffer.resize(nbRand*newSize.x*newSize.y*sizeof(float));
        randOffsetBuffer.resize(newSize.x*newSize.y*sizeof(uint32_t));
        hitCounterBuffer.resize(model->nbFacets_total*CORESPERMP*WARPSCHEDULERS*sizeof(CuFacetHitCounter));
        missCounterBuffer.resize(sizeof(uint32_t));

        // update the launch parameters that we'll pass to the optix
        // launch:
        launchParams.simConstants.nbRandNumbersPerThread = nbRand;
        launchParams.simConstants.scene_epsilon = SCENE_EPSILON;
        launchParams.simConstants.maxDepth  = MAX_DEPTH;
        launchParams.simConstants.size  = newSize;
        launchParams.simConstants.nbFacets  = model->nbFacets_total;
        // TODO: Total nb for indices and vertices
        launchParams.simConstants.nbIndices  = model->nbIndices_total;
        launchParams.simConstants.nbVertices  = model->nbVertices_total;
        launchParams.perThreadData.currentMoleculeData = (MolPRD*)moleculeBuffer.d_pointer();
        launchParams.perThreadData.randBufferOffset = (uint32_t*)randOffsetBuffer.d_pointer();

        crng::initializeRandHost(newSize.x * newSize.y, (float **) &randBuffer.d_ptr,  time(NULL));
        //crng::initializeRandHost(newSize.x * newSize.y, (float **) &randBuffer.d_ptr);
        launchParams.randomNumbers = (float*)randBuffer.d_pointer();
        launchParams.hitCounter = (CuFacetHitCounter*)hitCounterBuffer.d_pointer();

        launchParams.sharedData.missCounter = (uint32_t *)missCounterBuffer.d_pointer();

#ifdef DEBUGCOUNT
        memory_debug.detBuffer.resize(NCOUNTBINS*sizeof(uint32_t));
        memory_debug.uBuffer.resize(NCOUNTBINS*sizeof(uint32_t));
        memory_debug.vBuffer.resize(NCOUNTBINS*sizeof(uint32_t));

        launchParams.debugCounter.detCount  = (uint32_t*)memory_debug.detBuffer.d_pointer();
        launchParams.debugCounter.uCount  = (uint32_t*)memory_debug.uBuffer.d_pointer();
        launchParams.debugCounter.vCount  = (uint32_t*)memory_debug.vBuffer.d_pointer();
#endif

#ifdef DEBUGPOS
        memory_debug.posBuffer.resize(newSize.x * newSize.y*NBCOUNTS*sizeof(vec3f));
        memory_debug.posOffsetBuffer.resize(newSize.x * newSize.y*sizeof(uint32_t));
        launchParams.perThreadData.positionsBuffer_debug = (vec3f*)memory_debug.posBuffer.d_pointer();
        launchParams.perThreadData.posOffsetBuffer_debug = (uint32_t*)memory_debug.posOffsetBuffer.d_pointer();
#endif

#ifdef DEBUGMISS
        memory_debug.missBuffer.resize(NMISSES*newSize.x * newSize.y*sizeof(uint32_t));
        launchParams.perThreadData.missBuffer = (uint32_t*)memory_debug.missBuffer.d_pointer();
#endif
        //TODO: Only necessary if changes will be made after initialization
        launchParamsBuffer.upload(&launchParams,1);


        // -- one time init of device data --
        MolPRD *perThreadData = new MolPRD[newSize.x*newSize.y]();
        moleculeBuffer.upload(perThreadData, newSize.x * newSize.y);
        delete[] perThreadData;

    }
    
    /*static vec3f dir_max(-10.e10);
    static vec3f dir_min(10.e10);*/
    static std::vector<unsigned int> counter2;
    static std::vector<unsigned int> absorb;
    static std::vector<unsigned int> desorb;

    /*! download the rendered color buffer and return the total amount of hits (= followed rays) */
    unsigned long long int OptixControl::downloadDataFromDevice(/*uint32_t *h_pixels*/)
    {
        //colorBuffer.download(h_pixels,launchParams.frame.size.x*launchParams.frame.size.y);
        // get the 'per thread' data
        MolPRD* hit = new MolPRD[launchParams.simConstants.size.x * launchParams.simConstants.size.y];
        CuFacetHitCounter* hitCounter = new CuFacetHitCounter[model->nbFacets_total * CORESPERMP*WARPSCHEDULERS];
        uint32_t * missCounter = new uint32_t;
        moleculeBuffer.download(hit, launchParams.simConstants.size.x * launchParams.simConstants.size.y);
        hitCounterBuffer.download(hitCounter, model->nbFacets_total * CORESPERMP*WARPSCHEDULERS);
        missCounterBuffer.download(missCounter, 1);


#ifdef DEBUGCOUNT
        uint32_t* detCounter = new uint32_t[NCOUNTBINS];
        uint32_t* uCounter = new uint32_t[NCOUNTBINS];
        uint32_t* vCounter = new uint32_t[NCOUNTBINS];

        memory_debug.detBuffer.download(detCounter, NCOUNTBINS);
        memory_debug.uBuffer.download(uCounter, NCOUNTBINS);
        memory_debug.vBuffer.download(vCounter, NCOUNTBINS);

        /*std::cout << "Determinant Distribution:"<<std::endl; for(int i=0;i<NBCOUNTS;i++) std::cout << "["<< ((float)i/NBCOUNTS)*(DETHIGH-DETLOW)+DETLOW << "] " <<detCounter[i] << std::endl;
        std::cout << "U Distribution:"<<std::endl; for(int i=0;i<NBCOUNTS;i++) std::cout << "["<< ((float)i/NBCOUNTS)*(UHIGH-ULOW)+ULOW << "] " <<uCounter[i] << std::endl;
        std::cout << "V Distribution:"<<std::endl; for(int i=0;i<NBCOUNTS;i++) std::cout << "["<< ((float)i/NBCOUNTS)*(VHIGH-VLOW)+VLOW << "] " <<vCounter[i] << std::endl;*/
        std::ofstream detfile,ufile,vfile;
        detfile.open ("det_counter.txt");ufile.open ("u_counter.txt");vfile.open ("v_counter.txt");

        for(int i=0;i<NCOUNTBINS;i++) detfile << "" << ((float)i/NBCOUNTS)*(DETHIGH - DETLOW) + DETLOW << " " << detCounter[i] << std::endl;
        for(int i=0;i<NCOUNTBINS;i++) ufile << "" << ((float)i/NBCOUNTS)*(UHIGH - ULOW) + ULOW << " " << uCounter[i] << std::endl;
        for(int i=0;i<NCOUNTBINS;i++) vfile << "" << ((float)i/NBCOUNTS)*(VHIGH - VLOW) + VLOW << " " << vCounter[i] << std::endl;
        detfile.close();
        ufile.close();
        vfile.close();

        delete[] detCounter;
        delete[] uCounter;
        delete[] vCounter;
#endif

#ifdef DEBUGPOS
        vec3f* positions = new vec3f[NBCOUNTS*LAUNCH_SIZE_X*LAUNCH_SIZE_Y];
        memory_debug.posBuffer.download(positions, NBCOUNTS*LAUNCH_SIZE_X*LAUNCH_SIZE_Y);
        uint32_t * posOffset = new uint32_t[LAUNCH_SIZE_X*LAUNCH_SIZE_Y];
        memory_debug.posOffsetBuffer.download(posOffset, LAUNCH_SIZE_X*LAUNCH_SIZE_Y);

        std::ofstream posFile;
        posFile.open ("debug_positions.txt");

        int nbPos = NBCOUNTS;

        std::cout << "LIMIT " <<NBCOUNTS<< std::endl;

        for(int i=0;i<NBCOUNTS*LAUNCH_SIZE_X*LAUNCH_SIZE_Y;){
            //posFile << i/(NBCOUNTS) << " " << posOffset[i/(NBCOUNTS)] << " ";
            posFile <<"{";
            for(int pos=0;pos<30;pos++){
                size_t index = i/(NBCOUNTS)*NBCOUNTS+pos;
                posFile <<"{"<<positions[index].x << "," << positions[index].y << "," << positions[index].z <<"}";
                if(pos!=30-1)posFile <<",";
                //posFile <<positions[index].x << "," << positions[index].y << "," << positions[index].z <<"   ";
            }
            i+=nbPos;
            posFile <<"},"<<std::endl;
        }
        posFile.close();
        std::cout << "POS closed " <<nbPos << std::endl;

        delete[] positions;
        delete[] posOffset;
#endif

        counter2.clear();
        counter2.resize(this->model->nbFacets_total);

        absorb.clear();
        absorb.resize(this->model->nbFacets_total);

        desorb.clear();
        desorb.resize(this->model->nbFacets_total);

        //std::ofstream myfile;
        //myfile.open ("gpu_counter_blockid.txt");

        for(unsigned int i = 0; i < launchParams.simConstants.nbFacets * CORESPERMP*WARPSCHEDULERS; i++) {
            unsigned int facIndex = i%launchParams.simConstants.nbFacets;
            counter2[facIndex] += hitCounter[i].nbMCHit; // let misses count as 0 (-1+1)
            absorb[facIndex] += hitCounter[i].nbAbsEquiv; // let misses count as 0 (-1+1)
            desorb[facIndex] += hitCounter[i].nbDesorbed; // let misses count as 0 (-1+1)

            /*if(hitCounter[i].nbMCHit > 0 || hitCounter[i].nbAbsEquiv > 0){
                std::cout << "["<<i/(launchParams.simConstants.nbFacets)<<"] on facet #"<<i%launchParams.simConstants.nbFacets<<" has total hits >>> "<< hitCounter[i].nbMCHit<< " / total abs >>> " << hitCounter[i].nbAbsEquiv<<" ---> "<< i<<std::endl;
            }*/
        }
        //myfile.close();
        unsigned long long int total_counter = 0;
        unsigned long long int total_abs = 0;
        unsigned long long int total_des = 0;
        for(unsigned int i = 0; i < this->model->nbFacets_total; i++){
            if(counter2[i]>0 || absorb[i]> 0 || desorb[i]>0) std::cout << i+1 << " " << counter2[i]<<" " << desorb[i]<<" " << absorb[i]<<std::endl;
            total_counter += counter2[i];
            total_abs += absorb[i];
            total_des += desorb[i];
        }
        std::cout << " total hits >>> "<< total_counter<<std::endl;
        std::cout << " total  abs >>> "<< total_abs<<std::endl;
        std::cout << " total  des >>> "<< total_des<<std::endl;
        std::cout << " total miss >>> "<< *missCounter<< " -- miss/hit ratio: "<<(double)(*missCounter) / total_counter <<std::endl;

        delete missCounter;
        delete[] hitCounter;
        delete[] hit;

        return total_counter;
    }

    void OptixControl::cleanup()
    {
        for (int meshID=0;meshID<model->triangle_meshes.size();meshID++) {
            tri_memory.vertexBuffer[meshID].free();
            tri_memory.indexBuffer[meshID].free();
            tri_memory.polyBuffer[meshID].free();
            tri_memory.cdfBuffer[meshID].free();
            tri_memory.facprobBuffer[meshID].free();
        }

        for (int meshID=0;meshID<model->poly_meshes.size();meshID++) {
            poly_memory.aabbBuffer[meshID].free();
            poly_memory.vertex2Buffer[meshID].free();
            poly_memory.vertexBuffer[meshID].free();
            poly_memory.indexBuffer[meshID].free();
            poly_memory.polyBuffer[meshID].free();
            poly_memory.cdfBuffer[meshID].free();
            poly_memory.facprobBuffer[meshID].free();
        }
        raygenRecordsBuffer.free();
        missRecordsBuffer.free();
        hitgroupRecordsBuffer.free();

        moleculeBuffer.free();
        crng::destroyRandHost((float**) &randBuffer.d_ptr),
        //randBuffer.free();
        randOffsetBuffer.free();
        hitCounterBuffer.free();
        missCounterBuffer.free();

#ifdef DEBUGCOUNT
        memory_debug.detBuffer.free();
        memory_debug.uBuffer.free();
        memory_debug.vBuffer.free();
#endif

    }

} // ::flowgpu
