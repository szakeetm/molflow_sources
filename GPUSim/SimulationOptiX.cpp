//
// Created by pbahr on 04/10/2019.
//

// common gdt helper tools
#include "gdt/gdt.h"
#include "optix7.h"
#include "SimulationOptiX.h"

SimulationOptiX::SimulationOptiX(){
    tracer = nullptr;
    model = nullptr;
}

SimulationOptiX::~SimulationOptiX(){
    if(tracer){
        delete tracer;
    }
    if(model){
        delete model;
    }
}


/**
 *
 * @return 1=could not load GPU Sim, 0=successfully loaded
 */
int SimulationOptiX::LoadSimulation(flowgpu::Model* loaded_model) {
    if(loaded_model == nullptr)
        return 1;

/*#ifdef DEBUG
    std::cout << "Vertices 3D:"<<std::endl;
    for(auto vert : loaded_model->poly_meshes[0]->vertices3d){
        std::cout << vert << "\t" <<std::endl;
    }
    std::cout<<std::endl;
    std::cout << "Vertices 2D:"<<std::endl;
    for(auto vert : loaded_model->poly_meshes[0]->vertices2d){
        std::cout << vert << "\t" <<std::endl;
    }
    std::cout<<std::endl;
    std::cout << "Offset 2D:"<<std::endl;
    for(auto poly : loaded_model->poly_meshes[0]->poly){
        std::cout << poly.vertOffset << "\t" <<std::endl;
    }
    std::cout<<std::endl;
#endif*/

    try {
        model = loaded_model;
    } catch (std::runtime_error& e) {
        std::cout << GDT_TERMINAL_RED << "FATAL ERROR: " << e.what()
                  << GDT_TERMINAL_DEFAULT << std::endl;
        std::cout << "Does GPUMolflow support this geometry yet?" << std::endl;
        exit(1);
    }

    tracer = new flowgpu::MolTracer(model);
    tracer->setup();
    return 0;
}

/**
 *
 * @return 1=could not load GPU Sim, 0=successfully loaded
 */
/*int SimulationOptiX::LoadSimulation(const std::vector<Vector3d> &geomVertices, const std::vector<SuperStructure> &structures, const WorkerParams& wp, const std::vector<std::vector<std::pair<double,double>>> CDFs) {
    try {
        model = flowgpu::loadFromMolflow(geomVertices, structures, wp, CDFs);
    } catch (std::runtime_error& e) {
        std::cout << GDT_TERMINAL_RED << "FATAL ERROR: " << e.what()
                  << GDT_TERMINAL_DEFAULT << std::endl;
        std::cout << "Does GPUMolflow support this geometry yet?" << std::endl;
        exit(1);
    }

    *//*std::ofstream file( "test_geom.xml" );
    cereal::XMLOutputArchive archive( file );
    archive(
            //CEREAL_NVP(model->poly_meshes[0]->poly) ,
            CEREAL_NVP(model->poly_meshes[0]->facetProbabilities) ,
            CEREAL_NVP(model->poly_meshes[0]->cdfs) ,
            CEREAL_NVP(model->poly_meshes[0]->vertices2d) ,
            CEREAL_NVP(model->poly_meshes[0]->vertices3d) ,
            CEREAL_NVP(model->poly_meshes[0]->indices) ,
            CEREAL_NVP(model->poly_meshes[0]->nbFacets) ,
            CEREAL_NVP(model->poly_meshes[0]->nbVertices) ,
            CEREAL_NVP(model->poly_meshes[0]->nbIndices) ,
            CEREAL_NVP(model->nbFacets_total) ,
            CEREAL_NVP(model->nbVertices_total) ,
            CEREAL_NVP(model->nbIndices_total) ,
            CEREAL_NVP(model->useMaxwell) ,
            CEREAL_NVP(model->bounds.lower) ,
            CEREAL_NVP(model->bounds.upper)
            );*//*


    tracer = new flowgpu::MolTracer(model);
    tracer->setup();
    return 0;
}*/

/**
 *
 * @return 1=could not load GPU Sim, 0=successfully loaded
 */
int SimulationOptiX::RunSimulation() {

    try {
        tracer->run();
    } catch (std::runtime_error& e) {
        std::cout << GDT_TERMINAL_RED << "FATAL ERROR: " << e.what()
                  << GDT_TERMINAL_DEFAULT << std::endl;
        exit(1);
    }
    return 0;
}

/**
 *
 * @return 1=could not load GPU Sim, 0=successfully loaded
 */
unsigned long long int SimulationOptiX::GetSimulationData() {

    try {
        return tracer->fetchData();
    } catch (std::runtime_error& e) {
        std::cout << GDT_TERMINAL_RED << "FATAL ERROR: " << e.what()
                  << GDT_TERMINAL_DEFAULT << std::endl;
        exit(1);
    }
    return 0;
}

/**
 *
 * @return 1=could not load GPU Sim, 0=successfully loaded
 */
int SimulationOptiX::CloseSimulation() {

    try {
        delete tracer;
    } catch (std::runtime_error& e) {
        std::cout << GDT_TERMINAL_RED << "FATAL ERROR: " << e.what()
                  << GDT_TERMINAL_DEFAULT << std::endl;
        exit(1);
    }
    return 0;
}

/*
SimulationOptiX::SimulationOptiX(){
    window = nullptr;
    model = nullptr;
}

SimulationOptiX::~SimulationOptiX(){
    if(window){
        glfwDestroyWindow(reinterpret_cast<GLFWwindow *>(window));
        delete window;
    }
    if(model){
        delete model;
    }
}

*/
/**
 *
 * @return 1=could not load GPU Sim, 0=successfully loaded
 *//*

int SimulationOptiX::LoadSimulation(const std::vector<Vector3d> &geomVertices, const std::vector<SuperStructure> &structures, const WorkerParams& wp, const std::vector<std::vector<std::pair<double,double>>> CDFs) {
    try {
        model = flowgpu::loadFromMolflow(geomVertices, structures, wp, CDFs);
        */
/*(model = flowgpu::loadOBJ(
#ifdef _WIN32
                // on windows, visual studio creates _two_ levels of build dir
                // (x86/Release)
                "../../models/sponza.obj"
#else
        // on linux, common practice is to have ONE level of build dir
      // (say, <project>/build/)...
      "../models/sponza.obj"
#endif
        );*//*

    } catch (std::runtime_error& e) {
        std::cout << GDT_TERMINAL_RED << "FATAL ERROR: " << e.what()
                  << GDT_TERMINAL_DEFAULT << std::endl;
        std::cout << "Did you forget to copy sponza.obj and sponza.mtl into your optix7course/models directory?" << std::endl;
        exit(1);
    }

    */
/*flowgpu::Camera camera = { *//*
*/
/*model->bounds.center()-*//*
*/
/*flowgpu::vec3f(7.637f, -2.58f, -5.03f),
                                                                                        model->bounds.center()-flowgpu::vec3f(0.f, 0.f, 0.0f),
                                                                                        flowgpu::vec3f(0.f,1.f,0.f) };*//*

    flowgpu::Camera camera = {flowgpu::vec3f(11.0f, -4.6f, -4.0f),
                              model->bounds.center(),
                              flowgpu::vec3f(0.f, 1.f, 0.f) };
    */
/*flowgpu::Camera camera = { model->bounds.center(),
                               model->bounds.center()-flowgpu::vec3f(0.1f, 0.1f, 0.1f),
                       flowgpu::vec3f(0.f,1.f,0.f) };*//*


    // something approximating the scale of the world, so the
    // camera knows how much to move for any given user interaction:
    const float worldScale = length(model->bounds.span());
    const std::string windowTitle = "Optix 7 OBJ Model";
    window = new flowgpu::SampleWindow(windowTitle, model, camera, worldScale);
    return 0;
}

*/
/**
 *
 * @return 1=could not load GPU Sim, 0=successfully loaded
 *//*

int SimulationOptiX::RunSimulation() {

    try {
        window->run();
    } catch (std::runtime_error& e) {
        std::cout << GDT_TERMINAL_RED << "FATAL ERROR: " << e.what()
                  << GDT_TERMINAL_DEFAULT << std::endl;
        exit(1);
    }
    return 0;
}

*/
/**
 *
 * @return 1=could not load GPU Sim, 0=successfully loaded
 *//*

int SimulationOptiX::CloseSimulation() {

    try {
        delete window;
    } catch (std::runtime_error& e) {
        std::cout << GDT_TERMINAL_RED << "FATAL ERROR: " << e.what()
                  << GDT_TERMINAL_DEFAULT << std::endl;
        exit(1);
    }
    return 0;
}*/
