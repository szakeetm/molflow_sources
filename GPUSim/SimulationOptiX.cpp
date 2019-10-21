//
// Created by pbahr on 04/10/2019.
//

// common gdt helper tools
#include "gdt/gdt.h"
#include "optix7.h"
#include "SimulationOptiX.h"

SimulationOptiX::SimulationOptiX(){

}

SimulationOptiX::~SimulationOptiX(){
    if(window){
        glfwDestroyWindow(reinterpret_cast<GLFWwindow *>(window));
        delete window;
    }
    if(model)
        delete model;
}

/**
 *
 * @return 1=could not load GPU Sim, 0=successfully loaded
 */
int SimulationOptiX::LoadSimulation(const std::vector<Vector3d> &geomVertices, const std::vector<SuperStructure> &structures) {
    try {
        model = osc::loadFromMolflow(geomVertices, structures);
        /*(model = osc::loadOBJ(
#ifdef _WIN32
                // on windows, visual studio creates _two_ levels of build dir
                // (x86/Release)
                "../../models/sponza.obj"
#else
        // on linux, common practice is to have ONE level of build dir
      // (say, <project>/build/)...
      "../models/sponza.obj"
#endif
        );*/
    } catch (std::runtime_error& e) {
        std::cout << GDT_TERMINAL_RED << "FATAL ERROR: " << e.what()
                  << GDT_TERMINAL_DEFAULT << std::endl;
        std::cout << "Did you forget to copy sponza.obj and sponza.mtl into your optix7course/models directory?" << std::endl;
        exit(1);
    }

    osc::Camera camera = { osc::vec3f(-2.0f, 0.1f, -2.0f),
            model->bounds.center()-osc::vec3f(-0.3f, 0.1f, 1.0f),
            osc::vec3f(0.f,1.f,0.f) };

    // something approximating the scale of the world, so the
    // camera knows how much to move for any given user interaction:
    const float worldScale = length(model->bounds.span());
    const std::string windowTitle = "Optix 7 OBJ Model";
    window = new osc::SampleWindow(windowTitle,model,camera,worldScale);
    return 0;
}

/**
 *
 * @return 1=could not load GPU Sim, 0=successfully loaded
 */
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

/**
 *
 * @return 1=could not load GPU Sim, 0=successfully loaded
 */
int SimulationOptiX::CloseSimulation() {

    try {
        delete window;
    } catch (std::runtime_error& e) {
        std::cout << GDT_TERMINAL_RED << "FATAL ERROR: " << e.what()
                  << GDT_TERMINAL_DEFAULT << std::endl;
        exit(1);
    }
    return 0;
}