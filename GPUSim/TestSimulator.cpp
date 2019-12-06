//
// Created by pbahr on 15/11/2019.
//

#include "SimulationOptiX.h"
#include "gdt/math/vec.h"
// debug output
#include <fstream>
#include <cereal/archives/xml.hpp>
#include <cereal/types/vector.hpp>

#include <chrono>
#include "GPUDefines.h"

void main() {
    flowgpu::Model* model = new flowgpu::Model();
    //model->poly_meshes.resize(1);
    //model->poly_meshes[0] = new flowgpu::PolygonMesh();
    model->poly_meshes.push_back(new flowgpu::PolygonMesh());

    // --- Initialise model with a Molflow-exported geometry ---
    std::ifstream file( "test_geom.xml" );
    cereal::XMLInputArchive archive( file );
    archive(
            cereal::make_nvp("poly", model->poly_meshes[0]->poly) ,
            cereal::make_nvp("facetProbabilities", model->poly_meshes[0]->facetProbabilities) ,
            cereal::make_nvp("cdfs", model->poly_meshes[0]->cdfs) ,
            cereal::make_nvp("vertices2d", model->poly_meshes[0]->vertices2d) ,
            cereal::make_nvp("vertices3d", model->poly_meshes[0]->vertices3d) ,
            cereal::make_nvp("indices", model->poly_meshes[0]->indices) ,
            cereal::make_nvp("nbFacets", model->poly_meshes[0]->nbFacets) ,
            cereal::make_nvp("nbVertices", model->poly_meshes[0]->nbVertices) ,
            cereal::make_nvp("nbIndices", model->poly_meshes[0]->nbIndices) ,
            cereal::make_nvp("nbFacetsTotal", model->nbFacets_total) ,
            cereal::make_nvp("nbVerticesTotal", model->nbVertices_total) ,
            cereal::make_nvp("nbIndicesTotal", model->nbIndices_total) ,
            cereal::make_nvp("useMaxwell", model->useMaxwell) ,
            cereal::make_nvp("bounds.lower", model->bounds.lower) ,
            cereal::make_nvp("bounds.upper", model->bounds.upper)
    );

    SimulationOptiX gpuSim;
    gpuSim.LoadSimulation(model);
    auto start_total = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < N_LOOPS; ++i){
        auto start = std::chrono::high_resolution_clock::now();
        gpuSim.RunSimulation();
        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double,std::milli> elapsed = finish - start;
        std::cout << "--- Run #"<<i<< " - Elapsed Time: " << elapsed.count() << " ms ---" << std::endl;
    }
    auto finish_total = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = finish_total - start_total;
    std::cout << "-- Total Elapsed Time: " << elapsed.count() << " s ---" << std::endl;
    double raysPerSecond = (double)(gpuSim.GetSimulationData())/(elapsed.count());
    std::cout << "-- Rays per second: " << raysPerSecond/1000.0/1000.0 << " MRay/s ---" << std::endl;
    /*std::chrono::duration<double,std::milli> elapsed = finish_total - start_total;
    std::cout << "-- Total Elapsed Time: " << elapsed.count() << " ms ---" << std::endl;
    double raysPerSecond = (double)(gpuSim.GetSimulationData())/(elapsed.count()/1000);
    std::cout << "-- Rays per second: " << raysPerSecond << " Ray/s ---" << std::endl;*/

    //gpuSim.CloseSimulation();

}