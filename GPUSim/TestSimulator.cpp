//
// Created by pbahr on 15/11/2019.
//

#include "SimulationOptiX.h"
#include "gdt/math/vec.h"
// debug output
#include <fstream>
#include <cereal/archives/xml.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/tuple.hpp>

#include <chrono>
#include "GPUDefines.h"

int main(int argc, char **argv) {
    flowgpu::Model* model = new flowgpu::Model();
    //model->poly_meshes.resize(1);
    //model->poly_meshes[0] = new flowgpu::PolygonMesh();

    // --- Initialise model with a Molflow-exported geometry ---

    std::cout << "You have entered " << argc
         << " arguments:" << "\n";

    for (int i = 0; i < argc; ++i)
        std::cout << argv[i] << "\n";

    std::string fileName = "test_geom.xml";
    if(argc>1)
        fileName = argv[1];

    size_t nbLoops = N_LOOPS;
    if(argc>2)
        nbLoops = atoi(argv[2]);

    std::ifstream file( fileName );
    cereal::XMLInputArchive archive( file );

#ifdef WITHTRIANGLES
    model->triangle_meshes.push_back(new flowgpu::TriangleMesh());
    archive(
            cereal::make_nvp("poly", model->triangle_meshes[0]->poly) ,
            cereal::make_nvp("facetProbabilities", model->triangle_meshes[0]->facetProbabilities) ,
            cereal::make_nvp("cdfs", model->triangle_meshes[0]->cdfs) ,
            //cereal::make_nvp("vertices2d", nullptr) ,
            cereal::make_nvp("vertices3d", model->triangle_meshes[0]->vertices3d) ,
            cereal::make_nvp("indices", model->triangle_meshes[0]->indices) ,
            cereal::make_nvp("nbFacets", model->triangle_meshes[0]->nbFacets) ,
            cereal::make_nvp("nbVertices", model->triangle_meshes[0]->nbVertices) ,
            cereal::make_nvp("nbIndices", model->triangle_meshes[0]->nbIndices) ,
            cereal::make_nvp("nbFacetsTotal", model->nbFacets_total) ,
            cereal::make_nvp("nbVerticesTotal", model->nbVertices_total) ,
            cereal::make_nvp("nbIndicesTotal", model->nbIndices_total) ,
            cereal::make_nvp("useMaxwell", model->useMaxwell) ,
            cereal::make_nvp("bounds.lower", model->bounds.lower) ,
            cereal::make_nvp("bounds.upper", model->bounds.upper)
    );
#else
    model->poly_meshes.push_back(new flowgpu::PolygonMesh());
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
#endif

    SimulationOptiX gpuSim;
    gpuSim.LoadSimulation(model);
    auto start_total = std::chrono::high_resolution_clock::now();
    for(size_t i = 0; i < nbLoops; ++i){
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

    return 0;
}