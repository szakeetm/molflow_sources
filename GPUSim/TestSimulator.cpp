//
// Created by pbahr on 15/11/2019.
//

#include "SimulationOptiX.h"
// debug output
#include <fstream>
#include <cereal/archives/xml.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/tuple.hpp>

#include <chrono>
#include "GPUDefines.h"

template<class Archive>
void serialize(Archive & archive,
               float2 & m)
{
    archive( m.x, m.y);
}

template<class Archive>
void serialize(Archive & archive,
               float3 & m)
{
    archive( m.x, m.y, m.z);
}

template<class Archive>
void serialize(Archive & archive,
               int3 & m)
{
    archive( m.x, m.y, m.z);
}

/*! --- Initialise model with a Molflow-exported geometry --- */
flowgpu::Model* initializeModel(std::string fileName){

    std::cout << "#GPUTestsuite: Loading input file: " << fileName << std::endl;

    flowgpu::Model* model = new flowgpu::Model();
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
            cereal::make_nvp("useMaxwell", model->useMaxwell)
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
            cereal::make_nvp("useMaxwell", model->useMaxwell)
    );
#endif

    std::cout << "#GPUTestsuite: Loading completed!" << std::endl;

    return model;
}

void printUsageAndExit( const char* argv0 )
{
    fprintf( stderr, "Usage  : %s [options]\n", argv0 );
    fprintf( stderr, "Options: --file  | -f <filename>        Specify file for model input\n" );
    fprintf( stderr, "         --help  | -h                   Print this usage message\n" );
    fprintf( stderr, "         --size  | -s <launchsize>      Set kernel launch size\n" );
    fprintf( stderr, "         --size=<width>x<height>x<depth>\n" );
    fprintf( stderr, "         --loop  | -l <nbLoops>         Set number of simulation loops\n" );
    fprintf( stderr, "         --nhit  | -n <nbHits>          Set approx. number of hits for the simulation\n" );
    //fprintf( stderr, "         --dim=<width>x<height>        Set image dimensions; defaults to 512x384\n" );
    exit(1);
}

void parseSize( const char* arg, size_t& kernelSize )
{

    // look for an 'x': <width>x<height>
    size_t width_end = strchr( arg, 'x' ) - arg;
    size_t height_begin = width_end + 1;

    if ( height_begin < strlen( arg ) )
    {
        // find the beginning of the height string/
        const char *height_arg = &arg[height_begin];

        // copy width to null-terminated string
        char width_arg[32];
        strncpy( width_arg, arg, width_end );
        width_arg[width_end] = '\0';

        // terminate the width string
        width_arg[width_end] = '\0';

        kernelSize = atoi(width_arg) * atoi(height_arg);

        size_t height_end = strchr( strchr( arg, 'x' )+1, 'x' ) - arg;
        size_t depth_begin = height_end + 1;

        if (depth_begin < strlen( arg ) )
        {
            const char *depth_arg = &arg[depth_begin];
            kernelSize *= atoi(depth_arg);
            if(kernelSize>0)
                return;
        }

        if(kernelSize>0)
            return;
    }

    std::cout << "#GPUTestsuite: Failed to parse width, height from string '" << std::string( arg ) << "'"  << std::endl;
    throw;
}

int main(int argc, char **argv) {

#ifdef WITHTRIANGLES
    std::string fileName = "test_geom_tri.xml"; // Input file
#else
    std::string fileName = "test_geom.xml"; // Input file
#endif

    size_t nbLoops = 1;               // Number of Simulation loops
    size_t launchSize = 1;                  // Kernel launch size

    for(int i = 1; i < argc; ++i ) {
        if( strcmp( argv[i], "--help" ) == 0 || strcmp( argv[i], "-h" ) == 0 ) {
            printUsageAndExit( argv[0] );
        } else if( strcmp( argv[i], "--file" ) == 0 || strcmp( argv[i], "-f" ) == 0 ) {
            if( i < argc-1 ) {
                fileName = argv[++i];
            } else {
                printUsageAndExit( argv[0] );
            }
        } else if ( strcmp( argv[i], "--size") == 0 || strcmp( argv[i], "-s" ) == 0 ) {
            if( i < argc-1 ) {
                launchSize = atoi(argv[++i]);
            } else {
                printUsageAndExit( argv[0] );
            }
        } else if ( strncmp( argv[i], "--size=", 7 ) == 0 ) {
                const char *size_arg = &argv[i][7];
                parseSize(size_arg, launchSize);

        } else if ( strcmp( argv[i], "--loop") == 0  || strcmp( argv[i], "-l" ) == 0 ) {
            if( i < argc-1 ) {
                nbLoops = atoi(argv[++i]);
            } else {
                printUsageAndExit( argv[0] );
            }
        } else if ( strcmp( argv[i], "--nhit") == 0  || strcmp( argv[i], "-n" ) == 0 ) {
            if( i < argc-1 ) {
                nbLoops = (size_t)(atof(argv[++i]) / launchSize);
            } else {
                printUsageAndExit( argv[0] );
            }
        } else {
            fprintf( stderr, "Unknown option '%s'\n", argv[i] );
            printUsageAndExit( argv[0] );
        }
    }

    SimulationOptiX gpuSim;
    flowgpu::Model* model = initializeModel(fileName);

    std::cout << "#GPUTestsuite: Loading simulation with kernel size: " << launchSize << std::endl;
    gpuSim.LoadSimulation(model, launchSize);

    std::cout << "#GPUTestsuite: Starting simulation with " << launchSize << " threads per launch => " << nbLoops << " runs "<<std::endl;

    const int printPerNRuns = nbLoops/100;

    auto start_total = std::chrono::steady_clock::now();
    auto t0 = start_total;

    for(size_t i = 0; i < nbLoops; ++i){
        //auto start = std::chrono::high_resolution_clock::now();

        gpuSim.RunSimulation();


        if((i+1)%printPerNRuns==0){
            auto t1 = std::chrono::steady_clock::now();
            std::chrono::duration<double,std::milli> elapsed = t1 - t0;
            t0 = t1;
            std::cout << "--- Run #"<<i+1<< " \t- Elapsed Time: " << elapsed.count() << " ms \t--- " << (double)launchSize * printPerNRuns / elapsed.count() / 1000.0 << " MRay/s ---" << std::endl;
        }
    }
    auto finish_total = std::chrono::steady_clock::now();

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