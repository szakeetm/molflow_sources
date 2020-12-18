//
// Created by pbahr on 15/11/2019.
//

#include "SimulationControllerGPU.h"
// debug output
#include <fstream>

#include <chrono>
#include <algorithm>
#include <iostream>
#include <filesystem>
#include "ModelReader.h"

void printUsageAndExit( const char* argv0 )
{
    fprintf( stderr, "Usage  : %s [options]\n", argv0 );
    fprintf( stderr, "Options: --file  | -f <filename>        Specify file for model input\n" );
    fprintf( stderr, "         --help  | -h                   Print this usage message\n" );
    fprintf( stderr, "         --size  | -s <launchsize>      Set kernel launch size\n" );
    fprintf( stderr, "         --size=<width>x<height>[x<depth>]\n" );
    fprintf( stderr, "         --loop  | -l <nbLoops>         Set number of simulation loops\n" );
    fprintf( stderr, "         --nhit  | -n <nbHits>          Set approx. number of hits for the simulation\n" );
    fprintf( stderr, "         --quiet | -q                   Set terminal output messages to a minimum\n" );
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
    std::string fileName = "minimalout.xml"; // Input file
#else
    std::string fileName = "test_geom.xml"; // Input file
#endif

    size_t nbLoops = 1;               // Number of Simulation loops
    size_t launchSize = 1;                  // Kernel launch size
    bool silentMode = false;

    for(int i = 1; i < argc; ++i ) {
        if( strcmp( argv[i], "--help" ) == 0 || strcmp( argv[i], "-h" ) == 0 ) {
            printUsageAndExit( argv[0] );
        } else if( strcmp( argv[i], "--file" ) == 0 || strcmp( argv[i], "-f" ) == 0 ) {
            if( i < argc-1 ) {
                fileName = argv[++i];
                if(!std::filesystem::exists(fileName)){
                    std::cout << "File "<<fileName<<" does not exist!"<< std::endl;
                    exit(0);
                }
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
        } else if ( strcmp( argv[i], "--quiet") == 0  || strcmp( argv[i], "-q" ) == 0 ) {
            silentMode = true;
        } else {
            fprintf( stderr, "Unknown option '%s'\n", argv[i] );
            printUsageAndExit( argv[0] );
        }
    }

    SimulationControllerGPU gpuSim;
    //flowgpu::Model* model = flowgeom::initializeModel(fileName);
    flowgpu::Model* model = flowgpu::loadFromExternalSerialization(fileName);
            //flowgpu::Model* test = flowgeom::loadFromExternalSerialization("minimalout.xml");
    //delete model;
    //model = test;

    std::cout << "#GPUTestsuite: Loading simulation with kernel size: " << launchSize << std::endl;
    gpuSim.LoadSimulation(model, launchSize);

    std::cout << "#GPUTestsuite: Starting simulation with " << launchSize << " threads per launch => " << nbLoops << " runs "<<std::endl;

    const uint64_t nPrints = 10;
    //const int printPerNRuns = std::max(1, static_cast<int>(nbLoops/nPrints)); // prevent n%0 operation
    uint64_t printPerNRuns = std::min(static_cast<uint64_t>(10000), static_cast<uint64_t>(nbLoops/nPrints)); // prevent n%0 operation
    printPerNRuns = std::max(printPerNRuns, static_cast<uint64_t>(1));

    auto start_total = std::chrono::steady_clock::now();
    auto t0 = start_total;

    double raysPerSecondMax = 0.0;
    /*double raysPerSecondSum = 0.0;
    uint64_t nRaysSum = 0.0;*/

    for(size_t i = 0; i < nbLoops; ++i){
        //auto start = std::chrono::high_resolution_clock::now();

        gpuSim.RunSimulation();

        if(!silentMode && (i+1)%printPerNRuns==0){
            auto t1 = std::chrono::steady_clock::now();
            std::chrono::duration<double,std::milli> elapsed = t1 - t0;
            t0 = t1;

            uint64_t nbHits = gpuSim.GetSimulationData();
            static const uint64_t nRays = launchSize * printPerNRuns;
            //double rpsRun = (double)nRays / elapsed.count() / 1000.0;
            double rpsRun = (double)(nbHits) / elapsed.count() / 1000.0;
            raysPerSecondMax = std::max(raysPerSecondMax,rpsRun);
            //raysPerSecondSum += rpsRun;
            std::cout << "--- Run #"<<i+1<< " \t- Elapsed Time: " << elapsed.count()/1000.0 << " s \t--- " << rpsRun << " MRay/s ---" << std::endl;

            /*auto tdat0 = std::chrono::steady_clock::now();
            gpuSim.GetSimulationData();
            std::chrono::duration<double,std::milli> tdiffdat = std::chrono::steady_clock::now() - tdat0;
            std::cout << "--- Data#"<<i+1<< " \t- Elapsed Time: " << tdiffdat.count() << " ms \t--- " << std::endl;
*/
            /*// 1. method
            double rpsRun = (double)launchSize * printPerNRuns / elapsed.count() / 1000.0;
            raysPerSecondMax = std::max(raysPerSecondMax,rpsRun);
            raysPerSecondSum += rpsRun;
            std::cout << "--- Run #"<<i+1<< " \t- Elapsed Time: " << elapsed.count() << " ms \t--- " << rpsRun << " MRay/s ---" << std::endl;

            // 2. method
            uint64_t nRays = gpuSim.GetSimulationData();
            nRaysSum += nRays;
            double raysPerSecond = (double)(nRays)/(elapsed.count())/1000.0;
            std::cout << "--- Run #"<<i+1<< " \t- Elapsed Time: " << elapsed.count() << " ms \t--- " << raysPerSecond << " MRay/s (fixed) ---" << std::endl;
            */
        }
    }
    auto finish_total = std::chrono::steady_clock::now();
    gpuSim.GetSimulationData(false);
    printf("Trans Prob: %lf\n",gpuSim.GetTransProb(1u));
    std::chrono::duration<double> elapsed = finish_total - start_total;
    std::cout << "--         Total Elapsed Time: " << elapsed.count() / 60.0 << " min ---" << std::endl;
    std::cout << "--         Avg. Rays per second: " << (double)launchSize*nbLoops/elapsed.count()/1.0e6 << " MRay/s ---" << std::endl;
    //std::cout << "--         Avg. Rays per second: " << raysPerSecondSum/(nbLoops/printPerNRuns) << " MRay/s ---" << std::endl;
    std::cout << "--         Max  Rays per second: " << raysPerSecondMax << " MRay/s ---" << std::endl;

    /*std::chrono::duration<double,std::milli> elapsed = finish_total - start_total;
    std::cout << "-- Total Elapsed Time: " << elapsed.count() << " ms ---" << std::endl;
    double raysPerSecond = (double)(gpuSim.GetSimulationData())/(elapsed.count()/1000);
    std::cout << "-- Rays per second: " << raysPerSecond << " Ray/s ---" << std::endl;*/

    //gpuSim.CloseSimulation();

    return 0;
}