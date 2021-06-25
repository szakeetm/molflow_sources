//
// Created by pbahr on 15/11/2019.
//

#include "SimulationControllerGPU.h"
// debug output
#include <fstream>

#include <list>
#include <chrono>
#include <algorithm>
#include <iostream>
#include <filesystem>
#include "ModelReader.h"

void printUsageAndExit( const char* argv0 )
{
    fprintf( stderr, "Usage  : %s [options]\n", argv0 );
    fprintf( stderr, "Options: --file    | -f <filename>        Specify file for model input\n" );
    fprintf( stderr, "         --help    | -h                   Print this usage message\n" );
    fprintf( stderr, "         --size    | -s <launchsize>      Set kernel launch size\n" );
    fprintf( stderr, "         --size=<width>x<height>[x<depth>]\n" );
    fprintf( stderr, "         --loop    | -l <nbLoops>         Set number of simulation loops\n" );
    fprintf( stderr, "         --ndes    | -d <nbDesorptions>   Set number of desorptions for sim. end\n" );
    fprintf( stderr, "         --ndes=<d1>x<d2>[x<d3>]\n" );
    fprintf( stderr, "         --nhit    | -n <nbHits>          Set approx. number of hits for the simulation (overwrites --loop)\n" );
    fprintf( stderr, "         --quiet   | -q                   Set terminal output messages to a minimum\n" );
    fprintf( stderr, "         --time    | -t                   Time limit for simulation in seconds (e.g. 0.5)\n" );
    fprintf( stderr, "         --nprints | -i                   Print runtime output n times based on loops {default 10}\n" );
    fprintf( stderr, "      --printevery | -j                   Print runtime output every n_th loop\n" );
    fprintf( stderr, "      --depth      |                      Unimplemented\n" );
    fprintf( stderr, "      --directRand |                      Unimplemented\n" );
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

void parseDesLimits( const char* arg, std::list<size_t>& limits )
{

    // look for an 'x': <width>x<height>
    size_t lim_end = strchr( arg, 'x' ) - arg;
    size_t lim_begin = 0;
    lim_end = std::min(lim_end,strlen(arg));

    while( lim_begin < strlen( arg ) )
    {
        // find the beginning of the height string/
        //const char *this_arg = &arg[next_begin];

        // copy width to null-terminated string
        char this_arg[32];
        strncpy( this_arg, &arg[lim_begin], lim_end-lim_begin );

        // terminate the width string
        this_arg[lim_end-lim_begin] = '\0';

        limits.emplace_back(strtod(this_arg,nullptr));
        lim_begin = lim_end + 1;
        lim_end = strchr( &arg[lim_end]+1, 'x' ) - arg;
        lim_end = std::min(lim_end,strlen(arg));
    }
    if(!limits.empty())
        return;
    std::cout << "#GPUTestsuite: Failed to parse desorption limits from string '" << std::string( arg ) << "'"  << std::endl;
    throw;
}

int main(int argc, char **argv) {

#ifdef WITHTRIANGLES
    std::string fileName = "minimalout.xml"; // Input file
#else
    std::string fileName = "test_geom.xml"; // Input file
#endif

    std::list<size_t> limits; // Number of desorptions: empty=use other end condition
    size_t nbLoops = 1;               // Number of Simulation loops
    size_t launchSize = 1;                  // Kernel launch size
    size_t nPrints = 10;
    size_t printPerN = 100000;
    double timeLimit = 0.0;
    bool silentMode = false;
    flowgpu::MolflowGlobal simParams{};

    for(int i = 1; i < argc; ++i ) {
        char* p;
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
                launchSize = strtoul(argv[++i],&p,10);
            } else {
                printUsageAndExit( argv[0] );
            }
        } else if ( strncmp( argv[i], "--size=", 7 ) == 0 ) {
            const char *size_arg = &argv[i][7];
            parseSize(size_arg, launchSize);
        } else if ( strcmp( argv[i], "--loop") == 0  || strcmp( argv[i], "-l" ) == 0 ) {
            if( i < argc-1 ) {
                nbLoops = strtoul(argv[++i],&p,10);
            } else {
                printUsageAndExit( argv[0] );
            }
        } else if ( strcmp( argv[i], "--nhit") == 0  || strcmp( argv[i], "-n" ) == 0 ) {
            if( i < argc-1 ) {
                nbLoops = (size_t)(strtoul(argv[++i],&p,10) / launchSize);
            } else {
                printUsageAndExit( argv[0] );
            }
        } else if ( strncmp( argv[i], "--ndes=",7) == 0  || strcmp( argv[i], "-d" ) == 0 ) {
            if ( strncmp( argv[i], "--ndes=",7) == 0) {
                const char *size_arg = &argv[i][7];
                parseDesLimits(size_arg, limits);
                //++i;
                for(auto& lim : limits) std::cout << " lim : "<< lim << std::endl;
            }
            else if( i < argc-1 ) {
                limits.emplace_back(strtoul(argv[++i],&p,10));
            } else {
                printUsageAndExit( argv[0] );
            }
        } else if ( strcmp( argv[i], "--nprints") == 0  || strcmp( argv[i], "-i" ) == 0 ) {
            if( i < argc-1 ) {
                nPrints = strtoul(argv[++i],&p,10);
            } else {
                printUsageAndExit( argv[0] );
            }
        } else if ( strcmp( argv[i], "--time") == 0  || strcmp( argv[i], "-t" ) == 0 ) {
            if( i < argc-1 ) {
                timeLimit = strtod(argv[++i],&p);
            } else {
                printUsageAndExit( argv[0] );
            }
        } else if ( strcmp( argv[i], "--printevery") == 0  || strcmp( argv[i], "-j" ) == 0 ) {
            if( i < argc-1 ) {
                printPerN = strtoul(argv[++i],&p,10);
            } else {
                printUsageAndExit( argv[0] );
            }
        } else if ( strcmp( argv[i], "--depth") == 0) {
            if( i < argc-1 ) {
                simParams.recursiveMaxDepth = strtoul(argv[++i],&p,10);
            } else {
                printUsageAndExit( argv[0] );
            }
        } else if ( strcmp( argv[i], "--directRand") == 0) {
            simParams.randomNumberMethod = true;
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

    // Set desorption limit if used
    if(!limits.empty()) {
        model->ontheflyParams.desorptionLimit = limits.front();
        limits.pop_front();
    }
    else{
        model->ontheflyParams.desorptionLimit = 0;
    }
    // Get parsed methods etc
    model->parametersGlobal = simParams;

    std::cout << "#GPUTestsuite: Loading simulation with kernel size: " << launchSize << std::endl;
    gpuSim.LoadSimulation(model, launchSize);

    std::cout << "#GPUTestsuite: Starting simulation with " << launchSize << " threads per launch => " << nbLoops << " runs "<<std::endl;

    //const int printPerNRuns = std::max(1, static_cast<int>(nbLoops/nPrints)); // prevent n%0 operation
    uint64_t printPerNRuns = std::min(static_cast<uint64_t>(printPerN), static_cast<uint64_t>(nbLoops/nPrints)); // prevent n%0 operation
    printPerNRuns = std::max(printPerNRuns, static_cast<uint64_t>(1));

    auto start_total = std::chrono::steady_clock::now();
    auto t0 = start_total;

    double raysPerSecondMax = 0.0;
    /*double raysPerSecondSum = 0.0;
    uint64_t nRaysSum = 0.0;*/

    size_t refreshForStop = std::numeric_limits<size_t>::max();
    size_t loopN;
    for(loopN = 0; loopN < nbLoops && !gpuSim.hasEnded; ++loopN){
        //auto start = std::chrono::high_resolution_clock::now();

        gpuSim.RunSimulation();

        auto t1 = std::chrono::steady_clock::now();
        std::chrono::duration<double,std::ratio<60,1>> elapsedMinutes = t1 - start_total;
        if(timeLimit >= 1e-6 && elapsedMinutes.count() >= timeLimit)
            gpuSim.hasEnded = true;

        if((!silentMode && (loopN + 1) % printPerNRuns == 0) || refreshForStop <= loopN || gpuSim.hasEnded){
            auto t1 = std::chrono::steady_clock::now();
            std::chrono::duration<double,std::milli> elapsed = t1 - t0;
            t0 = t1;

            uint64_t nbHits = gpuSim.GetSimulationData();

            if(model->ontheflyParams.desorptionLimit != 0) {
                // add remaining steps to current loop count, this is the new approx. stop until desorption limit is reached
               refreshForStop = loopN + gpuSim.RemainingStepsUntilStop();
               std::cout << " Stopping at " << loopN << " / " << refreshForStop << std::endl;
            }
            if(gpuSim.hasEnded){
                // if there is a next des limit, handle that
                if(!limits.empty()) {
                    model->ontheflyParams.desorptionLimit = limits.front();
                    limits.pop_front();
                    gpuSim.hasEnded = false;
                    gpuSim.endCalled = false;
                    gpuSim.AllowNewParticles();
                    std::cout << " Handling next des limit " << model->ontheflyParams.desorptionLimit << std::endl;
                }
            }
            static const uint64_t nRays = launchSize * printPerNRuns;
            //double rpsRun = (double)nRays / elapsed.count() / 1000.0;
            double rpsRun = (double)(nbHits) / elapsed.count() / 1000.0;
            raysPerSecondMax = std::max(raysPerSecondMax,rpsRun);
            //raysPerSecondSum += rpsRun;
            std::cout << "--- Run #" << loopN + 1 << " \t- Elapsed Time: " << elapsed.count() / 1000.0 << " s \t--- " << rpsRun << " MRay/s ---" << std::endl;
            printf("--- Trans Prob: %e\n",gpuSim.GetTransProb(1u));
        }
    }
    auto finish_total = std::chrono::steady_clock::now();
    gpuSim.GetSimulationData(false);
    std::chrono::duration<double> elapsed = finish_total - start_total;
    if(elapsed.count() / 60.0 >= 1.0)
        std::cout << "--         Total Elapsed Time: " << elapsed.count() / 60.0 << " min ---" << std::endl;
    else
        std::cout << "--         Total Elapsed Time: " << elapsed.count() << " sec ---" << std::endl;
    std::cout << "--         Avg. Rays per second: " << (double)gpuSim.figures.total_counter*loopN/elapsed.count()/1.0e6 << " MRay/s ---" << std::endl;
    //std::cout << "--         Avg. Rays per second: " << raysPerSecondSum/(nbLoops/printPerNRuns) << " MRay/s ---" << std::endl;
    std::cout << "--         Max  Rays per second: " << raysPerSecondMax << " MRay/s ---" << std::endl;

    /*std::chrono::duration<double,std::milli> elapsed = finish_total - start_total;
    std::cout << "-- Total Elapsed Time: " << elapsed.count() << " ms ---" << std::endl;
    double raysPerSecond = (double)(gpuSim.GetSimulationData())/(elapsed.count()/1000);
    std::cout << "-- Rays per second: " << raysPerSecond << " Ray/s ---" << std::endl;*/

    //gpuSim.CloseSimulation();

    return 0;
}