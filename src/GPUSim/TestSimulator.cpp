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
#include "SimulationManager.h"
#include "Initializer.h"
#include "IO/WriterXML.h"
#include "fmt/core.h"
#include <SettingsIO.h>

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
    fprintf( stderr, "         --nprints | -i                   Print runtime output n times based on -l or -t {default 10}\n" );
    fprintf( stderr, "      --printevery | -j                   Print runtime output every n_th loop\n" );
    fprintf( stderr, "  --printEveryNMin | -k                   Print runtime output every k minutes\n" );
    fprintf( stderr, "      --depth      |                      Recursive max depth for secondary rays (reflections)\n" );
    fprintf( stderr, "    --cyclesForRNG | -r                   Number of cycles the RNG should be buffered for {default 1}\n" );
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

// helper class for flexible argument initialization
class CharPVec {
public:
    CharPVec() {};

    CharPVec(int size) : vec(size, nullptr) {};

    CharPVec(const std::vector<std::string> &svec) : vec(svec.size(), nullptr) {
        cpy(svec);
    };

    ~CharPVec() { clear(); };

    char **data() { return vec.data(); }

    void cpy(const std::vector<std::string> &svec) {
        if (vec.size() != svec.size()) vec.resize(svec.size());
        for (int i = 0; i < vec.size(); ++i) {
            vec[i] = new char[svec[i].size() + 1];
            std::snprintf(vec[i], std::size(svec[i])+1, "%s", svec[i].c_str());
        }
    }

    void clear() {
        for (auto &c: vec) delete[] c;
    }

    // member
    std::vector<char *> vec;
};

// TODO: model->tri_facetOffset all 0
int main(int argc, char **argv) {

#ifdef WITHTRIANGLES
    std::string fileName = "minimalout.xml"; // Input file
#else
    std::string fileName = "test_geom.xml"; // Input file
#endif

    std::list<size_t> limits; // Number of desorptions: empty=use other end condition
    size_t nbLoops = 0;               // Number of Simulation loops
    size_t launchSize = 1;                  // Kernel launch size
    size_t nPrints = 10;
    size_t printPerN = 100000;
    double printEveryNMinutes = 0.0; // timeLimit / -i or direct by -k
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
        } else if ( strcmp( argv[i], "--printEveryNMin") == 0  || strcmp( argv[i], "-k" ) == 0 ) {
            if( i < argc-1 ) {
                printEveryNMinutes = strtod(argv[++i],&p);
            } else {
                printUsageAndExit( argv[0] );
            }
        } else if ( strcmp( argv[i], "--depth") == 0) {
            if( i < argc-1 ) {
                simParams.recursiveMaxDepth = strtoul(argv[++i],&p,10);
            } else {
                printUsageAndExit( argv[0] );
            }
        }
#ifdef RNG_BULKED
        else if ( strcmp( argv[i], "--cyclesForRNG") == 0 || strcmp( argv[i], "-r") == 0) {
            if( i < argc-1 ) {
                simParams.cyclesRNG = strtoul(argv[++i],&p,10);
            } else {
                printUsageAndExit( argv[0] );
            }
        }
#endif
        else if ( strcmp( argv[i], "--directRand") == 0) {
            simParams.randomNumberMethod = true;
        } else if ( strcmp( argv[i], "--quiet") == 0  || strcmp( argv[i], "-q" ) == 0 ) {
            silentMode = true;
        } else {
            fprintf( stderr, "Unknown option '%s'\n", argv[i] );
            printUsageAndExit( argv[0] );
        }
    }

    SimulationControllerGPU gpuSim;
    std::shared_ptr<flowgpu::Model> model;

    // for export
    GlobalSimuState globState{};
    std::shared_ptr<SimulationModel> simModel = std::make_shared<SimulationModel>();

    //flowgpu::Model* model = flowgeom::initializeModel(fileName);
    {
        SimulationManager simManager{};
        simManager.interactiveMode = true;

        std::vector<std::string> argv = {"tester", "--reset", "-t", "60",
                                         "--file", fileName,
                                         "--outputPath", fileName+"out.xml"};

        CharPVec argc_v(argv);
        char **args = argc_v.data();
       /* Initializer::initFromArgv(argv.size(), (args), &simManager, model);
        Initializer::initFromFile(&simManager, model, &globState);
*/
        if(-1 < Initializer::initFromArgv(argv.size(), (args), &simManager, simModel)){
#if defined(USE_MPI)
            MPI_Finalize();
#endif
            return 41;
        }

#if defined(USE_MPI)
        MFMPI::mpi_transfer_simu();
#endif

        if(Initializer::initFromFile(&simManager, simModel, &globState)){
#if defined(USE_MPI)
            MPI_Finalize();
#endif
            return 42;
        }


        model = flowgpu::loadFromSimModel(*simModel.get());
    }
    //flowgpu::Model* model = flowgpu::loadFromExternalSerialization(fileName);
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
    gpuSim.LoadSimulation(model.get(), launchSize);

    std::cout << "#GPUTestsuite: Starting simulation with " << launchSize << " threads per launch => " << nbLoops << " runs "<<std::endl;

    //const int printPerNRuns = std::max(1, static_cast<int>(nbLoops/nPrints)); // prevent n%0 operation
    uint64_t printPerNRuns = std::min(static_cast<uint64_t>(printPerN), static_cast<uint64_t>(nbLoops/nPrints)); // prevent n%0 operation
    printPerNRuns = std::max(printPerNRuns, static_cast<uint64_t>(1));

    // -i
    if(nbLoops==0 && printEveryNMinutes <= 0.0 && timeLimit > 1e-6){
            printEveryNMinutes = timeLimit / nPrints;
    }

    auto start_total = std::chrono::steady_clock::now();
    auto t0 = start_total;
    double nextPrintAtMin = printEveryNMinutes;

    double raysPerSecondMax = 0.0;
    double desPerSecondMax = 0.0;
    /*double raysPerSecondSum = 0.0;
    uint64_t nRaysSum = 0.0;*/

    size_t refreshForStop = std::numeric_limits<size_t>::max();
    size_t loopN = 0;
    do{
    //for(loopN = 0; loopN < nbLoops && !gpuSim.hasEnded; ++loopN){
        //auto start = std::chrono::high_resolution_clock::now();

        gpuSim.RunSimulation();

        auto t1 = std::chrono::steady_clock::now();
        std::chrono::duration<double,std::ratio<60,1>> elapsedMinutes = t1 - start_total;

        // Fetch end conditions
        if(nbLoops > 0 && loopN >= nbLoops)
            gpuSim.hasEnded = true;
        else if(timeLimit >= 1e-6 && elapsedMinutes.count() >= timeLimit)
            gpuSim.hasEnded = true;

        if((!silentMode && ((nbLoops > 0 && (loopN + 1) % printPerNRuns == 0) || (printEveryNMinutes > 0.0 && elapsedMinutes.count() > nextPrintAtMin)) || refreshForStop <= loopN || gpuSim.hasEnded)){
            if(printEveryNMinutes > 0.0 && elapsedMinutes.count() > nextPrintAtMin)
                nextPrintAtMin += printEveryNMinutes;

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
            printf("--- Trans Prob: %e\n",gpuSim.GetTransProb());
        }

        ++loopN;
    } while(!gpuSim.hasEnded);
    auto finish_total = std::chrono::steady_clock::now();
    gpuSim.GetSimulationData(false);
    std::chrono::duration<double> elapsed = finish_total - start_total;
    if(elapsed.count() / 60.0 >= 1.0)
        fmt::print("--         Total Elapsed Time: {} min ---\n",elapsed.count() / 60.0);
    else
        fmt::print("--         Total Elapsed Time: {} sec ---\n",elapsed.count());
    std::cout << "--         Avg. Rays per second: " << (double)gpuSim.figures.total_counter/elapsed.count()/1.0e6 << " MRay/s ---" << std::endl;
    //std::cout << "--         Avg. Rays per second: " << raysPerSecondSum/(nbLoops/printPerNRuns) << " MRay/s ---" << std::endl;
    std::cout << "--         Max  Rays per second: " << raysPerSecondMax << " MRay/s ---" << std::endl;

    // Export results
    // 1st convert from GPU types to CPU types
    // 2nd save with XML
    gpuSim.ConvertSimulationData(globState);
    FlowIO::WriterXML writer(false, false);
    pugi::xml_document newDoc;
    //newDoc.load_file(SettingsIO::workFile.c_str());
    std::string fullOutFile = "./testout.xml";

    // Copy full file description first, in case outputFile is different
    /*std::filesystem::copy_file(SettingsIO::workFile, fullOutFile,
                               std::filesystem::copy_options::overwrite_existing);*/

    writer.SaveGeometry(newDoc, simModel);
    writer.SaveSimulationState(newDoc, simModel, globState);
    if (!newDoc.save_file(fullOutFile.c_str())) {
        fmt::print(stderr,"Error writing XML file {}\n", fullOutFile);
        return 42;
    }
    /*std::chrono::duration<double,std::milli> elapsed = finish_total - start_total;
    std::cout << "-- Total Elapsed Time: " << elapsed.count() << " ms ---" << std::endl;
    double raysPerSecond = (double)(gpuSim.GetSimulationData())/(elapsed.count()/1000);
    std::cout << "-- Rays per second: " << raysPerSecond << " Ray/s ---" << std::endl;*/

    //gpuSim.CloseSimulation();

    return 0;
}