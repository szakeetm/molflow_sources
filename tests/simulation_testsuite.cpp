//
// Created by pascal on 8/8/19.
//

#include <omp.h>
#include "../src_shared/SimulationManager.h"
#include "gtest/gtest.h"
#include "../src/Initializer.h"
#include "../src/ParameterParser.h"

//#define MOLFLOW_PATH ""

#include <filesystem>
#include <fstream>

// hash time to create random file name
#include <ctime>
#include <functional>

namespace {

// Tests factorial of positive numbers.
    TEST(SubProcessInit, Zero) {

        {
            SimulationManager simMan("molflow","MFLW");
            EXPECT_EQ(0, simMan.InitSimUnits());
        }

        {
            SimulationManager simMan("molflow","MFLW");
            simMan.useCPU = true;
            simMan.nbThreads = 0;
            simMan.InitSimUnits();
            EXPECT_NE(0, simMan.nbThreads);
        }

        {
            SimulationManager simMan("molflow","MFLW");
            simMan.useCPU = true;
            simMan.nbThreads = 1; // more not possible,
            simMan.InitSimUnits();
            EXPECT_EQ(1, simMan.nbThreads);
        }
    }

    TEST(SubProcessCreateAndKill, CPU) {
        {
            SimulationManager simMan("molflow","MFLW");
            simMan.useCPU = true;
            simMan.nbThreads = 1;
            simMan.InitSimUnits();
            EXPECT_EQ(1, simMan.nbThreads);
            simMan.KillAllSimUnits();
            EXPECT_EQ(0, simMan.nbThreads);
        }
    }

    /*TEST(RunningSimu, Run) {
        SimulationManager simManager("molflow","MFLW");
        SimulationModel model{};
        GlobalSimuState globState{};

        std::vector<char*> argv = {"tester", "--config", "simulation.cfg", "--reset"};
        char** args = argv.data();

        Initializer::init(argv.size(), (args), &simManager, &model, &globState);
        size_t oldHitsNb = globState.globalHits.globalHits.hit.nbMCHit;
        size_t oldDesNb = globState.globalHits.globalHits.hit.nbDesorbed;

        // Skip desorptions if limit was already reached
        if(!Settings::desLimit.empty())
        {
            size_t listSize = Settings::desLimit.size();
            for(size_t l = 0; l < listSize; ++l) {
                if (oldDesNb > Settings::desLimit.front()){
                    printf("Skipping desorption limit: %lu\n",Settings::desLimit.front());
                    Settings::desLimit.pop_front();
                }
                else{
                    printf("Starting with desorption limit: %lu from %zu\n",Settings::desLimit.front(), oldDesNb);
                    model.otfParams.desorptionLimit = Settings::desLimit.front();
                    simManager.ForwardOtfParams(&model.otfParams);
                    break;
                }
            }
            if(Settings::desLimit.empty()){
                exit(0);
            }
        }

        try {
            simManager.StartSimulation();
        }
        catch (std::runtime_error& e) {
            exit(0);
        }

        double timeNow = omp_get_wtime();
        double timeStart = omp_get_wtime();
        double timeEnd = (Settings::simDuration > 0) ? timeStart + 1.0 * Settings::simDuration : std::numeric_limits<double>::max();

        bool endCondition = false;
        do {
            ProcessSleep(1000);
            timeNow = omp_get_wtime();
            if(model.otfParams.desorptionLimit != 0)
                endCondition = globState.globalHits.globalHits.hit.nbDesorbed>= model.otfParams.desorptionLimit;
        } while(timeNow < timeEnd && !endCondition);

        // Stop and copy results
        simManager.StopSimulation();
        simManager.KillAllSimUnits();

        EXPECT_EQ(0, oldDesNb);
        EXPECT_EQ(0, oldHitsNb);
        EXPECT_LT(0, globState.globalHits.globalHits.hit.nbDesorbed);
        EXPECT_LT(0, globState.globalHits.globalHits.hit.nbMCHit);

    }*/

    TEST(Performance, Pipe100) {
        SimulationManager simManager("molflow","MFLW");
        SimulationModel model{};
        GlobalSimuState globState{};

        std::vector<char*> argv = {"tester", "--config", "simulation.cfg", "--reset"};
        char** args = argv.data();

        Initializer::init(argv.size(), (args), &simManager, &model, &globState);
        size_t oldHitsNb = globState.globalHits.globalHits.hit.nbMCHit;
        size_t oldDesNb = globState.globalHits.globalHits.hit.nbDesorbed;

        // Skip desorptions if limit was already reached
        if(!Settings::desLimit.empty()) {
            size_t listSize = Settings::desLimit.size();
            for(size_t l = 0; l < listSize; ++l) {
                if (oldDesNb > Settings::desLimit.front()){
                    printf("Skipping desorption limit: %lu\n",Settings::desLimit.front());
                    Settings::desLimit.pop_front();
                }
                else{
                    printf("Starting with desorption limit: %lu from %zu\n",Settings::desLimit.front(), oldDesNb);
                    model.otfParams.desorptionLimit = Settings::desLimit.front();
                    simManager.ForwardOtfParams(&model.otfParams);
                    break;
                }
            }
            if(Settings::desLimit.empty()){
                exit(0);
            }
        }

        try {
            simManager.StartSimulation();
        }
        catch (std::runtime_error& e) {
            exit(0);
        }

        double timeNow = omp_get_wtime();
        double timeStart = omp_get_wtime();
        double timeEnd = (Settings::simDuration > 0) ? timeStart + 1.0 * Settings::simDuration : std::numeric_limits<double>::max();

        bool endCondition = false;
        do {
            ProcessSleep(1000);
            timeNow = omp_get_wtime();
            if(model.otfParams.desorptionLimit != 0)
                endCondition = globState.globalHits.globalHits.hit.nbDesorbed>= model.otfParams.desorptionLimit;
        } while(timeNow < timeEnd && !endCondition);

        // Stop and copy results
        simManager.StopSimulation();
        simManager.KillAllSimUnits();

        printf("Hit/s => %e\n", (double)(globState.globalHits.globalHits.hit.nbMCHit - oldHitsNb) / (timeNow-timeStart));
        EXPECT_EQ(0, oldDesNb);
        EXPECT_EQ(0, oldHitsNb);
        EXPECT_LT(0, globState.globalHits.globalHits.hit.nbDesorbed);
        EXPECT_LT(0, globState.globalHits.globalHits.hit.nbMCHit);
        EXPECT_GT((double)(globState.globalHits.globalHits.hit.nbMCHit - oldHitsNb) / (timeNow-timeStart), 0.1e7);
    }

    TEST(ParameterParsing, Sweep) {

        // generate hash name for tmp working file
        std::string paramFile = std::to_string(std::hash<time_t>()(time(nullptr))) + ".cfg";
        std::ofstream outfile (paramFile);

        outfile << "facet.42.opacity=0.5\n"
                   "facet.3.sticking=10.01\n"
                   "facet.50-90.outgassing=42e5\n"
                   "facet.100.temperature=290.92\n"
                   "simulation.mass=42.42";
        outfile.close();
        ParameterParser::Parse(paramFile);

        WorkerParams wp;
        ASSERT_FALSE(std::abs(wp.gasMass - 42.42) < 1e-5);
        ParameterParser::ChangeSimuParams(wp);
        ASSERT_TRUE(std::abs(wp.gasMass - 42.42) < 1e-5);


        std::vector<SubprocessFacet> facets(200);
        ASSERT_FALSE(std::abs(facets[41].sh.opacity - 0.5) < 1e-5);
        ASSERT_FALSE(std::abs(facets[2].sh.sticking - 10.01) < 1e-5);
        ASSERT_FALSE(std::abs(facets[49].sh.outgassing - 42e5) < 1e-5); // first
        ASSERT_FALSE(std::abs(facets[69].sh.outgassing - 42e5) < 1e-5); // mid
        ASSERT_FALSE(std::abs(facets[89].sh.outgassing - 42e5) < 1e-5); // last
        ASSERT_FALSE(std::abs(facets[99].sh.temperature - 290.92) < 1e-5);
        ParameterParser::ChangeFacetParams(facets);
        ASSERT_DOUBLE_EQ(facets[41].sh.opacity, 0.5);
        ASSERT_DOUBLE_EQ(facets[2].sh.sticking, 10.01);
        ASSERT_DOUBLE_EQ(facets[49].sh.outgassing, 42e5); // first
        ASSERT_DOUBLE_EQ(facets[69].sh.outgassing, 42e5); // mid
        ASSERT_DOUBLE_EQ(facets[89].sh.outgassing, 42e5); // last
        ASSERT_DOUBLE_EQ(facets[99].sh.temperature, 290.92);

        std::filesystem::remove(paramFile);
    }
}  // namespace

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);



    return RUN_ALL_TESTS();
}