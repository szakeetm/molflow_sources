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
#include <Helper/Chronometer.h>

#ifndef GIT_COMMIT_HASH
#define GIT_COMMIT_HASH "?"
#endif

namespace {

    class SimulationFixture : public ::testing::TestWithParam<std::string> {
    protected:
        SimulationFixture() {
            // You can do set-up work for each test here.
        }

        ~SimulationFixture() override {
            // You can do clean-up work that doesn't throw exceptions here.
        }

        void SetUp() override {
            // Code here will be called immediately after the constructor (right
            // before each test).
        }

        void TearDown() override {
            // Code here will be called immediately after each test (right
            // before the destructor).
        }

        // Objects declared here can be used by all tests in the test suite for Foo.
    };

    INSTANTIATE_TEST_CASE_P(
            Performance,
            SimulationFixture,
            ::testing::Values(
                    "test_lr1000_pipe.xml", "test_lr10_pipe_tex.xml"
            ));

    TEST_P(SimulationFixture, PerformanceOkay) {
        std::string testFile = GetParam();
        printf("Filename: %s\n",testFile.c_str());
        size_t nbFails = 0;
        bool fastEnough = false;
        do {
            SimulationManager simManager("molflow", "MFLW");
            SimulationModel model{};
            GlobalSimuState globState{};

            std::vector<char *> argv = {"tester", "--config", "simulation.cfg", "--reset", "--file"};
            char * fileName_c = new char[testFile.size() + 1];
            std::copy(testFile.begin(), testFile.end(), fileName_c);
            fileName_c[testFile.size()] = '\0';
            argv.push_back(fileName_c);
            char **args = argv.data();
            Initializer::init(argv.size(), (args), &simManager, &model, &globState);
            delete[] fileName_c;

            size_t oldHitsNb = globState.globalHits.globalHits.hit.nbMCHit;
            size_t oldDesNb = globState.globalHits.globalHits.hit.nbDesorbed;

            EXPECT_NO_THROW(simManager.StartSimulation());

            Chronometer simTimer;
            simTimer.Start();
            double elapsedTime;

            bool endCondition = false;
            do {
                ProcessSleep(1000);
                elapsedTime = simTimer.Elapsed();
                if (model.otfParams.desorptionLimit != 0)
                    endCondition = globState.globalHits.globalHits.hit.nbDesorbed >= model.otfParams.desorptionLimit;
                // Check for potential time end
                if (Settings::simDuration > 0) {
                    endCondition |= elapsedTime >= Settings::simDuration;
                }
            } while (!endCondition);
            simTimer.Stop();

            // Stop and copy results
            simManager.StopSimulation();
            simManager.KillAllSimUnits();

            double hitPS = (double) (globState.globalHits.globalHits.hit.nbMCHit - oldHitsNb) / (elapsedTime);
            EXPECT_EQ(0, oldDesNb);
            EXPECT_EQ(0, oldHitsNb);
            EXPECT_LT(0, globState.globalHits.globalHits.hit.nbDesorbed);
            EXPECT_LT(0, globState.globalHits.globalHits.hit.nbMCHit);

            double prev = -1;
            std::string testName(::testing::UnitTest::GetInstance()->current_test_info()->name());
            std::string timeRecFile = "./time_record_" + testFile.substr(0,testFile.size() - 4) + ".txt";
            {
                std::ifstream ifs(timeRecFile);
                std::string fromCommit;
                ifs >> fromCommit;
                ifs >> prev;
            }
            if (prev < 0) prev = 1e5;

            //EXPECT_GT(hitPS, 0.9 * prev);
            fastEnough = hitPS > 0.9 * prev;
            if(!fastEnough) {
                ++nbFails;
            }

            printf("Current Hit/s: %e vs. Prev Hit/s: %e\n", hitPS, prev);

            if (hitPS > prev) {
                std::string hash = GIT_COMMIT_HASH;
                std::ofstream ofs(timeRecFile);
                ofs << hash.substr(0, 8) << ' ' << hitPS << std::endl;
            }
        } while (!fastEnough && nbFails < 3);
        EXPECT_LT(nbFails, 3);
    }

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