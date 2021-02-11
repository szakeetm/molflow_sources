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
#include <numeric>

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

    INSTANTIATE_TEST_SUITE_P(
            Performance,
            SimulationFixture,
            ::testing::Values(
                    "test_lr1000_pipe.xml", "test_lr10_pipe_tex.xml"
            ));


    struct Stats {
        std::string commitHash;
        double min{-1.0}, max{-1.0}, avg{-1.0}, med{-1.0};
        friend std::istream& operator>>(std::istream& is, Stats& r)       {
            // Get current position
            int len = is.tellg();
            char line[256];
            is.getline(line, 256);
            // Return to position before "Read line".
            is.seekg(len ,std::ios_base::beg);

            size_t countDelimiter = 0;
            size_t pos = 0;
            while((line[pos] != '\n' && line[pos] != '\0') && pos < 256){
                if(line[pos] == ' ' || line[pos] == '\t') {
                    ++countDelimiter;
                }
                ++pos;
            }

            r = Stats();
            if (countDelimiter == 4){
                is >> r.commitHash >> r.max >> r.min >> r.med >> r.avg;
                is.getline(line, 256);
                return is;
            }
            else if(countDelimiter == 1) {
                is >> r.commitHash >> r.max;
                is.getline(line, 256);
                //is.getline(line, 256);
                return is;
            }
            else { // unknown format, skip line
                is.getline(line, 256);
                return is;
            }
        }
        friend std::ostream& operator<<(std::ostream& os, Stats const& r) {
            return os << r.commitHash << "\t" << std::scientific << std::setprecision(8) << r.max << "\t" << r.min << "\t" << r.med << "\t" << r.avg;
        }

        bool operator<(Stats const& other) const {
            return other.max < max;
        }

        /*Stats& operator=(const OldStats& src){
            commitHash = src.commitHash;
            max = src.max;
            min = src.min;
            avg = src.avg;
            med = src.med;

            return *this;
        }*/
    };

    /*struct OldStats {
        std::string commitHash;
        double min{-1.0}, max{-1.0}, avg{-1.0}, med{-1.0};
        friend std::istream& operator>>(std::istream& is, OldStats& r)       {

            // Get current position
            int len = is.tellg();
            char line[256];
            is.getline(line, 256);
            // Return to position before "Read line".
            is.seekg(len ,std::ios_base::beg);

            size_t countDelimiter = 0;
            size_t pos = 0;
            while((line[pos] != '\n' && line[pos] != '\0') && pos < 256 && countDelimiter < 2){
                if(line[pos] == ' ' || line[pos] == '\t') {
                    ++countDelimiter;
                }
                ++pos;
            }
            if(countDelimiter >= 2) {
                is >> r.commitHash >> r.max;
                is.getline(line, 256);
                return is;
            }
            else
                return is >> r.commitHash >> r.max >> r.min >> r.med >> r.avg;
        }

        bool operator<(OldStats const& other) const {
            return other.max < max;
        }

        //implicit conversion
        operator Stats() const {
            Stats newStat;
            newStat.commitHash = commitHash;
            newStat.max = max;
            newStat.min = min;
            newStat.avg = avg;
            newStat.med = med;
        return newStat; }

    };*/

    TEST_P(SimulationFixture, PerformanceOkay) {
        std::string testFile = GetParam();
        printf("Filename: %s\n",testFile.c_str());
        size_t nbFails = 0;
        bool fastEnough = false;
        const size_t nRuns = 10;
        const size_t keepNEntries = 20;
        const size_t runForTSec = 20;
        std::vector<double> perfTimes;
        for(size_t runNb = 0; runNb < nRuns; ++runNb){
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

            perfTimes.emplace_back((double) (globState.globalHits.globalHits.hit.nbMCHit - oldHitsNb) / (elapsedTime));
            EXPECT_EQ(0, oldDesNb);
            EXPECT_EQ(0, oldHitsNb);
            EXPECT_LT(0, globState.globalHits.globalHits.hit.nbDesorbed);
            EXPECT_LT(0, globState.globalHits.globalHits.hit.nbMCHit);

            printf("[Run %zu/%zu] Current Hit/s: %e\n", runNb, nRuns, perfTimes.back());
        };

        // Compare to old performance values here

        if(!perfTimes.empty()) {
            std::sort(perfTimes.begin(), perfTimes.end());

            Stats currentRun;
            currentRun.commitHash = GIT_COMMIT_HASH;
            currentRun.commitHash = currentRun.commitHash.substr(0,8); // only first 8 hex digits
            currentRun.min = perfTimes.front();
            currentRun.max = perfTimes.back();
            double sum = std::accumulate(perfTimes.begin(), perfTimes.end(), 0.0);
            currentRun.avg = sum / perfTimes.size();
            currentRun.med = (perfTimes.size() % 2 == 0)
                             ? (perfTimes[perfTimes.size() / 2 - 1] + perfTimes[perfTimes.size() / 2]) / 2
                             : perfTimes[perfTimes.size() / 2];

            std::vector<Stats> prevRun;
            std::string testName(::testing::UnitTest::GetInstance()->current_test_info()->name());
            std::string timeRecFile = "./time_record_" + testFile.substr(0, testFile.size() - 4) + ".txt";
            {
                std::ifstream ifs(timeRecFile);
                //prevRun.insert( prevRun.begin(), std::istream_iterator<OldStats>(ifs), std::istream_iterator<OldStats>() );
                prevRun.insert( prevRun.begin(), std::istream_iterator<Stats>(ifs), std::istream_iterator<Stats>() );
                //ifs >> prevRun;
                /*ifs >> fromCommit;
                ifs >> prevRun.max;
                ifs >> prevRun.min;
                ifs >> prevRun.med;
                ifs >> prevRun.avg;*/
            }

            // Check either
            fastEnough = currentRun.max >= 0.95 * prevRun.front().max;
            if(fastEnough)
                EXPECT_GE(currentRun.max, 0.95 * prevRun.front().max);
            else
                EXPECT_GE(currentRun.med, 0.95 * prevRun.front().med);
            {
                // Keep top 20 in list
                prevRun.push_back(currentRun);
                std::sort(prevRun.begin(), prevRun.end());
                for(auto iter_o = prevRun.begin(); iter_o != prevRun.end(); ++iter_o) {
                    for(auto iter_i = iter_o+1; iter_i != prevRun.end(); ++iter_i) { // "slower entry"
                        if(iter_i->commitHash == iter_o->commitHash){
                            prevRun.erase(iter_i--); // remove outer entry, as it is slower
                        }
                    }
                }
                std::ofstream ofs(timeRecFile);
                for(size_t lineN = 0; lineN < std::min(prevRun.size(), keepNEntries); ++lineN) {
                    if(!prevRun[lineN].commitHash.empty())
                        ofs << prevRun[lineN] << std::endl;
                }
            }

        }
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