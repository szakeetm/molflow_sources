//
// Created by pascal on 8/8/19.
//

#include "../src_shared/SimulationManager.h"
#include "../src/Initializer.h"
#include "../src/ParameterParser.h"

#include <filesystem>
#include <fstream>

// hash time to create random file name
#include <ctime>
#include <Helper/Chronometer.h>
#include <memory>
#include <cmath>
#include <IO/CSVExporter.h>
#include <SettingsIO.h>
#include <fmt/core.h>

// helper class for flexible argument initialization
class CharPVec {
public:
    CharPVec() = default;

    [[maybe_unused]] explicit CharPVec(size_t size) : vec(size, nullptr) {};

    explicit CharPVec(const std::vector<std::string> &str_vec) : vec(str_vec.size(), nullptr) {
        cpy(str_vec);
    };

    ~CharPVec() { clear(); };

    char **data() { return vec.data(); }

    void cpy(const std::vector<std::string> &str_vec) {
        if (vec.size() != str_vec.size()) vec.resize(str_vec.size());
        for (size_t i = 0; i < vec.size(); ++i) {
            vec[i] = new char[str_vec[i].size() + 1];
            std::snprintf(vec[i], std::size(str_vec[i]) + 1, "%s", str_vec[i].c_str());
        }
    }

    void clear() {
        for (auto &c: vec) delete[] c;
    }

    // member
    std::vector<char *> vec;
};

enum class BenchAlgo {
    ALGO_BVH_SAH,
    ALGO_BVH_Prob,
    ALGO_KD_SAH,
    ALGO_KD_SAH_ROPE,
    ALGO_KD_SAH_ROPERESTART,
    ALGO_KD_Prob,
    ALGO_KD_Prob_ROPE,
    ALGO_KD_Prob_ROPERESTART
};

int main(int argc, char **argv) {

    static std::unordered_map<BenchAlgo, std::string> const tableDetail = {
            {BenchAlgo::ALGO_BVH_SAH,            "BVHxSAH"},
            {BenchAlgo::ALGO_BVH_Prob,            "BVHxProb"},
            {BenchAlgo::ALGO_KD_SAH,      "KDxSAH"},
            {BenchAlgo::ALGO_KD_SAH_ROPE,       "KDxSAHxRope"},
            {BenchAlgo::ALGO_KD_SAH_ROPERESTART,     "KDxSAHxRopeRestart"},
            {BenchAlgo::ALGO_KD_Prob,      "KDxProb"},
            {BenchAlgo::ALGO_KD_Prob_ROPE,       "KDxProbxRope"},
            {BenchAlgo::ALGO_KD_Prob_ROPERESTART,     "KDxProbxRopeRestart"}
    };

    std::string outPath = "AlgoPath_" + std::to_string(std::hash<time_t>()(time(nullptr)));

    size_t nbFails = 0;
    bool fastEnough = false;
    const size_t nRuns = 5;
    const size_t keepNEntries = 20;
    const size_t runForTSec = 40;
    const int n_algo = 4;

    std::map<std::string, std::vector<std::pair<int,double>>> perfTimes;
    std::filesystem::create_directory(outPath);
    std::string timeRecFile = std::filesystem::path(outPath).append("time_bench.txt").string();
    //std::string timeRecFile = std::filesystem::path("time_bench.txt").string();
    fmt::print("Writing to {}\n", timeRecFile);
    std::ofstream ofs(timeRecFile);
    ofs << fmt::format("[filename] algorithm time\n");
    ofs.close();

    for(auto const& dir_entry : std::filesystem::directory_iterator{"./AlgoCases"}){
        if(!(dir_entry.path().extension() == ".zip" || dir_entry.path().extension() == ".xml")){
            continue;
        }

        std::string testFile = dir_entry.path().string();
        fmt::print("Filename: {}\n", testFile.c_str());
        perfTimes.emplace(testFile, std::vector<std::pair<int,double>>());
        std::vector<BenchAlgo> run_algos {
            BenchAlgo::ALGO_BVH_SAH,
            BenchAlgo::ALGO_KD_SAH,
            BenchAlgo::ALGO_KD_SAH_ROPE,
            BenchAlgo::ALGO_KD_SAH_ROPERESTART
            /*BenchAlgo::ALGO_BVH_Prob,
            BenchAlgo::ALGO_KD_Prob,
            BenchAlgo::ALGO_KD_Prob_ROPE,
            BenchAlgo::ALGO_KD_Prob_ROPERESTART*/
        };

        for(auto current_algo : run_algos){
            std::shared_ptr<SimulationModel> model = std::make_shared<SimulationModel>();
            if(current_algo == BenchAlgo::ALGO_BVH_SAH){
                model->wp.accel_type = 0;
                model->wp.splitMethod = static_cast<int>(BVHAccel::SplitMethod::SAH);
            }
            else if(current_algo == BenchAlgo::ALGO_KD_SAH){
                model->wp.accel_type = 1;
                model->wp.splitMethod = static_cast<int>(KdTreeAccel::SplitMethod::SAH);
                model->wp.kd_with_ropes = false;
            }
            else if(current_algo == BenchAlgo::ALGO_KD_SAH_ROPE){
                model->wp.accel_type = 1;
                model->wp.splitMethod = static_cast<int>(KdTreeAccel::SplitMethod::SAH);
                model->wp.kd_with_ropes = true;
                model->wp.kd_restart_ropes = false;
            }
            else if(current_algo == BenchAlgo::ALGO_KD_SAH_ROPERESTART){
                model->wp.accel_type = 1;
                model->wp.splitMethod = static_cast<int>(KdTreeAccel::SplitMethod::SAH);
                model->wp.kd_with_ropes = true;
                model->wp.kd_restart_ropes = true;
            }
            else if(current_algo == BenchAlgo::ALGO_BVH_Prob){
                model->wp.accel_type = 0;
                model->wp.splitMethod = static_cast<int>(BVHAccel::SplitMethod::TestSplit);
            }
            else if(current_algo == BenchAlgo::ALGO_KD_Prob){
                model->wp.accel_type = 1;
                model->wp.splitMethod = static_cast<int>(KdTreeAccel::SplitMethod::HybridSplit);
                model->wp.kd_with_ropes = false;
            }
            else if(current_algo == BenchAlgo::ALGO_KD_Prob_ROPE){
                model->wp.accel_type = 1;
                model->wp.splitMethod = static_cast<int>(KdTreeAccel::SplitMethod::HybridSplit);
                model->wp.kd_with_ropes = true;
                model->wp.kd_restart_ropes = false;
            }
            else if(current_algo == BenchAlgo::ALGO_KD_Prob_ROPERESTART){
                model->wp.accel_type = 1;
                model->wp.splitMethod = static_cast<int>(KdTreeAccel::SplitMethod::HybridSplit);
                model->wp.kd_with_ropes = true;
                model->wp.kd_restart_ropes = true;
            }
            else{
                fmt::print("Unavailable algorithm {}\n", static_cast<int>(current_algo));
                exit(0);
            }

            SimulationManager simManager;
            //std::shared_ptr<SimulationModel> model = std::make_shared<SimulationModel>();
            GlobalSimuState globState{};

            {
                std::vector<std::string> argvec = {"algobench", "-t", fmt::format("{}", runForTSec), "--reset",
                                                 "--file", testFile,
                                                 "--outputPath", outPath};

                CharPVec argc_v(argvec);
                char **args = argc_v.data();
                Initializer::initFromArgv(static_cast<int>(argvec.size()), (args), &simManager, model);
                Initializer::initFromFile(&simManager, model, &globState);
            }

            simManager.StartSimulation();
            simManager.StopSimulation();
            simManager.StartSimulation();

            Chronometer simTimer(false);
            simTimer.Start();
            double elapsedTime;

            bool endCondition = false;
            do {
                ProcessSleep(1000);
                elapsedTime = simTimer.Elapsed();
                if (model->otfParams.desorptionLimit != 0)
                    endCondition = globState.globalHits.globalHits.nbDesorbed >= model->otfParams.desorptionLimit;
                // Check for potential time end
                if (Settings::simDuration > 0) {
                    endCondition |= elapsedTime >= static_cast<double>(Settings::simDuration);
                }
            } while (!endCondition);
            simTimer.Stop();

            // Stop and copy results
            simManager.StopSimulation();
            simManager.KillAllSimUnits();

            perfTimes[testFile].emplace_back(std::make_pair((int)current_algo, (double) globState.globalHits.globalHits.nbHitEquiv / (elapsedTime)));
            std::string out_str = fmt::format("\"{}\" {} \"{}\" {:.4e}\n", testFile, perfTimes[testFile].back().first, tableDetail.at((BenchAlgo)perfTimes[testFile].back().first), perfTimes[testFile].back().second);
            fmt::print(out_str);

            ofs.open(timeRecFile, std::ios_base::app);
            ofs << out_str;
            ofs.close();
        }
    };

    // Cleanup
    SettingsIO::cleanup_files();

    return 0;
}