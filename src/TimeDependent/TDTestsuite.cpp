//
// Created by Pascal Baehr on 28.04.20.
//

#include <fstream>
#include <vector>
#include <iostream>
#include <iterator>
#include "../../include/CLI11/CLI11.hpp"
#include "../../src_shared/Helper/Chronometer.h"
#include "../TimeMoments.h"

typedef std::pair<std::string,double> UserMoment;
typedef std::pair<double,double> Moment;

namespace Settings {
    std::string inputFile;
    std::string outputFile;
    std::vector<Moment> intervals;
    std::vector<double> time_points;
    std::vector<std::vector<double>> time_points_wbreak; // use when utilizing startIndex
}

void parsePoT(const std::string& fileName, std::vector<double>& numbers){
    std::ifstream timeStream(fileName);
    std::istream_iterator<double> start(timeStream), end;
    numbers.assign(start, end);
    std::cout << "Read " << numbers.size() << " numbers" << std::endl;
}

void parsePoTWStop(const std::string& fileName, std::vector<std::vector<double>>& numbers){
    std::ifstream timeStream(fileName);
    std::string line;
    while (std::getline(timeStream, line))
    {
        if(line.empty()) continue;
        numbers.emplace_back();
        std::istringstream iss(line);
        double val = 0.0;
        while(iss >> val){
            numbers.back().push_back(val);
        }
    }

    size_t totalNumbers = 0;
    for(auto& numb : numbers){
        totalNumbers += numb.size();
    }
    std::cout << "Read " << totalNumbers << " numbers" << std::endl;
}

int parseCommands(int argc, char** argv) {
    CLI::App app{"Time dependent algorithm test suite"};

    app.add_option("-f,--file", Settings::inputFile, "Required input file (XML only)")
            ->required()
            ->check(CLI::ExistingFile);
    app.add_option("-o,--output", Settings::outputFile, "Output file if different from input file");
    CLI11_PARSE(app, argc, argv);

    // Save to inputfile
    if(Settings::outputFile.empty())
        Settings::outputFile = "out.txt";

    return 0;
}

/*!
 * @brief Lookup the index of the interval related to a given key and a start position for accelerated lookup
 * @param key specific moment
 * @param moments vector of time intervals
 * @param startIndex offset to only look in a subset of moments
 * @return -1 if moment doesnt relate to an interval, else index of moment (+1 to account for [0]== steady state)
 * TODO: Removed +1 for test runs
 */
int LookupMomentIndex(double key, const std::vector<std::pair<double, double>>& moments, size_t startIndex){

    if(!moments.empty()) {
        auto lowerBound = std::lower_bound(moments.begin() + startIndex, moments.end(), std::make_pair(key, key));
        if(lowerBound == moments.begin())
            return -1;
        --lowerBound; //even moments.end() can be a bound

        if (lowerBound->first <= key && key < lowerBound->second) {
            return static_cast<int>(std::distance(moments.begin() + startIndex, lowerBound) + startIndex);
        }
    }
    return -1;
}

/**
* \brief Parses a user input and returns a vector of time moments
* \param userInput string of the form "%lf,%lf,%lf" for beginning, interval step, ending of the moment series
* \return message Type of the source (button)
*/
std::vector<Moment> ParseMoment(const std::string& userInput, double timeWindow) {
    std::vector<Moment> parsedResult;
    double begin, interval, end;

    int nb = sscanf(userInput.c_str(), "%lf,%lf,%lf", &begin, &interval, &end);
    if (nb == 1 && (begin >= 0.0)) {
        //One moment
        parsedResult.emplace_back(begin,timeWindow);
        //} else if (nb==3 && (begin>0.0) && (end>begin) && (interval<(end-begin)) && ((end-begin)/interval<300.0)) {
    }
    else if (nb == 3 && (begin >= 0.0) && (end > begin) && (interval < (end - begin))) {
        //Range
        // First check for potential overlap due to interval<timeWindow
        if(interval >= timeWindow){
            for (double time = begin; time <= end; time += interval)
                parsedResult.emplace_back(time, timeWindow);
        }
    }
    return parsedResult;
}


/**
* \brief Adds a time series to moments and returns the number of elements
* \param newMoments vector of new moments that should be inserted
* \return amount of new moments
*/
int AddMoment(std::vector<Moment> newMoments) {
    std::vector<Moment>& moments = Settings::intervals;

    int nb = newMoments.size();
    moments.insert(moments.end(),newMoments.begin(),newMoments.end());
    std::sort(moments.begin(),moments.end());
    return nb;
}

void momentsReader(const std::vector<UserMoment>& userMoments){
    std::vector<std::vector<Moment>> parsedMoments;
    for (size_t u = 0; u != userMoments.size(); u++) {
        parsedMoments.emplace_back(ParseMoment(userMoments[u].first, userMoments[u].second));
    }

    auto overlapPair = TimeMoments::CheckIntervalOverlap(parsedMoments);
    if(overlapPair.first != 0 || overlapPair.second != 0){
        char tmp[128];
        if(overlapPair.second < 0)
            sprintf(tmp, "Interval length and time window would create overlap! Check line %d.", overlapPair.first+1);
        else
            sprintf(tmp, "Overlapping time window detected! Check lines %d and %d.", overlapPair.first+1,overlapPair.second+1);
        std::cerr << tmp << std::endl;

        return;
    }

    Settings::intervals.clear();
    for(auto& newMoment : parsedMoments)
        AddMoment(newMoment);
}

void parseMoments(const std::vector<Moment>& userIntervals, std::vector<Moment>& outputIntervals){
    outputIntervals.reserve(userIntervals.size());
    for(auto& moment : userIntervals){
        outputIntervals.emplace_back(std::make_pair(moment.first - (0.5 * moment.second), moment.first + (0.5 * moment.second)));
    }

    std::cout << "Parsed " << outputIntervals.size() << " time bins from " << userIntervals.size() << " intervals " << std::endl;
}

int main(int argc, char** argv) {
    parseCommands(argc, argv);
    parsePoT(Settings::inputFile, Settings::time_points);
    std::vector<UserMoment> uMoments{
            {"1e-13,0.0000001,0.001",0.0000001}
    };
    momentsReader(uMoments);
    std::vector<Moment> intervalMoments;
    parseMoments(Settings::intervals, intervalMoments);

    std::vector<size_t> timeBins(intervalMoments.size(),0);
    // 1. Vector binary search
    Chronometer time;
    time.Start();
    printf("[%.4lfms] Vector -- Binary search [START]\n", time.Elapsed());
    for(auto moment : Settings::time_points){
        int ind = LookupMomentIndex(moment, intervalMoments, 0);
        if(ind >=0)
            ++timeBins[ind];
    }
    time.Stop();
    printf("[%.4lfms] Vector -- Binary search [ END ]\n", time.Elapsed());
    size_t i = 0;
    for(auto& bin : timeBins){
        printf("%zu<>%zu ", i,timeBins[i]);
        i++;
        if((i%10)==0) printf("\n");
        if(i>=20) break;
    }
    printf("\n");

    //Free some memory
    Settings::time_points.clear();
    std::vector<size_t>(intervalMoments.size(),0).swap(timeBins);
    parsePoTWStop(Settings::inputFile, Settings::time_points_wbreak);

    time.ReInit();
    time.Start();
    printf("[%.4lfms] Vector -- Binary search -- Index [START]\n", time.Elapsed());
    size_t lastIndex = 0;
    for(const auto& particleTrace : Settings::time_points_wbreak) {
        for (const auto moment : particleTrace) {
            int ind = LookupMomentIndex(moment, intervalMoments, lastIndex);
            if (ind >= 0) {
                ++timeBins[ind];
                lastIndex = ind;
            }
        }
        lastIndex = 0;
    }
    time.Stop();
    printf("[%.4lfms] Vector -- Binary search -- Index [ END ]\n", time.Elapsed());
    i = 0;
    for(auto& bin : timeBins){
        printf("%zu<>%zu ", i,timeBins[i]);
        i++;
        if((i%10)==0) printf("\n");
        if(i>=20) break;
    }
    printf("\n");
    return 0;
}