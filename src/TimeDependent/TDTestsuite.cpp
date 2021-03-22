//
// Created by Pascal Baehr on 28.04.20.
//

#include <fstream>
#include <vector>
#include <iostream>
#include <iterator>
#include <random>

#include "../../include/CLI11/CLI11.hpp"
#include "../../src_shared/Helper/Chronometer.h"
#include "../TimeMoments.h"

typedef std::pair<std::string,double> UserMoment;
typedef std::pair<double,double> Moment;
//typedef std::tuple<double,double,double,double> MomentInterval; // begin, interval, end, timewindow

struct MomentInterval {
    double start;
    double interval;
    double end;
    double timeWindow;
    size_t startIndex;
};
namespace Settings {
    std::string inputFile;
    std::string outputFile;
    static std::vector<Moment> intervals;
    static std::vector<MomentInterval> uIntervals;
    static std::vector<double> time_points;
    static std::vector<std::vector<double>> time_points_wbreak; // use when utilizing startIndex
    static std::vector<size_t> breakPoints;
}

void generatePoT(size_t nTimes, double fromTime, double toTime, double deltaMin, double deltaMax, std::vector<double>& numbers, size_t seed = 0){

    if(seed == 0) {
        std::random_device rd;  //Will be used to obtain a seed for the random number engine
        seed = rd();
    }
    std::mt19937 gen(seed); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(deltaMin, deltaMax);

    numbers.clear();
    numbers.reserve(nTimes);
    double num = fromTime;
    for(size_t n = 0 ; n < nTimes; ++n){
        num += dis(gen);
        numbers.push_back(num);

        if(num >= toTime)
            num = fromTime;
    }

    std::cout << "Generated " << numbers.size() << " numbers" << std::endl;
}

void generatePoTW(size_t nTimes, double fromTime, double toTime, double deltaMin, double deltaMax, std::vector<std::vector<double>>& numbers, size_t seed = 0){

    if(seed == 0) {
        std::random_device rd;  //Will be used to obtain a seed for the random number engine
        seed = rd();
    }
    std::mt19937 gen(seed); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(deltaMin, deltaMax);

    numbers.clear();
    double num = fromTime;
    numbers.emplace_back();
    for(size_t n = 0 ; n < nTimes; ++n){
        num += dis(gen);
        numbers.back().push_back(num);

        if(num >= toTime) {
            numbers.emplace_back();
            num = fromTime;
        }
    }

    size_t totalNumbers = 0;
    for(auto& numb : numbers){
        totalNumbers += numb.size();
    }
    std::cout << "Generated " << totalNumbers << " numbers with " <<  numbers.size() << " desorptions " << std::endl;
}

void generatePoTW2(size_t nTimes, double fromTime, double toTime, double deltaMin, double deltaMax, std::vector<double>& numbers, std::vector<size_t>& resetPoints, size_t seed = 0){

    if(seed == 0) {
        std::random_device rd;  //Will be used to obtain a seed for the random number engine
        seed = rd();
    }
    std::mt19937 gen(seed); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(deltaMin, deltaMax);

    numbers.clear();
    double num = fromTime;
    for(size_t n = 0 ; n < nTimes; ++n){
        num += dis(gen);
        numbers.push_back(num);

        if(num >= toTime) {
            resetPoints.push_back(numbers.size());
            num = fromTime;
        }
    }

    std::cout << "Generated " << numbers.size() << " numbers with " <<  resetPoints.size()+1 << " desorptions "<<std::endl;

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

void parsePoTWStop2(const std::string& fileName, std::vector<double>& numbers, std::vector<size_t>& resetPoints){
    std::ifstream timeStream(fileName);
    std::string line;
    while (std::getline(timeStream, line))
    {
        if(line.empty()) continue;
        std::istringstream iss(line);
        double val = 0.0;
        while(iss >> val){
            numbers.push_back(val);
        }
        resetPoints.push_back(numbers.size());
    }

    std::cout << "Read " << numbers.size() << " numbers with " <<  resetPoints.size() << " desorptions "<<std::endl;
}

int parseCommands(int argc, char** argv) {
    CLI::App app{"Time dependent algorithm test suite"};

    app.add_option("-f,--file", Settings::inputFile, "Required input file (XML only)")
            /*->required()*/
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

        if (lowerBound->first <= key && key <= lowerBound->second) {
            return static_cast<int>(std::distance(moments.begin() + startIndex, lowerBound) + startIndex);
        }
    }
    return -1;
}

int LookupMomentIndex(double key, const std::vector<Moment>& moments){

    if(!moments.empty()) {
        auto lowerBound = std::lower_bound(moments.begin(), moments.end(), std::make_pair(key, key));
        if(lowerBound == moments.begin())
            return -1;
        --lowerBound; //even moments.end() can be a bound

        if (lowerBound->first <= key && key <= lowerBound->second) {
            return static_cast<int>(std::distance(moments.begin(), lowerBound));
        }
    }
    return -1;
}

int quadraticSearch(double key, const std::vector<Moment>& moments){

    if(!moments.empty()) {
        int p1=0;
        int p2=0;
        int mid=0;
        int first = 0;
        int last = moments.size()-1;

        while(first <= last){
            mid = (first+last) * 0.5;
            p1 = first + (last-first) * 0.25;
            p2 = first + (last-first) * 0.75;

            // TODO: Out of the box needs 6 comparisons to check for a hit, should be reducable by nesting
            /*if(key >= moments[p1].first){ // found
                if(key <= moments[p1].second)
                    return p1;
                else if(key >= moments[mid].first){ // found
                    if(key <= moments[mid].second)
                        return mid;
                    else if(key >= moments[p2].first && key <= moments[p2].second){ // found
                        return p2;
                    }
                }
            }*/
            if(p1 > mid || p2 < mid){
                printf("Indices messed up %d < %d < %d\n", p1, mid, p2);
            }
            if(key >= moments[mid].first && key <= moments[mid].second){ // found
                return mid;
            }
            else if(key >= moments[p1].first && key <= moments[p1].second){ // found
                return p1;
            }
            else if(key >= moments[p2].first && key <= moments[p2].second){ // found
                return p2;
            }
            else if(key < moments[mid].second && key < moments[p1].first){
                last = p1 - 1;
            }
            else if(key < moments[mid].second && key > moments[p1].first){
                first = p1 + 1;
                last = mid - 1;
            }
            else if(key > moments[mid].first && key > moments[p2].second){
                first = p2 + 1;
            }
            else if(key > moments[mid].first && key < moments[p2].second){
                first = mid + 1;
                last = p2 - 1;
            }
        }
        return -1;
    }
    return -1;
}

int interpolationSearch(double key, const std::vector<Moment>& moments){

    // Find indexes of two corners
    int lo = 0, hi = (moments.size() - 1);

    // Since array is sorted, an element present
    // in array must be in range defined by corner
    while (lo <= hi && key >= moments[lo].first && key <= moments[hi].second)
    {
        if (lo == hi)
        {
            if (moments[lo].first <= key && key <= moments[lo].second) return lo;
            return -1;
        }
        // Probing the position with keeping
        // uniform distribution in mind.
        int pos = lo + (((double)(hi - lo) /
                         (moments[hi].second - moments[lo].first)) * (key - moments[lo].first));
        //printf("%d -> %d <- %d\n", lo, pos, hi);
        // Condition of target found
        if (moments[pos].first <= key && key <= moments[pos].second)
            return pos;

        // If x is larger, x is in upper part
        if (moments[pos].second < key)
            lo = pos + 1;

            // If x is smaller, x is in the lower part
        else
            hi = pos - 1;
    }
    return -1;
}

int jumpSearchProg(const std::vector<Moment>& arr, double noToSearch, int ArrayLim)
{
    int previous = 0;
    const int stepSize = std::sqrt(ArrayLim);
    int step = stepSize;
    //Step to skip elements for jumping

    while (arr[std::min(step,ArrayLim)-1].first < noToSearch)
    {
        if(step >= ArrayLim) return -1;
        previous = step;
        step += stepSize;

    }

    /*Applying linear Search and starting from the previous elements*/
    while (arr[previous].second < noToSearch)
    {
        previous++;
        /*If element has not found yet then it means element is not present in the array*/
        if (previous == std::min(step, ArrayLim)) return -1;
    }
    // if we found the element then
    if (arr[previous].first < noToSearch && noToSearch < arr[previous].second) {
        return previous;
    }

    return -1;
}

int jumpSearchProg(const std::vector<Moment>& arr, double noToSearch, int ArrayLim, int startIndex)
{
    int previous = startIndex;
    const int stepSize = std::sqrt(ArrayLim-startIndex);
    int step = startIndex + stepSize;
    //Step to skip elements for jumping

    while (arr[std::min(step,ArrayLim)-1].first < noToSearch)
    {
        if(step >= ArrayLim) return -1;
        previous = step;
        step += stepSize;

    }

    /*Applying linear Search and starting from the previous elements*/
    while (arr[previous].second < noToSearch)
    {
        previous++;
        /*If element has not found yet then it means element is not present in the array*/
        if (previous == std::min(step, ArrayLim)) return -1;
    }
    // if we found the element then
    if (arr[previous].first <= noToSearch && noToSearch <= arr[previous].second) {
        return previous;
    }

    return -1;
}

int calcSearch(double key, const std::vector<Moment>& moments, const std::vector<MomentInterval>& userMoments){

    int start = -1; // TODO: ....
    //size_t indexOffset = 0;

#ifdef DEBUG
    int controlIndex = LookupMomentIndex(key, moments);
#endif
    for(auto& uMom : userMoments){
        //printf("Parsed %e , %e , %e\n", uMom.start, uMom.interval, uMom.end);
        const double halfTimeWindow = uMom.timeWindow * 0.5;
        if(key <= uMom.end + halfTimeWindow){
            //printf("Found %e <= %e (??? %e < ???)\n", key, uMom.end + uMom.interval * 0.5, uMom.start - uMom.timeWindow * 0.5);

            if(key >= uMom.start - halfTimeWindow){
                // found it
                const double nbMoments = std::floor((uMom.end - uMom.start + 1e-3) / uMom.interval);
                start = std::min(((key - uMom.start + halfTimeWindow) / uMom.interval), nbMoments); // can go above limits on edge values
                //printf("[%e] Potential find at start %d (%d) [%e , %e] [[%e , %e]] [[[%e , %e]]]\n", key, start, nbMoments, moments[start].first, moments[start].second, uMom.start, uMom.end, uMom.start - uMom.timeWindow * 0.5 , uMom.end + uMom.timeWindow * 0.5);
                if(key <= uMom.start + start * uMom.interval - halfTimeWindow
                && key >= uMom.start + start * uMom.interval + halfTimeWindow) {
                    //printf("[%e] Calc not in window [%e , %e] , %d * %lf (%lf)\n", key, uMom.start + start * uMom.interval - uMom.timeWindow, uMom.start + start * uMom.interval + uMom.timeWindow, start, nbMoments, uMom.interval);
                    //return -1;
                    start = -1;
                }
                else{
                    start += uMom.startIndex;
                }

                if(start != -1 && !(key >= moments[start].first && key <= moments[start].second)) {
                    //printf("[%e] Calc passed [%e , %e] , %lu * %lf (%lf)\n", key, uMom.start + (start-uMom.startIndex) * uMom.interval - halfTimeWindow, uMom.start + (start-uMom.startIndex) * uMom.interval + halfTimeWindow, (start-uMom.startIndex), nbMoments, uMom.interval);
                    //printf("[%e] Moments not in window [%e , %e] [[%e , %e]]\n", key, moments[start].first, moments[start].second, uMom.start - halfTimeWindow, uMom.end + halfTimeWindow);
                    start = -1;
                }
#ifdef DEBUG
                if(start != controlIndex) {
                    printf("[%e] %d (%lf - %d) vs %d [%e , %e] [[%e , %e]]\n", key, start, std::floor((key - uMom.start + halfTimeWindow) / uMom.interval), (int)((key - uMom.start + halfTimeWindow) / uMom.interval), controlIndex, moments[start].first,
                           moments[start].second, moments[controlIndex].first, moments[controlIndex].second);
                    printf("[%e] pre   [%e , %e]\n", key, moments[start-1].first, moments[start-1].second);
                    printf("[%e] post  [%e , %e]\n", key, moments[start+1].first, moments[start+1].second);
                    printf("[%e] first [%e , %e]\n", key, moments.front().first, moments.front().second);
                    printf("[%e] last  [%e , %e]\n", key, moments.back().first, moments.back().second);
                }
#endif
                return start;
            }
            else {
                /*if(controlIndex != -1) {
                    printf("[%e] start but not end %d [%e , %e] [[%e , %e]]\n", key, controlIndex,
                           moments[controlIndex].first, moments[controlIndex].second);
                    printf("[%e] first [%e , %e]\n", key, moments.front().first, moments.front().second);
                    printf("[%e] last [%e , %e]\n", key, moments.back().first, moments.back().second);
                }*/
                return -1;
            }
        }
        //indexOffset += ((uMom.end - uMom.start) / uMom.interval + 1);
        //printf("Indexoffset %zu\n", indexOffset);

    }
    /*if(controlIndex != -1) {
        printf("[%e] no start no end %d [%e , %e] [[%e , %e]]\n", key, controlIndex, moments[controlIndex].first,
               moments[controlIndex].second);
        printf("[%e] first [%e , %e]\n", key, moments.front().first, moments.front().second);
        printf("[%e] last [%e , %e]\n", key, moments.back().first, moments.back().second);
    }*/
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
            double time = begin;
            while(time <= end + 0.1 * interval){ // TODO: Otherwise it's up to float error if we stay below or go above limit
                parsedResult.emplace_back(time, timeWindow);
                time += interval;
            }
        }
    }
    return parsedResult;
}

/**
* \brief Parses a user input and returns a vector of time moments
* \param userInput string of the form "%lf,%lf,%lf" for beginning, interval step, ending of the moment series
* \return message Type of the source (button)
*/
MomentInterval ParseUserMomentToTuple(const std::string& userInput, double timeWindow) {
    double begin, interval, end;
    MomentInterval ret{};
    ret.start = -1.0;
    ret.interval = -1.0;
    ret.end = -1.0;
    ret.timeWindow = -1.0;
    ret.startIndex = 0;

    int nb = sscanf(userInput.c_str(), "%lf,%lf,%lf", &begin, &interval, &end);
    if (nb == 3 && (begin >= 0.0) && (end > begin) && (interval < (end - begin))) {
        // First check for potential overlap due to interval<timeWindow
        if(interval >= timeWindow){
            ret.start = begin;
            ret.interval = interval;
            ret.end = end;
            ret.timeWindow = timeWindow;
            /*size_t index = 0;
            for (double time = begin; time <= end - 0.1 * interval; time += interval)
                index++;
            double time = begin;
            size_t index2 = 0;
            while(time <= end - 0.1 * interval){
                time += interval;
                index2++;
            }
            printf("Index calc from %lf to %lf [%lf]\n", begin, time, end);
            printf("Index calc: %lf vs %zu vs %zu\n", std::ceil((end - begin) / interval), index, index2);
*/
            ret.startIndex = std::floor((end - begin) / interval + 1e-3) + 1; // add one time window on interval end
        }
    }
    else if(nb == 1){
        ret.start = begin;
        ret.interval = timeWindow;
        ret.end = begin;
        ret.timeWindow = timeWindow;
        ret.startIndex = 1;
        //printf("Index calc: %lf\n", 1.0);

    }

    return ret;
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

void momentIntervalReader(const std::vector<UserMoment>& userMoments){
    std::vector<MomentInterval> parsedMoments;
    for (size_t u = 0; u != userMoments.size(); u++) {
        parsedMoments.emplace_back(ParseUserMomentToTuple(userMoments[u].first, userMoments[u].second));
    }
    size_t cummulativeIndex = 0;
    for (size_t u = 0; u != parsedMoments.size(); u++) {
        size_t index = cummulativeIndex;
        auto& moment = parsedMoments[u];
        cummulativeIndex += moment.startIndex;
        moment.startIndex = index;
        printf("[%zu] Parsed %e , %e , %e [-> %zu]\n", u, moment.start, moment.interval, moment.end, moment.startIndex);
    }

    /*auto overlapPair = TimeMoments::CheckIntervalOverlap(parsedMoments);
    if(overlapPair.first != 0 || overlapPair.second != 0){
        char tmp[128];
        if(overlapPair.second < 0)
            sprintf(tmp, "Interval length and time window would create overlap! Check line %d.", overlapPair.first+1);
        else
            sprintf(tmp, "Overlapping time window detected! Check lines %d and %d.", overlapPair.first+1,overlapPair.second+1);
        std::cerr << tmp << std::endl;

        return;
    }*/

    Settings::uIntervals.clear();
    Settings::uIntervals = parsedMoments;
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
    for(auto& newMoment : parsedMoments) {
        printf("new size: %lu\n", Settings::intervals.size());
        size_t oldSize = Settings::intervals.size();
        AddMoment(newMoment);
        printf("first: %lf , %lf\n", Settings::intervals[oldSize].first, Settings::intervals[oldSize].second);
        printf("last: %lf , %lf\n", Settings::intervals.back().first, Settings::intervals.back().second);
    }
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

    std::vector<UserMoment> uMoments;

    if(Settings::inputFile.empty()){
        uMoments = std::vector<UserMoment>{
                {"0.001,0.1,100.0", 0.1}
        };
        uMoments = std::vector<UserMoment>{
                {"0.001,0.001,100.0", 0.001}
        };
        uMoments = std::vector<UserMoment>{
                {"0.001,0.001,1.0", 0.001},
                {"1.01,0.001,10.0", 0.001},
                {"10.01,0.001,20.0", 0.001},
                {"20.01,0.001,30.0", 0.001},
                {"30.01,0.001,40.0", 0.001},
                {"40.01,0.001,50.0", 0.001},
                {"50.01,0.001,60.0", 0.001},
                {"60.01,0.001,70.0", 0.001},
                {"70.01,0.001,80.0", 0.001},
                {"80.01,0.001,90.0", 0.001},
                {"90.01,0.001,100.0", 0.001}
        };
        uMoments = std::vector<UserMoment>{
                {"0.001, 0.01,1.0", 0.01},
                {"1.1, 0.01,10.0", 0.01},
                {"10.1, 0.01,20.0", 0.01},
                {"20.1, 0.01,30.0", 0.01},
                {"30.1, 0.01,40.0", 0.01},
                {"40.1, 0.01,50.0", 0.01},
                {"50.1, 0.01,60.0", 0.01},
                {"60.1, 0.01,70.0", 0.01},
                {"70.1, 0.01,80.0", 0.01},
                {"80.1, 0.01,90.0", 0.01},
                {"90.1, 0.01,100.0", 0.01}
        };
        uMoments = std::vector<UserMoment>{
                {"0.001, 0.1,1.0", 0.01},
                {"1.1, 0.1,10.0", 0.01},
                {"10.1, 0.1,20.0", 0.01},
                {"20.1, 0.1,30.0", 0.01},
                {"30.1, 0.1,40.0", 0.01},
                {"40.1, 0.1,50.0", 0.01},
                {"50.1, 0.1,60.0", 0.01},
                {"60.1, 0.1,70.0", 0.01},
                {"70.1, 0.1,80.0", 0.01},
                {"80.1, 0.1,90.0", 0.01},
                {"90.1, 0.1,100.0", 0.01}
        };
        uMoments = std::vector<UserMoment>{
                {"0.001,0.01,100.0", 0.01}
        };
        uMoments = std::vector<UserMoment>{
                {"0.001, 0.1,1.0", 0.01},
                {"1.1, 0.001,10.0", 0.001},
                {"10.1, 0.1,20.0", 0.00001},
                {"20.1, 0.1,30.0", 0.0001},
                {"30.1, 0.1,40.0", 0.001},
                {"40.1, 0.00001,50.0", 0.00001},
                {"50.1, 0.05,60.0", 0.001},
                {"60.1, 0.1,70.0", 0.01},
                {"70.1, 0.7,80.0", 0.07},
                {"80.1, 0.1,90.0", 0.01},
                {"90.1, 0.09,100.0", 0.001}
        };
        uMoments = std::vector<UserMoment>{
                {"0.001,0.01,100.0", 0.00001}
        };
    }
    /*std::vector<UserMoment> uMoments{
            {"1.0e-13,1.0e-5,0.001", 1.0e-5}
    };*/
    /*std::vector<UserMoment> uMoments{
            {"0.1", 0.5},
            {"1,1,60", 0.5},
            {"70,10,3600", 0.5},
            {"3660,60,30000", 0.5}
    };*/
    /*std::vector<UserMoment> uMoments{
            {"0.001,0.0001,0.08", 0.0001},
            {"0.1", 0.01},
            {"1,0.001,60", 0.001},
            {"70,0.01,3600", 0.01},
            {"3660,0.06,30000", 0.06}
    };*/
    momentsReader(uMoments);
    momentIntervalReader(uMoments);

    std::vector<Moment> intervalMoments;
    parseMoments(Settings::intervals, intervalMoments);
    if(intervalMoments.empty()) exit(0);
    std::vector<size_t> timeBins(intervalMoments.size(), 0);

    //exit(0);
    if(!Settings::inputFile.empty()){
        parsePoT(Settings::inputFile, Settings::time_points);
    }
    else{
        generatePoT(1.0e6, 0.0, 100.0, 0.001, 10.0, Settings::time_points, 42424242);
    }

    // 1. Vector binary search
    Chronometer time;
    const int totalRuns = 1;
    for(int runNb = 0; runNb < totalRuns; ++runNb) {
        time.ReInit();
        time.Start();
        printf("[%.4lfms] Vector -- Binary search [START]\n", time.ElapsedMs());
        for (auto moment : Settings::time_points) {
            int ind = LookupMomentIndex(moment, intervalMoments);
            if (ind >= 0)
                ++timeBins[ind];
        }
        time.Stop();
        printf("[%.4lfms] Vector -- Binary search [ END ]\n", time.ElapsedMs());

        size_t i = 0;
        for (auto &bin : timeBins) {
            printf("%zu<>%zu ", i, timeBins[i]);
            i++;
            if ((i % 10) == 0) printf("\n");
            if (i >= 20) break;
        }
        printf("\n");
        std::vector<size_t>(intervalMoments.size(), 0).swap(timeBins);

    }

    for(int runNb = 0; runNb < totalRuns; ++runNb) {
        time.ReInit();
        time.Start();
        printf("[%.4lfms] Vector -- Quad search [START]\n", time.ElapsedMs());
        for (auto moment : Settings::time_points) {
            int ind = quadraticSearch(moment, intervalMoments);
            if (ind >= 0)
                ++timeBins[ind];
        }
        time.Stop();
        printf("[%.4lfms] Vector -- Quad search [ END ]\n", time.ElapsedMs());

        size_t i = 0;
        for (auto &bin : timeBins) {
            printf("%zu<>%zu ", i, timeBins[i]);
            i++;
            if ((i % 10) == 0) printf("\n");
            if (i >= 20) break;
        }
        printf("\n");
        std::vector<size_t>(intervalMoments.size(), 0).swap(timeBins);
    }

    for(int runNb = 0; runNb < totalRuns; ++runNb) {
        time.ReInit();
        time.Start();
        printf("[%.4lfms] Vector -- Interp search [START]\n", time.ElapsedMs());
        for (auto moment : Settings::time_points) {
            int ind = interpolationSearch(moment, intervalMoments);
            if (ind >= 0)
                ++timeBins[ind];
        }
        time.Stop();
        printf("[%.4lfms] Vector -- Interp search [ END ]\n", time.ElapsedMs());

        size_t i = 0;
        for (auto &bin : timeBins) {
            printf("%zu<>%zu ", i, timeBins[i]);
            i++;
            if ((i % 10) == 0) printf("\n");
            if (i >= 20) break;
        }
        printf("\n");
        std::vector<size_t>(intervalMoments.size(), 0).swap(timeBins);
    }

    // jump search
    for(int runNb = 0; runNb < totalRuns; ++runNb) {
        time.ReInit();
        time.Start();
        printf("[%.4lfms] Vector -- Jump search [START]\n", time.ElapsedMs());
        for (auto moment : Settings::time_points) {
            int ind = jumpSearchProg(intervalMoments, moment, intervalMoments.size());
            if (ind >= 0)
                ++timeBins[ind];
        }
        time.Stop();
        printf("[%.4lfms] Vector -- Jump search [ END ]\n", time.ElapsedMs());

        size_t i = 0;
        for (auto &bin : timeBins) {
            printf("%zu<>%zu ", i, timeBins[i]);
            i++;
            if ((i % 10) == 0) printf("\n");
            if (i >= 20) break;
        }
        printf("\n");
        std::vector<size_t>(intervalMoments.size(), 0).swap(timeBins);

    }


    for(int runNb = 0; runNb < totalRuns; ++runNb) {
        time.ReInit();
        time.Start();
        printf("[%.4lfms] Vector -- Calc search [START]\n", time.ElapsedMs());
        size_t lastIndex = 0;
        size_t u = 0;
        for (auto moment : Settings::time_points) {
            int ind = calcSearch(moment, intervalMoments, Settings::uIntervals);
            if(ind >= (int)timeBins.size()){
                printf("Calc search error: %d >= %zu for %lf (should be %d)\n", ind, timeBins.size(), moment, LookupMomentIndex(moment, intervalMoments));
            }
            else if (ind >= 0) {
                ++timeBins[ind];
                lastIndex = ind;
            }
            ++u;
            //if(u >= 10000) exit(0);
        }

        time.Stop();
        printf("[%.4lfms] Vector -- Calc search [ END ]\n", time.ElapsedMs());
        size_t i = 0;
        for (auto &bin : timeBins) {
            printf("%zu<>%zu ", i, timeBins[i]);
            i++;
            if ((i % 10) == 0) printf("\n");
            if (i >= 20) break;
        }
        printf("\n");
        std::vector<size_t>(intervalMoments.size(), 0).swap(timeBins);
    }


    //Free some memory
    Settings::time_points.clear();
    Settings::time_points_wbreak.clear();

    if(!Settings::inputFile.empty()){
        parsePoTWStop(Settings::inputFile, Settings::time_points_wbreak);
    }
    else{
        generatePoTW(1.0e6, 0.0, 100.0, 0.001, 10.0, Settings::time_points_wbreak, 42424242);
    }

    for(int runNb = 0; runNb < totalRuns; ++runNb) {
        time.ReInit();
        time.Start();
        printf("[%.4lfms] Vector -- Binary search -- Index [START]\n", time.ElapsedMs());
        size_t lastIndex = 0;
        for (const auto &particleTrace : Settings::time_points_wbreak) {
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
        printf("[%.4lfms] Vector -- Binary search -- Index [ END ]\n", time.ElapsedMs());
        size_t i = 0;
        for (auto &bin : timeBins) {
            printf("%zu<>%zu ", i, timeBins[i]);
            i++;
            if ((i % 10) == 0) printf("\n");
            if (i >= 20) break;
        }
        printf("\n");
        std::vector<size_t>(intervalMoments.size(), 0).swap(timeBins);
    }

    for(int runNb = 0; runNb < totalRuns; ++runNb) {
        time.ReInit();
        time.Start();
        printf("[%.4lfms] Vector -- Jump search -- Index [START]\n", time.ElapsedMs());
        size_t lastIndex = 0;
        for (const auto &particleTrace : Settings::time_points_wbreak) {
            for (const auto moment : particleTrace) {
                int ind = jumpSearchProg(intervalMoments, moment, intervalMoments.size(), lastIndex);
                if (ind >= 0) {
                    ++timeBins[ind];
                    lastIndex = ind;
                }
            }
            lastIndex = 0;
        }
        time.Stop();
        printf("[%.4lfms] Vector -- Jump search -- Index [ END ]\n", time.ElapsedMs());
        size_t i = 0;
        for (auto &bin : timeBins) {
            printf("%zu<>%zu ", i, timeBins[i]);
            i++;
            if ((i % 10) == 0) printf("\n");
            if (i >= 20) break;
        }
        printf("\n");
        std::vector<size_t>(intervalMoments.size(), 0).swap(timeBins);
    }

    //Free some memory
    Settings::time_points_wbreak.clear();
    //std::vector<size_t>(intervalMoments.size(), 0).swap(timeBins);
    if(!Settings::inputFile.empty()){
        parsePoTWStop2(Settings::inputFile, Settings::time_points, Settings::breakPoints);
    }
    else{
        generatePoTW2(1.0e6, 0.0, 100.0, 0.001, 10.0, Settings::time_points, Settings::breakPoints, 42424242);
    }

    for(int i = 0; i<20; ++i)
        printf("Break %d : %zu : %e\n",i,Settings::breakPoints[i], Settings::time_points[i]);
    for(int runNb = 0; runNb < totalRuns; ++runNb) {
        time.ReInit();
        time.Start();
        printf("[%.4lfms] Vector -- Binary search -- Index2 [START]\n", time.ElapsedMs());
        size_t lastIndex = 0;
        size_t desNb = 0;
        const size_t total_events = Settings::time_points.size();
        size_t breakAt = Settings::breakPoints[0];
        for (size_t eventNb = 0; eventNb < total_events; ++eventNb) {
            if(eventNb >= breakAt) {
                lastIndex = 0;
                ++desNb;
                breakAt = Settings::breakPoints[desNb];
            }

            int ind = LookupMomentIndex(Settings::time_points[eventNb], intervalMoments, lastIndex);
            if (ind >= 0) {
                ++timeBins[ind];
                lastIndex = ind;
            }

        }
        lastIndex = 0;
        time.Stop();
        printf("[%.4lfms] Vector -- Binary search -- Index2 [ END ]\n", time.ElapsedMs());
        size_t i = 0;
        for (auto &bin : timeBins) {
            printf("%zu<>%zu ", i, timeBins[i]);
            i++;
            if ((i % 10) == 0) printf("\n");
            if (i >= 20) break;
        }
        printf("\n");
        std::vector<size_t>(intervalMoments.size(), 0).swap(timeBins);

    }

    Settings::time_points.clear();
    Settings::breakPoints.clear();

    return 0;
}