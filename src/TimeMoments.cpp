//
// Helper functions to parse and check time related structures: Moment, UserMoment
//

#include "TimeMoments.h"
#include <vector>
#include <cfloat>
#include <iostream>
#include <algorithm>

/*!
 * @brief Check for 2 unsorted interval vectors (a and b), if any of the contained intervals (a_i and b_j) overlap
 * @param vecA first Vector [a_low,a_high[
 * @param vecB second Vector [b_low,b_high[
 * @return 0=no overlap, 1=overlap
 */
int TimeMoments::CheckIntervalOverlap(const std::vector<Moment>& vecA, const std::vector<Moment>& vecB) {
    if(vecA.empty() || vecB.empty())
        return 0;

    // Get min_max values for largest vector and compare with single elements of smaller vector
    if(vecA.size()>=vecB.size()) {
        double a_low = std::numeric_limits<double>::max();
        double a_high = std::numeric_limits<double>::lowest();

        for (auto &a_i : vecA) {
            a_low = std::min(a_low, a_i.first - 0.5 * a_i.second);
            a_high = std::max(a_high, a_i.first + 0.5 * a_i.second);
        }
        for (auto &b_j : vecB) {
            const double bj_low = b_j.first - 0.5 * b_j.second;
            const double bj_high = b_j.first + 0.5 * b_j.second;

            if (bj_low + DBL_EPSILON < a_high &&
                a_low + DBL_EPSILON <= bj_high) { //b_low < a_high && a_low <= b_high
                return 1; // overlap
            } else
                return 0; // no overlap
        }
    }
    else {
        double b_low = std::numeric_limits<double>::max();
        double b_high = std::numeric_limits<double>::lowest();
        for (auto &b_j : vecB) {
            b_low = std::min(b_low, b_j.first - 0.5 * b_j.second);
            b_high = std::max(b_high, b_j.first + 0.5 * b_j.second);
        }
        for (auto &a_i : vecA) {
            const double ai_low = a_i.first - 0.5 * a_i.second;
            const double ai_high = a_i.first + 0.5 * a_i.second;

            if (ai_low + DBL_EPSILON < b_high &&
                b_low + DBL_EPSILON <= ai_high) { //b_low < a_high && a_low <= b_high
                return 1; // overlap
            } else
                return 0; // no overlap
        }
    }

    return 0;
}

/*!
 * @brief Check for 2 unsorted interval vectors (a and b), if any of the contained intervals (a_i and b_j) overlap
 * @param vecA first Vector [a_low,a_high[
 * @param vecB second Vector [b_low,b_high[
 * @return 0=no overlap, 1=overlap
 */
std::pair<int, int> TimeMoments::CheckIntervalOverlap(const std::vector<std::vector<Moment>>& vecParsedMoments) {
    if(vecParsedMoments.empty())
        return std::make_pair<int,int>(0,0);

    // Overlap when parsedMoment is empty
    for(auto vec = vecParsedMoments.begin(); vec != vecParsedMoments.end(); ++vec) {
        if(vec->empty())
            return std::make_pair<int,int>(vec - vecParsedMoments.begin(),-1);
    }

    std::vector<std::pair<double,double>> intervalBoundaries;
    for(auto& vec : vecParsedMoments){
        double a_low = std::numeric_limits<double>::max();
        double a_high = std::numeric_limits<double>::lowest();

        for (auto &a_i : vec) {
            a_low = std::min(a_low, a_i.first - 0.5 * a_i.second);
            a_high = std::max(a_high, a_i.first + 0.5 * a_i.second);
        }

        intervalBoundaries.emplace_back(std::make_pair(a_low,a_high));
    }

    for(auto vecOuter = intervalBoundaries.begin(); vecOuter != intervalBoundaries.end(); vecOuter++){
        for(auto vecInner = vecOuter + 1; vecInner != intervalBoundaries.end() && vecInner != vecOuter; vecInner++){
            if ((*vecOuter).first + DBL_EPSILON < (*vecInner).second &&
                (*vecInner).first + DBL_EPSILON <= (*vecOuter).second) { //b_low < a_high && a_low <= b_high
                return std::make_pair<int,int>
                        (vecOuter-intervalBoundaries.begin(),vecInner-intervalBoundaries.begin()); // overlap
            }
        }
    }
    return std::make_pair<int,int>(0,0);
}

/**
* \brief Parses a user input and returns a vector of time moments
* \param userInput string of format "%lf,%lf,%lf" describing start, interval and end for a list of new moments
* \return vector containing parsed moments
*/
std::vector<Moment> TimeMoments::ParseMoment(const std::string& userInput, double timeWindow) {
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
        if(!(interval<timeWindow)){
            for (double time = begin; time <= end; time += interval)
                parsedResult.emplace_back(time, timeWindow);
        }
    }
    return parsedResult;
}

int TimeMoments::ParseAndCheckUserMoments(std::vector<Moment> *moments, const std::vector<UserMoment> &userMoments) {
    std::vector<std::vector<Moment>> parsedMoments;
    for (size_t u = 0; u != userMoments.size(); u++) {
        parsedMoments.emplace_back(ParseMoment(userMoments[u].first, userMoments[u].second));
    }

    auto overlapPair = CheckIntervalOverlap(parsedMoments);
    if (overlapPair.first != 0 || overlapPair.second != 0) {
        moments->clear();
        std::cerr << "Overlap in time moments detected! Check in Moments Editor (GUI)!" << std::endl;
        return 1;
        //GLMessageBox::Display("Overlap in time moments detected! Check in Moments Editor!", "Warning", GLDLG_OK, GLDLG_ICONWARNING);
    }
    else{
        for (auto &newMoment : parsedMoments) {
            AddMoment(moments, newMoment);
        }
        std::sort(moments->begin(),moments->end());
    }

    return 0;
}

/**
* \brief Adds a time serie to moments and returns the number of elements
* \param newMoments vector containing a list of new moments that should be added
* \return number of new moments that got added
*/
int TimeMoments::AddMoment(std::vector<Moment> *moments, std::vector<Moment> newMoments) {

    if(CheckIntervalOverlap(*moments, newMoments)){
        return -1; // error
    }
    int nb = newMoments.size();
    moments->insert(moments->end(),newMoments.begin(),newMoments.end());
    return nb;
}