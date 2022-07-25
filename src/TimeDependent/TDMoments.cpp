//
// Created by pascal on 7/19/22.
//

#include "TDMoments.h"
#include "fmt/core.h"
#include "TimeMoments.h"
#include "fmt/color.h"
#include <cmath>
#include <c++/8/algorithm>

namespace MFTD {
/**
 * \brief Parses a user input and returns a vector of time moments
 * \param userInput string of the form "{},{},{}" for beginning, interval step,
 * ending of the moment series \return message Type of the source (button)
 */
/*std::vector<Moment> ParseMoment(const std::string &userInput,
                                double timeWindow) {
    std::vector<Moment> parsedResult;
    double begin, interval, end;

    int nb = sscanf(userInput.c_str(), "%lf,%lf,%lf", &begin, &interval, &end);
    if (nb == 1 && (begin >= 0.0)) {
        // One moment
        parsedResult.emplace_back(begin, timeWindow);
        //} else if (nb==3 && (begin>0.0) && (end>begin) && (interval<(end-begin))
        //&& ((end-begin)/interval<300.0)) {
    } else if (nb == 3 && (begin >= 0.0) && (end > begin) &&
               (interval < (end - begin))) {
        // Range
        // First check for potential overlap due to interval<timeWindow
        if (interval >= timeWindow) {
            double time = begin;
            while (time <=
                   end + 0.1 * interval) { // TODO: Otherwise it's up to float error
                // if we stay below or go above limit
                parsedResult.emplace_back(time, timeWindow);
                time += interval;
            }
        }
    }
    return parsedResult;
}*/

// new parser to compare for calculative results, starts at start - ... instead of start
    std::vector<Moment> ParseMoment(const std::string &userInput,
                                    double timeWindow) {
        std::vector<Moment> parsedResult;
        double begin, interval, end;

        int nb = sscanf(userInput.c_str(), "%lf,%lf,%lf", &begin, &interval, &end);
        if (nb == 1 && (begin >= 0.0)) {
            // One moment
            parsedResult.emplace_back(begin, timeWindow);
            //} else if (nb==3 && (begin>0.0) && (end>begin) && (interval<(end-begin))
            //&& ((end-begin)/interval<300.0)) {
        } else if (nb == 3 && (begin >= 0.0) && (end > begin) &&
                   (interval < (end - begin))) {
            // Range
            // First check for potential overlap due to interval<timeWindow
            if (interval >= timeWindow) {
                double time = begin;
                int n_bins = std::ceil((end - begin) / interval);
                while (time <
                       end + 0.499 * timeWindow) { // TODO: Otherwise it's up to float error
                    // if we stay below or go above limit
                    parsedResult.emplace_back(time/*std::min(time,end)*/, timeWindow);
                    time += interval;
                }
            }
        }
        return parsedResult;
    }

/**
 * \brief Parses a user input and returns a vector of time moments
 * \param userInput string of the form "{},{},{}" for beginning, interval step,
 * ending of the moment series \return message Type of the source (button)
 */
    MomentInterval ParseUserMomentToTuple(const std::string &userInput,
                                          double timeWindow) {
        double begin, interval, end;
        MomentInterval ret{};
        ret.start = -1.0;
        ret.interval = -1.0;
        ret.end = -1.0;
        ret.timeWindow = -1.0;
        ret.startIndex = 0;

        int nb = sscanf(userInput.c_str(), "%lf,%lf,%lf", &begin, &interval, &end);
        if (nb == 3 && (begin >= 0.0) && (end > begin) &&
            (interval < (end - begin))) {
            // First check for potential overlap due to interval<timeWindow
            if (interval >= timeWindow) {
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
                fmt::print("Index calc from {} to {} [{}]\n", begin, time, end);
                fmt::print("Index calc: {} vs {} vs {}\n", std::ceil((end - begin) /
                interval), index, index2);
          */
                ret.startIndex = std::floor((end - begin) / interval + 1e-3) + 1;
                // add one time window on interval end
            }
        } else if (nb == 1) {
            ret.start = begin;
            ret.interval = timeWindow;
            ret.end = begin;
            ret.timeWindow = timeWindow;
            ret.startIndex = 1;
            // fmt::print("Index calc: {}\n", 1.0);
        }

        return ret;
    }

/**
 * \brief Adds a time series to moments and returns the number of elements
 * \param newMoments vector of new moments that should be inserted
 * \return amount of new moments
 */
    int AddMoment(const std::vector<Moment> &newMoments, std::vector<Moment> &intervals) {
        std::vector<Moment> &moments = intervals;

        int nb = newMoments.size();
        moments.insert(moments.end(), newMoments.begin(), newMoments.end());
        std::sort(moments.begin(), moments.end());
        return nb;
    }

    void momentIntervalReader(const std::vector<UserMoment> &userMoments, std::vector<MomentInterval> &uIntervals) {
        std::vector<MomentInterval> parsedMoments;
        for (size_t u = 0; u != userMoments.size(); u++) {
            parsedMoments.emplace_back(
                    ParseUserMomentToTuple(userMoments[u].first, userMoments[u].second));
        }
        size_t cummulativeIndex = 0;
        for (size_t u = 0; u != parsedMoments.size(); u++) {
            size_t index = cummulativeIndex;
            auto &moment = parsedMoments[u];
            cummulativeIndex += moment.startIndex;
            moment.startIndex = index;
            fmt::print("[{}] Parsed {:e} , {:e} , {:e} [{} -> {}]\n", u, moment.start,
                       moment.interval, moment.end, moment.startIndex, cummulativeIndex);
        }

        /*auto overlapPair = TimeMoments::CheckIntervalOverlap(parsedMoments);
        if(overlapPair.first != 0 || overlapPair.second != 0){
            char tmp[128];
            if(overlapPair.second < 0)
                sprintf(tmp, "Interval length and time window would create overlap!
        Check line {}.", overlapPair.first+1); else sprintf(tmp, "Overlapping time
        window detected! Check lines {} and {}.",
        overlapPair.first+1,overlapPair.second+1); std::cerr << tmp << std::endl;
    
            return;
        }*/

        uIntervals.clear();
        uIntervals = parsedMoments;
    }

    void momentsReader(const std::vector<UserMoment> &userMoments, std::vector<Moment> &intervals) {
        std::vector<std::vector<Moment>> parsedMoments;
        parsedMoments.reserve(userMoments.size());
        for (const auto &userMoment: userMoments) {
            parsedMoments.emplace_back(
                    ParseMoment(userMoment.first, userMoment.second));
        }

        auto overlapPair = TimeMoments::CheckIntervalOverlap(parsedMoments);
        if (overlapPair.first != 0 || overlapPair.second != 0) {
            std::string tmp;
            if (overlapPair.second < 0)
                tmp = fmt::format(fg(fmt::color::medium_violet_red),
                                  "Interval length and time window would create overlap! "
                                  "Check line {}.",
                                  overlapPair.first + 1);
            else
                tmp = fmt::format(
                        fg(fmt::color::medium_violet_red),
                        "Overlapping time window detected! Check lines {} and {}.",
                        overlapPair.first + 1, overlapPair.second + 1);
            fmt::print(stderr, tmp);

            return;
        }

        intervals.clear();
        for (auto &newMoment: parsedMoments) {
            // fmt::print("new size: {}\n", intervals.size());
            size_t oldSize = intervals.size();
            AddMoment(newMoment, intervals);
            fmt::print(fg(fmt::color::teal), "Split intervals: [{} , {}] --",
                       intervals[oldSize].first,
                       intervals[oldSize].second);
            fmt::print(fg(fmt::color::teal), " [{} , {}]\n",
                       intervals.back().first,
                       intervals.back().second);
        }
    }

    void parseMoments(const std::vector<Moment> &userIntervals,
                      std::vector<Moment> &outputIntervals) {
        outputIntervals.reserve(userIntervals.size());
        for (auto &moment: userIntervals) {
            outputIntervals.emplace_back(
                    std::make_pair(moment.first - (0.5 * moment.second),
                                   moment.first + (0.5 * moment.second)));
        }

        fmt::print("Parsed {} time bins from {} intervals\n", outputIntervals.size(),
                   userIntervals.size());
    }
}