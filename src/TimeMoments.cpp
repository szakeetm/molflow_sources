//
// Helper functions to parse and check time related structures: Moment, UserMoment
//

#include "TimeMoments.h"
#include <vector>
#include <cfloat>
#include <iostream>
#include <algorithm>
#include <limits>


/*!
 * @brief Check for 2 unsorted interval vectors (a and b), if any of the contained intervals (singleMoment and b_j) overlap
 * return value is the offending index (on its own or with its pair)
 */
std::optional<std::pair<int,int>> TimeMoments::HasIntervalOverlap(const std::vector<std::vector<Moment>>& moments) {
    if(moments.empty())
        return std::nullopt;

    // Overlap when parsedMoment is empty
    for(auto vec = moments.begin(); vec != moments.end(); ++vec) {
        if (vec->empty())
            return { { vec - moments.begin(),-1 } };
    }

    //Convert moments to intervals
    std::vector<Interval> intervals;
    for(auto& vec : moments){
        double startTime = std::numeric_limits<double>::max();
        double endTime = std::numeric_limits<double>::lowest();

        for (auto &singleMoment : vec) {
            startTime = std::min(startTime, singleMoment.time - 0.5 * singleMoment.window);
            endTime = std::max(endTime, singleMoment.time + 0.5 * singleMoment.window);
        }

        intervals.push_back({ startTime,endTime });
    }

    for(auto vecOuter = intervals.begin(); vecOuter != intervals.end(); vecOuter++){
        for(auto vecInner = vecOuter + 1; vecInner != intervals.end() && vecInner != vecOuter; vecInner++){
            if (vecOuter->startTime + DBL_EPSILON < vecInner->endTime &&
                vecInner->startTime + DBL_EPSILON <= vecOuter->endTime) {
                return { { vecOuter - intervals.begin(),vecInner - intervals.begin() } }; // overlap
            }
        }
    }
    return std::nullopt;
}

/**
* \brief Parses a user input and returns a vector of time moments
* \param userInput string of format "%lf,%lf,%lf" describing start, interval and end for a list of new moments
* \return vector containing parsed moments
*/
std::vector<Moment> TimeMoments::ParseUserMoment(const UserMoment& input) {
    std::vector<Moment> parsedResult;
    double begin, spacing, end;

    int nb = sscanf(input.content.c_str(), "%lf,%lf,%lf", &begin, &spacing, &end);
    if (nb == 1 && (begin >= 0.0)) {
        //One moment
        Moment m;
        m.time=begin;
        m.window=input.timeWindow;
        parsedResult.emplace_back(m);
    }
    else if (nb == 3 && (begin >= 0.0) && (end > begin) && (spacing < (end - begin))) {
        //Range
        // First check for potential overlap due to spacing<timeWindow
        if(!(spacing<input.timeWindow)){
            for (double time = begin; time <= end; time += spacing) {
                Moment m;
                m.time=time;
                m.window=input.timeWindow;
                parsedResult.emplace_back(m);
            }
        }
    }
    return parsedResult;
}

void TimeMoments::ParseAndCheckUserMoments(std::vector<Moment>& moments, const std::vector<UserMoment>& userMoments,
	GLProgress_Abstract& prg) {
	std::vector<std::vector<Moment>> parsedMoments;
	prg.SetMessage("Parsing moments...", false);
	for (size_t u = 0; u != userMoments.size(); u++) {
		parsedMoments.emplace_back(ParseUserMoment(userMoments[u]));
		prg.SetProgress(0.5 * (double)u / (double)userMoments.size());
	}

	auto overlapPair = HasIntervalOverlap(parsedMoments);
	if (overlapPair.has_value()) {
		moments.clear();
		throw Error("Overlap in time moments detected! Check in Moments Editor (GUI)!");
	}
	int m = 0;
	for (const auto& vec : parsedMoments) {
		for (const auto& newMoment : vec) {
			moments.push_back(newMoment);
		}
		prg.SetProgress(0.5 + (double)m++ / (double)parsedMoments.size());
	}
	std::sort(moments.begin(), moments.end(), [](const Moment& a, const Moment& b) {
		return a.time - .5 * a.window < b.time - .5 * b.window; // This will sort in ascending order based on the 'startTime' member
		});

}