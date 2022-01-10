//
// Created by Pascal Baehr on 20.07.20.
//

#ifndef MOLFLOW_PROJ_TIMEMOMENTS_H
#define MOLFLOW_PROJ_TIMEMOMENTS_H

#include <vector>
#include "MolflowTypes.h"

class TimeMoments {
public:
    static int CheckIntervalOverlap(const std::vector<Moment>& vecA, const std::vector<Moment>& vecB);
    static std::pair<int, int> CheckIntervalOverlap(const std::vector<std::vector<Moment>>& vecParsedMoments);
    static std::vector<Moment> ParseMoment(const std::string& userInput, double timeWindow);
    static int
    ParseAndCheckUserMoments(std::vector<Moment> *moments, std::vector<UserMoment> *userMoments,
                             double *progress);

    static int AddMoment(std::vector<Moment> *moments, std::vector<Moment> newMoments); //Adds a time serie to moments and returns the number of elements
};


#endif //MOLFLOW_PROJ_TIMEMOMENTS_H
