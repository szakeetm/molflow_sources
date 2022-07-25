//
// Created by pascal on 7/19/22.
//

#ifndef MOLFLOW_PROJ_TDMOMENTS_H
#define MOLFLOW_PROJ_TDMOMENTS_H


#include "TDTypes.h"
#include <vector>

namespace MFTD {
    int AddMoment(const std::vector<Moment> &newMoments, std::vector<Moment> &intervals);
    void
    momentIntervalReader(const std::vector<UserMoment> &userMoments, std::vector<MomentInterval> &uIntervals);
    void momentsReader(const std::vector<UserMoment> &userMoments, std::vector<Moment> &intervals);
    void parseMoments(const std::vector<Moment> &userIntervals,
                             std::vector<Moment> &outputIntervals);
};


#endif //MOLFLOW_PROJ_TDMOMENTS_H
