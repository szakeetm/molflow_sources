//
// Created by pascal on 7/19/22.
//

#ifndef MOLFLOW_PROJ_TDTYPES_H
#define MOLFLOW_PROJ_TDTYPES_H

#include <string>

namespace MFTD {
    typedef std::pair<double, double> Moment;
    typedef std::pair<std::string, double> UserMoment;
    struct MomentInterval {
        double start;
        double interval;
        double end;
        double timeWindow;
        size_t startIndex;
    };
}

#endif //MOLFLOW_PROJ_TDTYPES_H
