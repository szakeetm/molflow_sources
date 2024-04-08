

#ifndef MOLFLOW_PROJ_TIMEMOMENTS_H
#define MOLFLOW_PROJ_TIMEMOMENTS_H

#include <vector>
#include "MolflowTypes.h"
#include <optional>
class GLProgress_Abstract;

class TimeMoments {
public:
    static std::optional<std::pair<int, int>> HasIntervalOverlap(const std::vector<std::vector<Moment>>& vecParsedMoments);
    static std::vector<Moment> ParseUserMoment(const UserMoment& input);
    static void ParseAndCheckUserMoments(std::vector<Moment>& moments, const std::vector<UserMoment>& userMoments,
        GLProgress_Abstract& prg);
};


#endif //MOLFLOW_PROJ_TIMEMOMENTS_H
