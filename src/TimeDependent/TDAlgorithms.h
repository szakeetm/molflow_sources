//
// Created by pascal on 7/9/22.
//

#ifndef MOLFLOW_PROJ_TDALGORITHMS_H
#define MOLFLOW_PROJ_TDALGORITHMS_H

#include "TDTypes.h"
#include <vector>



namespace MFTD {
    int LookupMomentIndex(double key, const std::vector<std::pair<double, double>>& moments, size_t startIndex);
    int LookupMomentIndex(double key, const std::vector<Moment>& moments);
    int quadraticSearch(double key, const std::vector<Moment>& moments);
    int interpolationSearch(double key, const std::vector<Moment>& moments);
    int interpolationSearch(double key, const std::vector<Moment>& moments, int startIndex);
    int jumpSearchProg(const std::vector<Moment>& arr, double noToSearch, int ArrayLim);
    int jumpSearchProg(const std::vector<Moment>& arr, double noToSearch, int ArrayLim, int startIndex);
    int calcSearch(double key, const std::vector<Moment>& moments, const std::vector<MomentInterval>& userMoments);
} // MFTD

#endif //MOLFLOW_PROJ_TDALGORITHMS_H
