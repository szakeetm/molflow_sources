//
// Created by pascal on 2/9/21.
//

#ifndef MOLFLOW_PROJ_CDFGENERATION_H
#define MOLFLOW_PROJ_CDFGENERATION_H

#include <list>
#include "GeometrySimu.h"

namespace CDFGeneration {
    // CDF
    int GetCDFId(const std::vector<double> &temperatureList, double temperature);
    std::pair<int,std::vector<CDF_p>>GenerateNewCDF(std::vector<double> &temperatureList, double temperature, double gasMass);
    std::vector<CDF_p> Generate_CDF(double gasTempKelvins, double gasMassGramsPerMol, size_t size);
};


#endif //MOLFLOW_PROJ_CDFGENERATION_H
