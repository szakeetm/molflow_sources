

#ifndef MOLFLOW_PROJ_CDFGENERATION_H
#define MOLFLOW_PROJ_CDFGENERATION_H

#include <list>
#include "MolflowSimGeom.h"

namespace CDFGeneration {
    // CDF
    int GetCDFId(const std::vector<double> &temperatureList, double temperature);
    std::pair<int,std::vector<IntegratedVelocityEntry>>GenerateNewCDF(std::vector<double> &temperatureList, double temperature, double gasMass);
    std::vector<IntegratedVelocityEntry> Generate_CDF(double gasTempKelvins, double gasMassGramsPerMol);
};


#endif //MOLFLOW_PROJ_CDFGENERATION_H
