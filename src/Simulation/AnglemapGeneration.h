

#ifndef MOLFLOW_PROJ_ANGLEMAPGENERATION_H
#define MOLFLOW_PROJ_ANGLEMAPGENERATION_H


#include <MolflowTypes.h>
#include "MolflowSimFacet.h" // Anglemap

namespace AnglemapGeneration {
    double GetTheta(const double thetaIndex, const AnglemapParams &anglemapParams);

    double GetPhi(const double phiIndex, const AnglemapParams &anglemapParams);

    double GetPhiNormalizedPdfValue(const double thetaIndex, const int phiLowerIndex,
                                 const AnglemapParams &anglemapParams, const Anglemap & anglemap);

    double GetPhiCDFValue(const double thetaIndex, const int phiLowerIndex,
                                 const AnglemapParams &anglemapParams, const Anglemap &anglemap);

    double GetPhiCDFSum(const double thetaIndex, const AnglemapParams &anglemapParams,
                               const Anglemap &anglemap);

    std::tuple<double, int, double>
    GenerateThetaFromAngleMap(const AnglemapParams &anglemapParams, Anglemap &anglemap,
                              double lookupValue);

    double GeneratePhiFromAngleMap(const int thetaLowerIndex, const double thetaOvershoot,
                                   const AnglemapParams &anglemapParams,
                                   Anglemap &anglemap, double lookupValue);
};


#endif //MOLFLOW_PROJ_ANGLEMAPGENERATION_H
